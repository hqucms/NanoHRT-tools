import os
import tempfile
import shutil
import itertools
import logging
import numpy as np
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

from .utils import polarP4, p4, configLogger
from .jetSmearingHelper import jetSmearer, find_and_extract_tarball

logger = logging.getLogger('jme')
configLogger('jme', loglevel=logging.INFO)


def rndSeed(event, jets, extra=0):
    seed = (event.run << 20) + (event.luminosityBlock << 10) + event.event + extra
    if len(jets) > 0:
        seed += int(jets[0].eta / 0.01)
    return seed


def _sf(vals, syst='nominal'):
    if syst == 'nominal':
        return vals[0]
    elif syst == 'up':
        return vals[1]
    elif syst == 'down':
        return vals[2]
    else:
        raise ValueError('Invalid syst type: %s' % str(syst))


class JetCorrector(object):

    def __init__(self, globalTag, jetType, jecPath, applyResidual=True):
        self.jecLevels = ['L2Relative', 'L3Absolute'] if 'Puppi' in jetType else [
            'L1FastJet', 'L2Relative', 'L3Absolute']
        if applyResidual:
            self.jecLevels += ['L2L3Residual']
        self.vPar = ROOT.vector(ROOT.JetCorrectorParameters)()
        logger.info('Init JetCorrector: %s, %s, %s', globalTag, jetType, str(self.jecLevels))
        for level in self.jecLevels:
            self.vPar.push_back(ROOT.JetCorrectorParameters(os.path.join(
                jecPath, "%s_%s_%s.txt" % (globalTag, level, jetType)), ""))
        self.corrector = ROOT.FactorizedJetCorrector(self.vPar)

    def getCorrection(self, jet, rho, level=None):
        try:
            raw_pt = jet.rawP4.pt()
        except RuntimeError:
            raw_pt = jet.pt * (1. - jet.rawFactor)
        self.corrector.setJetPt(raw_pt)
        self.corrector.setJetPhi(jet.phi)
        self.corrector.setJetEta(jet.eta)
        self.corrector.setRho(rho)
        try:
            self.corrector.setJetA(jet.area)
        except RuntimeError:
            pass
        if level is None:
            return self.corrector.getCorrection()
        else:
            idx = self.jecLevels.index(level)
            return self.corrector.getSubCorrections()[idx]


class JetMETCorrector(object):

    def __init__(
            self, year, jetType="AK4PFchs", jec=False, jes=None, jes_source=None, jes_uncertainty_file_prefix=None,
            jer='nominal', jmr=None, met_unclustered=None, smearMET=True, applyHEMUnc=False):
        '''
        jec: re-apply jet energy correction (True|False)
        jes: Jet energy scale options
            - None: do nothing
            - 'nominal': update JES central values
            - 'up', 'down': up/down variation, using the total uncertainty
            - '[UncertaintySource]_(up|down)': up/down variation, using per source unc.
        jes_uncertainty_file_prefix: Prefix of the JES uncertainty file, use `None` for the full set, 'Regrouped(V2)_' for the reduced set.
        jer: Jet energy resolution options
            - None: do nothing
            - 'nominal': apply nominal smearing
            - 'up', 'down': up/down variation
        jmr: Jet mass resolution options
            - None: do nothing
            - 'nominal': apply nominal smearing
            - 'up', 'down': up/down variation
        met_unclustered: MET unclustered energy options
            - None: do nothing
            - 'up', 'down': up/down variation of the unclustered energy
        '''

        self.year = year
        self.jetType = jetType
        self.jec = jec
        self.jes = jes
        self.jes_source = '' if jes_source is None else jes_source
        self.jes_uncertainty_file_prefix = '' if jes_uncertainty_file_prefix is None else jes_uncertainty_file_prefix
        self.jer = jer
        self.jmr = jmr
        self.met_unclustered = met_unclustered
        self.correctMET = (jetType == 'AK4PFchs')  # FIXME
        self.smearMET = smearMET
        self.applyHEMUnc = applyHEMUnc

        self.excludeJetsForMET = None

        # set up tags for each year
        if self.year == 2015:
            # hack, actually UL2016 preVFP (APV)
            self.globalTag = 'Summer19UL16APV_V7_MC'
            self.jerTag = 'Summer20UL16APV_JRV3_MC'
            self.dataTags = (
                # set the name of the tarball with a dummy run number
                (0, 'Summer19UL16APV_V7_DATA'),
                # (start run number (inclusive), 'tag name')
                (272007, 'Summer19UL16APV_RunBCD_V7_DATA'),
                (276831, 'Summer19UL16APV_RunEF_V7_DATA'),
            )
        elif self.year == 2016:
            # hack, actually UL2016 postVFP
            self.globalTag = 'Summer19UL16_V7_MC'
            self.jerTag = 'Summer20UL16_JRV3_MC'
            self.dataTags = (
                # set the name of the tarball with a dummy run number
                (0, 'Summer19UL16_V7_DATA'),
                # (start run number (inclusive), 'tag name')
                (277772, 'Summer19UL16_RunFGH_V7_DATA'),
            )
        elif self.year == 2017:
            self.globalTag = 'Summer19UL17_V6_MC'
            self.jerTag = 'Summer19UL17_JRV2_MC'
            self.dataTags = (
                # set the name of the tarball with a dummy run number
                (0, 'Summer19UL17_V6_DATA'),
                # (start run number (inclusive), 'tag name')
                (297020, 'Summer19UL17_RunB_V6_DATA'),
                (299337, 'Summer19UL17_RunC_V6_DATA'),
                (302030, 'Summer19UL17_RunD_V6_DATA'),
                (303435, 'Summer19UL17_RunE_V6_DATA'),
                (304911, 'Summer19UL17_RunF_V6_DATA'),
            )
        elif self.year == 2018:
            self.globalTag = 'Summer19UL18_V5_MC'
            self.jerTag = 'Summer19UL18_JRV2_MC'
            self.dataTags = (
                # set the name of the tarball with a dummy run number
                (0, 'Summer19UL18_V5_DATA'),
                # (start run number (inclusive), 'tag name')
                (315252, 'Summer19UL18_RunA_V5_DATA'),
                (316998, 'Summer19UL18_RunB_V5_DATA'),
                (319313, 'Summer19UL18_RunC_V5_DATA'),
                (320394, 'Summer19UL18_RunD_V5_DATA'),
            )
        else:
            raise RuntimeError('Invalid year: %s' % (str(self.year)))

    def beginJob(self):
        # set up JEC
        if self.jec or self.jes in ['up', 'down'] or self.correctMET:
            for library in ["libCondFormatsJetMETObjects", "libPhysicsToolsNanoAODTools"]:
                if library not in ROOT.gSystem.GetLibraries():
                    logger.info("Load Library '%s'" % library.replace("lib", ""))
                    ROOT.gSystem.Load(library)

            self.jesInputFilePath = tempfile.mkdtemp()
            # extract the MC and unc files
            find_and_extract_tarball(self.globalTag, self.jesInputFilePath,
                                     copy_txt_with_prefix=self.jes_uncertainty_file_prefix)

            # updating JEC/re-correct MET
            self.jetCorrectorMC = JetCorrector(globalTag=self.globalTag,
                                               jetType=self.jetType,
                                               jecPath=self.jesInputFilePath,
                                               applyResidual=False)
            self.jetCorrectorsDATA = {}
            for iov, tag in self.dataTags:
                find_and_extract_tarball(tag, self.jesInputFilePath)
                if iov > 0:
                    self.jetCorrectorsDATA[tag] = JetCorrector(globalTag=tag,
                                                               jetType=self.jetType,
                                                               jecPath=self.jesInputFilePath,
                                                               applyResidual=True)

        # JES uncertainty
        if self.jes in ['up', 'down']:
            if not self.jes_source:
                # total unc.
                self.jesUncertaintyInputFileName = self.globalTag + "_Uncertainty_" + self.jetType + ".txt"
            else:
                # unc. by source
                self.jesUncertaintyInputFileName = self.jes_uncertainty_file_prefix + self.globalTag + "_UncertaintySources_" + self.jetType + ".txt"

            pars = ROOT.JetCorrectorParameters(
                os.path.join(self.jesInputFilePath, self.jesUncertaintyInputFileName),
                self.jes_source)
            self.jesUncertainty = ROOT.JetCorrectionUncertainty(pars)

        # set up JER
        self.jetSmearer = None
        if self.jer is not None or self.jmr is not None:
            self.jetSmearer = jetSmearer(self.jerTag, jetType=self.jetType)
            self.jetSmearer.beginJob()

    def endJob(self):
        shutil.rmtree(self.jerInputFilePath)

    def setSeed(self, seed):
        if self.jetSmearer is not None:
            self.jetSmearer.setSeed(seed)

    def calcT1CorrEEFix(self, jet):
        zero = np.zeros(2, dtype='float')
        if self.excludeJetsForMET is None or not self.excludeJetsForMET(jet):
            return zero
        if jet.neEmEF + jet.chEmEF > 0.9:
            return zero
        rawP4 = jet.rawP4 * (1 - jet.muonSubtrFactor)
        corrP4 = rawP4 * jet._jecFactor  # FIXME: _jecFactor and _jecFactorL1 here should be the JEC used in the NanoAOD production
        if corrP4.pt() < 15:
            return zero
        delta = rawP4 * jet._jecFactorL1 - corrP4
        return np.array([delta.px(), delta.py()])

    def calcT1Corr(self, jet):
        zero = np.zeros(2, dtype='float')
        if self.excludeJetsForMET is not None and self.excludeJetsForMET(jet):
            return zero
        if jet.neEmEF + jet.chEmEF > 0.9:
            return zero
        rawP4 = jet.rawP4 * (jet._smearFactorNominal if self.jer and self.smearMET else 1 - jet.muonSubtrFactor)
        corrP4 = rawP4 * jet._jecFactor
        if corrP4.pt() < 15:
            return zero
        delta = rawP4 * jet._jecFactorL1 - corrP4
        if self.jer in ['up', 'down'] or self.jes in ['up', 'down'] or self.applyHEMUnc:
            nominalP4 = jet.rawP4 * jet._jecFactor * (jet._smearFactorNominal if self.jer and self.smearMET else 1)
            if self.jer in ['up', 'down']:
                delta -= nominalP4 * (jet._smearFactor - jet._smearFactorNominal)
            if self.jes in ['up', 'down']:
                delta -= nominalP4 * (jet._jesUncFactor - 1)
            if self.applyHEMUnc:
                delta -= nominalP4 * (jet._HEMUncFactor - 1)
        return np.array([delta.px(), delta.py()])

    def correctJetAndMET(self, jets, lowPtJets=None, met=None, rawMET=None, defaultMET=None,
                         rho=None, genjets=[], isMC=True, runNumber=None):
        # for MET correction, use 'Jet' (corr_pt>15) and 'CorrT1METJet' (corr_pt<15) collections
        # Type-1 MET correction: https://github.com/cms-sw/cmssw/blob/master/JetMETCorrections/Type1MET/interface/PFJetMETcorrInputProducerT.h
        if met is None:
            lowPtJets = []
        else:
            for j in lowPtJets:
                j.pt = j.rawPt
                j.mass = 0
                j.rawFactor = 0
                j.neEmEF = j.chEmEF = 0

        for j in itertools.chain(jets, lowPtJets):
            # set JEC factor ( = corrPt / rawPt)
            j.rawP4 = polarP4(j) * (1. - j.rawFactor)
            j._jecFactor = None
            j._jecFactorL1 = None
            if self.jec or (isMC and met is not None):
                if isMC:
                    jetCorrector = self.jetCorrectorMC
                else:
                    tag = next(run for iov, run in reversed(self.dataTags) if iov <= runNumber)
                    jetCorrector = self.jetCorrectorsDATA[tag]
                j._jecFactor = jetCorrector.getCorrection(j, rho)
                if self.jec:
                    j.pt = j.rawP4.pt() * j._jecFactor
                    j.mass = j.rawP4.mass() * j._jecFactor
                if met is not None:
                    j._jecFactorL1 = jetCorrector.getCorrection(j, rho, 'L1FastJet')

            # set JER factor
            j._smearFactorNominal = 1
            j._smearFactor = 1
            if isMC and self.jer is not None:
                jerFactors = self.jetSmearer.getSmearValsPt(j, genjets, rho)
                j._smearFactorNominal = _sf(jerFactors)
                j._smearFactor = _sf(jerFactors, self.jer)
                j.pt *= j._smearFactor
                j.mass *= j._smearFactor

            # set JES uncertainty ( = varied-Pt / Pt)
            j._jesUncFactor = 1
            if isMC and self.jes in ['up', 'down']:
                self.jesUncertainty.setJetPt(j.pt)  # corrected(+smeared) pt
                self.jesUncertainty.setJetEta(j.eta)
                delta = self.jesUncertainty.getUncertainty(True)
                j._jesUncFactor = 1 + delta if self.jes == 'up' else 1 - delta
                j.pt *= j._jesUncFactor
                j.mass *= j._jesUncFactor

            # set uncertainty due to HEM15/16 issue
            j._HEMUncFactor = 1
            if isMC and self.applyHEMUnc:
                if j.pt > 15 and j.phi > -1.57 and j.phi < -0.87:
                    try:
                        tightId = j.jetId & 2
                    except RuntimeError:
                        tightId = True
                    if tightId:
                        if j.eta > -2.5 and j.eta < -1.3:
                            j._HEMUncFactor = 0.8
                        elif j.eta > -3 and j.eta <= -2.5:
                            j._HEMUncFactor = 0.65
                        j.pt *= j._HEMUncFactor
                        j.mass *= j._HEMUncFactor

            # last thing: calc MET type-1 correction
            j._t1MetDelta = None
            if met is not None:
                j._t1MetDelta = self.calcT1Corr(j) + self.calcT1CorrEEFix(j)

        # correct MET
        if met is not None:
            met_shift = sum([j._t1MetDelta for j in itertools.chain(jets, lowPtJets)])
            # MET unclustered energy
            if isMC and self.met_unclustered:
                delta = np.array([met.MetUnclustEnUpDeltaX, met.MetUnclustEnUpDeltaY])
            if self.met_unclustered == 'up':
                met_shift += delta
            elif self.met_unclustered == 'down':
                met_shift -= delta
            rawMetP4 = p4(rawMET, eta=None, mass=None)
            newMET = rawMetP4 + ROOT.Math.XYZTVector(met_shift[0], met_shift[1], 0, 0)
            if self.excludeJetsForMET is not None:
                newMET += p4(met, eta=None, mass=None) - p4(defaultMET, eta=None, mass=None)
            met.pt, met.phi = newMET.pt(), newMET.phi()

    def smearJetMass(self, jets, gensubjets=[], isMC=True, runNumber=None):
        # jmr smearing (mass resolution)
        if isMC and self.jmr is not None:
            for j in jets:
                jmrsf = _sf(self.jetSmearer.getSmearValsM(j, gensubjets), self.jmr)
                j.mass *= jmrsf
                j.msoftdrop = j.mass
