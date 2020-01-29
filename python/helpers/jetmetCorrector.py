import os, tarfile, tempfile, shutil
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoHRTTools.helpers.jetSmearingHelper import jetSmearer, find_and_extract_tarball
from PhysicsTools.NanoAODTools.postprocessing.modules.jme.JetReCalibrator import JetReCalibrator


def _sf(vals, syst='nominal'):
    if syst == 'nominal':
        return vals[0]
    elif syst == 'up':
        return vals[1]
    elif syst == 'down':
        return vals[2]
    else:
        raise ValueError('Invalid syst type: %s' % str(syst))


def selectJetsForMET(jet):
    return jet.pt > 15 and jet.chEmEF + jet.neEmEF < 0.9


class JetMETCorrector(object):

    def __init__(self, year, jetType="AK4PFchs", jec=False, jes=None, jes_source=None, jer='nominal', jmr=None, met_unclustered=None):
        '''
        jec: re-apply jet energy correction (True|False)
        jes: Jet energy scale options
            - None: do nothing
            - 'nominal': update JES central values
            - 'up', 'down': up/down variation, using the total uncertainty
            - '[UncertaintySource]_(up|down)': up/down variation, using per source unc.
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
        self.jer = jer
        self.jmr = jmr
        self.met_unclustered = met_unclustered

        # set up tags for each year
        if self.year == 2016:
            self.globalTag = 'Summer16_07Aug2017_V11_MC'
            self.jerTag = 'Summer16_25nsV1_MC'
            self.dataTags = (
                # set the name of the tarball with a dummy run number
                (0, 'Summer16_07Aug2017_V11_DATA'),
                # (start run number (inclusive), 'tag name')
                (272007, 'Summer16_07Aug2017BCD_V11_DATA'),
                (276831, 'Summer16_07Aug2017EF_V11_DATA'),
                (278820, 'Summer16_07Aug2017GH_V11_DATA'),
            )
        elif self.year == 2017:
            self.globalTag = 'Fall17_17Nov2017_V32_MC'
            self.jerTag = 'Fall17_V3_MC'
            self.dataTags = (
                # set the name of the tarball with a dummy run number
                (0, 'Fall17_17Nov2017_V32_DATA'),
                # (start run number (inclusive), 'tag name')
                (297020, 'Fall17_17Nov2017B_V32_DATA'),
                (299337, 'Fall17_17Nov2017C_V32_DATA'),
                (302030, 'Fall17_17Nov2017DE_V32_DATA'),
                (304911, 'Fall17_17Nov2017F_V32_DATA'),
            )
        elif self.year == 2018:
            self.globalTag = 'Autumn18_V19_MC'
            self.jerTag = 'Autumn18_V7b_MC'
            self.dataTags = (
                # set the name of the tarball with a dummy run number
                (0, 'Autumn18_V19_DATA'),
                # (start run number (inclusive), 'tag name')
                (315252, 'Autumn18_RunA_V19_DATA'),
                (316998, 'Autumn18_RunB_V19_DATA'),
                (319313, 'Autumn18_RunC_V19_DATA'),
                (320394, 'Autumn18_RunD_V19_DATA'),
            )
        else:
            raise RuntimeError('Invalid year: %s' % (str(self.year)))

    def beginJob(self):
        # set up JEC
        if self.jec or self.jes in ['up', 'down']:
            for library in ["libCondFormatsJetMETObjects", "libPhysicsToolsNanoAODTools"]:
                if library not in ROOT.gSystem.GetLibraries():
                    print("Load Library '%s'" % library.replace("lib", ""))
                    ROOT.gSystem.Load(library)

            self.jesInputFilePath = tempfile.mkdtemp()
            # extract the MC and unc files
            find_and_extract_tarball(self.globalTag, self.jesInputFilePath)

        # updating JEC
        if self.jec:
            self.jetReCalibratorMC = JetReCalibrator(globalTag=self.globalTag,
                                                     jetFlavour=self.jetType,
                                                     doResidualJECs=False,
                                                     jecPath=self.jesInputFilePath,
                                                     calculateSeparateCorrections=False,
                                                     calculateType1METCorrection=False)
            self.jetReCalibratorsDATA = {}
            for iov, tag in self.dataTags:
                find_and_extract_tarball(tag, self.jesInputFilePath)
                if iov > 0:
                    self.jetReCalibratorsDATA[tag] = JetReCalibrator(globalTag=tag,
                                                                     jetFlavour=self.jetType,
                                                                     doResidualJECs=True,
                                                                     jecPath=self.jesInputFilePath,
                                                                     calculateSeparateCorrections=False,
                                                                     calculateType1METCorrection=False)

        # JES uncertainty
        if self.jes in ['up', 'down']:
            if not self.jes_source:
                # total unc.
                self.jesUncertaintyInputFileName = self.globalTag + "_Uncertainty_" + self.jetType + ".txt"
            else:
                # unc. by source
                self.jesUncertaintyInputFileName = self.globalTag + "_UncertaintySources_" + self.jetType + ".txt"

            pars = ROOT.JetCorrectorParameters(os.path.join(self.jesInputFilePath, self.jesUncertaintyInputFileName), self.jes_source)
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

    def correctJetAndMET(self, jets, met=None, rho=None, genjets=[], isMC=True, runNumber=None):
        # for MET correction, use 'Jet' (corr_pt>15) and 'CorrT1METJet' (corr_pt<15) collections
        # Type-1 MET correction: https://github.com/cms-sw/cmssw/blob/master/JetMETCorrections/Type1MET/interface/PFJetMETcorrInputProducerT.h
        # FIXME: FIX MET correction
        # prepare to propogate corrections to MET
        if met is not None:
            sumJetsP4Orig = ROOT.TLorentzVector()
            for j in jets:
                j._forMET = selectJetsForMET(j)
                if j._forMET:
                    sumJetsP4Orig += j.p4()

        # updating JEC
        if self.jec:
            if isMC:
                jetReCalibrator = self.jetReCalibratorMC
            else:
                tag = next(run for iov, run in reversed(self.dataTags) if iov <= runNumber)
                jetReCalibrator = self.jetReCalibratorsDATA[tag]
            for j in jets:
                j.pt, j.mass = jetReCalibrator.correct(j, rho)

        # JES uncertainty
        if isMC and self.jes in ['up', 'down']:
            for j in jets:
                self.jesUncertainty.setJetPt(j.pt)
                self.jesUncertainty.setJetEta(j.eta)
                delta = self.jesUncertainty.getUncertainty(True)
                sf = 1 + delta if self.jes == 'up' else 1 - delta
                j.pt *= sf
                j.mass *= sf

        # jer smearing
        if isMC and self.jer is not None:
            for j in jets:
                jersf = _sf(self.jetSmearer.getSmearValsPt(j, genjets, rho), self.jer)
                j.pt *= jersf
                j.mass *= jersf

        # propogate to MET
        if met is not None:
            sumJetsP4New = sum([j.p4() for j in jets if j._forMET], ROOT.TLorentzVector())
            newMET = met.p4() + (sumJetsP4Orig - sumJetsP4New)
            met.pt, met.phi = newMET.Pt(), newMET.Phi()

        # MET unclustered energy
        if isMC and met is not None and self.met_unclustered:
            delta = ROOT.TLorentzVector(met.MetUnclustEnUpDeltaX, met.MetUnclustEnUpDeltaY, 0, 0)
            if self.met_unclustered == 'up':
                newMET = met.p4() + delta
            elif self.met_unclustered == 'down':
                newMET = met.p4() - delta
            met.pt, met.phi = newMET.Pt(), newMET.Phi()

    def smearJetMass(self, jets, gensubjets=[], isMC=True, runNumber=None):
        # jmr smearing (mass resolution)
        if isMC and self.jmr is not None:
            for j in jets:
                jmrsf = _sf(self.jetSmearer.getSmearValsM(j, gensubjets), self.jmr)
                j.mass *= jmrsf
                j.msoftdrop = j.mass
