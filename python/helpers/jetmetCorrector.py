import os
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoHRTTools.helpers.jetSmearingHelper import jetSmearer
from PhysicsTools.NanoAODTools.postprocessing.modules.jme.JetReCalibrator import JetReCalibrator

import logging
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] %(levelname)s: %(message)s')


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

    def __init__(self, jetType="AK4PFchs", jec=False, jes=None, jes_source=None, jer='nominal', jmr=None, met_unclustered=None):
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

        self.jetType = jetType
        self.jec = jec
        self.jes = jes
        self.jes_source = '' if jes_source is None else jes_source
        self.jer = jer
        self.jmr = jmr
        self.met_unclustered = met_unclustered

    def beginJob(self):
        # set up JEC
        if self.jec or self.jes in ['up', 'down']:
            for library in ["libCondFormatsJetMETObjects", "libPhysicsToolsNanoAODTools"]:
                if library not in ROOT.gSystem.GetLibraries():
                    logging.info("Load Library '%s'" % library.replace("lib", ""))
                    ROOT.gSystem.Load(library)

            self.jesInputFilePath = os.environ['CMSSW_BASE'] + "/src/PhysicsTools/NanoHRTTools/data/2016/jme/"

        # updating JEC
        if self.jec:
            logging.info('Loading JEC parameters...')
            self.jetReCalibratorMC = None
#             self.jetReCalibratorMC = JetReCalibrator(globalTag='Summer16_23Sep2016V4_MC',
#                                                    jetFlavour=self.jetType,
#                                                    doResidualJECs=False,
#                                                    jecPath=self.jesInputFilePath + 'Summer16_23Sep2016V4_MC/',
#                                                    calculateSeparateCorrections=False,
#                                                    calculateType1METCorrection=False)
            self.jetReCalibratorDATA_BCD = JetReCalibrator(globalTag='Summer16_23Sep2016BCDV4_DATA',
                                                   jetFlavour=self.jetType,
                                                   doResidualJECs=True,
                                                   jecPath=self.jesInputFilePath + 'Summer16_23Sep2016BCDV4_DATA/',
                                                   calculateSeparateCorrections=False,
                                                   calculateType1METCorrection=False)
            self.jetReCalibratorDATA_EF = JetReCalibrator(globalTag='Summer16_23Sep2016EFV4_DATA',
                                                   jetFlavour=self.jetType,
                                                   doResidualJECs=True,
                                                   jecPath=self.jesInputFilePath + 'Summer16_23Sep2016EFV4_DATA/',
                                                   calculateSeparateCorrections=False,
                                                   calculateType1METCorrection=False)
            self.jetReCalibratorDATA_G = JetReCalibrator(globalTag='Summer16_23Sep2016GV4_DATA',
                                                   jetFlavour=self.jetType,
                                                   doResidualJECs=True,
                                                   jecPath=self.jesInputFilePath + 'Summer16_23Sep2016GV4_DATA/',
                                                   calculateSeparateCorrections=False,
                                                   calculateType1METCorrection=False)
            self.jetReCalibratorDATA_H = JetReCalibrator(globalTag='Summer16_23Sep2016HV4_DATA',
                                                   jetFlavour=self.jetType,
                                                   doResidualJECs=True,
                                                   jecPath=self.jesInputFilePath + 'Summer16_23Sep2016HV4_DATA/',
                                                   calculateSeparateCorrections=False,
                                                   calculateType1METCorrection=False)

        # JES uncertainty
        if self.jes in ['up', 'down']:
            if not self.jes_source:
                # total unc.
                self.jesUncertaintyInputFileName = "Summer16_23Sep2016V4_MC_Uncertainty_" + self.jetType + ".txt"
            else:
                # unc. by source
                self.jesUncertaintyInputFileName = "Summer16_23Sep2016V4_MC_UncertaintySources_" + self.jetType + ".txt"

            pars = ROOT.JetCorrectorParameters(os.path.join(self.jesInputFilePath, self.jesUncertaintyInputFileName), self.jes_source)
            self.jesUncertainty = ROOT.JetCorrectionUncertainty(pars)

        # set up JER
        self.jetSmearer = None
        if self.jer is not None or self.jmr is not None:
            self.jetSmearer = jetSmearer(globalTag='2016',
                                         jetType=self.jetType,
                                         jerInputFileName="Spring16_25nsV10a_MC_PtResolution_%s.txt" % self.jetType,
                                         jerUncertaintyInputFileName="Spring16_25nsV10a_MC_SF_%s.txt" % self.jetType
                                         )
            self.jetSmearer.jerInputFilePath = os.environ['CMSSW_BASE'] + "/src/PhysicsTools/NanoHRTTools/data/2016/jme/"
            self.jetSmearer.beginJob()

    def setSeed(self, seed):
        if self.jetSmearer is not None:
            self.jetSmearer.setSeed(seed)

    def correctJetAndMET(self, jets, met=None, rho=None, genjets=[], isMC=True, runNumber=None):

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
                if runNumber >= 272007 and runNumber <= 276811:
                    jetReCalibrator = self.jetReCalibratorDATA_BCD
                elif runNumber >= 276831 and runNumber <= 278801:
                    jetReCalibrator = self.jetReCalibratorDATA_EF
                elif runNumber >= 278802 and runNumber <= 280385:
                    jetReCalibrator = self.jetReCalibratorDATA_G
                elif runNumber >= 280919 and runNumber <= 284044:
                    jetReCalibrator = self.jetReCalibratorDATA_H
                else:
                    raise RuntimeError("Run %d out of range" % runNumber)
            for j in jets:
                newpt = jetReCalibrator.correct(j, rho)
                j.mass = newpt / j.pt * j.mass
                j.pt = newpt

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
