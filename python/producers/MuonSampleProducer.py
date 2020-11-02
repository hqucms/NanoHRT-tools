import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.tools import deltaPhi, deltaR, closest

from PhysicsTools.NanoHRTTools.producers.HRTBaseProducer import HRTBaseProducer
from PhysicsTools.NanoHRTTools.helpers.triggerHelper import passTrigger

import logging
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] %(levelname)s: %(message)s')


class MuonSampleProducer(HRTBaseProducer):

    def __init__(self, **kwargs):
        super(MuonSampleProducer, self).__init__(channel='muon', **kwargs)

    def beginJob(self):
        super(MuonSampleProducer, self).beginJob()

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        super(MuonSampleProducer, self).beginFile(inputFile, outputFile, inputTree, wrappedOutputTree)

        ## trigger variables
        self.out.branch("passMuTrig", "O")

        ## event variables
        self.out.branch("muon_pt", "F")
        self.out.branch("muon_eta", "F")
        self.out.branch("muon_miniIso", "F")
        self.out.branch("leptonicW_pt", "F")

    def prepareEvent(self, event):

        logging.debug('processing event %d' % event.event)

        ## met selection
        if event.met.pt < 50.0:
            return False

        ## muon selection
        event._allMuons = Collection(event, "Muon")
        event.muons = []
        for muon in event._allMuons:
            if muon.pt > 55 and abs(muon.eta) < 2.4 and muon.tightId and abs(muon.dxy) < 0.2 and abs(muon.dz) < 0.5:
                if muon.miniPFRelIso_all < 0.10:
                    event.muons.append(muon)
        if len(event.muons) != 1:
            return False

        # #leptonic W pt cut
        event.mu = event.muons[0]
        event.leptonicW = event.mu.p4() + event.met.p4()
        if event.leptonicW.Pt() < 100.0:
            return False

        ## b-tag AK4 jet selection
        event.bjets = []
        for j in event._allJets:
            if not (j.pt > 25.0 and abs(j.eta) < 2.4 and (j.jetId & 2)):
                continue
            if j.btagDeepB > self.DeepCSV_WP_M and abs(deltaPhi(j, event.muons[0])) < 2.0:
                event.bjets.append(j)

        if len(event.bjets) < 1:
            return False

        # # selection on AK8 jets
        event.ak8jets = []
        for fj in event._allAK8jets:
            if not (fj.pt > 150 and abs(fj.eta) < 2.4 and (fj.jetId & 2)):
                continue
            if abs(deltaPhi(fj, event.muons[0])) > 2.0:
                event.ak8jets.append(fj)

        if len(event.ak8jets) < 1:
            return False

        ## return True if passes selection
        return True

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""

        self.correctJetsAndMET(event)

        if self.prepareEvent(event) is False:
            return False

        self.loadGenHistory(event)

        # fill
        self.out.fillBranch("passMuTrig", passTrigger(event, ['HLT_Mu50', 'HLT_TkMu50']))
        self.out.fillBranch("muon_pt", event.mu.pt)
        self.out.fillBranch("muon_eta", event.mu.eta)
        self.out.fillBranch("muon_miniIso", event.mu.miniPFRelIso_all)
        self.out.fillBranch("leptonicW_pt", event.leptonicW.Pt())

        self.fillBaseEventInfo(event)
        self.fillFatJetInfo(event, fillGenMatching=True)

        return True


# define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed
MuonTree_2016 = lambda: MuonSampleProducer(year=2016)
MuonTree_2017 = lambda: MuonSampleProducer(year=2017)
MuonTree_2018 = lambda: MuonSampleProducer(year=2018)
