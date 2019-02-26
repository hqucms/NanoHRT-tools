import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.tools import deltaPhi, deltaR, closest

from PhysicsTools.NanoHRTTools.producers.HRTBaseProducer import HRTBaseProducer

import logging
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] %(levelname)s: %(message)s')


class MuonSampleProducer(HRTBaseProducer):

    def __init__(self, **kwargs):
        #super(MuonSampleProducer, self).__init__(channel='muon', **kwargs)
        super(MuonSampleProducer, self).__init__(channel='muon')

    def beginJob(self):
        super(MuonSampleProducer, self).beginJob()

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        super(MuonSampleProducer, self).beginFile(inputFile, outputFile, inputTree, wrappedOutputTree)

        ## trigger variables
        self.out.branch("passMuTrig", "O")

        ## event variables
        self.out.branch("muon_pt", "F")
        self.out.branch("muon_eta", "F")
        self.out.branch("muon_pTrel", "F")
        self.out.branch("muon_drMuJet", "F")
        self.out.branch("leptonicW_pt", "F")

    def prepareEvent(self, event):

        logging.debug('processing event %d' % event.event)

        ## met selection
        if event.met.pt < 50.0:
            return False

        ## muon selection
        muSubJets = []
        for j in Collection(event, "CustomAK4CHS"):
            if j.pt > 30 and abs(j.eta) < 2.4:
                muSubJets.append(j)

        event._allMuons = Collection(event, "Muon")
        event.muons = []
        for muon in event._allMuons:
            if muon.pt > 55 and abs(muon.eta) < 2.4 and muon.tightId and abs(muon.dxy) < 0.2 and abs(muon.dz) < 0.5:
                j, muon.drjet = closest(muon, muSubJets)
                muon.ptrel = muon.p4().Perp(j.p4().Vect()) if j else 0
#                 if muon.pfRelIso04_all < 0.15:
                if muon.drjet > 0.4 or muon.ptrel > 25:
                    event.muons.append(muon)
        if len(event.muons) != 1:
            return False

        # #leptonic W pt cut
        event.mu = event.muons[0]
        event.leptonicW = event.mu.p4() + event.met.p4()
        if event.leptonicW.Pt() < 250.0:
            return False

        ## b-tag AK4 jet selection
        event.ak4jets = []
        for j in event._allJets:
            if not (j.pt > 25.0 and abs(j.eta) < 2.4 and (j.jetId & 2)):
                continue
            if j.btagCSVV2 > 0.8484 and\
               abs(deltaPhi(j, event.muons[0])) < 2.0:
                event.ak4jets.append(j)

        if len(event.ak4jets) < 1:
            return False

        # # selection on AK8 jets
        event.ak8jets = []
        for fj in event._allAK8jets:
            if not (fj.pt > 200 and abs(fj.eta) < 2.4 and (fj.jetId & 2)):
                continue
            if abs(deltaPhi(fj, event.muons[0])) > 2.0:
                event.ak8jets.append(fj)

        if len(event.ak8jets) < 1:
            return False

        # # selection on CA15 jets
        event.ca15jets = []
        for fj in event._allCA15jets:
            if not (fj.pt > 200 and abs(fj.eta) < 2.4 and (fj.jetId & 2)):
                continue
            if abs(deltaPhi(fj, event.muons[0])) > 2.0:
                event.ca15jets.append(fj)

        if len(event.ca15jets) < 1:
            return False

        ## require the leading ak8 & ca15 jets overlap
        if deltaR(event.ak8jets[0], event.ca15jets[0]) > 0.8:
            return False

        # # selection on HOTVR jets
        event.hotvrjets = []
        for fj in event._allHOTVRjets:
            if not (fj.pt > 200 and abs(fj.eta) < 2.4):
                continue
            if abs(deltaPhi(fj, event.muons[0])) > 2.0:
                event.hotvrjets.append(fj)

        ## return True if passes selection
        return True

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""

        self.correctJetsAndMET(event)

        if self.prepareEvent(event) is False:
            return False

        self.loadGenHistory(event)

        # fill
        passMuTrig_ = False
        try:
            if event.HLT_Mu50:
                passMuTrig_ = True
        except:
            passMuTrig_ = False
        try:
            if event.HLT_TkMu50:
                passMuTrig_ = True
        except:
            passMuTrig_ = False

        self.out.fillBranch("passMuTrig", passMuTrig_)
        self.out.fillBranch("muon_pt", event.mu.pt)
        self.out.fillBranch("muon_eta", event.mu.eta)
        self.out.fillBranch("muon_pTrel", event.mu.ptrel)
        self.out.fillBranch("muon_drMuJet", event.mu.drjet)
        self.out.fillBranch("leptonicW_pt", event.leptonicW.Pt())

        self.fillBaseEventInfo(event)
        self.fillFatJetInfo(event, fillGenMatching=True)

        return True


# define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed
MuonTree = lambda: MuonSampleProducer()
