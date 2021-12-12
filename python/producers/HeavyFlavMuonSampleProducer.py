from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection

from ..helpers.utils import deltaPhi, polarP4
from ..helpers.triggerHelper import passTrigger
from .HeavyFlavBaseProducer import HeavyFlavBaseProducer


class MuonSampleProducer(HeavyFlavBaseProducer):

    def __init__(self, **kwargs):
        super(MuonSampleProducer, self).__init__(channel='muon', **kwargs)

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        super(MuonSampleProducer, self).beginFile(inputFile, outputFile, inputTree, wrappedOutputTree)

        # trigger variables
        self.out.branch("passMuTrig", "O")

        # event variables
        self.out.branch("muon_pt", "F")
        self.out.branch("muon_eta", "F")
        self.out.branch("muon_miniIso", "F")
        self.out.branch("leptonicW_pt", "F")

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""

        # muon selection
        event._allMuons = Collection(event, "Muon")
        event.muons = [mu for mu in event._allMuons if mu.pt > 55 and abs(mu.eta) < 2.4 and abs(
            mu.dxy) < 0.2 and abs(mu.dz) < 0.5 and mu.tightId and mu.miniPFRelIso_all < 0.10]
        if len(event.muons) != 1:
            return False

        self.selectLeptons(event)
        self.correctJetsAndMET(event)

        # met selection
        if event.met.pt < 50.0:
            return False

        # leptonic W pt cut
        event.mu = event.muons[0]
        event.leptonicW = polarP4(event.mu) + event.met.p4()
        if event.leptonicW.Pt() < 100.0:
            return False

        # at least one b-jet, in the same hemisphere of the muon
        event.bjets = [j for j in event.ak4jets if j.btagDeepFlavB > self.DeepJet_WP_M and
                       abs(deltaPhi(j, event.mu)) < 2]
        if len(event.bjets) == 0:
            return False

        # require fatjet away from the muon
        probe_jets = [fj for fj in event.fatjets if abs(deltaPhi(fj, event.mu)) > 2]
        if len(probe_jets) == 0:
            return False

        probe_jets = probe_jets[:1]
        self.loadGenHistory(event, probe_jets)
        self.evalTagger(event, probe_jets)
        self.evalMassRegression(event, probe_jets)

        # fill output branches
        self.fillBaseEventInfo(event)
        self.fillFatJetInfo(event, probe_jets)

        # fill
        self.out.fillBranch("passMuTrig", passTrigger(event, ['HLT_Mu50', 'HLT_TkMu50']))
        self.out.fillBranch("muon_pt", event.mu.pt)
        self.out.fillBranch("muon_eta", event.mu.eta)
        self.out.fillBranch("muon_miniIso", event.mu.miniPFRelIso_all)
        self.out.fillBranch("leptonicW_pt", event.leptonicW.Pt())

        return True


# define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed
def MuonTree_2016(): return MuonSampleProducer(year=2016)
def MuonTree_2017(): return MuonSampleProducer(year=2017)
def MuonTree_2018(): return MuonSampleProducer(year=2018)
