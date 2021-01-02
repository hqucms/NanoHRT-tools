from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection

from ..helpers.utils import deltaPhi
from .HeavyFlavBaseProducer import HeavyFlavBaseProducer


class PhotonSampleProducer(HeavyFlavBaseProducer):

    def __init__(self, **kwargs):
        super(PhotonSampleProducer, self).__init__(channel='photon', **kwargs)

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        super(PhotonSampleProducer, self).beginFile(inputFile, outputFile, inputTree, wrappedOutputTree)

        # trigger variables
        self.out.branch("passTrigPhoton", "O")

        # photons
        self.out.branch("nphotons", "I")
        self.out.branch("pho_1_pt", "F")
        self.out.branch("pho_1_eta", "F")
        self.out.branch("pho_1_phi", "F")

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""

        # select leading photon
        event._allPhotons = Collection(event, "Photon")
        event.photons = [pho for pho in event._allPhotons if pho.pt > 200 and abs(
            pho.eta) < 2.5 and pho.cutBased >= 2 and pho.electronVeto]  # medium ID
        if len(event.photons) < 1:
            return False

        self.selectLeptons(event)
        self.correctJetsAndMET(event)

        # require jet and photon to be back-to-back
        probe_jets = [fj for fj in event.fatjets if abs(deltaPhi(event.photons[0], fj)) > 2]
        if len(probe_jets) == 0:
            return False
        probe_jets = probe_jets[:1]  # only consider the leading fatjet

        if self._opts['sfbdt_threshold'] > -99:
            # require: 2 subjets, 50 < sdmass < 200
            probe_jets = [fj for fj in probe_jets if len(fj.subjets) == 2 and fj.msoftdrop > 50 and fj.msoftdrop < 200]
            if len(probe_jets) == 0:
                return False

            self.selectSV(event)
            if len(event.secondary_vertices) < 2:
                return False
            self.matchSVToFatJets(event, probe_jets)

            # require: sfBDT > 0.5 (and >=1 SV matched to each subjet -- already included in the sfBDT requirement)
            probe_jets = [fj for fj in probe_jets if fj.sfBDT >= self._opts['sfbdt_threshold']]
            if len(probe_jets) == 0:
                return False

        self.loadGenHistory(event, probe_jets)
        self.evalTagger(event, probe_jets)
        self.evalMassRegression(event, probe_jets)

        # fill output branches
        self.fillBaseEventInfo(event)
        self.fillFatJetInfo(event, probe_jets)

        if self.year == 2016:
            self.out.fillBranch("passTrigPhoton", event.HLT_Photon175)
        else:
            self.out.fillBranch("passTrigPhoton", event.HLT_Photon200)

        # photon variables
        self.out.fillBranch("nphotons", len(event.photons))
        self.out.fillBranch("pho_1_pt", event.photons[0].pt)
        self.out.fillBranch("pho_1_eta", event.photons[0].eta)
        self.out.fillBranch("pho_1_phi", event.photons[0].phi)

        return True


# define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed
def PhotonTree_2016(): return PhotonSampleProducer(year=2016)
def PhotonTree_2017(): return PhotonSampleProducer(year=2017)
def PhotonTree_2018(): return PhotonSampleProducer(year=2018)
