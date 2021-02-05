from ..helpers.utils import deltaPhi, sumP4, deltaR, minValue
from ..helpers.triggerHelper import passTrigger
from .HeavyFlavBaseProducer import HeavyFlavBaseProducer


class DibosonSampleProducer(HeavyFlavBaseProducer):

    def __init__(self, **kwargs):
        super(DibosonSampleProducer, self).__init__(channel='diboson', **kwargs)

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        super(DibosonSampleProducer, self).beginFile(inputFile, outputFile, inputTree, wrappedOutputTree)

        # trigger variables
        self.out.branch("passTrigEl", "O")
        self.out.branch("passTrigMu", "O")

        # V boson
        self.out.branch("v_pt", "F")
        self.out.branch("v_eta", "F")
        self.out.branch("v_phi", "F")
        self.out.branch("v_mass", "F")

        # leptons
        self.out.branch("lep1_pt", "F")
        self.out.branch("lep1_eta", "F")
        self.out.branch("lep1_phi", "F")
        self.out.branch("lep1_mass", "F")
        self.out.branch("lep1_pdgId", "I")
        self.out.branch("lep2_pt", "F")
        self.out.branch("lep2_eta", "F")
        self.out.branch("lep2_phi", "F")
        self.out.branch("lep2_mass", "F")
        self.out.branch("lep2_pdgId", "I")
        self.out.branch("deltaR_ll", "F")

        # AK4 jets, cleaned vs FatJet
        self.out.branch("n_ak4_cleaned", "I")

        # event level
        self.out.branch("deta_v_fj", "F")
        self.out.branch("dphi_v_fj", "F")
        self.out.branch("min_deta_v_ak4", "F")
        self.out.branch("min_deta_fj_ak4", "F")

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""

        # Z->ll selection
        self.selectLeptons(event)  # select looseLeptons
        if len(event.looseLeptons) < 2:
            return False
        event.selectedLeptons = []  # used for reconstructing the vector boson
        for lep in event.looseLeptons:
            if lep.pt < 20:
                continue
            event.selectedLeptons.append(lep)
            # if abs(lep.pdgId) == 13:
            #     # mu
            #     if lep.pfRelIso04_all < 0.25:
            #         event.selectedLeptons.append(lep)
            # else:
            #     # el
            #     if lep.pfRelIso03_all < 0.15:
            #         event.selectedLeptons.append(lep)
        if len(event.selectedLeptons) < 2:
            return False
        if event.selectedLeptons[0].pdgId + event.selectedLeptons[1].pdgId != 0:
            return False  # OSSF
        vboson = sumP4(event.selectedLeptons[0], event.selectedLeptons[1])
        if vboson.M() < 70 or vboson.M() > 110:
            return False
        if vboson.Pt() < 150:
            return False

        # correct jets before making jet related selections
        self.correctJetsAndMET(event)

        if len(event.fatjets) == 0:
            return False
        fj = event.fatjets[0]
        if abs(deltaPhi(fj.phi, vboson.Phi())) < 2.5:
            return False

        probe_jets = [fj]
        self.loadGenHistory(event, probe_jets)
        self.evalTagger(event, probe_jets)
        self.evalMassRegression(event, probe_jets)

        # fill output branches
        self.fillBaseEventInfo(event)
        self.fillFatJetInfo(event, probe_jets)

        # fill producer-specific branches
        if self.year == 2016:
            self.out.fillBranch("passTrigEl", passTrigger(event, 'HLT_Ele27_WPTight_Gsf'))
            self.out.fillBranch("passTrigMu", passTrigger(event, ['HLT_IsoMu24', 'HLT_IsoTkMu24']))
        elif self.year == 2017:
            self.out.fillBranch("passTrigEl", passTrigger(event, 'HLT_Ele32_WPTight_Gsf_L1DoubleEG'))
            self.out.fillBranch("passTrigMu", passTrigger(event, 'HLT_IsoMu27'))
        elif self.year == 2018:
            self.out.fillBranch("passTrigEl", passTrigger(event, 'HLT_Ele32_WPTight_Gsf'))
            self.out.fillBranch("passTrigMu", passTrigger(event, 'HLT_IsoMu24'))

        # V boson
        self.out.fillBranch("v_pt", vboson.Pt())
        self.out.fillBranch("v_eta", vboson.Eta())
        self.out.fillBranch("v_phi", vboson.Phi())
        self.out.fillBranch("v_mass", vboson.M())

        # leptons
        self.out.fillBranch("lep1_pt", event.selectedLeptons[0].pt)
        self.out.fillBranch("lep1_eta", event.selectedLeptons[0].eta)
        self.out.fillBranch("lep1_phi", event.selectedLeptons[0].phi)
        self.out.fillBranch("lep1_mass", event.selectedLeptons[0].mass)
        self.out.fillBranch("lep1_pdgId", event.selectedLeptons[0].pdgId)
        self.out.fillBranch("lep2_pt", event.selectedLeptons[1].pt)
        self.out.fillBranch("lep2_eta", event.selectedLeptons[1].eta)
        self.out.fillBranch("lep2_phi", event.selectedLeptons[1].phi)
        self.out.fillBranch("lep2_mass", event.selectedLeptons[1].mass)
        self.out.fillBranch("lep2_pdgId", event.selectedLeptons[1].pdgId)
        self.out.fillBranch("deltaR_ll", deltaR(event.selectedLeptons[0], event.selectedLeptons[1]))

        # AK4 jets, cleaned vs FatJet
        cleaned_ak4jets = [j for j in event.ak4jets if deltaR(j, fj) >= self._jetConeSize]
        self.out.fillBranch("n_ak4_cleaned", len(cleaned_ak4jets))

        # event level
        self.out.fillBranch("deta_v_fj", abs(vboson.Eta() - fj.eta))
        self.out.fillBranch("dphi_v_fj", abs(deltaPhi(vboson.Phi(), fj.phi)))
        self.out.fillBranch("min_deta_v_ak4", minValue([abs(vboson.Eta() - j.eta) for j in cleaned_ak4jets]))
        self.out.fillBranch("min_deta_fj_ak4", minValue([abs(fj.eta - j.eta) for j in cleaned_ak4jets]))

        return True


# define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed
def DibosonTree_2016(): return DibosonSampleProducer(year=2016)
def DibosonTree_2017(): return DibosonSampleProducer(year=2017)
def DibosonTree_2018(): return DibosonSampleProducer(year=2018)
