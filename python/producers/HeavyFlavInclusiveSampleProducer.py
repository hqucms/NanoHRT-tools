from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection

from ..helpers.utils import deltaPhi, deltaR
from .HeavyFlavBaseProducer import HeavyFlavBaseProducer


class InclusiveSampleProducer(HeavyFlavBaseProducer):

    def __init__(self, **kwargs):
        super(InclusiveSampleProducer, self).__init__(channel='inclusive', **kwargs)

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        super(InclusiveSampleProducer, self).beginFile(inputFile, outputFile, inputTree, wrappedOutputTree)

        # event variables
        self.out.branch("htfwd", "F")
        self.out.branch("htveryfwd", "F")
        self.out.branch("nlb_fj_pihalf", "I")
        self.out.branch("nmb_fj_pihalf", "I")
        self.out.branch("ntb_fj_pihalf", "I")
        self.out.branch("nb_fj_pi", "I")
        self.out.branch("nb_away_fj", "I")

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""

        self.selectLeptons(event)
        self.correctJetsAndMET(event)

        if event.ht < 500:
            return False

        probe_jets = [fj for fj in event.fatjets if len(fj.subjets) == 2 and fj.msoftdrop > 30 and fj.msoftdrop < 250]
        if len(probe_jets) == 0:
            return False

        probe_jets = probe_jets[:1]
        self.loadGenHistory(event, probe_jets)
        self.evalTagger(event, probe_jets)
        self.evalMassRegression(event, probe_jets)

        # fill output branches
        self.fillBaseEventInfo(event)
        self.fillFatJetInfo(event, probe_jets)

        # event variables
        htfwd_ = 0.
        htveryfwd_ = 0.
        for j in event._allJets:
            if not (j.pt > 25 and (j.jetId & 2)):
                continue
            if (abs(j.eta) > 2.4):
                htfwd_ += j.pt
            if (abs(j.eta) > 3.5):
                htveryfwd_ += j.pt

        self.out.fillBranch("htfwd", htfwd_)
        self.out.fillBranch("htveryfwd", htveryfwd_)

        # b-tag AK4 jet selection
        event.bljets = []
        event.bmjets = []
        event.btjets = []
        for j in event._allJets:
            if not (j.pt > 20.0 and abs(j.eta) < 2.4 and (j.jetId & 2)):
                continue
            if j.btagDeepFlavB > self.DeepJet_WP_L:
                event.bljets.append(j)
            if j.btagDeepFlavB > self.DeepJet_WP_M:
                event.bmjets.append(j)
            if j.btagDeepFlavB > self.DeepJet_WP_T:
                event.btjets.append(j)

        # count bjets away from fatjet
        nlb_fj_pihalf_ = 0
        nmb_fj_pihalf_ = 0
        ntb_fj_pihalf_ = 0
        nb_fj_pi_ = 0
        for j in event.bljets:
            if abs(deltaPhi(j, event.fatjets[0])) > 3.14:
                nb_fj_pi_ += 1
            if abs(deltaPhi(j, event.fatjets[0])) > 3.14 / 2.:
                nlb_fj_pihalf_ += 1
        for j in event.bmjets:
            if abs(deltaPhi(j, event.fatjets[0])) > 3.14 / 2.:
                nmb_fj_pihalf_ += 1
        for j in event.btjets:
            if abs(deltaPhi(j, event.fatjets[0])) > 3.14 / 2.:
                ntb_fj_pihalf_ += 1

        nb_away_fj_ = 0
        for j in event.bljets:
            if deltaR(j, event.fatjets[0]) > 1.:
                nb_away_fj_ += 1

        self.out.fillBranch("nlb_fj_pihalf", nlb_fj_pihalf_)
        self.out.fillBranch("nmb_fj_pihalf", nmb_fj_pihalf_)
        self.out.fillBranch("ntb_fj_pihalf", ntb_fj_pihalf_)
        self.out.fillBranch("nb_fj_pi", nb_fj_pi_)
        self.out.fillBranch("nb_away_fj", nb_away_fj_)

        return True


# define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed
def InclusiveTree_2016(): return InclusiveSampleProducer(year=2016)
def InclusiveTree_2017(): return InclusiveSampleProducer(year=2017)
def InclusiveTree_2018(): return InclusiveSampleProducer(year=2018)
