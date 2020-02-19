import os
from collections import Counter
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.tools import deltaR, closest

from PhysicsTools.NanoHRTTools.helpers.ak8MassCorrectionHelper import get_corrected_sdmass
from PhysicsTools.NanoHRTTools.helpers.n2DDTHelper import N2DDTHelper
from PhysicsTools.NanoHRTTools.helpers.nnHelper import convert_prob


class _NullObject:
    '''An null object which does not store anything, and does not raise exception.'''
    def __bool__(self):
        return False

    def __nonzero__(self):
        return False

    def __getattr__(self, name):
        pass

    def __setattr__(self, name, value):
        pass


def get_subjets(jet, subjetCollection, idxNames=('subJetIdx1', 'subJetIdx2')):
    subjets = []
    for idxname in idxNames:
        idx = getattr(jet, idxname)
        if idx >= 0:
            subjets.append(subjetCollection[idx])
    return subjets


def get_sdmass(subjets):
    return sum([sj.p4() for sj in subjets], ROOT.TLorentzVector()).M()


class HRTMCTreeProducer(Module):

    def __init__(self):
        self._maxDeltaRJetParton = 0.6
        self._n2helper = N2DDTHelper(os.path.expandvars('$CMSSW_BASE/src/PhysicsTools/NanoHRTTools/data/N2DDT/OutputAK82016v13.root'))

    def beginJob(self):
        pass

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree

        self.out.branch("i_evt", "I")
        self.out.branch("i_parton", "I")
        self.out.branch("type", "I")
        self.out.branch("npv", "I")
        self.out.branch("genweight", "F")

        self.out.branch("gen_n_b_daus", "I")
        self.out.branch("gen_n_c_daus", "I")

        self.out.branch("gen_pt", "F")
        self.out.branch("gen_eta", "F")
        self.out.branch("gen_phi", "F")
        self.out.branch("gen_pdgid", "I")
        self.out.branch("gen_size", "F")
        self.out.branch("gentop_b_pt", "F")
        self.out.branch("gentop_b_eta", "F")
        self.out.branch("gentop_b_phi", "F")
        self.out.branch("gentop_w_pt", "F")
        self.out.branch("gentop_w_eta", "F")
        self.out.branch("gentop_w_phi", "F")
        self.out.branch("gentop_w_size", "F")

        # AK8
        self.out.branch("dR_gen_ak8", "F")
        self.out.branch("ak8_pt", "F")
        self.out.branch("ak8_eta", "F")
        self.out.branch("ak8_phi", "F")
        self.out.branch("ak8_sdmass", "F")
        self.out.branch("ak8_corr_sdmass", "F")
        self.out.branch("ak8_tau21", "F")
        self.out.branch("ak8_tau32", "F")
        self.out.branch("ak8_maxsubjetcsv", "F")
        self.out.branch("ak8_doubleb", "F")

        self.out.branch("ak8_DeepAK8_TvsQCD", "F")
        self.out.branch("ak8_DeepAK8_WvsQCD", "F")
        self.out.branch("ak8_DeepAK8_ZvsQCD", "F")
        self.out.branch("ak8_DeepAK8_HbbvsQCD", "F")

        self.out.branch("ak8_DeepAK8MD_TvsQCD", "F")
        self.out.branch("ak8_DeepAK8MD_WvsQCD", "F")
        self.out.branch("ak8_DeepAK8MD_ZvsQCD", "F")
        self.out.branch("ak8_DeepAK8MD_ZHbbvsQCD", "F")
        self.out.branch("ak8_DeepAK8MD_ZHccvsQCD", "F")

        self.out.branch("ak8_ecfN2", "F")
        self.out.branch("ak8_ecfN2DDT", "F")

    def _fillCommonInfo(self, event, i_parton, parton, daughters):
        self.out.fillBranch("i_evt", int(event.event))
        self.out.fillBranch("npv", event.PV_npvs)
        self.out.fillBranch("genweight", event.genWeight)

        pdgid2type = {24:1, 6:2, 23:3, 25:4,
                      1:0, 2:0, 3:0, 4:0, 5:0, 21:0}
        self.out.fillBranch("i_parton", i_parton)
        self.out.fillBranch("type", pdgid2type[abs(parton.pdgId)])
        self.out.fillBranch("gen_pt", parton.pt)
        self.out.fillBranch("gen_eta", parton.eta)
        self.out.fillBranch("gen_phi", parton.phi)
        self.out.fillBranch("gen_pdgid", parton.pdgId)
        self.out.fillBranch("gen_size", max([deltaR(parton, dau) for dau in daughters]) if len(daughters) else 0)

        c = Counter([abs(dau.pdgId) for dau in daughters])
        self.out.fillBranch("gen_n_b_daus", c[5])
        self.out.fillBranch("gen_n_c_daus", c[4])

        if abs(parton.pdgId) == 6:
            self.out.fillBranch("gentop_b_pt", parton.genB.pt)
            self.out.fillBranch("gentop_b_eta", parton.genB.eta)
            self.out.fillBranch("gentop_b_phi", parton.genB.phi)
            self.out.fillBranch("gentop_w_pt", parton.genW.pt)
            self.out.fillBranch("gentop_w_eta", parton.genW.eta)
            self.out.fillBranch("gentop_w_phi", parton.genW.phi)
            self.out.fillBranch("gentop_w_size", max([deltaR(parton.genW, dau) for dau in daughters[1:]]) if len(daughters) == 3 else 0)
        else:
            self.out.fillBranch("gentop_b_pt", 0)
            self.out.fillBranch("gentop_b_eta", 0)
            self.out.fillBranch("gentop_b_phi", 0)
            self.out.fillBranch("gentop_w_pt", 0)
            self.out.fillBranch("gentop_w_eta", 0)
            self.out.fillBranch("gentop_w_phi", 0)
            self.out.fillBranch("gentop_w_size", 0)

    def _do_matching(self, partons, fatjets):
        rlt = {}
        jets_to_match = [j for j in fatjets]
        for p in partons:
            fj, dR = closest(p, jets_to_match)
            if dR < self._maxDeltaRJetParton:
                rlt[p] = fj
                jets_to_match.remove(fj)
            else:
                rlt[p] = _NullObject()
        return rlt

    def _selectJets(self, event):
        event._allAK8jets = Collection(event, "FatJet")
        event.ak8Subjets = Collection(event, "SubJet")  # do not sort subjets after updating!!

        # link fatjet to subjets and recompute softdrop mass
        for fj in event._allAK8jets:
            fj.subjets = get_subjets(fj, event.ak8Subjets, ('subJetIdx1', 'subJetIdx2'))
            fj.msoftdrop = get_sdmass(fj.subjets)
            fj.corr_sdmass = get_corrected_sdmass(fj, fj.subjets)
            fj.n2b1ddt = self._n2helper.transform(fj.n2b1, pt=fj.pt, msd=fj.corr_sdmass)
        event._allAK8jets = sorted(event._allAK8jets, key=lambda x: x.pt, reverse=True)  # sort by pt

    def _get_filler(self, obj):
        def filler(branch, value, default=0):
            self.out.fillBranch(branch, value if obj else default)
        return filler

    def _fillAK8(self, parton, daughters, jet):
        self.out.fillBranch("dR_gen_ak8", deltaR(parton, jet) if jet else 999)

        fillBranch = self._get_filler(jet)

        fillBranch("ak8_pt", jet.pt, -1)
        fillBranch("ak8_eta", jet.eta)
        fillBranch("ak8_phi", jet.phi)
        fillBranch("ak8_sdmass", jet.msoftdrop)
        fillBranch("ak8_corr_sdmass", jet.corr_sdmass)
        fillBranch("ak8_tau21", jet.tau2 / jet.tau1 if jet.tau1 > 0 else 99, 99)
        fillBranch("ak8_tau32", jet.tau3 / jet.tau2 if jet.tau2 > 0 else 99, 99)
        fillBranch("ak8_maxsubjetcsv", max([sj.btagCSVV2 for sj in jet.subjets]) if jet and len(jet.subjets) else -99, -99)
        fillBranch("ak8_doubleb", jet.btagHbb, -99)

        fillBranch("ak8_DeepAK8_TvsQCD", jet.deepTag_TvsQCD)
        fillBranch("ak8_DeepAK8_WvsQCD", jet.deepTag_WvsQCD)
        fillBranch("ak8_DeepAK8_ZvsQCD", jet.deepTag_ZvsQCD)

        fillBranch("ak8_DeepAK8MD_TvsQCD", jet.deepTagMD_TvsQCD)
        fillBranch("ak8_DeepAK8MD_WvsQCD", jet.deepTagMD_WvsQCD)
        fillBranch("ak8_DeepAK8MD_ZvsQCD", jet.deepTagMD_ZvsQCD)
        fillBranch("ak8_DeepAK8MD_ZHbbvsQCD", jet.deepTagMD_ZHbbvsQCD)
        fillBranch("ak8_DeepAK8MD_ZHccvsQCD", jet.deepTagMD_ZHccvsQCD)

        fillBranch("ak8_ecfN2", jet.n2b1, 99)
        fillBranch("ak8_ecfN2DDT", jet.n2b1ddt, 99)

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        genparts = Collection(event, "GenPart")
        for idx, gp in enumerate(genparts):
            if not hasattr(gp, 'dauIdx'):
                gp.dauIdx = []
            if gp.genPartIdxMother >= 0:
                mom = genparts[gp.genPartIdxMother]
                if not hasattr(mom, 'dauIdx'):
                    mom.dauIdx = [idx]
                else:
                    mom.dauIdx.append(idx)

        def isHadronic(gp):
            if len(gp.dauIdx) == 0:
                raise ValueError('Particle has no daughters!')
            for idx in gp.dauIdx:
                if abs(genparts[idx].pdgId) < 6:
                    return True
            return False

        def getFinal(gp):
            for idx in gp.dauIdx:
                dau = genparts[idx]
                if dau.pdgId == gp.pdgId:
                    return getFinal(dau)
            return gp

        nGenTops = 0
        nGenWs = 0
        nGenZs = 0
        nGenHs = 0

        hadGenTops = []
        hadGenWs = []
        hadGenZs = []
        hadGenHs = []
        genQuarkGluons = []
        for gp in genparts:
            if gp.pt < 200 or abs(gp.eta) > 2.4:
                continue

            if (abs(gp.pdgId) <= 5 or abs(gp.pdgId) == 21) and gp.status == 23:
                genQuarkGluons.append(gp)

            if gp.statusFlags & (1 << 13) == 0:
                continue
            if abs(gp.pdgId) == 6:
                nGenTops += 1
                for idx in gp.dauIdx:
                    dau = genparts[idx]
                    if abs(dau.pdgId) == 24:
                        genW = getFinal(dau)
                        gp.genW = genW
                        if isHadronic(genW):
                            hadGenTops.append(gp)
                    elif abs(dau.pdgId) in (1, 3, 5):
                        gp.genB = dau
            elif abs(gp.pdgId) == 24:
                nGenWs += 1
                if isHadronic(gp):
                    hadGenWs.append(gp)
            elif abs(gp.pdgId) == 23:
                nGenZs += 1
                if isHadronic(gp):
                    hadGenZs.append(gp)
            elif abs(gp.pdgId) == 25:
                nGenHs += 1
                if isHadronic(gp):
                    hadGenHs.append(gp)

        genHadPartons = []
        if nGenTops > 0:
            # top sample
            genHadPartons = hadGenTops
        elif nGenHs > 0:
            # Higgs sample
            genHadPartons = hadGenHs
        elif nGenWs > 0:
            # W sample
            genHadPartons = hadGenWs
        elif nGenZs > 0:
            # Z sample
            genHadPartons = hadGenZs
        else:
            # QCD sample
            genHadPartons = genQuarkGluons
        genHadPartons.sort(key=lambda x: x.pt, reverse=True)  # sort by pt

        def get_daughters(parton):
            if abs(parton.pdgId) == 6:
                return (parton.genB, genparts[parton.genW.dauIdx[0]], genparts[parton.genW.dauIdx[1]])
            elif abs(parton.pdgId) in (23, 24, 25):
                return (genparts[parton.dauIdx[0]], genparts[parton.dauIdx[1]])
            elif abs(parton.pdgId) <= 5 or parton.pdgId == 21:
                return ()

        self._selectJets(event)
        # selection on AK8 jets
        event.ak8jets = []
        for fj in event._allAK8jets:
            if fj.pt > 200 and abs(fj.eta) < 2.4 and (fj.jetId & 2):
                event.ak8jets.append(fj)
        ak8match = self._do_matching(genHadPartons, event.ak8jets)

        for iparton, parton in enumerate(genHadPartons):
            daughters = get_daughters(parton)
            self._fillCommonInfo(event, iparton, parton, daughters)
            self._fillAK8(parton, daughters, ak8match[parton])
            self.out.fill()

        return False


# define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed
hrtMCTree = lambda: HRTMCTreeProducer()
