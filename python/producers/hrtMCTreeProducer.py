import os
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.tools import deltaR, closest

from PhysicsTools.NanoHRTTools.helpers.ak8MassCorrectionHelper import get_corrected_sdmass, get_sdmass_fromsubjets
from PhysicsTools.NanoHRTTools.helpers.deepAK8Helper import get_nominal_score, get_decorr_score
from PhysicsTools.NanoHRTTools.helpers.n2DDTHelper import N2DDTHelper


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


class HRTMCTreeProducer(Module):

    def __init__(self):
        self._maxDeltaRJetParton = 0.6
        self._deepAK8_nominal_scores = ('TvsQCD', 'WvsQCD', 'ZvsQCD', 'HbbvsQCD')
        self._deepAK8_decorr_scores = ('TvsQCD', 'WvsQCD', 'ZHbbvsQCD', 'bbvsLight')
        self._BEST_scores = ('bestT', 'bestW', 'bestZ', 'bestH', 'bestQCD', 'bestB')
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

        self.out.branch("gen_pt", "F")
        self.out.branch("gen_eta", "F")
        self.out.branch("gen_phi", "F")
        self.out.branch("gen_pdgid", "I")
        self.out.branch("gen_size", "F")

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

        for name in self._deepAK8_nominal_scores:
            self.out.branch("ak8_nn_%s" % name, "F")

        for name in self._deepAK8_decorr_scores:
            self.out.branch("ak8_decorr_nn_%s" % name, "F")

        for name in self._BEST_scores:
            self.out.branch("ak8_%s" % name, "F")

        self.out.branch("ak8_ecfN2", "F")
        self.out.branch("ak8_ecfN2DDT", "F")

        # HOTVR
        self.out.branch("dR_gen_hotvr", "F")
        self.out.branch("hotvr_pt", "F")
        self.out.branch("hotvr_eta", "F")
        self.out.branch("hotvr_phi", "F")
        self.out.branch("hotvr_mass", "F")
        self.out.branch("hotvr_tau32", "F")
        self.out.branch("hotvr_fpt", "F")
        self.out.branch("hotvr_mmin", "F")
        self.out.branch("hotvr_nsubjets", "F")
        self.out.branch("hotvr_mass_fromsubjets", "F")

        # CA15
        self.out.branch("dR_gen_ca15", "F")
        self.out.branch("ca15_pt", "F")
        self.out.branch("ca15_eta", "F")
        self.out.branch("ca15_phi", "F")
        self.out.branch("ca15_sdmass", "F")
        self.out.branch("ca15_ecfTopTagBDT", "F")
        self.out.branch("ca15_maxsubjetcsv", "F")

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

    def _fill_matching(self, parton, daughters, fatjetCollection, fjname):
        fj, dR = closest(parton, fatjetCollection)
        self.out.fillBranch("dR_gen_%s" % fjname, dR)
        return _NullObject() if dR > self._maxDeltaRJetParton else fj

    def _get_filler(self, obj):
        def filler(branch, value, default=0):
            self.out.fillBranch(branch, value if obj else default)
        return filler

    def _get_subjets(self, jet, subjetCollection, idxNames):
        subjets = []
        for idxname in idxNames:
            idx = getattr(jet, idxname)
            if idx >= 0:
                subjets.append(subjetCollection[idx])
        return subjets

    def _fillAK8(self, parton, daughters, fatjetCollection, subjetCollection):
        jet = self._fill_matching(parton, daughters, fatjetCollection, fjname='ak8')
        subjets = self._get_subjets(jet, subjetCollection, ('subJetIdx1', 'subJetIdx2')) if jet else []
        jet.corr_sdmass = get_corrected_sdmass(jet, subjets)

        fillBranch = self._get_filler(jet)

        fillBranch("ak8_pt", jet.pt, -1)
        fillBranch("ak8_eta", jet.eta)
        fillBranch("ak8_phi", jet.phi)
        fillBranch("ak8_sdmass", jet.msoftdrop)
        fillBranch("ak8_corr_sdmass", jet.corr_sdmass)
        fillBranch("ak8_tau21", jet.tau2 / jet.tau1 if jet.tau1 > 0 else 99, 99)
        fillBranch("ak8_tau32", jet.tau3 / jet.tau2 if jet.tau2 > 0 else 99, 99)
        fillBranch("ak8_maxsubjetcsv", max([sj.btagCSVV2 for sj in subjets]) if len(subjets) else -99, -99)
        fillBranch("ak8_doubleb", jet.btagHbb, -99)

        for name in self._deepAK8_nominal_scores:
            fillBranch("ak8_nn_%s" % name, get_nominal_score(jet, name), -1)

        for name in self._deepAK8_decorr_scores:
            fillBranch("ak8_decorr_nn_%s" % name, get_decorr_score(jet, name), -1)

        for name in self._BEST_scores:
            fillBranch("ak8_%s" % name, getattr(jet, name), -1)

        fillBranch("ak8_ecfN2", jet.n2b1, 99)
        jet.n2b1ddt = self._n2helper.transform(jet.n2b1, pt=jet.pt, msd=jet.corr_sdmass)
        fillBranch("ak8_ecfN2DDT", jet.n2b1ddt, 99)

    def _fillHOTVR(self, parton, daughters, fatjetCollection, subjetCollection):
        jet = self._fill_matching(parton, daughters, fatjetCollection, fjname='hotvr')
        subjets = self._get_subjets(jet, subjetCollection, ('subJetIdx1', 'subJetIdx2', 'subJetIdx3')) if jet else []
        jet.mass_from_subjets = get_sdmass_fromsubjets(jet, subjets)
        fillBranch = self._get_filler(jet)

        fillBranch("hotvr_pt", jet.pt, -1)
        fillBranch("hotvr_eta", jet.eta)
        fillBranch("hotvr_phi", jet.phi)
        fillBranch("hotvr_mass", jet.mass)
        fillBranch("hotvr_tau32", jet.tau3 / jet.tau2 if jet.tau2 > 0 else 99, 99)
        fillBranch("hotvr_fpt", jet.fpt)
        fillBranch("hotvr_mmin", jet.mmin)
        fillBranch("hotvr_nsubjets", jet.nsubjets)
        fillBranch("hotvr_mass_fromsubjets", jet.mass_from_subjets)

    def _fillCA15(self, parton, daughters, fatjetCollection, subjetCollection):
        jet = self._fill_matching(parton, daughters, fatjetCollection, fjname='ca15')
        subjets = self._get_subjets(jet, subjetCollection, ('subJetIdx1', 'subJetIdx2')) if jet else []
        fillBranch = self._get_filler(jet)

        fillBranch("ca15_pt", jet.pt, -1)
        fillBranch("ca15_eta", jet.eta)
        fillBranch("ca15_phi", jet.phi)
        fillBranch("ca15_sdmass", jet.msoftdrop)
        fillBranch("ca15_ecfTopTagBDT", jet.ecfTopTagBDT, -99)
        fillBranch("ca15_maxsubjetcsv", max([sj.btagCSVV2 for sj in subjets]) if len(subjets) else -99, -99)

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

        def get_daughters(parton):
            if abs(parton.pdgId) == 6:
                return (parton.genB, genparts[parton.genW.dauIdx[0]], genparts[parton.genW.dauIdx[1]])
            elif abs(parton.pdgId) in (23, 24, 25):
                return (genparts[parton.dauIdx[0]], genparts[parton.dauIdx[1]])
            elif abs(parton.pdgId) <= 5 or parton.pdgId == 21:
                return ()

        ak8jets = Collection(event, "CustomAK8Puppi")
        ak8subjets = Collection(event, "CustomAK8PuppiSubJet")

        hotvrjets = Collection(event, "HOTVRPuppi")
        hotvrsubjets = Collection(event, "HOTVRPuppiSubJet")

        ca15jets = Collection(event, "CA15Puppi")
        ca15subjets = Collection(event, "CA15PuppiSubJet")

        for iparton, parton in enumerate(genHadPartons):
            daughters = get_daughters(parton)
            self._fillCommonInfo(event, iparton, parton, daughters)
            self._fillAK8(parton, daughters, ak8jets, ak8subjets)
            self._fillHOTVR(parton, daughters, hotvrjets, hotvrsubjets)
            self._fillCA15(parton, daughters, ca15jets, ca15subjets)
            self.out.fill()

        return False


# define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed
hrtMCTree = lambda : HRTMCTreeProducer()
