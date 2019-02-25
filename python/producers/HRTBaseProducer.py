import os
import numpy as np
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.tools import deltaPhi, deltaR, closest
from PhysicsTools.NanoHRTTools.helpers.jetmetCorrector import JetMETCorrector

from PhysicsTools.NanoHRTTools.helpers.ak8MassCorrectionHelper import get_corrected_sdmass
from PhysicsTools.NanoHRTTools.helpers.deepAK8Helper import get_nominal_score, get_decorr_score
from PhysicsTools.NanoHRTTools.helpers.n2DDTHelper import N2DDTHelper

import logging
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] %(levelname)s: %(message)s')


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


class METObject(Object):

    def p4(self):
        ret = ROOT.TLorentzVector()
        ret.SetPtEtaPhiM(self.pt, 0, self.phi, 0)
        return ret


def get_subjets(jet, subjetCollection, idxNames=('subJetIdx1', 'subJetIdx2')):
    subjets = []
    for idxname in idxNames:
        idx = getattr(jet, idxname)
        if idx >= 0:
            subjets.append(subjetCollection[idx])
    return subjets


def get_sdmass(subjets):
    return sum([sj.p4() for sj in subjets], ROOT.TLorentzVector()).M()


def transverseMass(obj, met):
    cos_dphi = np.cos(deltaPhi(obj, met))
    return np.sqrt(2 * obj.pt * met.pt * (1 - cos_dphi))


def minValue(collection, fallback=99):
    if len(collection) == 0:
        return fallback
    else:
        return min(collection)


def maxValue(collection, fallback=0):
    if len(collection) == 0:
        return fallback
    else:
        return max(collection)


def rndSeed(event, jets, extra=0):
    seed = (event.run << 20) + (event.luminosityBlock << 10) + event.event + extra
    if len(jets) > 0:
        seed += int(jets[0].eta / 0.01)
    return seed


class HRTBaseProducer(Module, object):

    def __init__(self, channel, **kwargs):
        self._channel = channel
        self._systOpt = {'jec':False, 'jes':None, 'jes_source':'', 'jer':'nominal', 'jmr':None, 'met_unclustered':None}
        for k in kwargs:
            self._systOpt[k] = kwargs[k]

        logging.info('Running %s channel with systematics %s', self._channel, str(self._systOpt))

        self.jetmetCorr = JetMETCorrector(jetType="AK4PFchs",
                                          jec=self._systOpt['jec'],
                                          jes=self._systOpt['jes'],
                                          jes_source=self._systOpt['jes_source'],
                                          jer=self._systOpt['jer'],
                                          met_unclustered=self._systOpt['met_unclustered'])

        self.ak8Corr = JetMETCorrector(jetType="AK8PFPuppi",
                                          jec=self._systOpt['jec'],
                                          jes=self._systOpt['jes'],
                                          jes_source=self._systOpt['jes_source'],
                                          jer=self._systOpt['jer'],
                                          jmr=self._systOpt['jmr'],
                                          met_unclustered=self._systOpt['met_unclustered'])

        self.ak8SubjetCorr = JetMETCorrector(jetType="AK4PFPuppi",
                                          jec=self._systOpt['jec'],
                                          jes=self._systOpt['jes'],
                                          jes_source=self._systOpt['jes_source'],
                                          jer=self._systOpt['jer'],
                                          jmr=self._systOpt['jmr'],
                                          met_unclustered=self._systOpt['met_unclustered'])

        self.ca15Corr = JetMETCorrector(jetType="AK8PFPuppi",
                                          jec=self._systOpt['jec'],
                                          jes=self._systOpt['jes'],
                                          jes_source=self._systOpt['jes_source'],
                                          jer=self._systOpt['jer'],
                                          jmr=self._systOpt['jmr'],
                                          met_unclustered=self._systOpt['met_unclustered'])

        self.ca15SubjetCorr = JetMETCorrector(jetType="AK4PFPuppi",
                                          jec=self._systOpt['jec'],
                                          jes=self._systOpt['jes'],
                                          jes_source=self._systOpt['jes_source'],
                                          jer=self._systOpt['jer'],
                                          jmr=self._systOpt['jmr'],
                                          met_unclustered=self._systOpt['met_unclustered'])

        self.hotvrSubjetCorr = JetMETCorrector(jetType="AK4PFPuppi",
                                          jec=self._systOpt['jec'],
                                          jes=self._systOpt['jes'],
                                          jes_source=self._systOpt['jes_source'],
                                          jer=self._systOpt['jer'],
                                          met_unclustered=self._systOpt['met_unclustered'])

        self._n2helper = N2DDTHelper(os.path.expandvars('$CMSSW_BASE/src/PhysicsTools/NanoHRTTools/data/N2DDT/OutputAK82016v13.root'))

    def beginJob(self):
        self.jetmetCorr.beginJob()
        self.ak8Corr.beginJob()
        self.ak8SubjetCorr.beginJob()
        self.ca15Corr.beginJob()
        self.ca15SubjetCorr.beginJob()
        self.hotvrSubjetCorr.beginJob()

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.isMC = bool(inputTree.GetBranch('genWeight'))
        self.out = wrappedOutputTree

        self.out.branch("passmetfilters", "O")

        # Large-R jets
        self.out.branch("n_ak8", "I")
        self.out.branch("n_ca15", "I")
        self.out.branch("n_hotvr", "I")

        self.out.branch("ak8_1_pt", "F")
        self.out.branch("ak8_1_eta", "F")
        self.out.branch("ak8_1_phi", "F")
        self.out.branch("ak8_1_mass", "F")
        self.out.branch("ak8_1_corr_sdmass", "F")
        self.out.branch("ak8_1_n2b1ddt", "F")
        self.out.branch("ak8_1_n2b1", "F")
        self.out.branch("ak8_1_n3b1", "F")
        self.out.branch("ak8_1_tau1", "F")
        self.out.branch("ak8_1_tau2", "F")
        self.out.branch("ak8_1_tau3", "F")
        self.out.branch("ak8_1_DeepAK8_WvsQCD", "F")
        self.out.branch("ak8_1_DeepAK8_ZvsQCD", "F")
        self.out.branch("ak8_1_DeepAK8_HvsQCD", "F")
        self.out.branch("ak8_1_DeepAK8_TvsQCD", "F")
        self.out.branch("ak8_1_DeepAK8MD_WvsQCD", "F")
        self.out.branch("ak8_1_DeepAK8MD_ZvsQCD", "F")
        self.out.branch("ak8_1_DeepAK8MD_HvsQCD", "F")
        self.out.branch("ak8_1_DeepAK8MD_TvsQCD", "F")
        self.out.branch("ak8_1_best_WvsQCD", "F")
        self.out.branch("ak8_1_best_ZvsQCD", "F")
        self.out.branch("ak8_1_best_HvsQCD", "F")
        self.out.branch("ak8_1_best_TvsQCD", "F")
        self.out.branch("ak8_1_image_top", "F")
        self.out.branch("ak8_1_image_top_md", "F")
        self.out.branch("ak8_1_btagCSVV2", "F")
        self.out.branch("ak8_1_btagHbb", "F")
        self.out.branch("ak8_1_sj1_btagCSVV2", "F")
        self.out.branch("ak8_1_sj2_btagCSVV2", "F")

        self.out.branch("ca15_1_pt", "F")
        self.out.branch("ca15_1_eta", "F")
        self.out.branch("ca15_1_phi", "F")
        self.out.branch("ca15_1_mass", "F")
        self.out.branch("ca15_1_ecf0", "F")
        self.out.branch("ca15_1_ecfTopTagBDT", "F")
        self.out.branch("ca15_1_httFRec", "F")
        self.out.branch("ca15_1_tau32sd", "F")

        self.out.branch("hotvr_1_pt", "F")
        self.out.branch("hotvr_1_eta", "F")
        self.out.branch("hotvr_1_phi", "F")
        self.out.branch("hotvr_1_mass", "F")
        self.out.branch("hotvr_1_tau32", "F")
        self.out.branch("hotvr_1_fpt", "F")
        self.out.branch("hotvr_1_mmin", "F")
        self.out.branch("hotvr_1_nsubjets", "F")

        # matching variables
        if self.isMC:
            self.out.branch("ak8_1_dr_fj_top_b", "F")
            self.out.branch("ak8_1_dr_fj_top_wqmax", "F")
            self.out.branch("ak8_1_dr_fj_top_wqmin", "F")

            self.out.branch("ca15_1_dr_fj_top_b", "F")
            self.out.branch("ca15_1_dr_fj_top_wqmax", "F")
            self.out.branch("ca15_1_dr_fj_top_wqmin", "F")

            self.out.branch("hotvr_1_dr_fj_top_b", "F")
            self.out.branch("hotvr_1_dr_fj_top_wqmax", "F")
            self.out.branch("hotvr_1_dr_fj_top_wqmin", "F")

    def correctJetsAndMET(self, event):
        # correct Jets and MET
        event._allJets = Collection(event, "Jet")
        event.met = METObject(event, "MET")

        event._allAK8jets = Collection(event, "CustomAK8Puppi")
        event.ak8Subjets = Collection(event, "CustomAK8PuppiSubJet")  # do not sort subjets after updating!!

        event._allCA15jets = Collection(event, "CA15Puppi")
        event.ca15Subjets = Collection(event, "CA15PuppiSubJet")  # do not sort subjets after updating!!

        event._allHOTVRjets = Collection(event, "HOTVRPuppi")
        event.hotvrSubjets = Collection(event, "HOTVRPuppiSubJet")  # do not sort subjets after updating!!

        if self.isMC or self._systOpt['jec']:
            rho = event.fixedGridRhoFastjetAll
            # correct AK4 jets and MET
            self.jetmetCorr.setSeed(rndSeed(event, event._allJets))
            self.jetmetCorr.correctJetAndMET(jets=event._allJets, met=event.met, rho=rho,
                                             genjets=Collection(event, 'GenJet') if self.isMC else None,
                                             isMC=self.isMC, runNumber=event.run)
            event._allJets = sorted(event._allJets, key=lambda x : x.pt, reverse=True)  # sort by pt after updating

            # correct AK8 fatjets
            self.ak8Corr.setSeed(rndSeed(event, event._allAK8jets))
            self.ak8Corr.correctJetAndMET(jets=event._allAK8jets, met=None, rho=rho,
                                             genjets=Collection(event, 'CustomGenJetAK8') if self.isMC else None,
                                             isMC=self.isMC, runNumber=event.run)
            # correct AK8 subjets
            self.ak8SubjetCorr.setSeed(rndSeed(event, event.ak8Subjets))
            self.ak8SubjetCorr.correctJetAndMET(jets=event.ak8Subjets, met=None, rho=rho,
                                             genjets=Collection(event, 'CustomGenSubJetAK8') if self.isMC else None,
                                             isMC=self.isMC, runNumber=event.run)

            # correct AK8 fatjets
            self.ca15Corr.setSeed(rndSeed(event, event._allCA15jets))
            self.ca15Corr.correctJetAndMET(jets=event._allCA15jets, met=None, rho=rho,
                                             genjets=Collection(event, 'GenJetCA15') if self.isMC else None,
                                             isMC=self.isMC, runNumber=event.run)
            # correct AK8 subjets
            self.ca15SubjetCorr.setSeed(rndSeed(event, event.ca15Subjets))
            self.ca15SubjetCorr.correctJetAndMET(jets=event.ca15Subjets, met=None, rho=rho,
                                             genjets=Collection(event, 'GenSubJetCA15') if self.isMC else None,
                                             isMC=self.isMC, runNumber=event.run)

            # correct HOTVR subjets
            self.hotvrSubjetCorr.setSeed(rndSeed(event, event.hotvrSubjets))
            self.hotvrSubjetCorr.correctJetAndMET(jets=event.hotvrSubjets, met=None, rho=rho,
                                             genjets=Collection(event, 'GenJet') if self.isMC else None,
                                             isMC=self.isMC, runNumber=event.run)

        # jet mass resolution smearing
        if self.isMC and self._systOpt['jmr']:
            raise NotImplemented

        # link fatjet to subjets and recompute softdrop mass
        for fj in event._allAK8jets:
            fj.subjets = get_subjets(fj, event.ak8Subjets, ('subJetIdx1', 'subJetIdx2'))
            fj.msoftdrop = get_sdmass(fj.subjets)
            fj.corr_sdmass = get_corrected_sdmass(fj, fj.subjets)
            fj.n2b1ddt = self._n2helper.transform(fj.n2b1, pt=fj.pt, msd=fj.corr_sdmass)
        event._allAK8jets = sorted(event._allAK8jets, key=lambda x : x.pt, reverse=True)  # sort by pt

        for fj in event._allCA15jets:
            fj.subjets = get_subjets(fj, event.ca15Subjets, ('subJetIdx1', 'subJetIdx2'))
            fj.msoftdrop = get_sdmass(fj.subjets)
        event._allCA15jets = sorted(event._allCA15jets, key=lambda x : x.pt, reverse=True)  # sort by pt

        # for HOTVR we use the sum of subjets p4
        for fj in event._allHOTVRjets:
            fj.subjets = get_subjets(fj, event.hotvrSubjets, ('subJetIdx1', 'subJetIdx2', 'subJetIdx3'))
            newP4 = sum([sj.p4() for sj in fj.subjets], ROOT.TLorentzVector())
            fj.pt, fj.eta, fj.phi, fj.mass = newP4.Pt(), newP4.Eta(), newP4.Phi(), newP4.M()
        event._allHOTVRjets = sorted(event._allHOTVRjets, key=lambda x : x.pt, reverse=True)  # sort by pt

    def loadGenHistory(self, event):
        # gen matching
        if not self.isMC:
            return

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

        event.nGenTops = 0
        event.nGenWs = 0
        event.nGenZs = 0
        event.nGenHs = 0

        event.hadGenTops = []
        event.hadGenWs = []
        event.hadGenZs = []
        event.hadGenHs = []

        for gp in genparts:
            if gp.statusFlags & (1 << 13) == 0:
                continue
            if abs(gp.pdgId) == 6:
                event.nGenTops += 1
                for idx in gp.dauIdx:
                    dau = genparts[idx]
                    if abs(dau.pdgId) == 24:
                        genW = getFinal(dau)
                        gp.genW = genW
                        if isHadronic(genW):
                            event.hadGenTops.append(gp)
                    elif abs(dau.pdgId) in (1, 3, 5):
                        gp.genB = dau
            elif abs(gp.pdgId) == 24:
                event.nGenWs += 1
                if isHadronic(gp):
                    event.hadGenWs.append(gp)
            elif abs(gp.pdgId) == 23:
                event.nGenZs += 1
                if isHadronic(gp):
                    event.hadGenZs.append(gp)
            elif abs(gp.pdgId) == 25:
                event.nGenHs += 1
                if isHadronic(gp):
                    event.hadGenHs.append(gp)

        event.genparts = genparts

    def fillBaseEventInfo(self, event):
        met_filters = bool(
            event.Flag_goodVertices and
            event.Flag_globalSuperTightHalo2016Filter and
            event.Flag_HBHENoiseFilter and
            event.Flag_HBHENoiseIsoFilter and
            event.Flag_EcalDeadCellTriggerPrimitiveFilter and
            event.Flag_BadPFMuonFilter and
            event.Flag_BadChargedCandidateFilter
            )
        if not self.isMC:
            met_filters = met_filters and event.Flag_eeBadScFilter
        self.out.fillBranch("passmetfilters", met_filters)

    def _get_filler(self, obj):

        def filler(branch, value, default=0):
            self.out.fillBranch(branch, value if obj else default)

        return filler

    def fillFatJetInfo(self, event, fillGenMatching=False):
        # fill AK8
        self.out.fillBranch("n_ak8", len(event.ak8jets))
        ak8 = event.ak8jets[0] if len(event.ak8jets) > 0 else _NullObject()
        fillBranchAK8 = self._get_filler(ak8)  # wrapper, fill default value if ak8=None
        fillBranchAK8("ak8_1_pt", ak8.pt)
        fillBranchAK8("ak8_1_eta", ak8.eta)
        fillBranchAK8("ak8_1_phi", ak8.phi)
        fillBranchAK8("ak8_1_mass", ak8.msoftdrop)
        fillBranchAK8("ak8_1_corr_sdmass", ak8.corr_sdmass)
        fillBranchAK8("ak8_1_n2b1ddt", ak8.n2b1ddt)
        fillBranchAK8("ak8_1_n2b1", ak8.n2b1)
        fillBranchAK8("ak8_1_n3b1", ak8.n3b1)
        fillBranchAK8("ak8_1_tau1", ak8.tau1)
        fillBranchAK8("ak8_1_tau2", ak8.tau2)
        fillBranchAK8("ak8_1_tau3", ak8.tau3)
        fillBranchAK8("ak8_1_DeepAK8_TvsQCD", get_nominal_score(ak8, 'TvsQCD'))
        fillBranchAK8("ak8_1_DeepAK8_WvsQCD", get_nominal_score(ak8, 'WvsQCD'))
        fillBranchAK8("ak8_1_DeepAK8_ZvsQCD", get_nominal_score(ak8, 'ZvsQCD'))
        fillBranchAK8("ak8_1_DeepAK8_HvsQCD", get_nominal_score(ak8, 'HbbvsQCD'))
        fillBranchAK8("ak8_1_DeepAK8MD_TvsQCD", get_decorr_score(ak8, 'TvsQCD'))
        fillBranchAK8("ak8_1_DeepAK8MD_WvsQCD", get_decorr_score(ak8, 'WvsQCD'))
        fillBranchAK8("ak8_1_DeepAK8MD_ZvsQCD", get_decorr_score(ak8, 'ZvsQCD'))
        fillBranchAK8("ak8_1_DeepAK8MD_HvsQCD", get_decorr_score(ak8, 'ZHbbvsQCD'))
        fillBranchAK8("ak8_1_best_WvsQCD", ak8.bestW / (ak8.bestW + ak8.bestQCD + ak8.bestB))
        fillBranchAK8("ak8_1_best_ZvsQCD", ak8.bestZ / (ak8.bestZ + ak8.bestQCD + ak8.bestB))
        fillBranchAK8("ak8_1_best_HvsQCD", ak8.bestH / (ak8.bestH + ak8.bestQCD + ak8.bestB))
        fillBranchAK8("ak8_1_best_TvsQCD", ak8.bestT / (ak8.bestT + ak8.bestQCD + ak8.bestB))
        fillBranchAK8("ak8_1_image_top", ak8.itop)
        fillBranchAK8("ak8_1_image_top_md", ak8.iMDtop)
        fillBranchAK8("ak8_1_btagCSVV2", ak8.btagCSVV2)
        fillBranchAK8("ak8_1_btagHbb", ak8.btagHbb)
        fillBranchAK8("ak8_1_sj1_btagCSVV2", ak8.subjets[0].btagCSVV2 if len(ak8.subjets) > 0 else -9)
        fillBranchAK8("ak8_1_sj2_btagCSVV2", ak8.subjets[1].btagCSVV2 if len(ak8.subjets) > 1 else -9)

        # fill CA15
        self.out.fillBranch("n_ca15", len(event.ca15jets))
        ca15 = event.ca15jets[0] if len(event.ca15jets) > 0 else _NullObject()
        fillBranchCA15 = self._get_filler(ca15)  # wrapper, fill default value if ca15=None
        fillBranchCA15("ca15_1_pt", ca15.pt)
        fillBranchCA15("ca15_1_eta", ca15.eta)
        fillBranchCA15("ca15_1_phi", ca15.phi)
        fillBranchCA15("ca15_1_mass", ca15.msoftdrop)
        fillBranchCA15("ca15_1_ecfTopTagBDT", ca15.ecfTopTagBDT)
        fillBranchCA15("ca15_1_httFRec", ca15.httFRec)
        fillBranchCA15("ca15_1_tau32sd", ca15.tau32sd)

        # fill HOTVR
        self.out.fillBranch("n_hotvr", len(event.hotvrjets))
        hotvr = event.hotvrjets[0] if len(event.hotvrjets) > 0 else _NullObject()
        fillBranchHOTVR = self._get_filler(hotvr)  # wrapper, fill default value if hotvr=None
        fillBranchHOTVR("hotvr_1_pt", hotvr.pt)
        fillBranchHOTVR("hotvr_1_eta", hotvr.eta)
        fillBranchHOTVR("hotvr_1_phi", hotvr.phi)
        fillBranchHOTVR("hotvr_1_mass", hotvr.mass)
        fillBranchHOTVR("hotvr_1_tau32", hotvr.tau3 / hotvr.tau2 if hotvr.tau2 > 0 else 99, 99)
        fillBranchHOTVR("hotvr_1_fpt", hotvr.fpt)
        fillBranchHOTVR("hotvr_1_mmin", hotvr.mmin)
        fillBranchHOTVR("hotvr_1_nsubjets", hotvr.nsubjets)

        def drmatch(event, jet):
            dr_b, dr_wq1, dr_wq2 = 999, 999, 999
            if jet:
                top, _ = closest(jet, event.hadGenTops)
                if top:
                    dr_b = deltaR(jet, top.genB)
                    dr_wq1 = deltaR(jet, event.genparts[top.genW.dauIdx[0]])
                    dr_wq2 = deltaR(jet, event.genparts[top.genW.dauIdx[1]])
                else:
                    w, _ = closest(jet, event.hadGenWs)
                    if w:
                        dr_wq1 = deltaR(jet, event.genparts[w.dauIdx[0]])
                        dr_wq2 = deltaR(jet, event.genparts[w.dauIdx[1]])
            return dr_b, max(dr_wq1, dr_wq2), min(dr_wq1, dr_wq2)
        
        if self.isMC and fillGenMatching:
            drak8 = drmatch(event, ak8)
            self.out.fillBranch("ak8_1_dr_fj_top_b", drak8[0])
            self.out.fillBranch("ak8_1_dr_fj_top_wqmax", drak8[1])
            self.out.fillBranch("ak8_1_dr_fj_top_wqmin", drak8[2])

            drca15 = drmatch(event, ca15)
            self.out.fillBranch("ca15_1_dr_fj_top_b", drca15[0])
            self.out.fillBranch("ca15_1_dr_fj_top_wqmax", drca15[1])
            self.out.fillBranch("ca15_1_dr_fj_top_wqmin", drca15[2])

            drhotvr = drmatch(event, hotvr)
            self.out.fillBranch("hotvr_1_dr_fj_top_b", drhotvr[0])
            self.out.fillBranch("hotvr_1_dr_fj_top_wqmax", drhotvr[1])
            self.out.fillBranch("hotvr_1_dr_fj_top_wqmin", drhotvr[2])

