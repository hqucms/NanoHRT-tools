import os
import math
import itertools
import numpy as np
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.tools import deltaPhi, deltaR, closest

from PhysicsTools.NanoHRTTools.helpers.jetmetCorrector import JetMETCorrector
from PhysicsTools.NanoHRTTools.helpers.nnHelper import convert_prob
from PhysicsTools.NanoHRTTools.helpers.xgbHelper import XGBEnsemble

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


def corrected_svmass(sv):
    pproj = sv.p4().P() * math.sin(sv.pAngle)
    return math.sqrt(sv.mass * sv.mass + pproj * pproj) + pproj


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


class HeavyFlavBaseProducer(Module, object):

    def __init__(self, channel, **kwargs):
        self.year = int(kwargs['year'])
        self.jetType = kwargs.get('jetType', 'ak8').lower()
        if self.jetType == 'ak8':
            self._jetConeSize = 0.8
            self._fj_name = 'FatJet'
            self._sj_name = 'SubJet'
            self._fj_gen_name = 'GenJetAK8'
            self._sj_gen_name = 'SubGenJetAK8'
            # self._sfbdt_files = [os.path.expandvars('$CMSSW_BASE/src/PhysicsTools/NanoHRTTools/data/sfBDT/ak8/xgb_train_qcd.model.%d' % idx) for idx in range(10)]
            self._sfbdt_files = [os.path.expandvars('$CMSSW_BASE/src/PhysicsTools/NanoHRTTools/data/sfBDT/ak15/xgb_train_qcd.model.%d' % idx) for idx in range(10)] # FIXME: update to AK8 training
            self._sfbdt_vars = ['fj_2_tau21', 'fj_2_sj1_rawmass', 'fj_2_sj2_rawmass', 'fj_2_ntracks_sv12', 'fj_2_sj1_sv1_pt', 'fj_2_sj2_sv1_pt']
        elif self.jetType == 'ak15':
            self._jetConeSize = 1.5
            self._fj_name = 'AK15Puppi'
            self._sj_name = 'AK15PuppiSubJet'
            self._fj_gen_name = 'GenJetAK15'
            self._sj_gen_name = 'GenSubJetAK15'
            self._sfbdt_files = [os.path.expandvars('$CMSSW_BASE/src/PhysicsTools/NanoHRTTools/data/sfBDT/ak15/xgb_train_qcd.model.%d' % idx) for idx in range(10)]
            self._sfbdt_vars = ['fj_2_tau21', 'fj_2_sj1_rawmass', 'fj_2_sj2_rawmass', 'fj_2_ntracks_sv12', 'fj_2_sj1_sv1_pt', 'fj_2_sj2_sv1_pt']
        else:
            raise RuntimeError('Jet type %s is not recognized!' % self.jetType)
        print ('Running on %d DATA/MC for %s jets' % (self.year, self.jetType))

        #self.DeepCSV_WP_L = {2016: 0.2219, 2017: 0.2219, 2018: 0.2219}[self.year] ## 2017, 2018 placeholders
        #self.DeepCSV_WP_M = {2016: 0.6321, 2017: 0.4941, 2018: 0.4184}[self.year]
        #self.DeepCSV_WP_T = {2016: 0.8958, 2017: 0.8958, 2018: 0.8958}[self.year] ## 2017, 2018 placeholders
        
        self.DeepCSV_WP_L = {2016: 0.2219, 2017: 0.0521, 2018: 0.1241}[self.year]
        self.DeepCSV_WP_M = {2016: 0.6321, 2017: 0.4941, 2018: 0.4184}[self.year]
        self.DeepCSV_WP_T = {2016: 0.8958, 2017: 0.8001, 2018: 0.7527}[self.year]
        
        self._channel = channel
        self._systOpt = {'jec': False, 'jes': None, 'jes_source': '', 'jer': 'nominal', 'jmr': None, 'met_unclustered': None}
        for k in kwargs:
            self._systOpt[k] = kwargs[k]

        logging.info('Running %s channel with systematics %s', self._channel, str(self._systOpt))

        self.jetmetCorr = JetMETCorrector(self.year,
                                          jetType="AK4PFchs",
                                          jec=self._systOpt['jec'],
                                          jes=self._systOpt['jes'],
                                          jes_source=self._systOpt['jes_source'],
                                          jer=self._systOpt['jer'],
                                          met_unclustered=self._systOpt['met_unclustered'])

        self.fatjetCorr = JetMETCorrector(self.year,
                                          jetType="AK8PFPuppi",
                                          jec=self._systOpt['jec'],
                                          jes=self._systOpt['jes'],
                                          jes_source=self._systOpt['jes_source'],
                                          jer=self._systOpt['jer'],
                                          jmr=self._systOpt['jmr'],
                                          met_unclustered=self._systOpt['met_unclustered'])

        self.subjetCorr = JetMETCorrector(self.year,
                                          jetType="AK4PFPuppi",
                                          jec=self._systOpt['jec'],
                                          jes=self._systOpt['jes'],
                                          jes_source=self._systOpt['jes_source'],
                                          jer=self._systOpt['jer'],
                                          jmr=self._systOpt['jmr'],
                                          met_unclustered=self._systOpt['met_unclustered'])

    def beginJob(self):
        self.jetmetCorr.beginJob()
        self.fatjetCorr.beginJob()
        self.subjetCorr.beginJob()
        self.xgb = XGBEnsemble(self._sfbdt_files, self._sfbdt_vars)

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.isMC = bool(inputTree.GetBranch('genWeight'))
        self.isParticleNetV01 = bool(inputTree.GetBranch(self._fj_name + '_ParticleNetMD_probQCD'))
        self.out = wrappedOutputTree

        self.out.branch("jetR", "F")
        self.out.branch("passmetfilters", "O")
        self.out.branch("l1PreFiringWeight", "F")
        self.out.branch("l1PreFiringWeightUp", "F")
        self.out.branch("l1PreFiringWeightDown", "F")

        # Large-R jets
        self.out.branch("n_fatjet", "I")
        for idx in ([1, 2] if self._channel == 'qcd' else [1]):
            prefix = 'fj_%d_' % idx

            # tagger
            self.out.branch(prefix + "DeepAK8MD_ZHbbvsQCD", "F")
            self.out.branch(prefix + "DeepAK8MD_ZHccvsQCD", "F")
            self.out.branch(prefix + "DeepAK8MD_bbVsLight", "F")
            self.out.branch(prefix + "DeepAK8MD_bbVsTop", "F")
            self.out.branch(prefix + "DeepAK8_ZHbbvsQCD", "F")            
            self.out.branch(prefix + "ParticleNetMD_Xbb", "F")
            self.out.branch(prefix + "ParticleNetMD_Xcc", "F")
            self.out.branch(prefix + "ParticleNetMD_Xqq", "F")
            self.out.branch(prefix + "ParticleNetMD_QCD", "F")
            self.out.branch(prefix + "ParticleNetMD_bbVsLight", "F")
            self.out.branch(prefix + "ParticleNetMD_ccVsLight", "F")
            self.out.branch(prefix + "ParticleNetMD_XqqVsQCD", "F")
            self.out.branch(prefix + "ParticleNetMD_XbbVsQCD", "F")
            self.out.branch(prefix + "ParticleNetMD_XccVsQCD", "F")
            # fatjet
            self.out.branch(prefix + "isH", "F")
            self.out.branch(prefix + "isZ", "F")
            self.out.branch(prefix + "isW", "F")
            self.out.branch(prefix + "dr_H", "F")
            self.out.branch(prefix + "H_dau_pdgid", "I")
            self.out.branch(prefix + "dr_Z", "F")
            self.out.branch(prefix + "Z_dau_pdgid", "I")
            self.out.branch(prefix + "is_lep_overlap", "O")
            self.out.branch(prefix + "pt", "F")
            self.out.branch(prefix + "eta", "F")
            self.out.branch(prefix + "phi", "F")
            self.out.branch(prefix + "energy", "F")
            self.out.branch(prefix + "rawmass", "F")
            self.out.branch(prefix + "sdmass", "F")
            self.out.branch(prefix + "tau21", "F")
            self.out.branch(prefix + "btagcsvv2", "F")
            self.out.branch(prefix + "btagjp", "F")
            self.out.branch(prefix + "nsv", "I")
            self.out.branch(prefix + "nsv_ptgt25", "I")
            self.out.branch(prefix + "nsv_ptgt50", "I")
            self.out.branch(prefix + "ntracks", "I")
            self.out.branch(prefix + "ntracks_sv12", "I")
            self.out.branch(prefix + "deltaR_sj12", "F")
            # subjet #1
            self.out.branch(prefix + "sj1_pt", "F")
            self.out.branch(prefix + "sj1_eta", "F")
            self.out.branch(prefix + "sj1_phi", "F")
            self.out.branch(prefix + "sj1_rawmass", "F")
            self.out.branch(prefix + "sj1_energy", "F")
            self.out.branch(prefix + "sj1_btagdeepcsv", "F")
            self.out.branch(prefix + "sj1_btagcsvv2", "F")
            self.out.branch(prefix + "sj1_btagjp", "F")
            self.out.branch(prefix + "sj1_ntracks", "I")
            self.out.branch(prefix + "sj1_nsv", "I")
            self.out.branch(prefix + "sj1_sv1_pt", "F")
            self.out.branch(prefix + "sj1_sv1_mass", "F")
            self.out.branch(prefix + "sj1_sv1_masscor", "F")
            self.out.branch(prefix + "sj1_sv1_ntracks", "I")
            self.out.branch(prefix + "sj1_sv1_dxy", "F")
            self.out.branch(prefix + "sj1_sv1_dxysig", "F")
            self.out.branch(prefix + "sj1_sv1_dlen", "F")
            self.out.branch(prefix + "sj1_sv1_dlensig", "F")
            self.out.branch(prefix + "sj1_sv1_chi2ndof", "F")
            self.out.branch(prefix + "sj1_sv1_pangle", "F")
            # subjet #2
            self.out.branch(prefix + "sj2_pt", "F")
            self.out.branch(prefix + "sj2_eta", "F")
            self.out.branch(prefix + "sj2_phi", "F")
            self.out.branch(prefix + "sj2_rawmass", "F")
            self.out.branch(prefix + "sj2_energy", "F")
            self.out.branch(prefix + "sj2_btagdeepcsv", "F")
            self.out.branch(prefix + "sj2_btagcsvv2", "F")
            self.out.branch(prefix + "sj2_btagjp", "F")
            self.out.branch(prefix + "sj2_ntracks", "I")
            self.out.branch(prefix + "sj2_nsv", "I")
            self.out.branch(prefix + "sj2_sv1_pt", "F")
            self.out.branch(prefix + "sj2_sv1_mass", "F")
            self.out.branch(prefix + "sj2_sv1_masscor", "F")
            self.out.branch(prefix + "sj2_sv1_ntracks", "I")
            self.out.branch(prefix + "sj2_sv1_dxy", "F")
            self.out.branch(prefix + "sj2_sv1_dxysig", "F")
            self.out.branch(prefix + "sj2_sv1_dlen", "F")
            self.out.branch(prefix + "sj2_sv1_dlensig", "F")
            self.out.branch(prefix + "sj2_sv1_chi2ndof", "F")
            self.out.branch(prefix + "sj2_sv1_pangle", "F")
            self.out.branch(prefix + "sj12_masscor_dxysig", "F")
            # sfBDT
            self.out.branch(prefix + "sfBDT", "F")
            # matching variables
            if self.isMC:
                self.out.branch(prefix + "nbhadrons", "I")
                self.out.branch(prefix + "nchadrons", "I")
                self.out.branch(prefix + "partonflavour", "I")
                self.out.branch(prefix + "sj1_nbhadrons", "I")
                self.out.branch(prefix + "sj1_nchadrons", "I")
                self.out.branch(prefix + "sj1_partonflavour", "I")
                self.out.branch(prefix + "sj2_nbhadrons", "I")
                self.out.branch(prefix + "sj2_nchadrons", "I")
                self.out.branch(prefix + "sj2_partonflavour", "I")

    def correctJetsAndMET(self, event):
        # correct Jets and MET
        event._allJets = Collection(event, "Jet")
        event.met = METObject(event, "METFixEE2017") if self.year == 2017 else METObject(event, "MET")

        event._allFatJets = Collection(event, self._fj_name)
        event.subjets = Collection(event, self._sj_name)  # do not sort subjets after updating!!
        # prevent JetReCalibrator from crashing by setting a dummy jetArea -- this is never used for Puppi jets!
        for sj in event.subjets:
            sj.area = 0.5

        if self.isMC or self._systOpt['jec']:
            rho = event.fixedGridRhoFastjetAll
            # correct AK4 jets and MET
            self.jetmetCorr.setSeed(rndSeed(event, event._allJets))
            self.jetmetCorr.correctJetAndMET(jets=event._allJets, met=event.met, rho=rho,
                                             genjets=Collection(event, 'GenJet') if self.isMC else None,
                                             isMC=self.isMC, runNumber=event.run)
            event._allJets = sorted(event._allJets, key=lambda x: x.pt, reverse=True)  # sort by pt after updating

            # correct fatjets
            self.fatjetCorr.setSeed(rndSeed(event, event._allFatJets))
            self.fatjetCorr.correctJetAndMET(jets=event._allFatJets, met=None, rho=rho,
                                             genjets=Collection(event, self._fj_gen_name) if self.isMC else None,
                                             isMC=self.isMC, runNumber=event.run)
            # correct subjets
            self.subjetCorr.setSeed(rndSeed(event, event.subjets))
            self.subjetCorr.correctJetAndMET(jets=event.subjets, met=None, rho=rho,
                                             genjets=Collection(event, self._sj_gen_name) if self.isMC else None,
                                             isMC=self.isMC, runNumber=event.run)

        # jet mass resolution smearing
        if self.isMC and self._systOpt['jmr']:
            raise NotImplemented

        # link fatjet to subjets and recompute softdrop mass
        for fj in event._allFatJets:
            fj.subjets = get_subjets(fj, event.subjets, ('subJetIdx1', 'subJetIdx2'))
            fj.msoftdrop = get_sdmass(fj.subjets)
#             fj.corr_sdmass = get_corrected_sdmass(fj, fj.subjets)
        event._allFatJets = sorted(event._allFatJets, key=lambda x: x.pt, reverse=True)  # sort by pt

    def selectLeptons(self, event):
        # do lepton selection
        event.looseLeptons = []  # used for jet lepton cleaning and lepton counting

        electrons = Collection(event, "Electron")
        for el in electrons:
            el.etaSC = el.eta + el.deltaEtaSC
            if el.pt > 7 and abs(el.eta) < 2.4 and abs(el.dxy) < 0.05 and abs(el.dz) < 0.2 and el.pfRelIso03_all < 0.4:
                if el.mvaFall17V2noIso_WP90:
                    event.looseLeptons.append(el)

        muons = Collection(event, "Muon")
        for mu in muons:
            if mu.pt > 5 and abs(mu.eta) < 2.4 and abs(mu.dxy) < 0.5 and abs(mu.dz) < 1.0 and mu.pfRelIso04_all < 0.4:
                if mu.looseId:
                    event.looseLeptons.append(mu)

        event.looseLeptons.sort(key=lambda x: x.pt, reverse=True)

    def loadGenHistory(self, event):
        # gen matching
        if not self.isMC:
            return

        try:
            genparts = event.genparts
        except RuntimeError as e:
            genparts = Collection(event, "GenPart")
            for idx, gp in enumerate(genparts):
                if 'dauIdx' not in gp.__dict__:
                    gp.dauIdx = []
                if gp.genPartIdxMother >= 0:
                    mom = genparts[gp.genPartIdxMother]
                    if 'dauIdx' not in mom.__dict__:
                        mom.dauIdx = [idx]
                    else:
                        mom.dauIdx.append(idx)
            event.genparts = genparts

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

        def addDaughters(parton):
            if abs(parton.pdgId) == 6:
                parton.daughters = (parton.genB, genparts[parton.genW.dauIdx[0]], genparts[parton.genW.dauIdx[1]])
            elif abs(parton.pdgId) in (23, 24, 25):
                parton.daughters = (genparts[parton.dauIdx[0]], genparts[parton.dauIdx[1]])

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

        for gp in itertools.chain(event.hadGenTops, event.hadGenWs, event.hadGenZs, event.hadGenHs):
            addDaughters(gp)

    def fillBaseEventInfo(self, event):

        self.out.fillBranch("jetR", self._jetConeSize)

        met_filters = bool(
            event.Flag_goodVertices and
            event.Flag_globalSuperTightHalo2016Filter and
            event.Flag_HBHENoiseFilter and
            event.Flag_HBHENoiseIsoFilter and
            event.Flag_EcalDeadCellTriggerPrimitiveFilter and
            event.Flag_BadPFMuonFilter
#             event.Flag_BadChargedCandidateFilter
            )
        if self.year in (2017, 2018):
            met_filters = met_filters and event.Flag_ecalBadCalibFilterV2
        if not self.isMC:
            met_filters = met_filters and event.Flag_eeBadScFilter
        self.out.fillBranch("passmetfilters", met_filters)

        # L1 prefire weights
        if self.year == 2016 or self.year == 2017:
            self.out.fillBranch("l1PreFiringWeight", event.L1PreFiringWeight_Nom)
            self.out.fillBranch("l1PreFiringWeightUp", event.L1PreFiringWeight_Up)
            self.out.fillBranch("l1PreFiringWeightDown", event.L1PreFiringWeight_Dn)
        else:
            self.out.fillBranch("l1PreFiringWeight", 1.0)
            self.out.fillBranch("l1PreFiringWeightUp", 1.0)
            self.out.fillBranch("l1PreFiringWeightDown", 1.0)

    def _get_filler(self, obj):

        def filler(branch, value, default=0):
            self.out.fillBranch(branch, value if obj else default)

        return filler

    def matchSVToSubjets(self, event, fj):
        assert(len(fj.subjets) == 2)
        drcut = min(0.4, 0.5 * deltaR(*fj.subjets))
        for sj in fj.subjets:
            sj.sv_list = []
            for sv in event.secondary_vertices:
                if deltaR(sv, sj) < drcut:
                    sj.sv_list.append(sv)

    def matchSVToJets(self, event, fj):
        drcut = 0.8
        for ifj in fj:
            ifj.sv_list = []
            for sv in event.secondary_vertices:
                        if deltaR(sv, ifj) < drcut:
                                        ifj.sv_list.append(sv)

    def _matchSVToFatjet(self, event, fj):
        if 'sv_list' in fj.__dict__:
            return
        fj.sv_list = []
        for sv in event.secondary_vertices:
            if deltaR(sv, fj) < self._jetConeSize:
                fj.sv_list.append(sv)

    def fillFatJetInfo(self, event):
        self.out.fillBranch("n_fatjet", len(event.fatjets))
        
        for idx in ([1, 2] if self._channel == 'qcd' else [1]):
            prefix = 'fj_%d_' % idx
            fj = event.fatjets[idx - 1]

            if self.isMC:
                h, dr_h = closest(fj, event.hadGenHs)
                z, dr_z = closest(fj, event.hadGenZs)
                self.out.fillBranch(prefix + "dr_H", dr_h)
                self.out.fillBranch(prefix + "H_dau_pdgid", abs(h.daughters[0].pdgId) if h else 0)
                self.out.fillBranch(prefix + "dr_Z", dr_z)
                self.out.fillBranch(prefix + "Z_dau_pdgid", abs(z.daughters[0].pdgId) if z else 0)

            try:
                self.out.fillBranch(prefix + "DeepAK8MD_ZHbbvsQCD", fj.deepTagMD_ZHbbvsQCD)
                self.out.fillBranch(prefix + "DeepAK8MD_ZHccvsQCD", fj.deepTagMD_ZHccvsQCD)
                self.out.fillBranch(prefix + "DeepAK8MD_bbVsLight", fj.deepTagMD_bbvsLight)
                self.out.fillBranch(prefix + "DeepAK8MD_bbVsTop", (1 / (1 + (fj.deepTagMD_TvsQCD / fj.deepTagMD_HbbvsQCD) * (1 - fj.deepTagMD_HbbvsQCD) / (1 - fj.deepTagMD_TvsQCD))))
                self.out.fillBranch(prefix + "DeepAK8_ZHbbvsQCD", ((fj.deepTag_probZbb + fj.deepTag_probHbb))/(fj.deepTag_probZbb + fj.deepTag_probHbb + fj.deepTag_probQCDbb + fj.deepTag_probQCDb + fj.deepTag_probQCDcc + fj.deepTag_probQCDc + fj.deepTag_probQCDothers) )
            except RuntimeError:
                self.out.fillBranch(prefix + "DeepAK8MD_ZHbbvsQCD", -1)
                self.out.fillBranch(prefix + "DeepAK8MD_ZHccvsQCD", -1)
                self.out.fillBranch(prefix + "DeepAK8MD_bbVsLight", -1)
                self.out.fillBranch(prefix + "DeepAK8MD_bbVsTop", -1)
                self.out.fillBranch(prefix + "DeepAK8_ZHbbvsQCD", -1)

            try:
                self.out.fillBranch(prefix + "ParticleNetMD_Xbb", fj.ParticleNetMD_probXbb)
                self.out.fillBranch(prefix + "ParticleNetMD_Xcc", fj.ParticleNetMD_probXcc)
                self.out.fillBranch(prefix + "ParticleNetMD_Xqq", fj.ParticleNetMD_probXqq)
                self.out.fillBranch(prefix + "ParticleNetMD_XqqVsQCD", convert_prob(fj, ['Xbb','Xcc','Xqq'], prefix='ParticleNetMD_prob'))
                self.out.fillBranch(prefix + "ParticleNetMD_bbVsLight" , convert_prob(fj, ['Xbb','QCDbb'] , ['QCDb','QCDcc','QCDc','QCDothers',] , prefix='ParticleNetMD_prob')) 
                self.out.fillBranch(prefix + "ParticleNetMD_ccVsLight" , convert_prob(fj, ['Xcc','QCDcc'] , ['QCDc','QCDbb','QCDb','QCDothers',] , prefix='ParticleNetMD_prob'))
                if self.isParticleNetV01:
                    self.out.fillBranch(prefix + "ParticleNetMD_QCD", fj.ParticleNetMD_probQCD)
                    self.out.fillBranch(prefix + "ParticleNetMD_XbbVsQCD", convert_prob(fj, ['Xbb'], ['QCD'], prefix='ParticleNetMD_prob'))
                    self.out.fillBranch(prefix + "ParticleNetMD_XccVsQCD", convert_prob(fj, ['Xcc'], ['QCD'], prefix='ParticleNetMD_prob'))
                else:
                    self.out.fillBranch(prefix + "ParticleNetMD_QCD", convert_prob(fj, None, prefix='ParticleNetMD_prob'))
                    self.out.fillBranch(prefix + "ParticleNetMD_XbbVsQCD", convert_prob(fj, ['Xbb'], prefix='ParticleNetMD_prob'))
                    self.out.fillBranch(prefix + "ParticleNetMD_XccVsQCD", convert_prob(fj, ['Xcc'], prefix='ParticleNetMD_prob'))
            except RuntimeError:
                self.out.fillBranch(prefix + "ParticleNetMD_Xbb", -1)
                self.out.fillBranch(prefix + "ParticleNetMD_Xcc", -1)
                self.out.fillBranch(prefix + "ParticleNetMD_Xqq", -1)
                self.out.fillBranch(prefix + "ParticleNetMD_QCD", -1)
                self.out.fillBranch(prefix + "ParticleNetMD_HbbVsQCD", -1)
                self.out.fillBranch(prefix + "ParticleNetMD_HccVsQCD", -1)

            #if self.isMC and isSignal:
            if self.isMC:    
                h, _ = closest(fj, event.hadGenHs)
                z, _ = closest(fj, event.hadGenZs)
                w, _ = closest(fj, event.hadGenWs)
                dr_h, dr_z, dr_w = 999., 999., 999.;
                if h:
                    dr_h = deltaR(fj, h)
                if z:
                    dr_z = deltaR(fj, z)
                if w:
                    dr_w = deltaR(fj, w)
                self.out.fillBranch(prefix + "isH", dr_h)
                self.out.fillBranch(prefix + "isZ", dr_z)
                self.out.fillBranch(prefix + "isW", dr_w)

            self.out.fillBranch(prefix + "is_lep_overlap", closest(fj, event.looseLeptons)[1] < self._jetConeSize)
            self.out.fillBranch(prefix + "pt", fj.pt)
            self.out.fillBranch(prefix + "eta", fj.eta)
            self.out.fillBranch(prefix + "phi", fj.phi)
            self.out.fillBranch(prefix + "energy", fj.p4().E())
            self.out.fillBranch(prefix + "rawmass", fj.mass)
            self.out.fillBranch(prefix + "sdmass", fj.msoftdrop)
            self.out.fillBranch(prefix + "tau21", fj.tau2 / fj.tau1 if fj.tau1 > 0 else 99)
            self.out.fillBranch(prefix + "btagcsvv2", fj.btagCSVV2)
            try:
                self.out.fillBranch(prefix + "btagjp", fj.btagJP)
            except RuntimeError:
                self.out.fillBranch(prefix + "btagjp", -1)

            self._matchSVToFatjet(event, fj)
            nsv_ptgt25_   = 0
            nsv_ptgt50_   = 0
            ntracks_      = 0
            ntracks_sv12_ = 0
            for isv, sv in enumerate(fj.sv_list):
                ntracks_ += sv.ntracks
                if isv<2:
                    ntracks_sv12_ += sv.ntracks
                if sv.pt>25.:
                    nsv_ptgt25_ += 1
                if sv.pt>50.:
                    nsv_ptgt50_ += 1 
            self.out.fillBranch(prefix + "nsv", len(fj.sv_list))
            self.out.fillBranch(prefix + "nsv_ptgt25"   , nsv_ptgt25_)
            self.out.fillBranch(prefix + "nsv_ptgt50"   , nsv_ptgt50_)
            self.out.fillBranch(prefix + "ntracks"      , ntracks_)
            self.out.fillBranch(prefix + "ntracks_sv12" , ntracks_sv12_)

            assert(len(fj.subjets) == 2)
            self.out.fillBranch(prefix + "deltaR_sj12", deltaR(*fj.subjets[:2]))
            for idx_sj, sj in enumerate(fj.subjets):
                prefix_sj = prefix + 'sj%d_' % (idx_sj + 1)
                self.out.fillBranch(prefix_sj + "pt", sj.pt)
                self.out.fillBranch(prefix_sj + "eta", sj.eta)
                self.out.fillBranch(prefix_sj + "phi", sj.phi)
                self.out.fillBranch(prefix_sj + "energy", sj.p4().E())
                self.out.fillBranch(prefix_sj + "rawmass", sj.mass)
                self.out.fillBranch(prefix_sj + "btagcsvv2", sj.btagCSVV2)
                try:
                    self.out.fillBranch(prefix_sj + "btagdeepcsv", sj.btagDeepB)
                except RuntimeError:
                    self.out.fillBranch(prefix_sj + "btagdeepcsv", -1)
                try:
                    self.out.fillBranch(prefix_sj + "btagjp", sj.btagJP)
                except RuntimeError:
                    self.out.fillBranch(prefix_sj + "btagjp", -1)

                self.out.fillBranch(prefix_sj + "ntracks" , sum([sv.ntracks for sv in sj.sv_list]))
                self.out.fillBranch(prefix_sj + "nsv", len(sj.sv_list))
                sv = sj.sv_list[0] if len(sj.sv_list) else _NullObject()
                fill_sv = self._get_filler(sv)  # wrapper, fill default value if sv=None
                fill_sv(prefix_sj + "sv1_pt", sv.pt)
                fill_sv(prefix_sj + "sv1_mass", sv.mass)
                fill_sv(prefix_sj + "sv1_masscor", corrected_svmass(sv) if sv else 0)
                fill_sv(prefix_sj + "sv1_ntracks", sv.ntracks)
                fill_sv(prefix_sj + "sv1_dxy", sv.dxy)
                fill_sv(prefix_sj + "sv1_dxysig", sv.dxySig)
                fill_sv(prefix_sj + "sv1_dlen", sv.dlen)
                fill_sv(prefix_sj + "sv1_dlensig", sv.dlenSig)
                fill_sv(prefix_sj + "sv1_chi2ndof", sv.chi2)
                fill_sv(prefix_sj + "sv1_pangle", sv.pAngle)

            sj1, sj2 = fj.subjets
            try:
                sv1, sv2 = sj1.sv_list[0], sj2.sv_list[0]
                sv = sv1 if sv1.dxySig > sv2.dxySig else sv2
                self.out.fillBranch(prefix + "sj12_masscor_dxysig", corrected_svmass(sv) if sv else 0)
            except IndexError:
                # if len(sv_list) == 0
                self.out.fillBranch(prefix + "sj12_masscor_dxysig", 0)

            # sfBDT
            sfbdt_inputs = {k: self.out._branches[k.replace('fj_2_', prefix)].buff[0] for k in self._sfbdt_vars}
            self.out.fillBranch(prefix + "sfBDT", self.xgb.eval(sfbdt_inputs, model_idx=(event.event % 10)))

            # matching variables
            if self.isMC:
                self.out.fillBranch(prefix + "nbhadrons", fj.nBHadrons)
                self.out.fillBranch(prefix + "nchadrons", fj.nCHadrons)
                self.out.fillBranch(prefix + "sj1_nbhadrons", sj1.nBHadrons)
                self.out.fillBranch(prefix + "sj1_nchadrons", sj1.nCHadrons)
                self.out.fillBranch(prefix + "sj2_nbhadrons", sj2.nBHadrons)
                self.out.fillBranch(prefix + "sj2_nchadrons", sj2.nCHadrons)
                try:
                    self.out.fillBranch(prefix + "partonflavour", fj.partonFlavour)
                    self.out.fillBranch(prefix + "sj1_partonflavour", sj1.partonFlavour)
                    self.out.fillBranch(prefix + "sj2_partonflavour", sj2.partonFlavour)
                except:
                    self.out.fillBranch(prefix + "partonflavour", -1)
                    self.out.fillBranch(prefix + "sj1_partonflavour", -1)
                    self.out.fillBranch(prefix + "sj2_partonflavour", -1)
