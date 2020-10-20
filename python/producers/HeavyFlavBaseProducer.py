import math
import numpy as np
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.tools import deltaPhi, deltaR, closest

from PhysicsTools.NanoHRTTools.helpers.jetmetCorrector import JetMETCorrector
from PhysicsTools.NanoHRTTools.helpers.nnHelper import convert_prob

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
        self.jetType = kwargs.get('jetType', 'ak15').lower()
        if self.jetType == 'ak8':
            self._jetConeSize = 0.8
            self._fj_name = 'FatJet'
            self._sj_name = 'SubJet'
            self._fj_gen_name = 'GenJetAK8'
            self._sj_gen_name = 'SubGenJetAK8'
        elif self.jetType == 'ak15':
            self._jetConeSize = 1.5
            self._fj_name = 'AK15Puppi'
            self._sj_name = 'AK15PuppiSubJet'
            self._fj_gen_name = 'GenJetAK15'
            self._sj_gen_name = 'GenSubJetAK15'
        else:
            raise RuntimeError('Jet type %s is not recognized!' % self.jetType)
        print ('Running on %d DATA/MC for %s jets' % (self.year, self.jetType))

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

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.isMC = bool(inputTree.GetBranch('genWeight'))
        self.out = wrappedOutputTree
        self.out.branch("jetR", "F")
        self.out.branch("passmetfilters", "O")
        
        
        # Large-R jets
        self.out.branch("n_fatjet", "I")
        for idx in ([1, 2] if self._channel == 'qcd' else [1]):
            prefix = 'fj_%d_' % idx
        
            # tagger
            self.out.branch(prefix + "DeepAK8MD_ZHbbvsQCD", "F")
            self.out.branch(prefix + "DeepAK8MD_ZHccvsQCD", "F")
            self.out.branch(prefix + "DeepAK8MD_bbVsLight", "F")
            self.out.branch(prefix + "DeepAK8MD_bbVsTop", "F")
            self.out.branch(prefix + "ParticleNetMD_Xbb", "F")
            self.out.branch(prefix + "ParticleNetMD_Xcc", "F")
            self.out.branch(prefix + "ParticleNetMD_Xqq", "F")
            self.out.branch(prefix + "ParticleNetMD_QCD", "F")
            self.out.branch(prefix + "ParticleNetMD_XbbVsQCD", "F")
            self.out.branch(prefix + "ParticleNetMD_XccVsQCD", "F")
            self.out.branch(prefix + "ParticleNetMD_bbVsLight", "F")
            self.out.branch(prefix + "ParticleNetMD_ccVsLight", "F")
            # fatjet
            self.out.branch(prefix + "isH", "F")
            self.out.branch(prefix + "isZ", "F")
            self.out.branch(prefix + "is_lep_overlap", "O")
            self.out.branch(prefix + "pt", "F")
            self.out.branch(prefix + "eta", "F")
            self.out.branch(prefix + "phi", "F")
            self.out.branch(prefix + "energy", "F")
            self.out.branch(prefix + "rawmass", "F")
            self.out.branch(prefix + "sdmass", "F")
            self.out.branch(prefix + "tau21", "F")
            self.out.branch(prefix + "btagcsvv2", "F")
            self.out.branch(prefix + "btagcmva", "F")
            self.out.branch(prefix + "btagjp", "F")
            self.out.branch(prefix + "nsv", "I")
            self.out.branch(prefix + "nsv_ptgt25", "I")
            self.out.branch(prefix + "nsv_ptgt50", "I")
            self.out.branch(prefix + "ntracks", "I")
            self.out.branch(prefix + "ntracks_sv12", "I")
            self.out.branch(prefix + "deltar_ak4", "F")  
            self.out.branch(prefix + "ak4count_maxdr0p6", "F")  
            self.out.branch(prefix + "ak4count_maxdr0p8", "F")  
            self.out.branch(prefix + "ak4count_maxdr1p5", "F")                  
            self.out.branch(prefix + "deltar_sj1_sj2", "F")                  
                      
            # subjet #1
            self.out.branch(prefix + "sj1_pt", "F")
            self.out.branch(prefix + "sj1_eta", "F")
            self.out.branch(prefix + "sj1_phi", "F")
            self.out.branch(prefix + "sj1_rawmass", "F")
            self.out.branch(prefix + "sj1_energy", "F")
            self.out.branch(prefix + "sj1_sdmass", "F")
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
            self.out.branch(prefix + "sj1_deltar_ak4", "F")            
            self.out.branch(prefix + "sj1_deltar", "F")            
            # subjet #2
            self.out.branch(prefix + "sj2_pt", "F")
            self.out.branch(prefix + "sj2_eta", "F")
            self.out.branch(prefix + "sj2_phi", "F")
            self.out.branch(prefix + "sj2_rawmass", "F")
            self.out.branch(prefix + "sj2_energy", "F")
            self.out.branch(prefix + "sj2_sdmass", "F")
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
            self.out.branch(prefix + "sj2_deltar_ak4", "F")            
            self.out.branch(prefix + "sj2_deltar", "F")            
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

                self.out.branch("ak8_1_dr_fj_b", "F")
                self.out.branch("ak8_1_dr_fj_wqmax", "F")
                self.out.branch("ak8_1_dr_fj_wqmin", "F")
                self.out.branch("ak8_1_dr_fj_tophad", "F")            
                self.out.branch("ak8_1_dr_fj_Hhad", "F")
                self.out.branch("ak8_1_dr_fj_Hcc", "F")
                self.out.branch("ak8_1_dr_fj_Zcc", "F")
                self.out.branch("ak8_1_dr_fj_Zhad", "F")
                self.out.branch("ak8_1_dr_fj_Whad", "F")
                self.out.branch("ak8_1_dr_fj_Wcx", "F")
                self.out.branch("ak8_1_dr_fj_Wux", "F")
                self.out.branch("ak8_1_dr_fj_WcxFromTop", "F")
                self.out.branch("ak8_1_dr_fj_WuxFromTop", "F")


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
        event.preselLeptons = []  # used for jet lepton cleaning
        event.looseLeptons = []  # used for lepton counting

        electrons = Collection(event, "Electron")
        for el in electrons:
            el.etaSC = el.eta + el.deltaEtaSC
            if el.pt > 7 and abs(el.eta) < 2.4 and abs(el.dxy) < 0.05 and abs(el.dz) < 0.2 and el.pfRelIso03_all < 0.4:
                event.preselLeptons.append(el)
                if el.mvaFall17V2noIso_WP90:
                    event.looseLeptons.append(el)

        muons = Collection(event, "Muon")
        for mu in muons:
            if mu.pt > 5 and abs(mu.eta) < 2.4 and abs(mu.dxy) < 0.5 and abs(mu.dz) < 1.0 and mu.pfRelIso04_all < 0.4:
                event.preselLeptons.append(mu)
                if mu.looseId:
                    event.looseLeptons.append(mu)

        event.preselLeptons.sort(key=lambda x: x.pt, reverse=True)
        event.looseLeptons.sort(key=lambda x: x.pt, reverse=True)

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
            
        def isCharming(gp):
            if len(gp.dauIdx) == 0:
                raise ValueError('Particle has no daughters!')
            for idx in gp.dauIdx:
                if abs(genparts[idx].pdgId) == 4:
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
        event.hadGenWsToCharm = []
        event.hadGenWsNoCharm = []
        event.hadGenTopsToCharm = []
        event.hadGenTopsNoCharm = []
        event.hadGenZs = []
        event.hadGenZsToCharm = []
        event.hadGenHs = []
        event.hadGenHsToCharm = []
        event.hadGenHsNoCharm = []

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
                            if isCharming(genW):
                                event.hadGenTopsToCharm.append(gp)
                            else:
                                event.hadGenTopsNoCharm.append(gp)
                    elif abs(dau.pdgId) in (1, 3, 5):
                        gp.genB = dau
            elif abs(gp.pdgId) == 24:
                event.nGenWs += 1
                if isHadronic(gp):
                    event.hadGenWs.append(gp)
                    if isCharming(gp):
                        event.hadGenWsToCharm.append(gp)
                    else:
                        event.hadGenWsNoCharm.append(gp)
            elif abs(gp.pdgId) == 23:
                event.nGenZs += 1
                if isHadronic(gp):
                    event.hadGenZs.append(gp)
                    if isCharming(gp):
                        event.hadGenZsToCharm.append(gp)
            elif abs(gp.pdgId) == 25:
                event.nGenHs += 1
                if isHadronic(gp):
                    event.hadGenHs.append(gp)
                    if isCharming(gp):
                        event.hadGenHsToCharm.append(gp)
                    else:
                        event.hadGenHsNoCharm.append(gp)

        event.genparts = genparts

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


    def fillFatJetInfo(self, event, isSignal=False):
        self.out.fillBranch("n_fatjet", len(event.fatjets))

        def drmatch(event, jet):
            dr_b, dr_wq1, dr_wq2 = 999, 999, 999
            if jet:
                higgs, _ = closest(jet, event.hadGenHs)
                if higgs:
                    dr_wq1 = deltaR(jet, event.genparts[higgs.dauIdx[0]])
                    dr_wq2 = deltaR(jet, event.genparts[higgs.dauIdx[1]])
                z, _ = closest(jet, event.hadGenZs)
                if z:
                    dr_wq1 = deltaR(jet, event.genparts[z.dauIdx[0]])
                    dr_wq2 = deltaR(jet, event.genparts[z.dauIdx[1]])
                else:
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

        def getmindr(event, jet, genParticles, pickGenW=False):
            dr = 999
            if jet:
                gp, _ = closest(jet, genParticles)
                if gp:
                    if pickGenW:
                        gp = gp.genW                
                    dr = deltaR(jet, gp)
            return dr

        def drcount(event, jet, countingObjects, maxdr):
            count = 0
            for countingObject in countingObjects:
                dr = deltaR(jet, countingObject)
                if dr < maxdr:
                    count +=1
            return count

                    
        for idx in ([1, 2] if self._channel == 'qcd' else [1]):
            prefix = 'fj_%d_' % idx
            fj = event.fatjets[idx - 1]
        
            try:
                self.out.fillBranch(prefix + "DeepAK8MD_ZHbbvsQCD", fj.deepTagMD_ZHbbvsQCD)
                self.out.fillBranch(prefix + "DeepAK8MD_ZHccvsQCD", fj.deepTagMD_ZHccvsQCD)
                self.out.fillBranch(prefix + "DeepAK8MD_bbVsLight", fj.deepTagMD_bbvsLight)
                self.out.fillBranch(prefix + "DeepAK8MD_bbVsTop", (1 / (1 + (fj.deepTagMD_TvsQCD / fj.deepTagMD_HbbvsQCD) * (1 - fj.deepTagMD_HbbvsQCD) / (1 - fj.deepTagMD_TvsQCD))))
            except RuntimeError:
                self.out.fillBranch(prefix + "DeepAK8MD_ZHbbvsQCD", -1)
                self.out.fillBranch(prefix + "DeepAK8MD_ZHccvsQCD", -1)
                self.out.fillBranch(prefix + "DeepAK8MD_bbVsLight", -1)
                self.out.fillBranch(prefix + "DeepAK8MD_bbVsTop", -1)

            try:
                self.out.fillBranch(prefix + "ParticleNetMD_Xbb", fj.ParticleNetMD_probXbb)
                self.out.fillBranch(prefix + "ParticleNetMD_Xcc", fj.ParticleNetMD_probXcc)
                self.out.fillBranch(prefix + "ParticleNetMD_Xqq", fj.ParticleNetMD_probXqq)
                self.out.fillBranch(prefix + "ParticleNetMD_QCD", convert_prob(fj, None, prefix='ParticleNetMD_prob'))
                self.out.fillBranch(prefix + "ParticleNetMD_XbbVsQCD", convert_prob(fj, ['Xbb'], prefix='ParticleNetMD_prob'))
                self.out.fillBranch(prefix + "ParticleNetMD_XccVsQCD", convert_prob(fj, ['Xcc'], prefix='ParticleNetMD_prob'))
                self.out.fillBranch(prefix + "ParticleNetMD_bbVsLight" , convert_prob(fj, ['Xbb','QCDbb'] , ['QCDb','QCDcc','QCDc','QCDothers',] , prefix='ParticleNetMD_prob')) 
                self.out.fillBranch(prefix + "ParticleNetMD_ccVsLight" , convert_prob(fj, ['Xcc','QCDcc'] , ['QCDc','QCDbb','QCDb','QCDothers',] , prefix='ParticleNetMD_prob'))
            except RuntimeError:
                self.out.fillBranch(prefix + "ParticleNetMD_HbbVsQCD", -1)
                self.out.fillBranch(prefix + "ParticleNetMD_HccVsQCD", -1)

            if self.isMC and isSignal:
                h, _ = closest(fj, event.hadGenHs)
                z, _ = closest(fj, event.hadGenZs)
                dr_h, dr_z = 999., 999.;
                if h:
                    dr_h = deltaR(fj, h)
                if z:
                    dr_z = deltaR(fj, z)
                self.out.fillBranch(prefix + "isH", dr_h)
                self.out.fillBranch(prefix + "isZ", dr_z)

            self.out.fillBranch(prefix + "is_lep_overlap", closest(fj, event.preselLeptons)[1] < self._jetConeSize)
            self.out.fillBranch(prefix + "pt", fj.pt)
            self.out.fillBranch(prefix + "eta", fj.eta)
            self.out.fillBranch(prefix + "phi", fj.phi)
            fj_theta = 2.*math.atan(math.exp(-1.*fj.eta))
            fj_p = (fj.pt)/(math.sin(fj_theta))
            fj_e = math.sqrt((fj.mass*fj.mass) + (fj_p*fj_p))
            self.out.fillBranch(prefix + "energy", fj_e)
            self.out.fillBranch(prefix + "rawmass", fj.mass)
            self.out.fillBranch(prefix + "sdmass", fj.msoftdrop)
            self.out.fillBranch(prefix + "tau21", fj.tau2 / fj.tau1 if fj.tau1 > 0 else 99)
            self.out.fillBranch(prefix + "btagcsvv2", fj.btagCSVV2)
            try:
                self.out.fillBranch(prefix + "btagcmva", fj.btagCMVA)
            except RuntimeError:
                self.out.fillBranch(prefix + "btagcmva", -1)
            try:
                self.out.fillBranch(prefix + "btagjp", fj.btagJP)
            except RuntimeError:
                self.out.fillBranch(prefix + "btagjp", -1)

            self.out.fillBranch(prefix + "nsv", len(fj.sv_list))

            nsv_ptgt25_   = 0
            nsv_ptgt50_   = 0
            ntracks_      = 0
            ntracks_sv12_ = 0
            for isv, sv in enumerate(fj.sv_list):
                ntracks_ += sv.ntracks
                if isv<3:
                    ntracks_sv12_ += sv.ntracks
                if sv.pt>25.:
                    nsv_ptgt25_ += 1
                if sv.pt>50.:
                    nsv_ptgt50_ += 1 

            self.out.fillBranch(prefix + "nsv_ptgt25"   , nsv_ptgt25_)
            self.out.fillBranch(prefix + "nsv_ptgt50"   , nsv_ptgt50_)
            self.out.fillBranch(prefix + "ntracks"      , ntracks_)
            self.out.fillBranch(prefix + "ntracks_sv12" , ntracks_sv12_)

            event.centraljets = []
            for j in event._allJets:
                if not (j.pt > 25 and abs(j.eta) < 2.4 and (j.jetId == 6 or j.jetId == 7) and (j.puId == 6 or j.puId == 7 or j.pt > 50.0)):
                    continue
                event.centraljets.append(j)

            self.out.fillBranch(prefix + "deltar_ak4" ,getmindr(event, fj, event.centraljets))

            self.out.fillBranch(prefix + "ak4count_maxdr0p6", drcount(event, fj, event.centraljets, 0.6))  
            self.out.fillBranch(prefix + "ak4count_maxdr0p8", drcount(event, fj, event.centraljets, 0.8))  
            self.out.fillBranch(prefix + "ak4count_maxdr1p5", drcount(event, fj, event.centraljets, 1.5))                  
                 
            assert(len(fj.subjets) == 2)
            for idx_sj, sj in enumerate(fj.subjets):
                prefix_sj = prefix + 'sj%d_' % (idx_sj + 1)
                self.out.fillBranch(prefix_sj + "pt", sj.pt)
                self.out.fillBranch(prefix_sj + "eta", sj.eta)
                self.out.fillBranch(prefix_sj + "phi", sj.phi)
                sj_theta = 2.*math.atan(math.exp(-1.*sj.eta))
                sj_p = (sj.pt)/(math.sin(sj_theta))
                sj_e = math.sqrt((sj.mass*sj.mass) + (sj_p*sj_p))
                self.out.fillBranch(prefix_sj + "energy", sj_e)
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

                self.out.fillBranch(prefix_sj + "deltar_ak4" ,getmindr(event, sj, event.centraljets))

                self.out.fillBranch(prefix_sj + "deltar", deltaR(fj, sj))


                ntracks_sj_ = 0
                for isjsv, sj_sv in enumerate(sj.sv_list):
                    ntracks_sj_ += sj_sv.ntracks

                self.out.fillBranch(prefix_sj + "ntracks" , ntracks_sj_)
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
            self.out.fillBranch(prefix + "deltar_sj1_sj2" ,deltaR(sj1, sj2))
            try:
                sv1, sv2 = sj1.sv_list[0], sj2.sv_list[0]
                sv = sv1 if sv1.dxySig > sv2.dxySig else sv2
                self.out.fillBranch(prefix + "sj12_masscor_dxysig", corrected_svmass(sv) if sv else 0)
            except IndexError:
                # if len(sv_list) == 0
                self.out.fillBranch(prefix + "sj12_masscor_dxysig", 0)

            # matching variables
            if self.isMC:
                self.out.fillBranch(prefix + "nbhadrons",      fj.nBHadrons)
                self.out.fillBranch(prefix + "nchadrons",      fj.nCHadrons)
                self.out.fillBranch(prefix + "sj1_nbhadrons",     sj1.nBHadrons)
                self.out.fillBranch(prefix + "sj1_nchadrons",     sj1.nCHadrons)
                self.out.fillBranch(prefix + "sj2_nbhadrons",     sj2.nBHadrons)
                self.out.fillBranch(prefix + "sj2_nchadrons",     sj2.nCHadrons)
                try:
                    self.out.fillBranch(prefix + "partonflavour", fj.partonFlavour)
                    self.out.fillBranch(prefix + "sj1_partonflavour", sj1.partonFlavour)
                    self.out.fillBranch(prefix + "sj2_partonflavour", sj2.partonFlavour)
                except RuntimeError:
                    self.out.fillBranch(prefix + "partonflavour", -1)
                    self.out.fillBranch(prefix + "sj1_partonflavour", -1)
                    self.out.fillBranch(prefix + "sj2_partonflavour", -1)

            if self.isMC:
            
                drak8 = drmatch(event, fj)
                self.out.fillBranch("ak8_1_dr_fj_b", drak8[0])
                self.out.fillBranch("ak8_1_dr_fj_wqmax", drak8[1])
                self.out.fillBranch("ak8_1_dr_fj_wqmin", drak8[2])

                self.out.fillBranch("ak8_1_dr_fj_tophad", getmindr(event, fj, event.hadGenTops))            
                self.out.fillBranch("ak8_1_dr_fj_Hhad", getmindr(event, fj, event.hadGenHs))
                self.out.fillBranch("ak8_1_dr_fj_Hcc", getmindr(event, fj, event.hadGenHsToCharm))
                self.out.fillBranch("ak8_1_dr_fj_Zcc", getmindr(event, fj, event.hadGenZsToCharm))
                self.out.fillBranch("ak8_1_dr_fj_Zhad", getmindr(event, fj, event.hadGenZs))
                self.out.fillBranch("ak8_1_dr_fj_Whad", getmindr(event, fj, event.hadGenWs))
                self.out.fillBranch("ak8_1_dr_fj_Wcx", getmindr(event, fj, event.hadGenWsToCharm))
                self.out.fillBranch("ak8_1_dr_fj_Wux", getmindr(event, fj, event.hadGenWsNoCharm))
                self.out.fillBranch("ak8_1_dr_fj_WcxFromTop", getmindr(event, fj, event.hadGenTopsToCharm, pickGenW=True))
                self.out.fillBranch("ak8_1_dr_fj_WuxFromTop", getmindr(event, fj, event.hadGenTopsNoCharm, pickGenW=True))
