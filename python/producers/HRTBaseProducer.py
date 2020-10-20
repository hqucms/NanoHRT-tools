# TESTING
import os
import numpy as np
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.tools import deltaPhi, deltaR, closest

from PhysicsTools.NanoHRTTools.helpers.jetmetCorrector import JetMETCorrector
from PhysicsTools.NanoHRTTools.helpers.nnHelper import convert_prob
from PhysicsTools.NanoHRTTools.helpers.ak8MassCorrectionHelper import get_corrected_sdmass
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
        self.year =  int(kwargs['year'])
        self.jetType = kwargs.get('jetType', 'ak8').lower()
#        self.jetType = kwargs.get('jetType', 'ak8').lower()
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
        print ('Running on %d DATA/MC' % (self.year))

        self._channel = channel
        self._systOpt = {'jec': False, 'jes': None, 'jes_source': '', 'jer': 'nominal', 'jmr': None, 'met_unclustered': None}
        for k in kwargs:
            self._systOpt[k] = kwargs[k]

        logging.info('Running %s channel with systematics %s', self._channel, str(self._systOpt))

        self.DeepCSV_WP_M = {2016: 0.6321, 2017: 0.4941, 2018: 0.4184}[self.year]

        self.jetmetCorr = JetMETCorrector(self.year,
                                          jetType="AK4PFchs",
                                          jec=self._systOpt['jec'],
                                          jes=self._systOpt['jes'],
                                          jes_source=self._systOpt['jes_source'],
                                          jer=self._systOpt['jer'],
                                          met_unclustered=self._systOpt['met_unclustered'])

        self.ak8Corr = JetMETCorrector(self.year,
                                       jetType="AK8PFPuppi",
                                       jec=self._systOpt['jec'],
                                       jes=self._systOpt['jes'],
                                       jes_source=self._systOpt['jes_source'],
                                       jer=self._systOpt['jer'],
                                       jmr=self._systOpt['jmr'],
                                       met_unclustered=self._systOpt['met_unclustered'])

        self.ak8SubjetCorr = JetMETCorrector(self.year,
                                             jetType="AK4PFPuppi",
                                             jec=self._systOpt['jec'],
                                             jes=self._systOpt['jes'],
                                             jes_source=self._systOpt['jes_source'],
                                             jer=self._systOpt['jer'],
                                             jmr=self._systOpt['jmr'],
                                             met_unclustered=self._systOpt['met_unclustered'])

        self._n2helper = N2DDTHelper(os.path.expandvars('$CMSSW_BASE/src/PhysicsTools/NanoHRTTools/data/N2DDT/OutputAK82016v13.root'))

    def beginJob(self):
        self.jetmetCorr.beginJob()
        self.ak8Corr.beginJob()
        self.ak8SubjetCorr.beginJob()

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.isMC = bool(inputTree.GetBranch('genWeight'))
        self.hasParticleNet = bool(inputTree.GetBranch('FatJet_ParticleNet_probQCDbb'))
        self.out = wrappedOutputTree

        self.out.branch("passmetfilters", "O")

        prefix = self.jetType+'_1_'
        
        # Large-R jets
        self.out.branch("n_"+self.jetType, "I")

        self.out.branch(prefix + "pt", "F")
        self.out.branch(prefix + "eta", "F")
        self.out.branch(prefix + "phi", "F")
        self.out.branch(prefix + "mass", "F")
        self.out.branch(prefix + "corr_sdmass", "F")
        if self.jetType=="ak8":
            self.out.branch(prefix + "n2b1ddt", "F")
            self.out.branch(prefix + "n2b1", "F")
            self.out.branch(prefix + "n3b1", "F")
        self.out.branch(prefix + "tau1", "F")
        self.out.branch(prefix + "tau2", "F")
        self.out.branch(prefix + "tau3", "F")
#        self.out.branch(prefix + "DeepAK8_TvsQCD", "F")
#        self.out.branch(prefix + "DeepAK8_WvsQCD", "F")
#        self.out.branch(prefix + "DeepAK8_ZvsQCD", "F")
        if self.jetType=="ak8":
            self.out.branch(prefix + "DeepAK8MD_TvsQCD", "F")
            self.out.branch(prefix + "DeepAK8MD_WvsQCD", "F")
            self.out.branch(prefix + "DeepAK8MD_ZvsQCD", "F")
            self.out.branch(prefix + "DeepAK8MD_ZHbbvsQCD", "F")
            self.out.branch(prefix + "DeepAK8MD_ZHccvsQCD", "F")
            self.out.branch(prefix + "DeepAK8MD_bbvsLight", "F")
            self.out.branch(prefix + "DeepAK8MD_ccvsLight", "F")

        if self.hasParticleNet:
#            self.out.branch(prefix + "ParticleNet_TvsQCD", "F")
#            self.out.branch(prefix + "ParticleNet_WvsQCD", "F")
#            self.out.branch(prefix + "ParticleNet_ZvsQCD", "F")
            self.out.branch(prefix + "ParticleNetMD_Xbb", "F")
            self.out.branch(prefix + "ParticleNetMD_Xcc", "F")
            self.out.branch(prefix + "ParticleNetMD_Xqq", "F")
            self.out.branch(prefix + "ParticleNetMD_QCD", "F")
            self.out.branch(prefix + "ParticleNetMD_XbbVsQCD", "F")
            self.out.branch(prefix + "ParticleNetMD_XccVsQCD", "F")
            self.out.branch(prefix + "ParticleNetMD_XccOrXqqVsQCD", "F")
            self.out.branch(prefix + "ParticleNetMD_bbVsNonbb", "F")
            self.out.branch(prefix + "ParticleNetMD_ccVsNoncc", "F")

        # matching variables
        if self.isMC:
            self.out.branch(prefix + "dr_fj_top_b", "F")
            self.out.branch(prefix + "dr_fj_top_wqmax", "F")
            self.out.branch(prefix + "dr_fj_top_wqmin", "F")
            self.out.branch(prefix + "dr_fj_tophad", "F")            
            self.out.branch(prefix + "dr_fj_Hhad", "F")
            self.out.branch(prefix + "dr_fj_Hcc", "F")
            self.out.branch(prefix + "dr_fj_Zhad", "F")
            self.out.branch(prefix + "dr_fj_Whad", "F")
            self.out.branch(prefix + "dr_fj_Wcx", "F")
            self.out.branch(prefix + "dr_fj_Wux", "F")
            self.out.branch(prefix + "dr_fj_WcxFromTop", "F")
            self.out.branch(prefix + "dr_fj_WuxFromTop", "F")
            
        # gen info about tops and W's
        if self.isMC:
            self.out.branch("n_hadGenWs", "F")        
            self.out.branch("n_hadGenWsToCharm", "F")        
            self.out.branch("n_hadGenWsNoCharm", "F")      
            self.out.branch("n_hadGenTops", "F")  
            self.out.branch("n_hadGenTopsToCharm", "F")        
            self.out.branch("n_hadGenTopsNoCharm", "F")        


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

            # correct AK8 fatjets
            if self.jetType == "ak8":
                self.ak8Corr.setSeed(rndSeed(event, event._allFatJets))
                self.ak8Corr.correctJetAndMET(jets=event._allFatJets, met=None, rho=rho,
                                          genjets=Collection(event, 'GenJetAK8') if self.isMC else None,
                                          isMC=self.isMC, runNumber=event.run)
            # correct AK8 subjets
            if self.jetType == "ak8":
                self.ak8SubjetCorr.setSeed(rndSeed(event, event.subjets))
                self.ak8SubjetCorr.correctJetAndMET(jets=event.subjets, met=None, rho=rho,
                                                genjets=Collection(event, 'SubGenJetAK8') if self.isMC else None,
                                                isMC=self.isMC, runNumber=event.run)

        # jet mass resolution smearing
        if self.isMC and self._systOpt['jmr']:
            raise NotImplemented

        # link fatjet to subjets and recompute softdrop mass
        for fj in event._allFatJets:
            fj.subjets = get_subjets(fj, event.subjets, ('subJetIdx1', 'subJetIdx2'))
            fj.msoftdrop = get_sdmass(fj.subjets)
            fj.corr_sdmass = get_corrected_sdmass(fj, fj.subjets)
            if self.jetType=="ak8":
                fj.n2b1ddt = self._n2helper.transform(fj.n2b1, pt=fj.pt, msd=fj.corr_sdmass)
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

    def fillFatJetInfo(self, event, fillGenMatching=False):
        # fill AK8 / AK15 info
        self.out.fillBranch("n_"+self.jetType, len(event.fatjets))
        fj = event.fatjets[0] if len(event.fatjets) > 0 else _NullObject()
        fillBranchFJ = self._get_filler(fj)  # wrapper, fill default value if ak8=None
        prefix = self.jetType+'_1_'
        fillBranchFJ(prefix + "pt", fj.pt)
        fillBranchFJ(prefix + "eta", fj.eta)
        fillBranchFJ(prefix + "phi", fj.phi)
        fillBranchFJ(prefix + "mass", fj.msoftdrop)
        fillBranchFJ(prefix + "corr_sdmass", fj.corr_sdmass)
        if self.jetType=="ak8":
            fillBranchFJ(prefix + "n2b1ddt", fj.n2b1ddt)
            fillBranchFJ(prefix + "n2b1", fj.n2b1)
            fillBranchFJ(prefix + "n3b1", fj.n3b1)
        fillBranchFJ(prefix + "tau1", fj.tau1)
        fillBranchFJ(prefix + "tau2", fj.tau2)
        fillBranchFJ(prefix + "tau3", fj.tau3)
        if self.jetType=="ak8":
#            fillBranchFJ(prefix + "DeepAK8_TvsQCD", fj.deepTag_TvsQCD)
#            fillBranchFJ(prefix + "DeepAK8_WvsQCD", fj.deepTag_WvsQCD)
#            fillBranchFJ(prefix + "DeepAK8_ZvsQCD", fj.deepTag_ZvsQCD)
            fillBranchFJ(prefix + "DeepAK8MD_TvsQCD", fj.deepTagMD_TvsQCD)
            fillBranchFJ(prefix + "DeepAK8MD_WvsQCD", fj.deepTagMD_WvsQCD)
            fillBranchFJ(prefix + "DeepAK8MD_ZvsQCD", fj.deepTagMD_ZvsQCD)
            fillBranchFJ(prefix + "DeepAK8MD_ZHbbvsQCD", fj.deepTagMD_ZHbbvsQCD)
            fillBranchFJ(prefix + "DeepAK8MD_ZHccvsQCD", fj.deepTagMD_ZHccvsQCD)
            fillBranchFJ(prefix + "DeepAK8MD_bbvsLight", fj.deepTagMD_bbvsLight)
            fillBranchFJ(prefix + "DeepAK8MD_ccvsLight", fj.deepTagMD_ccvsLight)


        if self.hasParticleNet:
#            fillBranchFJ(prefix + "ParticleNet_TvsQCD", convert_prob(ak8, ['Tbcq', 'Tbqq'], prefix='ParticleNet_prob'))
#            fillBranchFJ(prefix + "ParticleNet_WvsQCD", convert_prob(ak8, ['Wcq', 'Wqq'], prefix='ParticleNet_prob'))
#            fillBranchFJ(prefix + "ParticleNet_ZvsQCD", convert_prob(ak8, ['Zbb', 'Zcc', 'Zqq'], prefix='ParticleNet_prob'))
            fillBranchFJ(prefix + "ParticleNetMD_Xbb", fj.ParticleNetMD_probXbb)
            fillBranchFJ(prefix + "ParticleNetMD_Xcc", fj.ParticleNetMD_probXcc)
            fillBranchFJ(prefix + "ParticleNetMD_Xqq", fj.ParticleNetMD_probXqq)
            fillBranchFJ(prefix + "ParticleNetMD_QCD", convert_prob(fj, None, prefix='ParticleNetMD_prob'))
            fillBranchFJ(prefix + "ParticleNetMD_XbbVsQCD", convert_prob(fj, ['Xbb'], prefix='ParticleNetMD_prob'))
            fillBranchFJ(prefix + "ParticleNetMD_XccVsQCD", convert_prob(fj, ['Xcc'], prefix='ParticleNetMD_prob'))
            fillBranchFJ(prefix + "ParticleNetMD_XccOrXqqVsQCD", convert_prob(fj, ['Xcc', 'Xqq'], prefix='ParticleNetMD_prob')) # WvsQCD renamed
            fillBranchFJ(prefix + "ParticleNetMD_bbVsNonbb" , convert_prob(fj, ['Xbb','QCDbb'] , ['QCDb','QCDcc','QCDc','QCDothers',] , prefix='ParticleNetMD_prob')) # bbVsLight renamed
            fillBranchFJ(prefix + "ParticleNetMD_ccVsNoncc" , convert_prob(fj, ['Xcc','QCDcc'] , ['QCDc','QCDbb','QCDb','QCDothers',] , prefix='ParticleNetMD_prob')) # ccVsLight renamed

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

        def getmindr(event, jet, genParticles, pickGenW=False):
            dr = 999
            if jet:
                gp, _ = closest(jet, genParticles)
                if gp:
                    if pickGenW:
                        gp = gp.genW                
                    dr = deltaR(jet, gp)
            return dr

        if self.isMC and fillGenMatching:
            
            drfj = drmatch(event, fj)
            self.out.fillBranch(prefix + "dr_fj_top_b", drfj[0])
            self.out.fillBranch(prefix + "dr_fj_top_wqmax", drfj[1])
            self.out.fillBranch(prefix + "dr_fj_top_wqmin", drfj[2])

            self.out.fillBranch(prefix + "dr_fj_tophad", getmindr(event, fj, event.hadGenTops))            
            self.out.fillBranch(prefix + "dr_fj_Hhad", getmindr(event, fj, event.hadGenHs))
            self.out.fillBranch(prefix + "dr_fj_Hcc", getmindr(event, fj, event.hadGenHsToCharm))
            self.out.fillBranch(prefix + "dr_fj_Zhad", getmindr(event, fj, event.hadGenZs))
            self.out.fillBranch(prefix + "dr_fj_Whad", getmindr(event, fj, event.hadGenWs))
            self.out.fillBranch(prefix + "dr_fj_Wcx", getmindr(event, fj, event.hadGenWsToCharm))
            self.out.fillBranch(prefix + "dr_fj_Wux", getmindr(event, fj, event.hadGenWsNoCharm))
            self.out.fillBranch(prefix + "dr_fj_WcxFromTop", getmindr(event, fj, event.hadGenTopsToCharm, pickGenW=True))
            self.out.fillBranch(prefix + "dr_fj_WuxFromTop", getmindr(event, fj, event.hadGenTopsNoCharm, pickGenW=True))


    def fillHadGenWsAndTops(self, event, fillGenMatching=False):

        if self.isMC and fillGenMatching:
            self.out.fillBranch("n_hadGenWs", len(event.hadGenWs))
            self.out.fillBranch("n_hadGenWsToCharm", len(event.hadGenWsToCharm))
            self.out.fillBranch("n_hadGenWsNoCharm", len(event.hadGenWsNoCharm))
            self.out.fillBranch("n_hadGenTops", len(event.hadGenTops))
            self.out.fillBranch("n_hadGenTopsToCharm", len(event.hadGenTopsToCharm))
            self.out.fillBranch("n_hadGenTopsNoCharm", len(event.hadGenTopsNoCharm))
