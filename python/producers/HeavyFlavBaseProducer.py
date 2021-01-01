import os
import itertools
import numpy as np
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

from ..helpers.utils import deltaR, closest, polarP4, sumP4, get_subjets, corrected_svmass, configLogger
from ..helpers.xgbHelper import XGBEnsemble
from ..helpers.nnHelper import convert_prob, ensemble
from ..helpers.jetmetCorrector import JetMETCorrector, rndSeed

import logging
logger = logging.getLogger('nano')
configLogger('nano', loglevel=logging.INFO)

lumi_dict = {2016: 35.92, 2017: 41.53, 2018: 59.74}


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
        return polarP4(self, eta=None, mass=None)


class HeavyFlavBaseProducer(Module, object):

    def __init__(self, channel, **kwargs):
        self._channel = channel  # 'qcd', 'photon', 'inclusive', 'muon'
        self.year = int(kwargs['year'])
        self.jetType = kwargs.get('jetType', 'ak8').lower()
        self._jmeSysts = {'jec': False, 'jes': None, 'jes_source': '', 'jes_uncertainty_file_prefix': '',
                          'jer': None, 'jmr': None, 'met_unclustered': None, 'smearMET': True, 'applyHEMUnc': False}
        self._opts = {'sfbdt_threshold': -1,
                      'run_tagger': False, 'tagger_versions': ['V02b', 'V02c', 'V02d'],
                      'run_mass_regression': False, 'mass_regression_versions': ['V01a', 'V01b', 'V01c'],
                      'WRITE_CACHE_FILE': True}
        for k in kwargs:
            if k in self._jmeSysts:
                self._jmeSysts[k] = kwargs[k]
            else:
                self._opts[k] = kwargs[k]
        self._needsJMECorr = any([self._jmeSysts['jec'], self._jmeSysts['jes'],
                                  self._jmeSysts['jer'], self._jmeSysts['jmr'],
                                  self._jmeSysts['met_unclustered'], self._jmeSysts['applyHEMUnc']])

        logger.info('Running %s channel for %s jets with JME systematics %s, other options %s',
                    self._channel, self.jetType, str(self._jmeSysts), str(self._opts))

        if self.jetType == 'ak8':
            self._jetConeSize = 0.8
            self._fj_name = 'FatJet'
            self._sj_name = 'SubJet'
            self._fj_gen_name = 'GenJetAK8'
            self._sj_gen_name = 'SubGenJetAK8'
            self._sfbdt_files = [
                os.path.expandvars(
                    '$CMSSW_BASE/src/PhysicsTools/NanoHRTTools/data/sfBDT/ak15/xgb_train_qcd.model.%d' % idx)
                for idx in range(10)]  # FIXME: update to AK8 training
            self._sfbdt_vars = ['fj_2_tau21', 'fj_2_sj1_rawmass', 'fj_2_sj2_rawmass',
                                'fj_2_ntracks_sv12', 'fj_2_sj1_sv1_pt', 'fj_2_sj2_sv1_pt']
        elif self.jetType == 'ak15':
            self._jetConeSize = 1.5
            self._fj_name = 'AK15Puppi'
            self._sj_name = 'AK15PuppiSubJet'
            self._fj_gen_name = 'GenJetAK15'
            self._sj_gen_name = 'GenSubJetAK15'
            self._sfbdt_files = [
                os.path.expandvars(
                    '$CMSSW_BASE/src/PhysicsTools/NanoHRTTools/data/sfBDT/ak15/xgb_train_qcd.model.%d' % idx)
                for idx in range(10)]
            self._sfbdt_vars = ['fj_2_tau21', 'fj_2_sj1_rawmass', 'fj_2_sj2_rawmass',
                                'fj_2_ntracks_sv12', 'fj_2_sj1_sv1_pt', 'fj_2_sj2_sv1_pt']
        else:
            raise RuntimeError('Jet type %s is not recognized!' % self.jetType)

        self._fill_sv = self._channel in ('qcd', 'photon', 'inclusive')

        if self._needsJMECorr:
            self.jetmetCorr = JetMETCorrector(year=self._year, jetType="AK4PFchs", **self._jmeSysts)
            self.fatjetCorr = JetMETCorrector(year=self._year, jetType="AK8PFPuppi", **self._jmeSysts)
            self.subjetCorr = JetMETCorrector(year=self._year, jetType="AK4PFPuppi", **self._jmeSysts)

        if self._opts['run_tagger'] or self._opts['run_mass_regression']:
            from ..helpers.makeInputs import ParticleNetTagInfoMaker
            from ..helpers.runPrediction import ParticleNetJetTagsProducer
            self.tagInfoMaker = ParticleNetTagInfoMaker(
                fatjet_branch=self._fj_name, pfcand_branch='PFCands', sv_branch='SV', jetR=self._jetConeSize)
            prefix = os.path.expandvars('$CMSSW_BASE/src/PhysicsTools/NanoHRTTools/data')
            if self._opts['run_tagger']:
                self.pnTaggers = [ParticleNetJetTagsProducer(
                    '%s/ParticleNet-MD/%s/{version}/particle-net.onnx' % (prefix, self.jetType),
                    '%s/ParticleNet-MD/%s/{version}/preprocess.json' % (prefix, self.jetType),
                    version=ver, cache_suffix='tagger') for ver in self._opts['tagger_versions']]
            if self._opts['run_mass_regression']:
                self.pnMassRegressions = [ParticleNetJetTagsProducer(
                    '%s/MassRegression/%s/{version}/particle_net_regression.onnx' % (prefix, self.jetType),
                    '%s/MassRegression/%s/{version}/preprocess.json' % (prefix, self.jetType),
                    version=ver, cache_suffix='mass') for ver in self._opts['mass_regression_versions']]

        # https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation
        self.DeepCSV_WP_L = {2016: 0.2217, 2017: 0.1522, 2018: 0.1241}[self.year]
        self.DeepCSV_WP_M = {2016: 0.6321, 2017: 0.4941, 2018: 0.4184}[self.year]
        self.DeepCSV_WP_T = {2016: 0.8953, 2017: 0.8001, 2018: 0.7527}[self.year]

    def beginJob(self):
        if self._needsJMECorr:
            self.jetmetCorr.beginJob()
            self.fatjetCorr.beginJob()
            self.subjetCorr.beginJob()
        self.xgb = XGBEnsemble(self._sfbdt_files, self._sfbdt_vars)

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.isMC = bool(inputTree.GetBranch('genWeight'))
        self.isParticleNetV01 = bool(inputTree.GetBranch(self._fj_name + '_ParticleNetMD_probQCD'))

        # remove all possible h5 cache files
        for f in os.listdir('.'):
            if f.endswith('.h5'):
                os.remove(f)

        if self._opts['run_tagger']:
            for p in self.pnTaggers:
                p.load_cache(inputFile)

        if self._opts['run_mass_regression']:
            for p in self.pnMassRegressions:
                p.load_cache(inputFile)

        if self._opts['run_tagger'] or self._opts['run_mass_regression']:
            self.tagInfoMaker.init_file(inputFile, fetch_step=1000)

        self.out = wrappedOutputTree

        # NOTE: branch names must start with a lower case letter
        # check keep_and_drop_output.txt
        self.out.branch("year", "I")
        self.out.branch("lumiwgt", "F")
        self.out.branch("jetR", "F")
        self.out.branch("passmetfilters", "O")
        self.out.branch("l1PreFiringWeight", "F")
        self.out.branch("l1PreFiringWeightUp", "F")
        self.out.branch("l1PreFiringWeightDown", "F")
        self.out.branch("nlep", "I")
        self.out.branch("ht", "F")
        self.out.branch("met", "F")
        self.out.branch("metphi", "F")

        # Large-R jets
        self.out.branch("n_fatjet", "I")
        for idx in ([1, 2] if self._channel == 'qcd' else [1]):
            prefix = 'fj_%d_' % idx

            # fatjet kinematics
            self.out.branch(prefix + "is_qualified", "O")
            self.out.branch(prefix + "pt", "F")
            self.out.branch(prefix + "eta", "F")
            self.out.branch(prefix + "phi", "F")
            self.out.branch(prefix + "rawmass", "F")
            self.out.branch(prefix + "sdmass", "F")
            self.out.branch(prefix + "regressed_mass", "F")
            self.out.branch(prefix + "tau21", "F")
            self.out.branch(prefix + "tau32", "F")
            self.out.branch(prefix + "btagcsvv2", "F")
            self.out.branch(prefix + "btagjp", "F")

            # subjets
            self.out.branch(prefix + "deltaR_sj12", "F")
            self.out.branch(prefix + "sj1_pt", "F")
            self.out.branch(prefix + "sj1_eta", "F")
            self.out.branch(prefix + "sj1_phi", "F")
            self.out.branch(prefix + "sj1_rawmass", "F")
            self.out.branch(prefix + "sj1_btagdeepcsv", "F")
            self.out.branch(prefix + "sj2_pt", "F")
            self.out.branch(prefix + "sj2_eta", "F")
            self.out.branch(prefix + "sj2_phi", "F")
            self.out.branch(prefix + "sj2_rawmass", "F")
            self.out.branch(prefix + "sj2_btagdeepcsv", "F")

            # taggers
            self.out.branch(prefix + "DeepAK8_TvsQCD", "F")
            self.out.branch(prefix + "DeepAK8_WvsQCD", "F")
            self.out.branch(prefix + "DeepAK8_ZvsQCD", "F")
            self.out.branch(prefix + "DeepAK8_ZHbbvsQCD", "F")
            self.out.branch(prefix + "DeepAK8MD_TvsQCD", "F")
            self.out.branch(prefix + "DeepAK8MD_WvsQCD", "F")
            self.out.branch(prefix + "DeepAK8MD_ZvsQCD", "F")
            self.out.branch(prefix + "DeepAK8MD_ZHbbvsQCD", "F")
            self.out.branch(prefix + "DeepAK8MD_ZHccvsQCD", "F")
            self.out.branch(prefix + "DeepAK8MD_bbVsLight", "F")
            self.out.branch(prefix + "DeepAK8MD_bbVsTop", "F")

            self.out.branch(prefix + "ParticleNet_TvsQCD", "F")
            self.out.branch(prefix + "ParticleNet_WvsQCD", "F")
            self.out.branch(prefix + "ParticleNet_ZvsQCD", "F")
            self.out.branch(prefix + "ParticleNetMD_Xbb", "F")
            self.out.branch(prefix + "ParticleNetMD_Xcc", "F")
            self.out.branch(prefix + "ParticleNetMD_Xqq", "F")
            self.out.branch(prefix + "ParticleNetMD_QCD", "F")
            self.out.branch(prefix + "ParticleNetMD_XbbVsQCD", "F")
            self.out.branch(prefix + "ParticleNetMD_XccVsQCD", "F")
            self.out.branch(prefix + "ParticleNetMD_XccOrXqqVsQCD", "F")

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

                # info of the closest hadGenH
                self.out.branch(prefix + "dr_H", "F")
                self.out.branch(prefix + "dr_H_daus", "F")
                self.out.branch(prefix + "H_pt", "F")
                self.out.branch(prefix + "H_decay", "I")

                # info of the closest hadGenZ
                self.out.branch(prefix + "dr_Z", "F")
                self.out.branch(prefix + "dr_Z_daus", "F")
                self.out.branch(prefix + "Z_pt", "F")
                self.out.branch(prefix + "Z_decay", "I")

                # info of the closest hadGenW
                self.out.branch(prefix + "dr_W", "F")
                self.out.branch(prefix + "dr_W_daus", "F")
                self.out.branch(prefix + "W_pt", "F")
                self.out.branch(prefix + "W_decay", "I")

                # info of the closest hadGenTop
                self.out.branch(prefix + "dr_T", "F")
                self.out.branch(prefix + "dr_T_b", "F")
                self.out.branch(prefix + "dr_T_Wq_max", "F")
                self.out.branch(prefix + "dr_T_Wq_min", "F")
                self.out.branch(prefix + "T_Wq_max_pdgId", "I")
                self.out.branch(prefix + "T_Wq_min_pdgId", "I")
                self.out.branch(prefix + "T_pt", "F")

            if self._fill_sv:
                # SV variables
                self.out.branch(prefix + "nsv", "I")
                self.out.branch(prefix + "nsv_ptgt25", "I")
                self.out.branch(prefix + "nsv_ptgt50", "I")
                self.out.branch(prefix + "ntracks", "I")
                self.out.branch(prefix + "ntracks_sv12", "I")

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

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        if self._opts['run_tagger'] and self._opts['WRITE_CACHE_FILE']:
            for p in self.pnTaggers:
                p.update_cache()

        if self._opts['run_mass_regression'] and self._opts['WRITE_CACHE_FILE']:
            for p in self.pnMassRegressions:
                p.update_cache()

        # remove all h5 cache files
        if self._opts['run_tagger'] or self._opts['run_mass_regression']:
            for f in os.listdir('.'):
                if f.endswith('.h5'):
                    os.remove(f)

    def selectLeptons(self, event):
        # do lepton selection
        event.looseLeptons = []  # used for jet lepton cleaning & lepton counting

        electrons = Collection(event, "Electron")
        for el in electrons:
            el.etaSC = el.eta + el.deltaEtaSC
            if el.pt > 10 and abs(el.eta) < 2.5 and abs(el.dxy) < 0.05 and abs(el.dz) < 0.2 \
                    and el.mvaFall17V2noIso_WP90 and el.miniPFRelIso_all < 0.4:
                event.looseLeptons.append(el)

        muons = Collection(event, "Muon")
        for mu in muons:
            if mu.pt > 10 and abs(mu.eta) < 2.4 and abs(mu.dxy) < 0.05 and abs(mu.dz) < 0.2 \
                    and mu.looseId and mu.miniPFRelIso_all < 0.4:
                event.looseLeptons.append(mu)

        event.looseLeptons.sort(key=lambda x: x.pt, reverse=True)

    def correctJetsAndMET(self, event):
        # correct Jets and MET
        event.idx = event._entry if event._tree._entrylist is None else event._tree._entrylist.GetEntry(event._entry)
        event._allJets = Collection(event, "Jet")
        event.met = METObject(event, "METFixEE2017") if self.year == 2017 else METObject(event, "MET")
        event._allFatJets = Collection(event, self._fj_name)
        event.subjets = Collection(event, self._sj_name)  # do not sort subjets after updating!!

        if self._needsJMECorr:
            rho = event.fixedGridRhoFastjetAll
            # correct AK4 jets and MET
            self.jetmetCorr.setSeed(rndSeed(event, event._allJets))
            self.jetmetCorr.correctJetAndMET(jets=event._allJets, lowPtJets=Collection(event, "CorrT1METJet"),
                                             met=event.met, rawMET=METObject(event, "RawMET"),
                                             defaultMET=METObject(event, "MET"),
                                             rho=rho, genjets=Collection(event, 'GenJet') if self.isMC else None,
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
        if self.isMC and self._jmeSysts['jmr']:
            raise NotImplementedError

        # link fatjet to subjets and recompute softdrop mass
        for fj in event._allFatJets:
            fj.subjets = get_subjets(fj, event.subjets, ('subJetIdx1', 'subJetIdx2'))
            fj.msoftdrop = sumP4(*fj.subjets).M()
        event._allFatJets = sorted(event._allFatJets, key=lambda x: x.pt, reverse=True)  # sort by pt

        # select lepton-cleaned jets
        event.fatjets = [fj for fj in event._allFatJets if fj.pt > 200 and abs(fj.eta) < 2.4 and (
            fj.jetId & 2) and closest(fj, event.looseLeptons)[1] >= self._jetConeSize]
        event.ak4jets = [j for j in event._allJets if j.pt > 25 and abs(j.eta) < 2.4 and (
            j.jetId & 4) and closest(j, event.looseLeptons)[1] >= 0.4]
        event.ht = sum([j.pt for j in event.ak4jets])

    def selectSV(self, event):
        event._allSV = Collection(event, "SV")
        event.secondary_vertices = []
        for sv in event._allSV:
            # if sv.ntracks > 2 and abs(sv.dxy) < 3. and sv.dlenSig > 4:
            # if sv.dlenSig > 4:
            if True:
                event.secondary_vertices.append(sv)
        event.secondary_vertices = sorted(event.secondary_vertices, key=lambda x: x.pt, reverse=True)  # sort by pt
        # event.secondary_vertices = sorted(event.secondary_vertices, key=lambda x : x.dxySig, reverse=True)  # sort by dxysig

    def matchSVToFatJets(self, event, fatjets):
        # match SV to fatjets
        for fj in fatjets:
            fj.sv_list = []
            for sv in event.secondary_vertices:
                if deltaR(sv, fj) < self._jetConeSize:
                    fj.sv_list.append(sv)
            # match SV to subjets
            assert(len(fj.subjets) == 2)
            drcut = min(0.4, 0.5 * deltaR(*fj.subjets))
            for sj in fj.subjets:
                sj.sv_list = []
                for sv in event.secondary_vertices:
                    if deltaR(sv, sj) < drcut:
                        sj.sv_list.append(sv)

            fj.nsv_ptgt25 = 0
            fj.nsv_ptgt50 = 0
            fj.ntracks = 0
            fj.ntracks_sv12 = 0
            for isv, sv in enumerate(fj.sv_list):
                fj.ntracks += sv.ntracks
                if isv < 2:
                    fj.ntracks_sv12 += sv.ntracks
                if sv.pt > 25:
                    fj.nsv_ptgt25 += 1
                if sv.pt > 50:
                    fj.nsv_ptgt50 += 1

            # sfBDT & sj12_masscor_dxysig
            sj1, sj2 = fj.subjets
            if len(sj1.sv_list) > 0 and len(sj2.sv_list) > 0:
                sj1_sv, sj2_sv = sj1.sv_list[0], sj2.sv_list[0]
                sfbdt_inputs = {
                    'fj_2_tau21': fj.tau2 / fj.tau1 if fj.tau1 > 0 else 99,
                    'fj_2_sj1_rawmass': sj1.mass,
                    'fj_2_sj2_rawmass': sj2.mass,
                    'fj_2_ntracks_sv12': fj.ntracks_sv12,
                    'fj_2_sj1_sv1_pt': sj1_sv.pt,
                    'fj_2_sj2_sv1_pt': sj2_sv.pt,
                }
                fj.sfBDT = self.xgb.eval(sfbdt_inputs, model_idx=(event.event % 10))
                fj.sj12_masscor_dxysig = corrected_svmass(sj1_sv if sj1_sv.dxySig > sj2_sv.dxySig else sj2_sv)
            else:
                fj.sfBDT = -1
                fj.sj12_masscor_dxysig = 0

    def loadGenHistory(self, event, fatjets):
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

        lepGenTops = []
        hadGenTops = []
        hadGenWs = []
        hadGenZs = []
        hadGenHs = []

        for gp in genparts:
            if gp.statusFlags & (1 << 13) == 0:
                continue
            if abs(gp.pdgId) == 6:
                for idx in gp.dauIdx:
                    dau = genparts[idx]
                    if abs(dau.pdgId) == 24:
                        genW = getFinal(dau)
                        gp.genW = genW
                        if isHadronic(genW):
                            hadGenTops.append(gp)
                        else:
                            lepGenTops.append(gp)
                    elif abs(dau.pdgId) in (1, 3, 5):
                        gp.genB = dau
            elif abs(gp.pdgId) == 24:
                if isHadronic(gp):
                    hadGenWs.append(gp)
            elif abs(gp.pdgId) == 23:
                if isHadronic(gp):
                    hadGenZs.append(gp)
            elif abs(gp.pdgId) == 25:
                if isHadronic(gp):
                    hadGenHs.append(gp)

        for parton in itertools.chain(lepGenTops, hadGenTops):
            parton.daus = (parton.genB, genparts[parton.genW.dauIdx[0]], genparts[parton.genW.dauIdx[1]])
            parton.genW.daus = parton.daus[1:]
        for parton in itertools.chain(hadGenWs, hadGenZs, hadGenHs):
            parton.daus = (genparts[parton.dauIdx[0]], genparts[parton.dauIdx[1]])

        for fj in fatjets:
            fj.genH, fj.dr_H = closest(fj, hadGenHs)
            fj.genZ, fj.dr_Z = closest(fj, hadGenZs)
            fj.genW, fj.dr_W = closest(fj, hadGenWs)
            fj.genT, fj.dr_T = closest(fj, hadGenTops)
            fj.genLepT, fj.dr_LepT = closest(fj, lepGenTops)

    def evalTagger(self, event, jets):
        for j in jets:
            if self._opts['run_tagger']:
                outputs = [p.predict_with_cache(self.tagInfoMaker, event.idx, j.idx, j) for p in self.pnTaggers]
                outputs = ensemble(outputs, np.mean)
                j.pn_Xbb = outputs['probXbb']
                j.pn_Xcc = outputs['probXcc']
                j.pn_Xqq = outputs['probXqq']
                j.pn_QCD = convert_prob(outputs, None, prefix='prob')
            else:
                j.pn_Xbb = j.ParticleNetMD_probXbb
                j.pn_Xcc = j.ParticleNetMD_probXcc
                j.pn_Xqq = j.ParticleNetMD_probXqq
                j.pn_QCD = j.ParticleNetMD_probQCD if self.isParticleNetV01 else convert_prob(
                    j, None, prefix='ParticleNetMD_prob')
            j.pn_XbbVsQCD = convert_prob(j, ['Xbb'], ['QCD'], prefix='pn_')
            j.pn_XccVsQCD = convert_prob(j, ['Xcc'], ['QCD'], prefix='pn_')
            j.pn_XccOrXqqVsQCD = convert_prob(j, ['Xcc', 'Xqq'], ['QCD'], prefix='pn_')

    def evalMassRegression(self, event, jets):
        for j in jets:
            if self._opts['run_mass_regression']:
                outputs = [p.predict_with_cache(self.tagInfoMaker, event.idx, j.idx, j) for p in self.pnMassRegressions]
                j.regressed_mass = ensemble(outputs, np.median)['mass']
            else:
                j.regressed_mass = 0

    def fillBaseEventInfo(self, event):
        self.out.fillBranch("jetR", self._jetConeSize)
        self.out.fillBranch("year", self._year)
        self.out.fillBranch("lumiwgt", lumi_dict[self._year])

        met_filters = bool(
            event.Flag_goodVertices and
            event.Flag_globalSuperTightHalo2016Filter and
            event.Flag_HBHENoiseFilter and
            event.Flag_HBHENoiseIsoFilter and
            event.Flag_EcalDeadCellTriggerPrimitiveFilter and
            event.Flag_BadPFMuonFilter
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

        self.out.fillBranch("nlep", len(event.looseLeptons))
        self.out.fillBranch("ht", event.ht)
        self.out.fillBranch("met", event.met.pt)
        self.out.fillBranch("metphi", event.met.phi)

    def _get_filler(self, obj):

        def filler(branch, value, default=0):
            self.out.fillBranch(branch, value if obj else default)

        return filler

    def fillFatJetInfo(self, event, fatjets):
        self.out.fillBranch("n_fatjet", len(fatjets))

        for idx in ([1, 2] if self._channel == 'qcd' else [1]):
            prefix = 'fj_%d_' % idx

            try:
                fj = fatjets[idx - 1]
            except IndexError:
                # fill zeros if `fatjets` does not contain enough jets
                for b in self.out._branches.keys():
                    if b.startswith(prefix):
                        self.out.fillBranch(b, 0)
                continue

            # fatjet kinematics
            self.out.fillBranch(prefix + "is_qualified", True)
            self.out.fillBranch(prefix + "pt", fj.pt)
            self.out.fillBranch(prefix + "eta", fj.eta)
            self.out.fillBranch(prefix + "phi", fj.phi)
            self.out.fillBranch(prefix + "rawmass", fj.mass)
            self.out.fillBranch(prefix + "sdmass", fj.msoftdrop)
            self.out.fillBranch(prefix + "regressed_mass", fj.regressed_mass)
            self.out.fillBranch(prefix + "tau21", fj.tau2 / fj.tau1 if fj.tau1 > 0 else 99)
            self.out.fillBranch(prefix + "tau32", fj.tau3 / fj.tau2 if fj.tau2 > 0 else 99)
            self.out.fillBranch(prefix + "btagcsvv2", fj.btagCSVV2)
            try:
                self.out.fillBranch(prefix + "btagjp", fj.btagJP)
            except RuntimeError:
                self.out.fillBranch(prefix + "btagjp", -1)

            # subjets
            self.out.fillBranch(prefix + "deltaR_sj12", deltaR(*fj.subjets) if len(fj.subjets) == 2 else 99)
            for idx_sj, sj in enumerate(fj.subjets):
                prefix_sj = prefix + 'sj%d_' % (idx_sj + 1)
                self.out.fillBranch(prefix_sj + "pt", sj.pt)
                self.out.fillBranch(prefix_sj + "eta", sj.eta)
                self.out.fillBranch(prefix_sj + "phi", sj.phi)
                self.out.fillBranch(prefix_sj + "rawmass", sj.mass)
                try:
                    self.out.fillBranch(prefix_sj + "btagdeepcsv", sj.btagDeepB)
                except RuntimeError:
                    self.out.fillBranch(prefix_sj + "btagdeepcsv", -1)

            # taggers
            try:
                # Full
                self.out.fillBranch(prefix + "DeepAK8_TvsQCD", fj.deepTag_TvsQCD)
                self.out.fillBranch(prefix + "DeepAK8_WvsQCD", fj.deepTag_WvsQCD)
                self.out.fillBranch(prefix + "DeepAK8_ZvsQCD", fj.deepTag_ZvsQCD)
                # MD
                self.out.fillBranch(prefix + "DeepAK8MD_TvsQCD", fj.deepTagMD_TvsQCD)
                self.out.fillBranch(prefix + "DeepAK8MD_WvsQCD", fj.deepTagMD_WvsQCD)
                self.out.fillBranch(prefix + "DeepAK8MD_ZvsQCD", fj.deepTagMD_ZvsQCD)
                self.out.fillBranch(prefix + "DeepAK8MD_ZHbbvsQCD", fj.deepTagMD_ZHbbvsQCD)
                self.out.fillBranch(prefix + "DeepAK8MD_ZHccvsQCD", fj.deepTagMD_ZHccvsQCD)
                self.out.fillBranch(prefix + "DeepAK8MD_bbVsLight", fj.deepTagMD_bbvsLight)
                self.out.fillBranch(prefix + "DeepAK8MD_bbVsTop",
                                    (1 / (1 + (fj.deepTagMD_TvsQCD / fj.deepTagMD_HbbvsQCD) * (1 - fj.deepTagMD_HbbvsQCD) / (1 - fj.deepTagMD_TvsQCD))))  # noqa
            except RuntimeError:
                # if no DeepAK8 branches
                self.out.fillBranch(prefix + "DeepAK8_TvsQCD", -1)
                self.out.fillBranch(prefix + "DeepAK8_WvsQCD", -1)
                self.out.fillBranch(prefix + "DeepAK8_ZvsQCD", -1)
                self.out.fillBranch(prefix + "DeepAK8MD_TvsQCD", -1)
                self.out.fillBranch(prefix + "DeepAK8MD_WvsQCD", -1)
                self.out.fillBranch(prefix + "DeepAK8MD_ZvsQCD", -1)
                self.out.fillBranch(prefix + "DeepAK8MD_ZHbbvsQCD", -1)
                self.out.fillBranch(prefix + "DeepAK8MD_ZHccvsQCD", -1)
                self.out.fillBranch(prefix + "DeepAK8MD_bbVsLight", -1)
                self.out.fillBranch(prefix + "DeepAK8MD_bbVsTop", -1)

            try:
                self.out.fillBranch(prefix + "DeepAK8_ZHbbvsQCD",
                                    convert_prob(fj, ['Zbb', 'Hbb'], prefix='deepTag_prob'))
            except RuntimeError:
                # if no DeepAK8 raw probs
                self.out.fillBranch(prefix + "DeepAK8_ZHbbvsQCD", -1)

            # ParticleNet
            try:
                self.out.fillBranch(prefix + "ParticleNet_TvsQCD",
                                    convert_prob(fj, ['Tbcq', 'Tbqq'], prefix='ParticleNet_prob'))
                self.out.fillBranch(prefix + "ParticleNet_WvsQCD",
                                    convert_prob(fj, ['Wcq', 'Wqq'], prefix='ParticleNet_prob'))
                self.out.fillBranch(prefix + "ParticleNet_ZvsQCD",
                                    convert_prob(fj, ['Zbb', 'Zcc', 'Zqq'], prefix='ParticleNet_prob'))
            except RuntimeError:
                # if no nominal ParticleNet
                self.out.fillBranch(prefix + "ParticleNet_TvsQCD", -1)
                self.out.fillBranch(prefix + "ParticleNet_WvsQCD", -1)
                self.out.fillBranch(prefix + "ParticleNet_ZvsQCD", -1)

            # ParticleNet-MD
            self.out.fillBranch(prefix + "ParticleNetMD_Xbb", fj.pn_Xbb)
            self.out.fillBranch(prefix + "ParticleNetMD_Xcc", fj.pn_Xcc)
            self.out.fillBranch(prefix + "ParticleNetMD_Xqq", fj.pn_Xqq)
            self.out.fillBranch(prefix + "ParticleNetMD_QCD", fj.pn_QCD)
            self.out.fillBranch(prefix + "ParticleNetMD_XbbVsQCD", fj.pn_XbbVsQCD)
            self.out.fillBranch(prefix + "ParticleNetMD_XccVsQCD", fj.pn_XccVsQCD)
            self.out.fillBranch(prefix + "ParticleNetMD_XccOrXqqVsQCD", fj.pn_XccOrXqqVsQCD)

            # matching variables
            if self.isMC:
                try:
                    sj1 = fj.subjets[0]
                except IndexError:
                    sj1 = None
                try:
                    sj2 = fj.subjets[1]
                except IndexError:
                    sj2 = None

                self.out.fillBranch(prefix + "nbhadrons", fj.nBHadrons)
                self.out.fillBranch(prefix + "nchadrons", fj.nCHadrons)
                self.out.fillBranch(prefix + "sj1_nbhadrons", sj1.nBHadrons if sj1 else -1)
                self.out.fillBranch(prefix + "sj1_nchadrons", sj1.nCHadrons if sj1 else -1)
                self.out.fillBranch(prefix + "sj2_nbhadrons", sj2.nBHadrons if sj2 else -1)
                self.out.fillBranch(prefix + "sj2_nchadrons", sj2.nCHadrons if sj2 else -1)
                try:
                    self.out.fillBranch(prefix + "partonflavour", fj.partonFlavour)
                    self.out.fillBranch(prefix + "sj1_partonflavour", sj1.partonFlavour if sj1 else -1)
                    self.out.fillBranch(prefix + "sj2_partonflavour", sj2.partonFlavour if sj2 else -1)
                except RuntimeError:
                    self.out.fillBranch(prefix + "partonflavour", -1)
                    self.out.fillBranch(prefix + "sj1_partonflavour", -1)
                    self.out.fillBranch(prefix + "sj2_partonflavour", -1)

                # info of the closest hadGenH
                self.out.fillBranch(prefix + "dr_H", fj.dr_H)
                self.out.fillBranch(prefix + "dr_H_daus",
                                    max([deltaR(fj, dau) for dau in fj.genH.daus]) if fj.genH else 99)
                self.out.fillBranch(prefix + "H_pt", fj.genH.pt if fj.genH else -1)
                self.out.fillBranch(prefix + "H_decay", abs(fj.genH.daus[0].pdgId) if fj.genH else 0)

                # info of the closest hadGenZ
                self.out.fillBranch(prefix + "dr_Z", fj.dr_Z)
                self.out.fillBranch(prefix + "dr_Z_daus",
                                    max([deltaR(fj, dau) for dau in fj.genZ.daus]) if fj.genZ else 99)
                self.out.fillBranch(prefix + "Z_pt", fj.genZ.pt if fj.genZ else -1)
                self.out.fillBranch(prefix + "Z_decay", abs(fj.genZ.daus[0].pdgId) if fj.genZ else 0)

                # info of the closest hadGenW
                self.out.fillBranch(prefix + "dr_W", fj.dr_W)
                self.out.fillBranch(prefix + "dr_W_daus",
                                    max([deltaR(fj, dau) for dau in fj.genW.daus]) if fj.genW else 99)
                self.out.fillBranch(prefix + "W_pt", fj.genW.pt if fj.genW else -1)
                self.out.fillBranch(prefix + "W_decay", abs(fj.genW.daus[0].pdgId) if fj.genW else 0)

                # info of the closest hadGenTop
                drwq1, drwq2 = [deltaR(fj, dau) for dau in fj.genT.genW.daus] if fj.genT else [99, 99]
                wq1_pdgId, wq2_pdgId = [dau.pdgId for dau in fj.genT.genW.daus] if fj.genT else [0, 0]
                if drwq1 < drwq2:
                    drwq1, drwq2 = drwq2, drwq1
                    wq1_pdgId, wq2_pdgId = wq2_pdgId, wq1_pdgId
                self.out.fillBranch(prefix + "dr_T", fj.dr_T)
                self.out.fillBranch(prefix + "dr_T_b", deltaR(fj, fj.genT.genB) if fj.genT else 99)
                self.out.fillBranch(prefix + "dr_T_Wq_max", drwq1)
                self.out.fillBranch(prefix + "dr_T_Wq_min", drwq2)
                self.out.fillBranch(prefix + "T_Wq_max_pdgId", wq1_pdgId)
                self.out.fillBranch(prefix + "T_Wq_min_pdgId", wq2_pdgId)
                self.out.fillBranch(prefix + "T_pt", fj.genT.pt if fj.genT else -1)

            if self._fill_sv:
                # SV variables
                self.out.fillBranch(prefix + "nsv", len(fj.sv_list))
                self.out.fillBranch(prefix + "nsv_ptgt25", fj.nsv_ptgt25)
                self.out.fillBranch(prefix + "nsv_ptgt50", fj.nsv_ptgt50)
                self.out.fillBranch(prefix + "ntracks", fj.ntracks)
                self.out.fillBranch(prefix + "ntracks_sv12", fj.ntracks_sv12)

                assert(len(fj.subjets) == 2)
                for idx_sj, sj in enumerate(fj.subjets):
                    prefix_sj = prefix + 'sj%d_' % (idx_sj + 1)
                    self.out.fillBranch(prefix_sj + "ntracks", sum([sv.ntracks for sv in sj.sv_list]))
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
                self.out.fillBranch(prefix + "sj12_masscor_dxysig", fj.sj12_masscor_dxysig)

                # sfBDT
                self.out.fillBranch(prefix + "sfBDT", fj.sfBDT)
