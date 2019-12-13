import os
import numpy as np
import ROOT, math
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.tools import deltaPhi, deltaR, closest

from PhysicsTools.NanoHRTTools.helpers.jetmetCorrector import JetMETCorrector

import logging
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] %(levelname)s: %(message)s')


class METObject(Object):
    def p4(self):
        ret = ROOT.TLorentzVector()
        ret.SetPtEtaPhiM(self.pt, 0, self.phi, 0)
        return ret

    
def correctedsvmass(self):
    svp4       = self.p4()
    svm2       = (self.mass) * (self.mass)
    svp2       = (svp4.P()) * (svp4.P())
    sin2theta   = (math.sin(self.pAngle)) * (math.sin(self.pAngle)) 
    svmasscorr = math.sqrt( svm2 + (svp2 * sin2theta) ) + (svp4.P())*(math.sqrt(sin2theta));
    return svmasscorr;

    
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


class HFQCDSFTreeProducer(Module, object):

    def __init__(self, year, **kwargs):
	self.year = year
        channel = 'qcd'
        
        print ('Channel: %s' %(channel))   
	print ('Running on %d DATA/MC' %(year))
  	
        #self._channel = channel  # '0L', '1L', '2L'
        self._systOpt = {'jec':False, 'jes':None, 'jes_source':'', 'jer':'nominal', 'met_unclustered':None}
        for k in kwargs:
            self._systOpt[k] = kwargs[k]
            
        #logging.info('Running %s channel with systematics %s', self._channel, str(self._systOpt))
        logging.info('Running systematics %s', str(self._systOpt))
        
        self.jetmetCorr = JetMETCorrector(year,
					  jetType="AK4PFchs",
                                          jec=self._systOpt['jec'],
                                          jes=self._systOpt['jes'],
                                          jes_source=self._systOpt['jes_source'],
                                          jer=self._systOpt['jer'],
                                          met_unclustered=self._systOpt['met_unclustered'])

        self.subjetCorr = JetMETCorrector(year,
					  jetType="AK4PFPuppi",
                                          jec=self._systOpt['jec'],
                                          jes=self._systOpt['jes'],
                                          jes_source=self._systOpt['jes_source'],
                                          jer=self._systOpt['jer'],
                                          met_unclustered=self._systOpt['met_unclustered'])

        
    def beginJob(self):
        self.jetmetCorr.beginJob()
        self.subjetCorr.beginJob()

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.isMC = bool(inputTree.GetBranch('genWeight'))

        self.out = wrappedOutputTree

        ## trigger variables
        self.out.branch("passHT800", "O")  ## 2016
        self.out.branch("passHT780", "O")  ## 2017/2018
        self.out.branch("passHT890", "O")  ## 2017/2018
        self.out.branch("passHT1050", "O") ## 2017/2018
        self.out.branch("passHTTrig", "O")

        ## event variables
        self.out.branch("ht", "F")
        self.out.branch("passmetfilters", "O")

        ## Large-R jets
        self.out.branch("n_fatjet", "I")
        self.out.branch("fj_1_nbhadrons", "F")
        self.out.branch("fj_1_nchadrons", "F")
        self.out.branch("fj_1_pt", "F")
        self.out.branch("fj_1_eta", "F")
        self.out.branch("fj_1_phi", "F")
        self.out.branch("fj_1_sdmass", "F")
        self.out.branch("fj_1_sjmass", "F")
        self.out.branch("fj_1_score_cc", "F")
        self.out.branch("fj_1_score_bb", "F")
        self.out.branch("fj_1_tau21", "F")
        self.out.branch("fj_1_btagcsvv2", "F")
        self.out.branch("fj_1_btagjp", "F")
        self.out.branch("fj_1_nsv", "I")
        self.out.branch("fj_1_sj1_nbhadrons", "F")
        self.out.branch("fj_1_sj1_nchadrons", "F")
        self.out.branch("fj_1_sj1_pt", "F")
        self.out.branch("fj_1_sj1_eta", "F")
        self.out.branch("fj_1_sj1_phi", "F")
        self.out.branch("fj_1_sj1_btagcsvv2", "F")
        self.out.branch("fj_1_sj1_btagjp", "F")
        self.out.branch("fj_1_sj1_nsv", "I")
        self.out.branch("fj_1_sj1_sv1_pt", "F")
        self.out.branch("fj_1_sj1_sv1_mass", "F")
        self.out.branch("fj_1_sj1_sv1_masscor", "F")
        self.out.branch("fj_1_sj1_sv1_ntracks", "I")
        self.out.branch("fj_1_sj1_sv1_dxy", "F")
        self.out.branch("fj_1_sj1_sv1_dxysig", "F")
        self.out.branch("fj_1_sj1_sv1_dlen", "F")
        self.out.branch("fj_1_sj1_sv1_dlensig", "F")
        self.out.branch("fj_1_sj1_sv1_chi2ndof", "F")
        self.out.branch("fj_1_sj1_sv1_pangle", "F")
        self.out.branch("fj_1_sj2_nbhadrons", "F")
        self.out.branch("fj_1_sj2_nchadrons", "F")
        self.out.branch("fj_1_sj2_pt", "F")
        self.out.branch("fj_1_sj2_eta", "F")
        self.out.branch("fj_1_sj2_phi", "F")
        self.out.branch("fj_1_sj2_btagcsvv2", "F")
        self.out.branch("fj_1_sj2_btagjp", "F")
        self.out.branch("fj_1_sj2_nsv", "I")
        self.out.branch("fj_1_sj2_sv1_pt", "F")
        self.out.branch("fj_1_sj2_sv1_mass", "F")
        self.out.branch("fj_1_sj2_sv1_masscor", "F")
        self.out.branch("fj_1_sj2_sv1_ntracks", "I")
        self.out.branch("fj_1_sj2_sv1_dxy", "F")
        self.out.branch("fj_1_sj2_sv1_dxysig", "F")
        self.out.branch("fj_1_sj2_sv1_dlen", "F")
        self.out.branch("fj_1_sj2_sv1_dlensig", "F")
        self.out.branch("fj_1_sj2_sv1_chi2ndof", "F")
        self.out.branch("fj_1_sj2_sv1_pangle", "F")
        self.out.branch("fj_1_sj12_masscor_dxysig", "F")

        self.out.branch("fj_2_nbhadrons", "F")
        self.out.branch("fj_2_nchadrons", "F")
        self.out.branch("fj_2_pt", "F")
        self.out.branch("fj_2_eta", "F")
        self.out.branch("fj_2_phi", "F")
        self.out.branch("fj_2_sdmass", "F")
        self.out.branch("fj_2_sjmass", "F")
        self.out.branch("fj_2_score_cc", "F")
        self.out.branch("fj_2_score_bb", "F")
        self.out.branch("fj_2_tau21", "F")
        self.out.branch("fj_2_btagcsvv2", "F")
        self.out.branch("fj_2_btagjp", "F")
        self.out.branch("fj_2_nsv", "I")
        self.out.branch("fj_2_sj1_nbhadrons", "F")
        self.out.branch("fj_2_sj1_nchadrons", "F")
        self.out.branch("fj_2_sj1_pt", "F")
        self.out.branch("fj_2_sj1_eta", "F")
        self.out.branch("fj_2_sj1_phi", "F")
        self.out.branch("fj_2_sj1_btagcsvv2", "F")
        self.out.branch("fj_2_sj1_btagjp", "F")
        self.out.branch("fj_2_sj1_nsv", "I")
        self.out.branch("fj_2_sj1_sv1_pt", "F")
        self.out.branch("fj_2_sj1_sv1_mass", "F")
        self.out.branch("fj_2_sj1_sv1_masscor", "F")
        self.out.branch("fj_2_sj1_sv1_ntracks", "I")
        self.out.branch("fj_2_sj1_sv1_dxy", "F")
        self.out.branch("fj_2_sj1_sv1_dxysig", "F")
        self.out.branch("fj_2_sj1_sv1_dlen", "F")
        self.out.branch("fj_2_sj1_sv1_dlensig", "F")
        self.out.branch("fj_2_sj1_sv1_chi2ndof", "F")
        self.out.branch("fj_2_sj1_sv1_pangle", "F")
        self.out.branch("fj_2_sj2_nbhadrons", "F")
        self.out.branch("fj_2_sj2_nchadrons", "F")
        self.out.branch("fj_2_sj2_pt", "F")
        self.out.branch("fj_2_sj2_eta", "F")
        self.out.branch("fj_2_sj2_phi", "F")
        self.out.branch("fj_2_sj2_btagcsvv2", "F")
        self.out.branch("fj_2_sj2_btagjp", "F")
        self.out.branch("fj_2_sj2_nsv", "I")
        self.out.branch("fj_2_sj2_sv1_pt", "F")
        self.out.branch("fj_2_sj2_sv1_mass", "F")
        self.out.branch("fj_2_sj2_sv1_masscor", "F")
        self.out.branch("fj_2_sj2_sv1_ntracks", "I")
        self.out.branch("fj_2_sj2_sv1_dxy", "F")
        self.out.branch("fj_2_sj2_sv1_dxysig", "F")
        self.out.branch("fj_2_sj2_sv1_dlen", "F")
        self.out.branch("fj_2_sj2_sv1_dlensig", "F")
        self.out.branch("fj_2_sj2_sv1_chi2ndof", "F")
        self.out.branch("fj_2_sj2_sv1_pangle", "F")
        self.out.branch("fj_2_sj12_masscor_dxysig", "F")

    def _correctJetAndMET(self, event):
        event._allJets = Collection(event, "Jet")
        event.ak8Subjets = Collection(event, "SubJet")  # do not sort after updating!!
        event.met = METObject(event, "MET")

        if self.isMC or self._systOpt['jec']:
            rho = event.fixedGridRhoFastjetAll
            ## correct AK4 jets and MET
            self.jetmetCorr.correctJetAndMET(jets=event._allJets, met=event.met, rho=rho,
                                             genjets=Collection(event, 'GenJet') if self.isMC else None,
                                             isMC=self.isMC, runNumber=event.run)
            event._allJets = sorted(event._allJets, key=lambda x : x.pt, reverse=True)  # sort by pt after updating

            ## correct AK15 subjets
            #self.subjetCorr.correctJetAndMET(jets=event.ak15Subjets, met=None, rho=rho,
            #                                 genjets=Collection(event, 'GenSubJetAK15') if self.isMC else None,
            #                                 isMC=self.isMC, runNumber=event.run)

        # construct AK8 p4 from (updated) subjets
        event._allAK8jets = Collection(event, "FatJet")
        for fj in event._allAK8jets:
            fj.subjets = get_subjets(fj, event.ak8Subjets, ('subJetIdx1', 'subJetIdx2'))
            newP4 = ROOT.TLorentzVector()
            if len(fj.subjets) == 2:
                newP4 = fj.subjets[0].p4() + fj.subjets[1].p4()
            fj.pt, fj.eta, fj.phi, fj.mass, fj.msoftdrop = newP4.Pt(), newP4.Eta(), newP4.Phi(), newP4.M(), newP4.M()
        event._allAK8jets = sorted(event._allAK8jets, key=lambda x : x.pt, reverse=True)  # sort by pt

    def _selectLeptons(self, event):
        # do lepton selection
        event.looseLeptons = []
        electrons = Collection(event, "Electron")
        for el in electrons:
            el.pt /= el.eCorr  # FIXME: turn off EGM correction
            el.etaSC = el.eta + el.deltaEtaSC
            if el.pt > 7 and abs(el.eta) < 2.4 and abs(el.dxy) < 0.05 and abs(el.dz) < 0.2 and el.pfRelIso03_all < 0.4:
                if el.mvaFall17V2noIso_WP90:
                    event.looseLeptons.append(el)

        muons = Collection(event, "Muon")
        for mu in muons:
            if mu.pt > 5 and abs(mu.eta) < 2.4 and abs(mu.dxy) < 0.5 and abs(mu.dz) < 1.0 and mu.pfRelIso04_all < 0.4:
                event.looseLeptons.append(mu)

        event.looseLeptons.sort(key=lambda x : x.pt, reverse=True)


        
    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

        
    def _prepareEvent(self, event):

        logging.debug('processing event %d' % event.event)

        ## veto leptons
        if len(event.looseLeptons) != 0:
            return False
        
        ## selection on fatjets           
	event.ak8jets = []
	for fj in event._allAK8jets:
	    if not (fj.pt > 200 and abs(fj.eta) < 2.4 and (fj.jetId & 2)):
		continue
	    event.ak8jets.append(fj)

	if (len(event.ak8jets)<2):
	    return False
	if (event.ak8jets[1].pt < 200):
	    return False
	if ( (event.ak8jets[1].msoftdrop < 50) or (event.ak8jets[1].msoftdrop > 200)):
	    return False

        ## selection on SV
        event._sv = Collection(event, "SV")        
        #event._sv = sorted(event._sv, key=lambda x : x.dxySig, reverse=True)  ## sort by dxysig
        if len(event._sv) < 2:
            return False

        event.sv = []
        for isv in event._sv:
            #if (isv.ntracks>2 and abs(isv.dxy)<3. and isv.dlenSig>4):
            #if (isv.dlenSig>4):
                event.sv.append(isv)

        drsj12_ = 0.5*(deltaR(event.ak8jets[1].subjets[0],event.ak8jets[1].subjets[1]))
        drcut_  = min(0.4,drsj12_);
        nsvsj1_ = 0 
        nsvsj2_ = 0

        ## sv mathced to subjets
        for isv in event.sv:
            if deltaR(isv,event.ak8jets[1].subjets[0]) < drcut_:
                nsvsj1_ += 1
            elif deltaR(isv,event.ak8jets[1].subjets[1]) < drcut_:
                nsvsj2_ += 1
        if nsvsj1_<1 or nsvsj2_<1:
            return False


        ## ht selection
        event.ak4jets = []
        for j in event._allJets:
            if not (j.pt > 25 and abs(j.eta) < 2.4 and (j.jetId & 2)):
                continue
            event.ak4jets.append(j)

        event.ht = 0
        for j in event.ak4jets:
            event.ht += j.pt
        if (event.ht<800.):
            return False
   
        ## return True if passes selection
        return True


    def _fillEventInfo(self, event):

        ## event cleaning
        self.out.fillBranch("passmetfilters", bool(
                event.Flag_goodVertices and
                event.Flag_globalSuperTightHalo2016Filter and
                event.Flag_HBHENoiseFilter and
                event.Flag_HBHENoiseIsoFilter and
                event.Flag_EcalDeadCellTriggerPrimitiveFilter and
                event.Flag_BadPFMuonFilter and
                event.Flag_BadChargedCandidateFilter)
                            )

        ## Triggers
        pass800_  = False
        pass780_  = False
        pass890_  = False
        pass1050_ = False
        passHTTrig_ = False
        if not self.isMC:
            try: 
                if event.HLT_PFHT800:
                    pass800_ = True
            except:
                pass800_ = False
            try: 
                if event.HLT_PFHT780:
                    pass780_ = True
            except:
                pass780_ = False
            try: 
                if event.HLT_PFHT890:
                    pass890_ = True
            except:
                pass890_ = False
            try: 
                if event.HLT_PFHT1050:
                    pass1050_ = True
            except:
                pass1050_ = False

        if pass780_ or pass800_ or pass890_ or pass1050_:
            passHTTrig_ = True
        
          
        self.out.fillBranch("passHT780" , pass780_)
        self.out.fillBranch("passHT800" , pass800_)
        self.out.fillBranch("passHT890" , pass890_)
        self.out.fillBranch("passHT1050", pass1050_)
        self.out.fillBranch("passHTTrig", passHTTrig_)

        ## event variables
        self.out.fillBranch("ht", event.ht)

        
        ## fj variables
        self.out.fillBranch("n_fatjet", len(event.ak8jets))

        if len(event.ak8jets)>0:
            self.out.fillBranch("fj_1_nbhadrons", event.ak8jets[0].nBHadrons)
            self.out.fillBranch("fj_1_nchadrons", event.ak8jets[0].nCHadrons)
            self.out.fillBranch("fj_1_pt", event.ak8jets[0].pt)
            self.out.fillBranch("fj_1_eta", event.ak8jets[0].eta)
            self.out.fillBranch("fj_1_phi", event.ak8jets[0].phi)
	    self.out.fillBranch("fj_1_sdmass", event.ak8jets[0].msoftdrop)
            self.out.fillBranch("fj_1_sjmass", get_sdmass(event.ak8jets[0].subjets))
	    self.out.fillBranch("fj_1_score_cc", event.ak8jets[0].deepTagMD_ccvsLight)
	    self.out.fillBranch("fj_1_score_bb", event.ak8jets[0].deepTagMD_bbvsLight)
    	    self.out.fillBranch("fj_1_tau21", (event.ak8jets[0].tau2/event.ak8jets[0].tau1) if event.ak8jets[0].tau1>0 else -1.)
            self.out.fillBranch("fj_1_btagcsvv2", event.ak8jets[0].btagCSVV2)
            #self.out.fillBranch("fj_1_btagjp", event.ak8jets[0].btagJP)

            ## go to sj
	    self.out.fillBranch("fj_1_sj1_pt", event.ak8jets[0].subjets[0].pt)
	    self.out.fillBranch("fj_1_sj1_eta", event.ak8jets[0].subjets[0].eta)
	    self.out.fillBranch("fj_1_sj1_phi", event.ak8jets[0].subjets[0].phi)
            self.out.fillBranch("fj_1_sj1_nbhadrons", event.ak8jets[0].subjets[0].nBHadrons)
	    self.out.fillBranch("fj_1_sj1_nchadrons", event.ak8jets[0].subjets[0].nCHadrons)
	    self.out.fillBranch("fj_1_sj1_btagcsvv2", event.ak8jets[0].subjets[0].btagCSVV2)
            #Self.out.fillBranch("fj_1_sj1_btagjp", event.ak8jets[0].subjets[0].btagJP)

	    self.out.fillBranch("fj_1_sj2_pt", event.ak8jets[0].subjets[1].pt)
	    self.out.fillBranch("fj_1_sj2_eta", event.ak8jets[0].subjets[1].eta)
	    self.out.fillBranch("fj_1_sj2_phi", event.ak8jets[0].subjets[1].phi)
            self.out.fillBranch("fj_1_sj2_nbhadrons", event.ak8jets[0].subjets[1].nBHadrons)
	    self.out.fillBranch("fj_1_sj2_nchadrons", event.ak8jets[0].subjets[1].nCHadrons)
            self.out.fillBranch("fj_1_sj2_btagcsvv2", event.ak8jets[0].subjets[1].btagCSVV2)
            #self.out.fillBranch("fj_1_sj2_btagjp", event.ak8jets[0].subjets[1].btagJP)

            ## sv vars related to each sj [for the leading in pt fatjet]
            sj_drsj12_ = 0.5*(deltaR(event.ak8jets[0].subjets[0],event.ak8jets[0].subjets[1]))
            sj_drcut_  = min(0.4,sj_drsj12_);
            sj1_nsv_ , sj2_nsv_ = 0 , 0
            sj1_masscor_ , sj1_dxysig_ = 0. , 0.
            sj2_masscor_ , sj2_dxysig_ = 0. , 0.
            fj_1_sj12_masscor_dxysig_ = -9.
            
            for isv in event.sv:
            
                if deltaR(isv,event.ak8jets[0].subjets[0]) < sj_drcut_:
                    sj1_nsv_ += 1
                    if sj1_nsv_ == 1:
                        self.out.fillBranch("fj_1_sj1_sv1_pt", isv.pt)
                        self.out.fillBranch("fj_1_sj1_sv1_mass", isv.mass)
                        self.out.fillBranch("fj_1_sj1_sv1_masscor", correctedsvmass(isv))
                        self.out.fillBranch("fj_1_sj1_sv1_ntracks", isv.ntracks)
                        self.out.fillBranch("fj_1_sj1_sv1_dxy", isv.dxy)
                        self.out.fillBranch("fj_1_sj1_sv1_dxysig", isv.dxySig)
                        self.out.fillBranch("fj_1_sj1_sv1_dlen", isv.dlen)
                        self.out.fillBranch("fj_1_sj1_sv1_dlensig", isv.dlenSig)
                        self.out.fillBranch("fj_1_sj1_sv1_chi2ndof", (isv.chi2/isv.ndof) if isv.ndof>0 else -1.)
                        self.out.fillBranch("fj_1_sj1_sv1_pangle", isv.pAngle)
                        sj1_masscor_ = correctedsvmass(isv)
                        sj1_dxysig_  = isv.dxySig 
                elif deltaR(isv,event.ak8jets[0].subjets[1]) < sj_drcut_:
                    sj2_nsv_ += 1
                    if sj2_nsv_ == 1:
                        self.out.fillBranch("fj_1_sj2_sv1_pt", isv.pt)
                        self.out.fillBranch("fj_1_sj2_sv1_mass", isv.mass)
                        self.out.fillBranch("fj_1_sj2_sv1_masscor", correctedsvmass(isv))
                        self.out.fillBranch("fj_1_sj2_sv1_ntracks", isv.ntracks)
                        self.out.fillBranch("fj_1_sj2_sv1_dxy", isv.dxy)
                        self.out.fillBranch("fj_1_sj2_sv1_dxysig", isv.dxySig)
                        self.out.fillBranch("fj_1_sj2_sv1_dlen", isv.dlen)
                        self.out.fillBranch("fj_1_sj2_sv1_dlensig", isv.dlenSig)
                        self.out.fillBranch("fj_1_sj2_sv1_chi2ndof", (isv.chi2/isv.ndof) if isv.ndof>0 else -1.)
                        self.out.fillBranch("fj_1_sj2_sv1_pangle", isv.pAngle)
                        sj2_masscor_ = correctedsvmass(isv)
                        sj2_dxysig_  = isv.dxySig 
            if (sj1_nsv_>0 and sj2_nsv_>0):
                if (sj1_dxysig_ > sj2_dxysig_):
                    fj_1_sj12_masscor_dxysig_ = sj1_masscor_
                else:
                    fj_1_sj12_masscor_dxysig_ = sj2_masscor_
	    self.out.fillBranch("fj_1_sj1_nsv", sj1_nsv_)
	    self.out.fillBranch("fj_1_sj2_nsv", sj2_nsv_)	    
            self.out.fillBranch("fj_1_sj12_masscor_dxysig", fj_1_sj12_masscor_dxysig_)



        if len(event.ak8jets)>1:
            
            self.out.fillBranch("fj_2_nbhadrons", event.ak8jets[1].nBHadrons)
            self.out.fillBranch("fj_2_nchadrons", event.ak8jets[1].nCHadrons)
            self.out.fillBranch("fj_2_pt", event.ak8jets[1].pt)
            self.out.fillBranch("fj_2_eta", event.ak8jets[1].eta)
            self.out.fillBranch("fj_2_phi", event.ak8jets[1].phi)
	    self.out.fillBranch("fj_2_sdmass", event.ak8jets[1].msoftdrop)
            self.out.fillBranch("fj_2_sjmass", get_sdmass(event.ak8jets[1].subjets))
            self.out.fillBranch("fj_2_score_cc", event.ak8jets[1].deepTagMD_ccvsLight)
            self.out.fillBranch("fj_2_score_bb", event.ak8jets[1].deepTagMD_bbvsLight)
    	    self.out.fillBranch("fj_2_tau21", (event.ak8jets[1].tau2/event.ak8jets[1].tau1) if event.ak8jets[1].tau1>0 else -1.)
            self.out.fillBranch("fj_2_btagcsvv2", event.ak8jets[1].btagCSVV2)
            #self.out.fillBranch("fj_2_btagjp", event.ak8jets[1].btagJP)

            ## go to sj
	    self.out.fillBranch("fj_2_sj1_pt", event.ak8jets[1].subjets[0].pt)
	    self.out.fillBranch("fj_2_sj1_eta", event.ak8jets[1].subjets[0].eta)
	    self.out.fillBranch("fj_2_sj1_phi", event.ak8jets[1].subjets[0].phi)
            self.out.fillBranch("fj_2_sj1_nbhadrons", event.ak8jets[1].subjets[0].nBHadrons)
	    self.out.fillBranch("fj_2_sj1_nchadrons", event.ak8jets[1].subjets[0].nCHadrons)
	    self.out.fillBranch("fj_2_sj1_btagcsvv2", event.ak8jets[1].subjets[0].btagCSVV2)
            #self.out.fillBranch("fj_2_sj1_btagjp", event.ak8jets[1].subjets[0].btagJP)

	    self.out.fillBranch("fj_2_sj2_pt", event.ak8jets[1].subjets[1].pt)
	    self.out.fillBranch("fj_2_sj2_eta", event.ak8jets[1].subjets[1].eta)
	    self.out.fillBranch("fj_2_sj2_phi", event.ak8jets[1].subjets[1].phi)
            self.out.fillBranch("fj_2_sj2_nbhadrons", event.ak8jets[1].subjets[1].nBHadrons)
	    self.out.fillBranch("fj_2_sj2_nchadrons", event.ak8jets[1].subjets[1].nCHadrons)
            self.out.fillBranch("fj_2_sj2_btagcsvv2", event.ak8jets[1].subjets[1].btagCSVV2)
            #self.out.fillBranch("fj_2_sj2_btagjp", event.ak8jets[1].subjets[1].btagJP)

            ## sv vars related to each sj
            sj_drsj12_ = 0.5*(deltaR(event.ak8jets[1].subjets[0],event.ak8jets[1].subjets[1]))
            sj_drcut_  = min(0.4,sj_drsj12_);
            sj1_nsv_ , sj2_nsv_ = 0 , 0
            sj1_masscor_ , sj1_dxysig_ = 0. , 0.
            sj2_masscor_ , sj2_dxysig_ = 0. , 0.
            fj_2_sj12_masscor_dxysig_ = -9.
            
            for isv in event.sv:

                if deltaR(isv,event.ak8jets[1].subjets[0]) < sj_drcut_:
                    sj1_nsv_ += 1
                    if sj1_nsv_ == 1:
                        self.out.fillBranch("fj_2_sj1_sv1_pt", isv.pt)
                        self.out.fillBranch("fj_2_sj1_sv1_mass", isv.mass)
                        self.out.fillBranch("fj_2_sj1_sv1_masscor", correctedsvmass(isv))
                        self.out.fillBranch("fj_2_sj1_sv1_ntracks", isv.ntracks)
                        self.out.fillBranch("fj_2_sj1_sv1_dxy", isv.dxy)
                        self.out.fillBranch("fj_2_sj1_sv1_dxysig", isv.dxySig)
                        self.out.fillBranch("fj_2_sj1_sv1_dlen", isv.dlen)
                        self.out.fillBranch("fj_2_sj1_sv1_dlensig", isv.dlenSig)
                        self.out.fillBranch("fj_2_sj1_sv1_chi2ndof", (isv.chi2/isv.ndof) if isv.ndof>0 else -1.)
                        self.out.fillBranch("fj_2_sj1_sv1_pangle", isv.pAngle)
                        sj1_masscor_ = correctedsvmass(isv)
                        sj1_dxysig_  = isv.dxySig
                elif deltaR(isv,event.ak8jets[1].subjets[1]) < sj_drcut_:
                    sj2_nsv_ += 1
                    if sj2_nsv_ == 1:
                        self.out.fillBranch("fj_2_sj2_sv1_pt", isv.pt)
                        self.out.fillBranch("fj_2_sj2_sv1_mass", isv.mass)
                        self.out.fillBranch("fj_2_sj2_sv1_masscor", correctedsvmass(isv))
                        self.out.fillBranch("fj_2_sj2_sv1_ntracks", isv.ntracks)
                        self.out.fillBranch("fj_2_sj2_sv1_dxy", isv.dxy)
                        self.out.fillBranch("fj_2_sj2_sv1_dxysig", isv.dxySig)
                        self.out.fillBranch("fj_2_sj2_sv1_dlen", isv.dlen)
                        self.out.fillBranch("fj_2_sj2_sv1_dlensig", isv.dlenSig)
                        self.out.fillBranch("fj_2_sj2_sv1_chi2ndof", (isv.chi2/isv.ndof) if isv.ndof>0 else -1.)
                        self.out.fillBranch("fj_2_sj2_sv1_pangle", isv.pAngle)
                        sj2_masscor_ = correctedsvmass(isv)
                        sj2_dxysig_  = isv.dxySig

            if (sj1_nsv_>0 and sj2_nsv_>0):
                if (sj1_dxysig_ > sj2_dxysig_):
                    fj_2_sj12_masscor_dxysig_ = sj1_masscor_
                else:
                    fj_2_sj12_masscor_dxysig_ = sj2_masscor_
	    self.out.fillBranch("fj_2_sj1_nsv", sj1_nsv_)
	    self.out.fillBranch("fj_2_sj2_nsv", sj2_nsv_)
	    self.out.fillBranch("fj_2_sj12_masscor_dxysig", fj_2_sj12_masscor_dxysig_)

        
    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""

        self._correctJetAndMET(event)
        self._selectLeptons(event)

        if self._prepareEvent(event) is False:
            return False

        # fill
        self._fillEventInfo(event)
        
        return True


# define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed
hfQCDTree = lambda : HFQCDSFTreeProducer()

