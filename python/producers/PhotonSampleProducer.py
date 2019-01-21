#!/usr/bin/env python
import os, sys
import ROOT, math
ROOT.PyConfig.IgnoreCommandLineOptions = True
from importlib import import_module
from PhysicsTools.NanoHRTTools.helpers.deepAK8Helper import get_nominal_score, get_decorr_score
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.tools import deltaPhi, deltaR, closest

import logging
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] %(levelname)s: %(message)s')


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

class PhotonSampleProducer(Module):

    def __init__(self, **kwargs):
        pass
        
    def beginJob(self):
        pass

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.isMC = bool(inputTree.GetBranch('genWeight'))

        self.out = wrappedOutputTree

        ## trigger variables
        self.out.branch("passPhoton165_HE10", "O")

        ## event variables
        self.out.branch("ht", "F")

        ## photons
        self.out.branch("nphotons", "I")
        self.out.branch("pho_1_pt", "F")
        self.out.branch("pho_1_eta", "F")

        ## Large-R jets
        self.out.branch("n_ak8" , "I")
        self.out.branch("n_ca15", "I")

        self.out.branch("ak8_1_pt" , "F")
        self.out.branch("ak8_1_eta" , "F")
        self.out.branch("ak8_1_phi" , "F")
        self.out.branch("ak8_1_mass" , "F")
        self.out.branch("ak8_1_n2b1" , "F")
        self.out.branch("ak8_1_n3b1" , "F")
        self.out.branch("ak8_1_tau1" , "F")
        self.out.branch("ak8_1_tau2" , "F")
        self.out.branch("ak8_1_tau3" , "F")
        self.out.branch("ak8_1_DeepAK8_WvsQCD" , "F")
        self.out.branch("ak8_1_DeepAK8_ZvsQCD" , "F")
        self.out.branch("ak8_1_DeepAK8_HvsQCD" , "F")
        self.out.branch("ak8_1_DeepAK8_TvsQCD" , "F")
        self.out.branch("ak8_1_DeepAK8MD_WvsQCD" , "F")
        self.out.branch("ak8_1_DeepAK8MD_ZvsQCD" , "F")
        self.out.branch("ak8_1_DeepAK8MD_HvsQCD" , "F")
        self.out.branch("ak8_1_DeepAK8MD_TvsQCD" , "F")
        self.out.branch("ak8_1_best_WvsQCD" , "F")
        self.out.branch("ak8_1_best_ZvsQCD" , "F")
        self.out.branch("ak8_1_best_HvsQCD" , "F")
        self.out.branch("ak8_1_best_TvsQCD" , "F")
        self.out.branch("ak8_1_btagCSVV2" , "F")
        self.out.branch("ak8_1_btagHbb" , "F")
        self.out.branch("ak8_1_nnHbb" , "F")
        self.out.branch("ak8_1_nnHcc" , "F")
        self.out.branch("ak8_1_sj1_btagCSVV2" , "F")
        self.out.branch("ak8_1_sj2_btagCSVV2" , "F")

        self.out.branch("ca15_1_pt"           , "F")
        self.out.branch("ca15_1_eta"          , "F")
        self.out.branch("ca15_1_phi"          , "F")
        self.out.branch("ca15_1_mass"         , "F")
        self.out.branch("ca15_1_ecf0"         , "F")
        self.out.branch("ca15_1_ecfTopTagBDT" , "F")
        self.out.branch("ca15_1_httFRec"      , "F")
        self.out.branch("ca15_1_tau32sd"      , "F")


    def _correctJetAndMET(self, event):
        event._allJets = Collection(event, "Jet")
        event.ak8Subjets = Collection(event, "CustomAK8PuppiSubJet")  # do not sort after updating!!
        event.ca15Subjets = Collection(event, "CA15PuppiSubJet")  # do not sort after updating!!
        
        if self.isMC or self._systOpt['jec']:
            rho = event.fixedGridRhoFastjetAll

        ## construct AK8 p4 from (updated) subjets
        event._allAK8jets = Collection(event, "CustomAK8Puppi")
        for fj in event._allAK8jets:
            fj.subjets = get_subjets(fj, event.ak8Subjets, ('subJetIdx1', 'subJetIdx2'))
            newP4 = ROOT.TLorentzVector()
            if len(fj.subjets) == 2:
                newP4 = fj.subjets[0].p4() + fj.subjets[1].p4()
            fj.pt, fj.eta, fj.phi, fj.mass, fj.msoftdrop = newP4.Pt(), newP4.Eta(), newP4.Phi(), newP4.M(), newP4.M()
        event._allAK8jets = sorted(event._allAK8jets, key=lambda x : x.pt, reverse=True)  # sort by pt

        ## construct CA15 p4 from (updated) subjets
        event._allCA15jets = Collection(event, "CA15Puppi")
        for fj in event._allCA15jets:
            fj.subjets = get_subjets(fj, event.ca15Subjets, ('subJetIdx1', 'subJetIdx2'))
            newP4 = ROOT.TLorentzVector()
            if len(fj.subjets) == 2:
                newP4 = fj.subjets[0].p4() + fj.subjets[1].p4()
            fj.pt, fj.eta, fj.phi, fj.mass, fj.msoftdrop = newP4.Pt(), newP4.Eta(), newP4.Phi(), newP4.M(), newP4.M()
        event._allCA15jets = sorted(event._allCA15jets, key=lambda x : x.pt, reverse=True)  # sort by pt


        
    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

        
    def _prepareEvent(self, event):

        logging.debug('processing event %d' % event.event)

        ## select leading photon
        event._allPhotons = Collection(event, "Photon")
        event.photons = []
        for ipho in event._allPhotons:
            if not (ipho.pt > 200 and abs(ipho.eta) < 2.4 and ipho.mvaID_WP90):
                continue
            event.photons.append(ipho)

        if (len(event.photons)<1):
            return False

        
        ## selection on AK8 jets / drop if overlaps with a photon
        event.ak8jets = []
        for fj in event._allAK8jets:
            if not (fj.pt > 200 and abs(fj.eta) < 2.4 and (fj.jetId & 2)):
                continue
            if deltaR(event.photons[0],fj)<0.8:
                continue
            event.ak8jets.append(fj)

        if (len(event.ak8jets)<1):
            return False

        ## selection on CA15 jets / drop if overlaps with a photon
        event.ca15jets = []
        for fj in event._allCA15jets:
            if not (fj.pt > 200 and abs(fj.eta) < 2.4 and (fj.jetId & 2)):
                continue
            if deltaR(event.photons[0],fj)<1.5:
                continue
            event.ca15jets.append(fj)
   
        if (len(event.ca15jets)<1):
            return False


        ## require the leading ak8 & ca15 jets overlap
        if deltaR(event.ak8jets[0],event.ca15jets[0])>0.8:
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
        if (event.ht<200.):
            return False
    
        ## return True if passes selection
        return True


    def _fillEventInfo(self, event):

        ## Triggers
        passPhoton165_HE10_ = False

        try: 
            if event.HLT_Photon165_HE10:
                passPhoton165_HE10_ = True
        except:
            passPhoton165_HE10_ = False

        self.out.fillBranch("passPhoton165_HE10", passPhoton165_HE10_)


        ## event variables
        self.out.fillBranch("ht", event.ht)

        ## photon variables
        self.out.fillBranch("nphotons" , len(event.photons))
        self.out.fillBranch("pho_1_pt" , event.photons[0].pt ) 
        self.out.fillBranch("pho_1_eta", event.photons[0].eta)        

        ## Large-R jets
        self.out.fillBranch("n_ak8", len(event.ak8jets))
        self.out.fillBranch("n_ca15", len(event.ca15jets))

        self.out.fillBranch("ak8_1_pt" , event.ak8jets[0].pt)
        self.out.fillBranch("ak8_1_eta" , event.ak8jets[0].eta)
        self.out.fillBranch("ak8_1_phi" , event.ak8jets[0].phi)
        self.out.fillBranch("ak8_1_mass" , event.ak8jets[0].msoftdrop)
        self.out.fillBranch("ak8_1_n2b1" , event.ak8jets[0].n2b1)
        self.out.fillBranch("ak8_1_n3b1" , event.ak8jets[0].n3b1)
        self.out.fillBranch("ak8_1_tau1" , event.ak8jets[0].tau1)
        self.out.fillBranch("ak8_1_tau2" , event.ak8jets[0].tau2)
        self.out.fillBranch("ak8_1_tau3" , event.ak8jets[0].tau3)
        self.out.fillBranch("ak8_1_DeepAK8_WvsQCD"    , get_nominal_score(event.ak8jets[0], 'WvsQCD'))
        self.out.fillBranch("ak8_1_DeepAK8_ZvsQCD"    , get_nominal_score(event.ak8jets[0], 'ZvsQCD'))
        self.out.fillBranch("ak8_1_DeepAK8_HvsQCD"    , get_nominal_score(event.ak8jets[0], 'HbbvsQCD'))
        self.out.fillBranch("ak8_1_DeepAK8_TvsQCD"    , get_nominal_score(event.ak8jets[0], 'TvsQCD'))
        self.out.fillBranch("ak8_1_DeepAK8MD_WvsQCD"  , get_decorr_score(event.ak8jets[0], 'WvsQCD'))
        self.out.fillBranch("ak8_1_DeepAK8MD_ZvsQCD"  , get_decorr_score(event.ak8jets[0], 'ZHbbvsQCD'))
        self.out.fillBranch("ak8_1_DeepAK8MD_HvsQCD"  , get_decorr_score(event.ak8jets[0], 'ZHbbvsQCD'))
        self.out.fillBranch("ak8_1_DeepAK8MD_TvsQCD"  , get_decorr_score(event.ak8jets[0], 'TvsQCD'))
        self.out.fillBranch("ak8_1_best_WvsQCD"       , event.ak8jets[0].bestW/(event.ak8jets[0].bestW + event.ak8jets[0].bestQCD + event.ak8jets[0].bestB))
        self.out.fillBranch("ak8_1_best_ZvsQCD"       , event.ak8jets[0].bestZ/(event.ak8jets[0].bestZ + event.ak8jets[0].bestQCD + event.ak8jets[0].bestB))
        self.out.fillBranch("ak8_1_best_HvsQCD"       , event.ak8jets[0].bestH/(event.ak8jets[0].bestH + event.ak8jets[0].bestQCD + event.ak8jets[0].bestB))
        self.out.fillBranch("ak8_1_best_TvsQCD"       , event.ak8jets[0].bestT/(event.ak8jets[0].bestT + event.ak8jets[0].bestQCD + event.ak8jets[0].bestB))
        self.out.fillBranch("ak8_1_btagCSVV2" , event.ak8jets[0].btagCSVV2)
        self.out.fillBranch("ak8_1_btagHbb" , event.ak8jets[0].btagHbb)
        self.out.fillBranch("ak8_1_nnHbb" , event.ak8jets[0].nnHbb)
        self.out.fillBranch("ak8_1_nnHcc" , event.ak8jets[0].nnHcc)
        self.out.fillBranch("ak8_1_sj1_btagCSVV2" , event.ak8jets[0].subjets[0].btagCSVV2)
        self.out.fillBranch("ak8_1_sj2_btagCSVV2" , event.ak8jets[0].subjets[1].btagCSVV2)

        self.out.fillBranch("ca15_1_pt"           , event.ca15jets[0].pt)
        self.out.fillBranch("ca15_1_eta"          , event.ca15jets[0].eta)
        self.out.fillBranch("ca15_1_phi"          , event.ca15jets[0].phi)
        self.out.fillBranch("ca15_1_mass"         , event.ca15jets[0].msoftdrop)
        self.out.fillBranch("ca15_1_ecf0"         , event.ca15jets[0].ecf0)
        self.out.fillBranch("ca15_1_ecfTopTagBDT" , event.ca15jets[0].ecfTopTagBDT)
        self.out.fillBranch("ca15_1_httFRec"      , event.ca15jets[0].httFRec)
        self.out.fillBranch("ca15_1_tau32sd"      , event.ca15jets[0].tau32sd)


        
    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""

        self._correctJetAndMET(event)
       
        if self._prepareEvent(event) is False:
            return False

        # fill
        self._fillEventInfo(event)
        
        return True


# define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed
PhotonTree = lambda : PhotonSampleProducer()

