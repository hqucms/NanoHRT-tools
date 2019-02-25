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


class topDecay:
    def __init__(self, topp4_, wp4_, bp4_, q1p4_, q2p4_):
        self.topp4 = topp4_
        self.wp4   = wp4_
        self.bp4   = bp4_
        self.q1p4  = q1p4_
        self.q2p4  = q2p4_


class genParticle:
    def __init__(self, p4_, idx_, pdgid_):
        self.p4    = p4_
        self.idx   = idx_
        self.pdgid = pdgid_


def get_subjets(jet, subjetCollection, idxNames=('subJetIdx1', 'subJetIdx2')):
    subjets = []
    for idxname in idxNames:
        idx = getattr(jet, idxname)
        if idx >= 0:
            subjets.append(subjetCollection[idx])
    return subjets


def get_sdmass(subjets):
    return sum([sj.p4() for sj in subjets], ROOT.TLorentzVector()).M()


def isQuark(gp):
    if ( (abs(gp.pdgId) == 1) or (abs(gp.pdgId) == 2) or (abs(gp.pdgId) == 3) or (abs(gp.pdgId) == 4) or (abs(gp.pdgId) == 5)):
        return True
    else:
        return False


class MuonSampleProducer(Module):

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
        self.out.branch("passMuTrig", "O")

        ## event variables

        ## Large-R jets
        self.out.branch("n_ak8", "I")
        self.out.branch("n_ca15", "I")

        self.out.branch("ak8_1_pt", "F")
        self.out.branch("ak8_1_eta", "F")
        self.out.branch("ak8_1_phi", "F")
        self.out.branch("ak8_1_mass", "F")
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
        self.out.branch("ak8_1_btagCSVV2", "F")
        self.out.branch("ak8_1_btagHbb", "F")
        self.out.branch("ak8_1_nnHbb", "F")
        self.out.branch("ak8_1_nnHcc", "F")
        self.out.branch("ak8_1_sj1_btagCSVV2", "F")
        self.out.branch("ak8_1_sj2_btagCSVV2", "F")

        self.out.branch("muon_pTrel", "F")
        self.out.branch("muon_drMuJet", "F")

        self.out.branch("ca15_1_pt", "F")
        self.out.branch("ca15_1_eta", "F")
        self.out.branch("ca15_1_phi", "F")
        self.out.branch("ca15_1_mass", "F")
        self.out.branch("ca15_1_ecf0", "F")
        self.out.branch("ca15_1_ecfTopTagBDT", "F")
        self.out.branch("ca15_1_httFRec", "F")
        self.out.branch("ca15_1_tau32sd", "F")

        self.out.branch("drgentop3q", "F")
        self.out.branch("drfj3q", "F")
        self.out.branch("drfj2q", "F")
        self.out.branch("drfjb", "F")


    def _correctJetAndMET(self, event):
        event._allJets = Collection(event, "Jet")
        event.ak8Subjets = Collection(event, "CustomAK8PuppiSubJet")  # do not sort after updating!!
        event.ca15Subjets = Collection(event, "CA15PuppiSubJet")  # do not sort after updating!!
        event.genParticles = Collection(event, "Jet")
        if self.isMC:
            event.genParticles = Collection(event, "GenPart")
        event.maxDR = -999.
        event.isFullyMerged = 0
        event.isSemiMerged = 0
        event.isUnMerged = 0
        event.muonpTrel = -999.9
        event.drMuJet = -999.9

        if self.isMC:
            rho = event.fixedGridRhoFastjetAll

        ## construct AK8 p4 from (updated) subjets
        event._allAK8jets = Collection(event, "CustomAK8Puppi")
        for fj in event._allAK8jets:
            fj.subjets = get_subjets(fj, event.ak8Subjets, ('subJetIdx1', 'subJetIdx2'))
            newP4 = ROOT.TLorentzVector()
            if len(fj.subjets) == 2:
                newP4 = fj.subjets[0].p4() + fj.subjets[1].p4()
            fj.pt, fj.eta, fj.phi, fj.mass, fj.msoftdrop = newP4.Pt(), newP4.Eta(), newP4.Phi(), newP4.M(), newP4.M()
        event._allAK8jets = sorted(event._allAK8jets, key=lambda x: x.pt, reverse=True)  # sort by pt

        ## construct CA15 p4 from (updated) subjets
        event._allCA15jets = Collection(event, "CA15Puppi")
        for fj in event._allCA15jets:
            fj.subjets = get_subjets(fj, event.ca15Subjets, ('subJetIdx1', 'subJetIdx2'))
            newP4 = ROOT.TLorentzVector()
            if len(fj.subjets) == 2:
                newP4 = fj.subjets[0].p4() + fj.subjets[1].p4()
            fj.pt, fj.eta, fj.phi, fj.mass, fj.msoftdrop = newP4.Pt(), newP4.Eta(), newP4.Phi(), newP4.M(), newP4.M()
        event._allCA15jets = sorted(event._allCA15jets, key=lambda x: x.pt, reverse=True)  # sort by pt



    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass


    def _prepareEvent(self, event):

        logging.debug('processing event %d' % event.event)

        ## met selection
        if event.MET_pt < 50.0:
            return False
        metLV = ROOT.TLorentzVector()
        metLV.SetPtEtaPhiM( event.MET_pt, 0.0, event.MET_phi, 0.0)


        ## muon selection
        event._allMuons = Collection(event, "Muon")
        event.muons = []
        for muon in event._allMuons:
            if abs(muon.eta) < 2.4 and\
               abs(muon.dxy) < 0.2 and\
               abs(muon.dz) < 0.5:# and\
               #muon.pfRelIso03_all < 0.15:
                drMuJet = deltaR(muon, event._allJets[muon.jetIdx]) if muon.jetIdx >= 0 else -999
                ptRel = muon.pt
                if drMuJet > -999:
                    jetforPtRel = event._allJets[muon.jetIdx]
                    jetLV = ROOT.TLorentzVector()
                    jetLV.SetPtEtaPhiM(jetforPtRel.pt, jetforPtRel.eta, jetforPtRel.phi, jetforPtRel.mass)
                    thismuonLV = ROOT.TLorentzVector()
                    thismuonLV.SetPtEtaPhiM(muon.pt, muon.eta, muon.phi, muon.mass)
                    ptRel = muon.pt * math.cos(thismuonLV.Angle(jetLV.Vect()))
                if ptRel > 55.0 or drMuJet > 0.4:
                    event.muons.append(muon)
                    event.muonpTrel = ptRel
                    event.drMuJet = drMuJet


        if (len(event.muons) != 1):
            return False

        muonLV = ROOT.TLorentzVector()
        muonLV.SetPtEtaPhiM( event.muons[0].pt, event.muons[0].eta, event.muons[0].phi, event.muons[0].mass)


        ## b-tag AK4 jet selection
        event.ak4jets = []
        for j in event._allJets:
            if not (j.pt > 25.0 and abs(j.eta) < 2.4 and (j.jetId & 2)):
                continue
            if j.btagCSVV2 > 0.8484 and\
               abs(deltaPhi(j, event.muons[0])) < 2.0:
                event.ak4jets.append(j)

        if (len(event.ak4jets) < 1):
            return False





        ## selection on AK8 jets / drop if overlaps with a photon
        event.ak8jets = []
        for fj in event._allAK8jets:
            if not (fj.pt > 200 and abs(fj.eta) < 2.4 and (fj.jetId & 2)):
                continue
            if abs(deltaPhi(fj, event.muons[0])) > 2.0:
                event.ak8jets.append(fj)

        if (len(event.ak8jets) < 1):
            return False

        ## selection on CA15 jets / drop if overlaps with a photon
        event.ca15jets = []
        for fj in event._allCA15jets:
            if not (fj.pt > 200 and abs(fj.eta) < 2.4 and (fj.jetId & 2)):
                continue
            if abs(deltaPhi(fj, event.muons[0])) > 2.0:
                event.ca15jets.append(fj)

        if (len(event.ca15jets) < 1):
            return False


        ## require the leading ak8 & ca15 jets overlap
        if deltaR(event.ak8jets[0], event.ca15jets[0]) > 0.8:
            return False

        ##leptonic W pt cut
        WLV = ROOT.TLorentzVector()
        WLV = muonLV + metLV
        if WLV.Pt() < 250.0:
            return False

        ##Find gen top and compute top size (-999. if unmatched)
        mindr       = 9999.
        event.drgentop3q_ = 999.
        event.drfj3q_     = 999.
        event.drfj2q_     = 999.
        event.drfjb_      = 999.
        if self.isMC:

            ## get the gen tops
            indexTop = -1
            event.genTop = []
            for p in event.genParticles:
                indexTop += 1
                if ( (abs(p.pdgId) == 6) and (p.status >= 60)):
                    pLV = ROOT.TLorentzVector()
                    #pLV = p.p4()
                    pLV.SetPtEtaPhiM(p.pt, p.eta, p.phi, p.mass)
                    tmpGenTop = genParticle(pLV, indexTop, p.pdgId)
                    event.genTop.append(tmpGenTop)
                else:
                    continue

            ## get the gen ws
            event.genW = []
            for top in event.genTop:
                indexW = -1
                for p in event.genParticles:
                    indexW += 1
                    if ( (abs(p.pdgId) == 24) and (p.status >= 50) and ((p.pdgId / abs(p.pdgId)) == (top.pdgid / abs(top.pdgid)))):
                        pLV = ROOT.TLorentzVector()
                        pLV.SetPtEtaPhiM(p.pt, p.eta, p.phi, p.mass)
                        tmpGenW = genParticle(pLV, indexW, p.pdgId)
                        event.genW.append(tmpGenW)
                    else:
                        continue

            ## get the W quarks and the b
            index = -1
            event.GenTopDecay = []
            for w in event.genW:
                index += 1
                indexQ = -1
                genQ = []
                genb = []
                for p in event.genParticles:
                    indexQ += 1
                    if (p.genPartIdxMother < 0):
                        continue
                    if (isQuark(p) and (abs(p.pdgId) != 5) and ((event.genParticles[p.genPartIdxMother].pdgId) / (w.pdgid) == 1)):
                        pLV = ROOT.TLorentzVector()
                        pLV.SetPtEtaPhiM(p.pt, p.eta, p.phi, p.mass)
                        tmpGenQ = genParticle(pLV, indexQ, p.pdgId)
                        genQ.append(tmpGenQ)
                    if ( (abs(p.pdgId) == 5) and (event.genParticles[p.genPartIdxMother].pdgId == top.pdgid)):
                        pLV = ROOT.TLorentzVector()
                        pLV.SetPtEtaPhiM(p.pt, p.eta, p.phi, p.mass)
                        tmpGenb = genParticle(pLV, indexQ, p.pdgId)
                        genb.append(tmpGenb)

                ## check if we have three quarks
                if ( (len(genQ) == 2) and (len(genb) == 1)):
                    tmpGenTopDecay = topDecay(event.genTop[index].p4, event.genW[index].p4, genb[0].p4, genQ[0].p4, genQ[1].p4)
                    event.GenTopDecay.append(tmpGenTopDecay)

            ## Define  categories
            for i in event.GenTopDecay:
                tmpdr = deltaR(event.ak8jets[0].eta, event.ak8jets[0].phi, i.topp4.Eta(), i.topp4.Phi())
                if (tmpdr < mindr):
                    mindr = tmpdr
                    drgentop3q_ = max(deltaR(i.topp4.Eta(), i.topp4.Phi(), i.bp4.Eta(), i.bp4.Phi()),\
                                          deltaR(i.topp4.Eta(), i.topp4.Phi(), i.q1p4.Eta(), i.q1p4.Phi()),\
                                          deltaR(i.topp4.Eta(), i.topp4.Phi(), i.q1p4.Eta(), i.q1p4.Phi()))

                    drfj3q_ = max(deltaR(event.ak8jets[0].eta, event.ak8jets[0].phi, i.bp4.Eta(), i.bp4.Phi()),\
                                      deltaR(event.ak8jets[0].eta, event.ak8jets[0].phi, i.q1p4.Eta(), i.q1p4.Phi()),\
                                      deltaR(event.ak8jets[0].eta, event.ak8jets[0].phi, i.q1p4.Eta(), i.q1p4.Phi()))

                    drfj2q_ = max(deltaR(event.ak8jets[0].eta, event.ak8jets[0].phi, i.q1p4.Eta(), i.q1p4.Phi()),\
                                      deltaR(event.ak8jets[0].eta, event.ak8jets[0].phi, i.q1p4.Eta(), i.q1p4.Phi()))

                    drfjb_ = deltaR(event.ak8jets[0].eta, event.ak8jets[0].phi, i.bp4.Eta(), i.bp4.Phi())

        ## return True if passes selection
        return True


    def _fillEventInfo(self, event):

        ## Triggers
        passMuTrig_ = False

        try:
            if event.HLT_Mu50:
                passMuTrig_ = True
        except:
            passMuTrig_ = False
        try:
            if event.HLT_TkMu50:
                passMuTrig_ = True
        except:
            passMuTrig_ = False

        self.out.fillBranch("passMuTrig", passMuTrig_)


        ## event variables

        ## Large-R jets
        self.out.fillBranch("n_ak8", len(event.ak8jets))
        self.out.fillBranch("n_ca15", len(event.ca15jets))

        self.out.fillBranch("ak8_1_pt", event.ak8jets[0].pt)
        self.out.fillBranch("ak8_1_eta", event.ak8jets[0].eta)
        self.out.fillBranch("ak8_1_phi", event.ak8jets[0].phi)
        self.out.fillBranch("ak8_1_mass", event.ak8jets[0].msoftdrop)
        self.out.fillBranch("ak8_1_n2b1", event.ak8jets[0].n2b1)
        self.out.fillBranch("ak8_1_n3b1", event.ak8jets[0].n3b1)
        self.out.fillBranch("ak8_1_tau1", event.ak8jets[0].tau1)
        self.out.fillBranch("ak8_1_tau2", event.ak8jets[0].tau2)
        self.out.fillBranch("ak8_1_tau3", event.ak8jets[0].tau3)
        self.out.fillBranch("ak8_1_DeepAK8_WvsQCD", get_nominal_score(event.ak8jets[0], 'WvsQCD'))
        self.out.fillBranch("ak8_1_DeepAK8_ZvsQCD", get_nominal_score(event.ak8jets[0], 'ZvsQCD'))
        self.out.fillBranch("ak8_1_DeepAK8_HvsQCD", get_nominal_score(event.ak8jets[0], 'HbbvsQCD'))
        self.out.fillBranch("ak8_1_DeepAK8_TvsQCD", get_nominal_score(event.ak8jets[0], 'TvsQCD'))
        self.out.fillBranch("ak8_1_DeepAK8MD_WvsQCD", get_decorr_score(event.ak8jets[0], 'WvsQCD'))
        self.out.fillBranch("ak8_1_DeepAK8MD_ZvsQCD", get_decorr_score(event.ak8jets[0], 'ZHbbvsQCD'))
        self.out.fillBranch("ak8_1_DeepAK8MD_HvsQCD", get_decorr_score(event.ak8jets[0], 'ZHbbvsQCD'))
        self.out.fillBranch("ak8_1_DeepAK8MD_TvsQCD", get_decorr_score(event.ak8jets[0], 'TvsQCD'))
        self.out.fillBranch("ak8_1_best_WvsQCD", event.ak8jets[0].bestW / (event.ak8jets[0].bestW + event.ak8jets[0].bestQCD + event.ak8jets[0].bestB))
        self.out.fillBranch("ak8_1_best_ZvsQCD", event.ak8jets[0].bestZ / (event.ak8jets[0].bestZ + event.ak8jets[0].bestQCD + event.ak8jets[0].bestB))
        self.out.fillBranch("ak8_1_best_HvsQCD", event.ak8jets[0].bestH / (event.ak8jets[0].bestH + event.ak8jets[0].bestQCD + event.ak8jets[0].bestB))
        self.out.fillBranch("ak8_1_best_TvsQCD", event.ak8jets[0].bestT / (event.ak8jets[0].bestT + event.ak8jets[0].bestQCD + event.ak8jets[0].bestB))
        self.out.fillBranch("ak8_1_btagCSVV2", event.ak8jets[0].btagCSVV2)
        self.out.fillBranch("ak8_1_btagHbb", event.ak8jets[0].btagHbb)
        self.out.fillBranch("ak8_1_nnHbb", event.ak8jets[0].nnHbb)
        self.out.fillBranch("ak8_1_nnHcc", event.ak8jets[0].nnHcc)
        self.out.fillBranch("ak8_1_sj1_btagCSVV2", event.ak8jets[0].subjets[0].btagCSVV2)
        self.out.fillBranch("ak8_1_sj2_btagCSVV2", event.ak8jets[0].subjets[1].btagCSVV2)
        self.out.fillBranch("muon_pTrel", event.muonpTrel)
        self.out.fillBranch("muon_drMuJet", event.drMuJet)

        self.out.fillBranch("ca15_1_pt", event.ca15jets[0].pt)
        self.out.fillBranch("ca15_1_eta", event.ca15jets[0].eta)
        self.out.fillBranch("ca15_1_phi", event.ca15jets[0].phi)
        self.out.fillBranch("ca15_1_mass", event.ca15jets[0].msoftdrop)
        self.out.fillBranch("ca15_1_ecf0", event.ca15jets[0].ecf0)
        self.out.fillBranch("ca15_1_ecfTopTagBDT", event.ca15jets[0].ecfTopTagBDT)
        self.out.fillBranch("ca15_1_httFRec", event.ca15jets[0].httFRec)
        self.out.fillBranch("ca15_1_tau32sd", event.ca15jets[0].tau32sd)

        self.out.fillBranch("drgentop3q", event.drgentop3q_)
        self.out.fillBranch("drfj3q", event.drfj3q_)
        self.out.fillBranch("drfj2q", event.drfj2q_)
        self.out.fillBranch("drfjb", event.drfjb_)


    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""

        self._correctJetAndMET(event)

        if self._prepareEvent(event) is False:
            return False

        # fill
        self._fillEventInfo(event)

        return True


# define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed
MuonTree = lambda: MuonSampleProducer()

