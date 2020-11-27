import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.tools import deltaPhi, deltaR, closest

from PhysicsTools.NanoHRTTools.producers.HeavyFlavBaseProducer import HeavyFlavBaseProducer

import logging
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] %(levelname)s: %(message)s')


class InclusiveSampleProducer(HeavyFlavBaseProducer):

    def __init__(self, **kwargs):
        super(InclusiveSampleProducer, self).__init__(channel='signal', **kwargs)

    def beginJob(self):
        super(InclusiveSampleProducer, self).beginJob()

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        super(InclusiveSampleProducer, self).beginFile(inputFile, outputFile, inputTree, wrappedOutputTree)

        ## event variables
        self.out.branch("ht", "F")
        self.out.branch("htfwd", "F")
        self.out.branch("htveryfwd", "F")
        self.out.branch("nlep", "I")
        self.out.branch("met", "F")
        self.out.branch("nlb_fj_pihalf" , "I")
        self.out.branch("nmb_fj_pihalf" , "I")
        self.out.branch("ntb_fj_pihalf" , "I")
        self.out.branch("nb_fj_pi" , "I")
        self.out.branch("nb_away_fj" , "I")

    def prepareEvent(self, event):

        logging.debug('processing event %d' % event.event)


        ## selection on AK8 jets
        event.fatjets = []
        for fj in event._allFatJets:
            if not (fj.pt > 200 and abs(fj.eta) < 2.4 and (fj.jetId & 2)):
                continue
            event.fatjets.append(fj)
        if len(event.fatjets) < 1:
            return False

        ## selection on SV
        event._allSV = Collection(event, "SV")
        event.secondary_vertices = []
        for sv in event._allSV:
#             if sv.ntracks > 2 and abs(sv.dxy) < 3. and sv.dlenSig > 4:
#             if sv.dlenSig > 4:
            if True:
                event.secondary_vertices.append(sv)
#        if len(event.secondary_vertices) < 2:
#            return False
        event.secondary_vertices = sorted(event.secondary_vertices, key=lambda x: x.pt, reverse=True)  # sort by pt
#         event.secondary_vertices = sorted(event.secondary_vertices, key=lambda x : x.dxySig, reverse=True)  # sort by dxysig

        self.matchSVToFatJets(event, event.fatjets)

        # selection on the probe jet (sub-leading in pT)
        probe_fj = event.fatjets[0]
        if not (probe_fj.pt > 200 and len(probe_fj.subjets) == 2 and probe_fj.msoftdrop > 30 and probe_fj.msoftdrop < 250):
            return False
        # require at least 1 SV matched to each subjet
        self.matchSVToSubjets(event, probe_fj)
        #if len(probe_fj.subjets[0].sv_list) == 0 or len(probe_fj.subjets[1].sv_list) == 0: ##LG put it back in
        #    return False        ## LG put it back in

        ## ht
        event.ak4jets = []
        for j in event._allJets:
            if not (j.pt > 25 and abs(j.eta) < 2.4 and (j.jetId & 2)):
                continue
            event.ak4jets.append(j)
        event.ht = sum([j.pt for j in event.ak4jets])
        if (event.ht<500):
            return False

        ## b-tag AK4 jet selection
        event.bljets = []
        event.bmjets = []
        event.btjets = []
        for j in event._allJets:
            if not (j.pt > 20.0 and abs(j.eta) < 2.4 and (j.jetId & 2)):
                continue
            if j.btagDeepB > self.DeepCSV_WP_L:
                event.bljets.append(j)
            if j.btagDeepB > self.DeepCSV_WP_M:
                event.bmjets.append(j)
            if j.btagDeepB > self.DeepCSV_WP_T:
                event.btjets.append(j)

        ## return True if passes selection
        return True

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""

        self.correctJetsAndMET(event)
        self.selectLeptons(event)

        if self.prepareEvent(event) is False:
            return False

        self.loadGenHistory(event)

        ## event variables
        htfwd_     = 0.;
        htveryfwd_ = 0.;
        for j in event._allJets:
            #if not (j.pt > 20 and abs(j.eta) < 2.4 and (j.jetId & 2)):
            if not (j.pt > 25 and (j.jetId & 2)):
                continue
            if ( abs(j.eta) > 2.4):
                htfwd_ += j.pt
            if ( abs(j.eta) > 3.5):
                htveryfwd_ += j.pt
    

        self.out.fillBranch("ht", event.ht)
        self.out.fillBranch("htfwd", htfwd_)
        self.out.fillBranch("htveryfwd", htveryfwd_)
        self.out.fillBranch("nlep", len(event.looseLeptons))
        self.out.fillBranch("met", event.met.pt)
        
        ## count bjets away from fatjet
        nlb_fj_pihalf_ = 0;
        nmb_fj_pihalf_ = 0;
        ntb_fj_pihalf_ = 0;
        nb_fj_pi_     = 0;
        for j in event.bljets:
            if abs(deltaPhi(j, event.fatjets[0])) > 3.14:
                nb_fj_pi_ += 1
            if abs(deltaPhi(j, event.fatjets[0])) > 3.14/2.:
                nlb_fj_pihalf_ += 1
        for j in event.bmjets:
            if abs(deltaPhi(j, event.fatjets[0])) > 3.14/2.:
                nmb_fj_pihalf_ += 1
        for j in event.btjets:
            if abs(deltaPhi(j, event.fatjets[0])) > 3.14/2.:
                ntb_fj_pihalf_ += 1
        
        nb_away_fj_ = 0;
        for j in event.bljets:
            if deltaR(j, event.fatjets[0])>1.:
                nb_away_fj_ += 1
            
        self.out.fillBranch("nlb_fj_pihalf" , nlb_fj_pihalf_)
        self.out.fillBranch("nmb_fj_pihalf" , nmb_fj_pihalf_)
        self.out.fillBranch("ntb_fj_pihalf" , ntb_fj_pihalf_)
        self.out.fillBranch("nb_fj_pi" , nb_fj_pi_)
        self.out.fillBranch("nb_away_fj" , nb_away_fj_)        

        self.fillBaseEventInfo(event)
        self.fillFatJetInfo(event)

        return True


# define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed
InclusiveTree_2016 = lambda: InclusiveSampleProducer(year=2016)
InclusiveTree_2017 = lambda: InclusiveSampleProducer(year=2017)
InclusiveTree_2018 = lambda: InclusiveSampleProducer(year=2018)

