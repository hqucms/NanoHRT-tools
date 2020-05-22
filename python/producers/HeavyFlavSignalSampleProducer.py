import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.tools import deltaPhi, deltaR, closest

from PhysicsTools.NanoHRTTools.producers.HeavyFlavBaseProducer import HeavyFlavBaseProducer

import logging
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] %(levelname)s: %(message)s')


class SignalSampleProducer(HeavyFlavBaseProducer):

    def __init__(self, **kwargs):
        super(SignalSampleProducer, self).__init__(channel='signal', **kwargs)

    def beginJob(self):
        super(SignalSampleProducer, self).beginJob()

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        super(SignalSampleProducer, self).beginFile(inputFile, outputFile, inputTree, wrappedOutputTree)

        ## event variables
        self.out.branch("ht", "F")
        self.out.branch("nlep", "I")

    def prepareEvent(self, event):

        logging.debug('processing event %d' % event.event)


        ## selection on AK8 jets / drop if overlaps with a photon
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

        self.matchSVToJets(event, event.fatjets)

        # selection on the probe jet (sub-leading in pT)
        probe_fj = event.fatjets[0]
        if not (probe_fj.pt > 200 and len(probe_fj.subjets) == 2 and probe_fj.msoftdrop > 50 and probe_fj.msoftdrop < 200):
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
        self.out.fillBranch("ht", event.ht)
        self.out.fillBranch("nlep", len(event.looseLeptons))
        self.fillBaseEventInfo(event)
        self.fillFatJetInfo(event,True)

        return True


# define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed
SignalTree_2016 = lambda: SignalSampleProducer(year=2016)
SignalTree_2017 = lambda: SignalSampleProducer(year=2017)
SignalTree_2018 = lambda: SignalSampleProducer(year=2018)

