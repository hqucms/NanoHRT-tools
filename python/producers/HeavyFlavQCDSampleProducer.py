import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.tools import deltaPhi, deltaR, closest

from PhysicsTools.NanoHRTTools.producers.HeavyFlavBaseProducer import HeavyFlavBaseProducer

import logging
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] %(levelname)s: %(message)s')


class QCDSampleProducer(HeavyFlavBaseProducer):

    def __init__(self, **kwargs):
        super(QCDSampleProducer, self).__init__(channel='qcd', **kwargs)

    def beginJob(self):
        super(QCDSampleProducer, self).beginJob()

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        super(QCDSampleProducer, self).beginFile(inputFile, outputFile, inputTree, wrappedOutputTree)

        ## trigger variables
        self.out.branch("passHTTrig", "O")

        ## event variables
        self.out.branch("ht", "F")
        self.out.branch("nlep", "I")

    def prepareEvent(self, event):

        logging.debug('processing event %d' % event.event)

        # # ht selection
        event.ak4jets = []
        for j in event._allJets:
            if not (j.pt > 25 and abs(j.eta) < 2.4 and (j.jetId & 2)):
                continue
            event.ak4jets.append(j)

        event.ht = sum([j.pt for j in event.ak4jets])
        if event.ht < 1000.:
            return False

        ## selection on AK8 jets
        event.fatjets = []
        for fj in event._allFatJets:
            if not (fj.pt > 200 and abs(fj.eta) < 2.4 and (fj.jetId & 2)):
                continue
            event.fatjets.append(fj)
        if len(event.fatjets) < 2:
            return False

        ## selection on SV
        event._allSV = Collection(event, "SV")
        event.secondary_vertices = []
        for sv in event._allSV:
#             if sv.ntracks > 2 and abs(sv.dxy) < 3. and sv.dlenSig > 4:
#             if sv.dlenSig > 4:
            if True:
                event.secondary_vertices.append(sv)
        if len(event.secondary_vertices) < 2:
            return False
        event.secondary_vertices = sorted(event.secondary_vertices, key=lambda x: x.pt, reverse=True)  # sort by pt
#         event.secondary_vertices = sorted(event.secondary_vertices, key=lambda x : x.dxySig, reverse=True)  # sort by dxysig

        # selection on the probe jet (sub-leading in pT)
        probe_fj = event.fatjets[1]
        if not (len(probe_fj.subjets) == 2 and probe_fj.msoftdrop > 50 and probe_fj.msoftdrop < 200):
            return False
        # require at least 1 SV matched to each subjet
        self.matchSVToSubjets(event, probe_fj)
        if len(probe_fj.subjets[0].sv_list) == 0 or len(probe_fj.subjets[1].sv_list) == 0:
            return False
        # match SV also to the leading jet
        if not (len(event.fatjets[0].subjets) == 2):
            return False
        self.matchSVToSubjets(event, event.fatjets[0])

        # load gen
        self.loadGenHistory(event)

        ## return True if passes selection
        return True

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""

        self.correctJetsAndMET(event)
        self.selectLeptons(event)

        if self.prepareEvent(event) is False:
            return False

        if self.year == 2016:
            self.out.fillBranch("passHTTrig", event.HLT_PFHT900)
        else:
            self.out.fillBranch("passHTTrig", event.HLT_PFHT1050)
        self.out.fillBranch("ht", event.ht)
        self.out.fillBranch("nlep", len(event.looseLeptons))

        self.fillBaseEventInfo(event)
        self.fillFatJetInfo(event)

        return True


# define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed
QCDTree_2016 = lambda: QCDSampleProducer(year=2016)
QCDTree_2017 = lambda: QCDSampleProducer(year=2017)
QCDTree_2018 = lambda: QCDSampleProducer(year=2018)
