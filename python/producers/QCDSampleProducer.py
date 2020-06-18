import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.tools import deltaPhi, deltaR, closest

from PhysicsTools.NanoHRTTools.producers.HRTBaseProducer import HRTBaseProducer

import logging
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] %(levelname)s: %(message)s')


class QCDSampleProducer(HRTBaseProducer):

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
        event.ak8jets = []
        for fj in event._allAK8jets:
            if not (fj.pt > 200 and abs(fj.eta) < 2.4 and (fj.jetId & 2)):
                continue
            event.ak8jets.append(fj)
        if len(event.ak8jets) < 2:
            return False

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

        if self._run_tagger:
            self.runParticleNet(event)

        return True


# define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed
QCDTree_2016 = lambda: QCDSampleProducer(year=2016)
QCDTree_2017 = lambda: QCDSampleProducer(year=2017)
QCDTree_2018 = lambda: QCDSampleProducer(year=2018)
