import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.tools import deltaPhi, deltaR, closest

from PhysicsTools.NanoHRTTools.producers.HRTBaseProducer import HRTBaseProducer

import logging
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] %(levelname)s: %(message)s')


class PhotonSampleProducer(HRTBaseProducer):

    def __init__(self, **kwargs):
        super(PhotonSampleProducer, self).__init__(channel='photon', **kwargs)

    def beginJob(self):
        super(PhotonSampleProducer, self).beginJob()

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        super(PhotonSampleProducer, self).beginFile(inputFile, outputFile, inputTree, wrappedOutputTree)

        ## trigger variables
        self.out.branch("passTrigPhoton", "O")

        ## event variables
        self.out.branch("ht", "F")

        ## photons
        self.out.branch("nphotons", "I")
        self.out.branch("pho_1_pt", "F")
        self.out.branch("pho_1_eta", "F")
        self.out.branch("pho_1_phi", "F")

    def prepareEvent(self, event):

        logging.debug('processing event %d' % event.event)

        ## select leading photon
        event._allPhotons = Collection(event, "Photon")
        event.photons = []
        for pho in event._allPhotons:
            if not (pho.pt > 200 and abs(pho.eta) < 2.4 and (pho.cutBased >= 1) and pho.electronVeto):  # loose ID
                continue
            event.photons.append(pho)

        if len(event.photons) < 1:
            return False


        ## selection on AK8 jets / drop if overlaps with a photon
        event.ak8jets = []
        for fj in event._allAK8jets:
            if not (fj.pt > 200 and abs(fj.eta) < 2.4 and (fj.jetId & 2)):
                continue
            if deltaR(event.photons[0], fj) < 0.8:
                continue
            event.ak8jets.append(fj)

        if len(event.ak8jets) < 1:
            return False

        ## ht selection
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

        if self.prepareEvent(event) is False:
            return False

        # fill
        if self.year == 2016:
            self.out.fillBranch("passTrigPhoton", event.HLT_Photon175)
        else:
            self.out.fillBranch("passTrigPhoton", event.HLT_Photon200)

        ## event variables
        self.out.fillBranch("ht", event.ht)

        ## photon variables
        self.out.fillBranch("nphotons", len(event.photons))
        self.out.fillBranch("pho_1_pt", event.photons[0].pt)
        self.out.fillBranch("pho_1_eta", event.photons[0].eta)
        self.out.fillBranch("pho_1_phi", event.photons[0].phi)

        self.fillBaseEventInfo(event)
        self.fillFatJetInfo(event)

        return True


# define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed
PhotonTree_2016 = lambda: PhotonSampleProducer(year=2016)
PhotonTree_2017 = lambda: PhotonSampleProducer(year=2017)
PhotonTree_2018 = lambda: PhotonSampleProducer(year=2018)

