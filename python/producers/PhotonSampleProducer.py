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
        self.out.branch("passPhoton165_HE10", "O")

        ## event variables
        self.out.branch("ht", "F")

        ## photons
        self.out.branch("nphotons", "I")
        self.out.branch("pho_1_pt", "F")
        self.out.branch("pho_1_eta", "F")

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

        ## selection on CA15 jets / drop if overlaps with a photon
        event.ca15jets = []
        for fj in event._allCA15jets:
            if not (fj.pt > 200 and abs(fj.eta) < 2.4 and (fj.jetId & 2)):
                continue
            if deltaR(event.photons[0], fj) < 1.5:
                continue
            event.ca15jets.append(fj)

        if len(event.ca15jets) < 1:
            return False


        ## require the leading ak8 & ca15 jets overlap
        if deltaR(event.ak8jets[0], event.ca15jets[0]) > 0.8:
            return False

        # # selection on HOTVR jets
        event.hotvrjets = []
        for fj in event._allHOTVRjets:
            if not (fj.pt > 200 and abs(fj.eta) < 2.4):
                continue
            event.hotvrjets.append(fj)

        ## ht selection
        event.ak4jets = []
        for j in event._allJets:
            if not (j.pt > 25 and abs(j.eta) < 2.4 and (j.jetId & 2)):
                continue
            event.ak4jets.append(j)

        event.ht = sum([j.pt for j in event.ak4jets])
        if event.ht < 200:
            return False

        ## return True if passes selection
        return True

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""

        self.correctJetsAndMET(event)

        if self.prepareEvent(event) is False:
            return False

        # fill
        self.out.fillBranch("passPhoton165_HE10", event.HLT_Photon165_HE10)

        ## event variables
        self.out.fillBranch("ht", event.ht)

        ## photon variables
        self.out.fillBranch("nphotons", len(event.photons))
        self.out.fillBranch("pho_1_pt", event.photons[0].pt)
        self.out.fillBranch("pho_1_eta", event.photons[0].eta)

        self.fillBaseEventInfo(event)
        self.fillFatJetInfo(event)

        return True


# define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed
PhotonTree = lambda: PhotonSampleProducer()

