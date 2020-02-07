import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.tools import deltaPhi, deltaR, closest

from PhysicsTools.NanoHRTTools.producers.HeavyFlavBaseProducer import HeavyFlavBaseProducer

import logging
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] %(levelname)s: %(message)s')


class PhotonSampleProducer(HeavyFlavBaseProducer):

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
        self.out.branch("nlep", "I")

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
            if not (pho.pt > 200 and abs(pho.eta) < 2.4 and (pho.cutBasedBitmap & 2) and pho.electronVeto):  # medium ID
                continue
            event.photons.append(pho)

        if len(event.photons) < 1:
            return False

        ## selection on AK8 jets / drop if overlaps with a photon
        event.fatjets = []
        for fj in event._allFatJets:
            if not (fj.pt > 200 and abs(fj.eta) < 2.4 and (fj.jetId & 2)):
                continue
            # require jet and photon to be back-to-back
            if deltaPhi(event.photons[0], fj) < 2:
                continue
#             if deltaR(event.photons[0], fj) < self._jetConeSize:
#                 continue
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
        if len(event.secondary_vertices) < 2:
            return False
        event.secondary_vertices = sorted(event.secondary_vertices, key=lambda x: x.pt, reverse=True)  # sort by pt
#         event.secondary_vertices = sorted(event.secondary_vertices, key=lambda x : x.dxySig, reverse=True)  # sort by dxysig

        # selection on the probe jet (sub-leading in pT)
        probe_fj = event.fatjets[0]
        if not (probe_fj.pt > 200 and probe_fj.msoftdrop > 50 and probe_fj.msoftdrop < 200):
            return False
        # require at least 1 SV matched to each subjet
        self.matchSVToSubjets(event, probe_fj)
        if len(probe_fj.subjets[0].sv_list) == 0 or len(probe_fj.subjets[1].sv_list) == 0:
            return False

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

        if self.year == 2016:
            self.out.fillBranch("passTrigPhoton", event.HLT_Photon175)
        else:
            self.out.fillBranch("passTrigPhoton", event.HLT_Photon200)

        ## event variables
        self.out.fillBranch("ht", event.ht)
        self.out.fillBranch("nlep", len(event.looseLeptons))

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

