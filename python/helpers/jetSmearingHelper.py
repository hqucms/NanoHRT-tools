import ROOT
import os, tarfile, tempfile, shutil
import numpy as np
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.tools import deltaR


def find_and_extract_tarball(name, destination):
    search_pathes = [os.environ['CMSSW_BASE'] + '/src/PhysicsTools/NanoHRTTools/data/jme/',
                     os.environ['CMSSW_BASE'] + '/src/PhysicsTools/NanoAODTools/data/jme/']
    for p in search_pathes:
        for ext in ['.tgz', '.tar.gz']:
            fullpath = os.path.join(p, name + ext)
            if os.path.exists(fullpath):
                with tarfile.open(fullpath, "r:gz") as tar:
                    tar.extractall(destination)
                print('... extracted %s to %s' % (fullpath, destination))
                return fullpath


def match(jet, genjets, resolution, coneSize=0.4):
    # Try to find a gen jet matching
    # dR < m_dR_max
    # dPt < m_dPt_max_factor * resolution
    minDR = 1e99
    matched = None
    for genj in genjets:
        dR = deltaR(genj, jet)
        if dR > minDR:
            continue
        if dR < 0.5 * coneSize:
            dPT = abs(genj.pt - jet.pt)
            if dPT > 3. * resolution:
                continue
        minDR = dR
        matched = genj
    return matched


class jetSmearer(Module):

    def __init__(self, jerTag, jetType="AK4PFchs"):

        self.jerTag = jerTag
        self.jetType = jetType
        if 'ak8' in jetType.lower():
            self.coneSize = 0.8
        elif 'ak4' in jetType.lower():
            self.coneSize = 0.4
        else:
            raise RuntimeError('Jet type %s is not recognized!' % jetType)

    def beginJob(self):
        # read jet energy resolution (JER) and JER scale factors and uncertainties
        # get latest version from: https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution
        self.jerInputFilePath = tempfile.mkdtemp()
        find_and_extract_tarball(self.jerTag, self.jerInputFilePath)
        self.jerInputFile = os.path.join(self.jerInputFilePath, '%s_PtResolution_%s.txt' % (self.jerTag, self.jetType))
        self.jerUncertaintyInputFile = os.path.join(self.jerInputFilePath, '%s_SF_%s.txt' % (self.jerTag, self.jetType))

        self.params_sf_and_uncertainty = ROOT.PyJetParametersWrapper()
        self.params_resolution = ROOT.PyJetParametersWrapper()

        # initialize random number generator
        # (needed for jet pT smearing)
        self.rnd = ROOT.TRandom3(12345)

        # load libraries for accessing JER scale factors and uncertainties from txt files
        for library in ["libCondFormatsJetMETObjects", "libPhysicsToolsNanoAODTools"]:
            if library not in ROOT.gSystem.GetLibraries():
                print("Load Library '%s'" % library.replace("lib", ""))
                ROOT.gSystem.Load(library)

        # initialize JER scale factors and uncertainties
        # (cf. PhysicsTools/PatUtils/interface/SmearedJetProducerT.h )
        print("Loading jet energy resolutions (JER) from file '%s'" % self.jerInputFile)
        self.jer = ROOT.PyJetResolutionWrapper(self.jerInputFile)
        print("Loading JER scale factors and uncertainties from file '%s'" % self.jerUncertaintyInputFile)
        self.jerSF_and_Uncertainty = ROOT.PyJetResolutionScaleFactorWrapper(self.jerUncertaintyInputFile)

    def endJob(self):
        shutil.rmtree(self.jerInputFilePath)

    def setSeed(self, seed):
        self.rnd.SetSeed(seed)

    def getSmearValsPt(self, jet, genjets, rho):

        # --------------------------------------------------------------------------------------------
        # CV: Smear jet pT to account for measured difference in JER between data and simulation.
        #     The function computes the nominal smeared jet pT simultaneously with the JER up and down shifts,
        #     in order to use the same random number to smear all three (for consistency reasons).
        #
        #     The implementation of this function follows PhysicsTools/PatUtils/interface/SmearedJetProducerT.h
        #
        # --------------------------------------------------------------------------------------------

        if not jet.pt > 0.:
            print("WARNING: jet pT = %1.1f !!" % jet.pt)
            return (1., 1., 1.)

        # --------------------------------------------------------------------------------------------
        # CV: define enums needed to access JER scale factors and uncertainties
        #    (cf. CondFormats/JetMETObjects/interface/JetResolutionObject.h)
        enum_nominal         = 0
        enum_shift_up        = 2
        enum_shift_down      = 1
        # --------------------------------------------------------------------------------------------

        self.params_resolution.setJetPt(jet.pt)
        self.params_resolution.setJetEta(jet.eta)
        self.params_resolution.setRho(rho)
        jet_pt_resolution = self.jer.getResolution(self.params_resolution)

        jet_pt_sf_and_uncertainty = {}
        for enum_central_or_shift in [enum_nominal, enum_shift_up, enum_shift_down]:
            self.params_sf_and_uncertainty.setJetEta(jet.eta)
            self.params_sf_and_uncertainty.setJetPt(jet.pt)
            jet_pt_sf_and_uncertainty[enum_central_or_shift] = self.jerSF_and_Uncertainty.getScaleFactor(self.params_sf_and_uncertainty, enum_central_or_shift)

        matched_genjet = match(jet, genjets, jet_pt_resolution * jet.pt, coneSize=self.coneSize)

        smear_vals = {}
        for central_or_shift in [enum_nominal, enum_shift_up, enum_shift_down]:

            smearFactor = None
            if matched_genjet:
                #
                # Case 1: we have a "good" generator level jet matched to the reconstructed jet
                #
                dPt = jet.pt - matched_genjet.pt
                smearFactor = 1. + (jet_pt_sf_and_uncertainty[central_or_shift] - 1.) * dPt / jet.pt
            elif jet_pt_sf_and_uncertainty[central_or_shift] > 1.:
                #
                # Case 2: we don't have a generator level jet. Smear jet pT using a random Gaussian variation
                #
                sigma = jet_pt_resolution * np.sqrt(jet_pt_sf_and_uncertainty[central_or_shift] ** 2 - 1.)
                smearFactor = self.rnd.Gaus(1., sigma)
            else:
                #
                # Case 3: we cannot smear this jet, as we don't have a generator level jet and the resolution in data is better than the resolution in the simulation,
                #         so we would need to randomly "unsmear" the jet, which is impossible
                #
                smearFactor = 1.

            # check that smeared jet energy remains positive,
            # as the direction of the jet would change ("flip") otherwise - and this is not what we want
            if (smearFactor * jet.pt) < 1.e-2:
                smearFactor = 1.e-2

            smear_vals[central_or_shift] = smearFactor

        return (smear_vals[enum_nominal], smear_vals[enum_shift_up], smear_vals[enum_shift_down])

    def getSmearValsM(self, jet, gensubjets):

        # --------------------------------------------------------------------------------------------
        # CV: Smear jet m to account for measured difference in JER between data and simulation.
        #     The function computes the nominal smeared jet m simultaneously with the JER up and down shifts,
        #     in order to use the same random number to smear all three (for consistency reasons).
        #
        #     The implementation of this function follows PhysicsTools/PatUtils/interface/SmearedJetProducerT.h
        #
        # --------------------------------------------------------------------------------------------

        if not hasattr(jet, 'subjets') or len(jet.subjets) != 2 or jet.mass <= 0:
        #    print("WARNING: jet does not have 2 subjets")
            return (1., 1., 1.)

        # --------------------------------------------------------------------------------------------
        # CV: define enums needed to access JER scale factors and uncertainties
        #    (cf. CondFormats/JetMETObjects/interface/JetResolutionObject.h)
        enum_nominal = 0
        enum_shift_up = 2
        enum_shift_down = 1
        # --------------------------------------------------------------------------------------------

        # from https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetWtagging
        jet_m_resolution = 10.1
        jet_m_sf_and_uncertainty = dict(zip([enum_nominal, enum_shift_up, enum_shift_down], [1.0, 1.2, 0.8]))

        matched_genjets = [match(sj, gensubjets, 1.e9) for sj in jet.subjets]
        gensdmass = None
        if matched_genjets[0] is not None and matched_genjets[1] is not None:
            gensdmass = (matched_genjets[0].p4() + matched_genjets[1].p4()).M()

        smear_vals = {}
        for central_or_shift in [enum_nominal, enum_shift_up, enum_shift_down]:

            smearFactor = None
            if gensdmass is not None and abs(jet.mass - gensdmass) < 3 * jet_m_resolution:
                #
                # Case 1: we have a "good" generator level jet matched to the reconstructed jet (and dM < 3 sigma)
                #
                dM = jet.mass - gensdmass
                smearFactor = 1. + (jet_m_sf_and_uncertainty[central_or_shift] - 1.) * dM / jet.mass
            elif jet_m_sf_and_uncertainty[central_or_shift] > 1.:
                #
                # Case 2: we don't have a generator level jet. Smear jet m using a random Gaussian variation
                #
                sigma = jet_m_resolution * np.sqrt(jet_m_sf_and_uncertainty[central_or_shift] ** 2 - 1.)
                smearFactor = self.rnd.Gaus(1., sigma)
            else:
                #
                # Case 3: we cannot smear this jet, as we don't have a generator level jet and the resolution in data is better than the resolution in the simulation,
                #         so we would need to randomly "unsmear" the jet, which is impossible
                #
                smearFactor = 1.

            # check that smeared jet energy remains positive,
            # as the direction of the jet would change ("flip") otherwise - and this is not what we want
            if (smearFactor * jet.mass) < 1.e-2:
                smearFactor = 1.e-2

            smear_vals[central_or_shift] = smearFactor

        return (smear_vals[enum_nominal], smear_vals[enum_shift_up], smear_vals[enum_shift_down])
