import itertools
import numpy as np
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.tools import deltaR


class TopCand(object):
    def __init__(self, j1, j2, j3):
        self.jets = [j1, j2, j3]
        self.jets.sort(key=lambda x: x.pt, reverse=True)
        self.candP4 = j1.p4() + j2.p4() + j3.p4()

    def doGenMatch(self, hadGenTops):
        def match(top, jet):
            if deltaR(top.genB, jet) < 0.4:
                return 'b'
            if deltaR(top.genW.daus[0], jet) < 0.4:
                return 'wq1'
            if deltaR(top.genW.daus[1], jet) < 0.4:
                return 'wq2'
            return None

        nMatchedJets = 0
        genMatchResult = set()
        matchedGenTop = None
        for genTop in hadGenTops:
            matches = set()
            for j in self.jets:
                rlt = match(genTop, j)
                if rlt:
                    matches.add(rlt)
            # print(matches, self.candP4.M(), max([j.btagDeepFlavB for j in self.jets]))
            if len(matches) > nMatchedJets:
                nMatchedJets = len(matches)
                genMatchResult = matches
                matchedGenTop = genTop

        if len(genMatchResult) == 3:
            self.genMatch = 'bqq'
            self.genMatchCode = 305
        elif len(genMatchResult) == 2:
            if 'b' in genMatchResult:
                self.genMatch = 'bq'
                self.genMatchCode = 25
            else:
                self.genMatch = 'qq'
                self.genMatchCode = 21
        elif len(genMatchResult) == 1:
            if 'b' in genMatchResult:
                self.genMatch = 'b'
                self.genMatchCode = 5
            else:
                self.genMatch = 'q'
                self.genMatchCode = 1
        else:
            self.genMatch = ''
            self.genMatchCode = 0

        self.matchedGenTop = matchedGenTop

    def makeFeatures(self):
        features = {}
        for idx in range(3):
            prefix = 'j%d_' % (idx + 1)
            features[prefix + "pt"] = self.jets[idx].pt
            features[prefix + "eta"] = self.jets[idx].eta
            features[prefix + "phi"] = self.jets[idx].phi
            features[prefix + "mass"] = self.jets[idx].mass
            features[prefix + "deepFlavB"] = self.jets[idx].btagDeepFlavB
            features[prefix + "deepFlavCvL"] = self.jets[idx].btagDeepFlavCvL
            features[prefix + "deepFlavCvB"] = self.jets[idx].btagDeepFlavCvB
            features[prefix + "deepFlavQG"] = self.jets[idx].btagDeepFlavQG

            for jdx in range(idx + 1, 3):
                prefix = 'j%d%d_' % (idx + 1, jdx + 1)
                features[prefix + "deltaR"] = deltaR(self.jets[idx], self.jets[jdx])
                features[prefix + "mass"] = (self.jets[idx].p4() + self.jets[jdx].p4()).M()

        features["topcand_mass"] = self.candP4.M()
        self.features = features


class ResTopTreeProducer(Module, object):

    def __init__(self):
        pass

    def beginJob(self):
        pass

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.isMC = bool(inputTree.GetBranch('genWeight'))
        self.hasParticleNetAK4 = bool(inputTree.GetBranch('Jet_particleNetAK4Tau_b'))

        self.out = wrappedOutputTree

        self.out.branch("event", "I")

        self.out.branch("isHadTop", "O")
        self.out.branch("genMatch", "I")
        self.out.branch("gentop_pt", "F")
        self.out.branch("gentop_eta", "F")

        for idx in range(3):
            prefix = 'j%d_' % (idx + 1)
            self.out.branch(prefix + "pt", "F")
            self.out.branch(prefix + "eta", "F")
            self.out.branch(prefix + "phi", "F")
            self.out.branch(prefix + "mass", "F")
            self.out.branch(prefix + "deepFlavB", "F")
            self.out.branch(prefix + "deepFlavCvL", "F")
            self.out.branch(prefix + "deepFlavCvB", "F")
            self.out.branch(prefix + "deepFlavQG", "F")

        for idx in range(3):
            for jdx in range(idx + 1, 3):
                prefix = 'j%d%d_' % (idx + 1, jdx + 1)
                self.out.branch(prefix + "deltaR", "F")
                self.out.branch(prefix + "mass", "F")

        self.out.branch("topcand_mass", "F")

    def loadGenHistory(self, event):
        # gen matching
        if not self.isMC:
            return

        try:
            genparts = event.genparts
        except RuntimeError as e:
            genparts = Collection(event, "GenPart")
            for idx, gp in enumerate(genparts):
                if 'dauIdx' not in gp.__dict__:
                    gp.dauIdx = []
                if gp.genPartIdxMother >= 0:
                    mom = genparts[gp.genPartIdxMother]
                    if 'dauIdx' not in mom.__dict__:
                        mom.dauIdx = [idx]
                    else:
                        mom.dauIdx.append(idx)
            event.genparts = genparts

        def isHadronic(gp):
            if len(gp.dauIdx) == 0:
                raise ValueError('Particle has no daughters!')
            for idx in gp.dauIdx:
                if abs(genparts[idx].pdgId) < 6:
                    return True
            return False

        def getFinal(gp):
            for idx in gp.dauIdx:
                dau = genparts[idx]
                if dau.pdgId == gp.pdgId:
                    return getFinal(dau)
            return gp

        lepGenTops = []
        hadGenTops = []
        hadGenWs = []
        hadGenZs = []
        hadGenHs = []

        for gp in genparts:
            if gp.statusFlags & (1 << 13) == 0:
                continue
            if abs(gp.pdgId) == 6:
                for idx in gp.dauIdx:
                    dau = genparts[idx]
                    if abs(dau.pdgId) == 24:
                        genW = getFinal(dau)
                        gp.genW = genW
                        if isHadronic(genW):
                            hadGenTops.append(gp)
                        else:
                            lepGenTops.append(gp)
                    elif abs(dau.pdgId) in (1, 3, 5):
                        gp.genB = dau
            elif abs(gp.pdgId) == 24:
                if isHadronic(gp):
                    hadGenWs.append(gp)
            elif abs(gp.pdgId) == 23:
                if isHadronic(gp):
                    hadGenZs.append(gp)
            elif abs(gp.pdgId) == 25:
                if isHadronic(gp):
                    hadGenHs.append(gp)

        for parton in itertools.chain(lepGenTops, hadGenTops):
            parton.daus = (parton.genB, genparts[parton.genW.dauIdx[0]], genparts[parton.genW.dauIdx[1]])
            parton.genW.daus = parton.daus[1:]
        for parton in itertools.chain(hadGenWs, hadGenZs, hadGenHs):
            parton.daus = (genparts[parton.dauIdx[0]], genparts[parton.dauIdx[1]])

        return hadGenTops

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""

        event._allJets = Collection(event, "Jet")
        event.ak4jets = []
        for idx, j in enumerate(event._allJets):
            j.idx = idx
            if not (j.pt > 25 and abs(j.eta) < 2.4 and (j.jetId & 4)):
                # pt, eta, tightIdLepVeto
                continue
            event.ak4jets.append(j)

        if len(event.ak4jets) < 3:
            return False

        # NOTE: use only up to 10 jets
        event.ak4jets = event.ak4jets[:10]

        hadGenTops = self.loadGenHistory(event)

        sigCands, bkgCands = [], []
        for idx in range(len(event.ak4jets)):
            for jdx in range(idx + 1, len(event.ak4jets)):
                for kdx in range(jdx + 1, len(event.ak4jets)):
                    topcand = TopCand(event.ak4jets[idx], event.ak4jets[jdx], event.ak4jets[kdx])
                    if abs(topcand.candP4.M() - 175) > 80:
                        continue
                    topcand.doGenMatch(hadGenTops)
                    if topcand.genMatch == 'bqq':
                        sigCands.append(topcand)
                    else:
                        bkgCands.append(topcand)
        # keep less bkg top candidates (up to 3x of signal candidates)
        np.random.shuffle(bkgCands)
        bkgCands = bkgCands[:3 * len(sigCands)]

        for topcand in itertools.chain(sigCands, bkgCands):
            self.out.fillBranch("event", event.event)

            self.out.fillBranch("isHadTop", topcand.genMatch == 'bqq')
            self.out.fillBranch("genMatch", topcand.genMatchCode)
            if topcand.genMatch:
                self.out.fillBranch("gentop_pt", topcand.matchedGenTop.pt)
                self.out.fillBranch("gentop_eta", topcand.matchedGenTop.eta)
            else:
                self.out.fillBranch("gentop_pt", -1)
                self.out.fillBranch("gentop_eta", -9)

            topcand.makeFeatures()
            for k in topcand.features:
                self.out.fillBranch(k, topcand.features[k])

            self.out.fill()

        return False


# define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed
def resTopTree(): return ResTopTreeProducer()
