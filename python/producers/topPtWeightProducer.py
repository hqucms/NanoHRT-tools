import ROOT
import math
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

from ..helpers.utils import clip


class TopPtWeightProducer(Module):

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.isMC = bool(inputTree.GetBranch('genWeight'))
        if self.isMC:
            self.out = wrappedOutputTree

            self.out.branch('topptWeight', "F")
            # accumulation hists
            self.h_genwgt = ROOT.TH1D('genweight', 'genweight', 1, 0, 1)
            self.h_ttbar_topptwgt = ROOT.TH1D('TT_toppt_weight', 'TT_toppt_weight', 1, 0, 1)

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        if self.isMC:
            cwd = ROOT.gDirectory
            outputFile.cd()
            self.h_genwgt.Write()
            self.h_ttbar_topptwgt.Write()
            cwd.cd()

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""

        if not self.isMC:
            return True

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

        genTops = []
        for gp in genparts:
            if gp.statusFlags & (1 << 13) == 0:
                # 13: isLastCopy
                continue
            if abs(gp.pdgId) == 6:
                genTops.append(gp)

        topptWeight = 1.

        if len(genTops) == 2:
            # ttbar (+X ?)

            def wgt(pt):
                return math.exp(0.0615 - 0.0005 * clip(pt, 0, 800))

            def wgt_nnlo(pt):
                x = clip(pt, 0, 2000)
                return 0.103 * math.exp(-0.0118 * x) - 0.000134 * x + 0.973

            # topptWeight = np.sqrt(wgt(genTops[0].pt) * wgt(genTops[1].pt))
            topptWeight = math.sqrt(wgt_nnlo(genTops[0].pt) * wgt_nnlo(genTops[1].pt))

        self.out.fillBranch('topptWeight', topptWeight)
        # accumulation hists
        self.h_genwgt.Fill(0.5, event.genWeight)
        self.h_ttbar_topptwgt.Fill(0.5, topptWeight * event.genWeight)

        return True


# define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed
def topPtWeight():
    return TopPtWeightProducer()
