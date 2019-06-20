import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
import numpy as np


class N2DDTHelper:

    def __init__(self, mapfile, mapname='Rho2D'):
        f_h2ddt = ROOT.TFile.Open(mapfile)
        self._trans_h2ddt = f_h2ddt.Get("Rho2D")
        self._trans_h2ddt.SetDirectory(0)
        f_h2ddt.Close()

    def transform(self, n2, pt, msd):
        # https://github.com/DAZSLE/ZPrimePlusJet/blob/master/analysis/sampleContainer.py#L939-L947
        if not n2 or not pt or not msd:
            return 99
        rho = 2 * np.log(max(msd, 0.01) / pt)
        cur_rho_index = np.clip(self._trans_h2ddt.GetXaxis().FindFixBin(rho), 1, self._trans_h2ddt.GetXaxis().GetNbins())
        cur_pt_index = np.clip(self._trans_h2ddt.GetYaxis().FindFixBin(pt), 1, self._trans_h2ddt.GetYaxis().GetNbins())
        return n2 - self._trans_h2ddt.GetBinContent(cur_rho_index, cur_pt_index)
