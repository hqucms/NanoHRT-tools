#!/usr/bin/env python

from __future__ import print_function
import ROOT
import argparse


class CovMatrix():

    def __init__(self, args):
        self.args = args

    def run(self):
        POIs = self.args.POIs.split(',')
        f_in = self.args.input.split(':')
        f = ROOT.TFile(f_in[0])
        fitres = f.Get(f_in[1])
        fitres_cov = ROOT.TMatrixDSym(len(POIs))
        fitres_cov_src = fitres.covarianceMatrix()
        fitres_cor = ROOT.TMatrixDSym(len(POIs))
        fitres_cor_src = fitres.correlationMatrix()
        ipos = []
        for p in POIs:
            ipos.append(fitres.floatParsFinal().index(p))
        for i, ip in enumerate(POIs):
            for j, jp in enumerate(POIs):
                fitres_cor[i][j] = ROOT.Double_t(
                    fitres_cor_src[ipos[i]][ipos[j]])
                fitres_cov[i][j] = ROOT.Double_t(
                    fitres_cov_src[ipos[i]][ipos[j]])
        print('RooFitResult correlation matrix:')
        fitres_cor.Print()
        print('RooFitResult covariance matrix:')
        fitres_cov.Print()
        if self.args.output is not None:
            out = self.args.output.split(':')
            fout = ROOT.TFile(out[0], 'RECREATE')
            prefix = out[1]
            fout.WriteTObject(fitres_cor, prefix + '_comp_cor')
            h_cor_compare = self.fix_TH2(ROOT.TH2D(fitres_cor), POIs)
            fout.WriteTObject(h_cor_compare, prefix + '_comp_h_cor')
            ROOT.gStyle.SetOptStat(0)
            ROOT.gStyle.SetPaintTextFormat(".2f")
            c1 = ROOT.TCanvas("c1", "c1", 800, 800)
            h_cor_compare.SetTitle('Correlation')
            h_cor_compare.GetZaxis().SetRangeUser(-1, 1)
            h_cor_compare.Draw("colztext")
            c1.Print(out[0].replace('.root', '') + '_' + prefix + '_comp_cor.pdf')
            fout.WriteTObject(fitres_cov, prefix + '_comp_cov')
            h_cov_compare = self.fix_TH2(ROOT.TH2D(fitres_cov), POIs)
            fout.WriteTObject(h_cov_compare, prefix + '_comp_h_cov')
            c2 = ROOT.TCanvas("c2", "c2", 800, 800)
            h_cov_compare.SetTitle('Covariance')
            h_cov_compare.Draw("colztext")
            c2.Print(out[0].replace('.root', '') + '_' + prefix + '_comp_cov.pdf')

    def fix_TH2(self, h, labels):
        h_fix = h.Clone()
        for y in range(1, h.GetNbinsY() + 1):
            for x in range(1, h.GetNbinsX() + 1):
                h_fix.SetBinContent(
                    x, y, h.GetBinContent(x, h.GetNbinsY() + 1 - y))
        for x in range(1, h_fix.GetNbinsX() + 1):
            h_fix.GetXaxis().SetBinLabel(x, labels[x - 1])
        for y in range(1, h_fix.GetNbinsY() + 1):
            h_fix.GetYaxis().SetBinLabel(y, labels[-y])
        return h_fix


if __name__ == '__main__':
    parser = argparse.ArgumentParser('Print norm SFs')
    parser.add_argument(
        '-i', '--input', help='The input file containing RooFitResult, in the format file:prefix')
    parser.add_argument(
        '-o', '--output', help='The output name in the format file:prefix')
    parser.add_argument(
        '-P', '--POIs', help='The params that were scanned (in scan order)')
    args = parser.parse_args()

    c = CovMatrix(args)
    c.run()
