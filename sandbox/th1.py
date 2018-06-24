#!/usr/bin/env python3
import matplotlib.pyplot as plt
import numpy as np
import ROOT
from ROOT import TH1D, TCanvas
import waveLibs as wl

def main():
    """ Make sure that the numpy histograms I'm using match the ROOT histogram class output.
    Examples:
        - Matching wl.GetHisto to a TH1D
        - Matching wl.npTH1D to a TH1D
        - Filling the bins of a TH1D manually
    """
    # start with 100 random values, from 0 to 10
    r = 10 * np.random.rand(100)

    # bin from 0 to 10 in intervals of 1 (11 bins.  nB = 10 b/c it is 0-indexed)
    xLo, xHi, xpb = 0, 10, 1
    nB = int((xHi-xLo)/xpb)

    # make the TH1D purely w/ ROOT stuff
    c = TCanvas()
    th1 = TH1D("th","",nB,xLo,xHi)
    for v in r: th1.Fill(v)

    # 1. getHisto and npTH1D
    x1, y1 = wl.GetHisto(r, xLo, xHi, xpb, shift=False)
    plt.plot(x1, y1, ls='steps', c='b', lw=6)
    plt.axvline(1, c='g', lw=1)

    x2, y2, xpb2 = wl.npTH1D(th1)
    plt.plot(x2, y2, ls='steps', c='r', lw=4)

    # 2. manually fill each histogram bin to match some known curve (used for pdfs)
    th2 = TH1D("t2","t2",nB,xLo,xHi)

    # the +1 is very important, it covers the last bin
    for iB in range(nB+1):
        ctr = (iB + 0.5)*xpb + xLo
        bLo, bHi = ctr-xpb/2, ctr+xpb/2
        idx = np.where((x2 >= bLo) & (x2 < bHi))
        hCts = np.sum(y2[idx])
        th2.SetBinContent(iB, hCts)

    th1.SetLineColor(ROOT.kBlue)
    th2.SetLineColor(ROOT.kRed)
    th1.SetLineWidth(2)
    th1.Draw("hist")
    th2.Draw("hist same")
    c.Print("../plots/htmp.pdf")

    x3, y3, xpb3 = wl.npTH1D(th2)
    plt.plot(x3, y3, ls='steps', c='g',lw=2)
    plt.xlim(xLo, xHi)

    plt.tight_layout()
    plt.savefig("../plots/hnp.pdf")



if __name__=="__main__":
    main()