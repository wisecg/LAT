#!/usr/bin/env python
""" DEPRECATED:
    We've been asked by Steve & others not to "inflate" the real data
    but to account for efficiencies in the spectrum fits.
"""

def main():
    getEff

def getEff():
    """ ./job-panda.py -getEff

    METHOD:
    open up the latskim file for each channel.
    loop over the good run ranges.
    for each good range, make an energy histogram.
    then calculate the efficiency curve based on the sigma value
    and convolve it with the histogram points.
    """
    import numpy as np
    import waveLibs as wl
    import scipy.special as spec
    import matplotlib.pyplot as plt
    from ROOT import TFile, TTree, TH1D, TF1, TCanvas, gROOT
    import ROOT, random

    gROOT.ProcessLine(".x ~/env/MJDClintPlotStyle.C")
    # gROOT.ProcessLine("gErrorIgnoreLevel = 3001;") # suppress ROOT messages

    bins, xlo, xhi = 50,0,15  # set it just high enough that the first bin center isn't negative

    hSumCorr = TH1D("hSumCorr","hSumCorr",bins,xlo,xhi)
    hSumUncr = TH1D("hSumUncr","hSumUncr",bins,xlo,xhi)

    dsNum = 1
    # ch = 578
    for ch in ds.GetGoodChanList(dsNum):

        inFile = TFile(homePath+"/project/latskim/latSkimDS%d_ch%d.root" % (dsNum,ch))
        tree = inFile.Get("skimTree")
        fileCut = inFile.Get("theCut").GetTitle()

        _,_,goodRunErfs = ds.GetThreshDicts(dsNum)

        hUnc = wl.H1D(tree,bins,xlo,xhi,"trapENFCal",fileCut)
        hSumUncr.Add(hUnc)

        for erfs in goodRunErfs[ch]:

            runCut = " && run >= %d && run <= %d" % (erfs[0],erfs[1])
            theCut = fileCut + runCut

            h1 = wl.H1D(tree,bins,xlo,xhi,"trapENFCal",theCut)
            h1x,h1y = wl.npTH1D(h1)

            # calculate efficiency curve
            thisErf = TF1("thisErf","0.5*(1+TMath::Erf((x-[0])/(TMath::Sqrt(2)*[1]) ))")
            thisErf.SetParameter(0,erfs[2]) # mu
            thisErf.SetParameter(1,erfs[3]) # sigma
            h1.Divide(thisErf)

            # thisErf = 0.5 * (1 + spec.erf( (h1x - mu) / (np.sqrt(2) * sig) ))
            # h1yScaled = h1y / thisErf
            # nameStr = str(random.uniform(1.,2.))
            # h2 = TH1D(nameStr,nameStr,bins,xlo,xhi)
            # for i in range(bins):
            #     h2.SetBinContent(i,h1yScaled[i])

            hSumCorr.Add(h1)


    # eff-corrected spectrum.
    c = TCanvas("c","c",800,600)
    c.SetLogy(1)

    hSumCorr.SetLineColor(ROOT.kBlue)
    hSumCorr.Draw("hist")

    hSumUncr.SetLineColor(ROOT.kRed)
    hSumUncr.Draw("hist same")

    # l1 = TLegend

    c.Print("./plots/effWeight/eff_DS%d.pdf" % dsNum)


if __name__=="__main__":
    main()