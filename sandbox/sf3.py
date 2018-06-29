#!/usr/bin/env python3
import sys, warnings
import numpy as np
import waveLibs as wl
import dsi
import matplotlib.pyplot as plt
plt.style.use('../pltReports.mplstyle')

sys.argv.append("-b")
import ROOT
from ROOT import RooFit as RF
# from ROOT import gROOT
# gROOT.ProcessLine("gErrorIgnoreLevel = 3001;") # suppress ROOT messages
# gROOT.ProcessLine("RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);")
# gROOT.ProcessLine("RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);")
# ROOT.Math.MinimizerOptions.SetDefaultMinimizer("Minuit") # the PLC doesn't work w/ minuit2

eLo, eHi, epb = 0, 50, 0.1
nB = int((eHi-eLo)/epb)

dsList = [0,1,2,3,4,"5A","5B","5C"]

enr = False

def main(argv):

    getUnscaledPDFs(makePlots=False)
    runFit()


def getUnscaledPDFs(ma=0, makePlots=False):
    """ Generate a set of TH1D's to be turned into RooDataHist objects.
    Be careful they have the same axis limits and binning as the RooDataSet.
    Takes axion mass (in keV) as a parameter.
    """
    from ROOT import TFile, TH1D, gROOT

    # output files
    rOut = "../data/specPDFs.root"
    tf = TFile(rOut,"RECREATE")
    td = gROOT.CurrentDirectory()

    # print("Generating unscaled PDFs, eLo %.1f  eHi %.1f  epb %.2f: %s" % (eLo, eHi, epb, rOut))


    # === 1. axion flux

    # axion flux scale.
    # NOTE: to do the fit and set a new limit, we set g_ae=1.
    # To plot an expected flux, we would use a real value.
    # Redondo's note: I calculated the flux using gae = 0.511*10^-10
    # for other values of gae use: FLUX = Table*[gae/(0.511*10^-10)]^2
    gae = 1
    gRat = (gae / 5.11e-11)
    redondoScale = 1e19 * gRat**2 # convert table to [flux / (keV cm^2 d)]

    axData = []
    with open("../data/redondoFlux.txt") as f1: # 23577 entries
        lines = f1.readlines()[11:]
        for line in lines:
            data = line.split()
            axData.append([float(data[0]),float(data[1])])
    axData = np.array(axData)

    def sig_ae(E,m):
        """ E, m are in units of keV.  must multiply result by sig_pe """
        beta = (1 - m**2./E**2.)**(1./2)
        return (1 - (1./3.)*beta**(2./3.)) * (3. * E**2.) / (16. * np.pi * (1./137.) * 511.**2. * beta)

    # === 2. ge photoelectric xs
    phoData = []
    with open("../data/ge76peXS.txt") as f2: # 2499 entries, 0.01 kev intervals
        lines = f2.readlines()
        for line in lines:
            data = line.split()
            phoData.append([float(data[0]),float(data[1])])
    phoData = np.array(phoData)

    # === 3. tritium
    tritData = []
    with open("../data/TritiumSpectrum.txt") as f3: # 20000 entries
        lines = f3.readlines()[1:]
        for line in lines:
            data = line.split()
            conv = float(data[2]) # raw spectrum convolved w/ ge cross section
            if conv < 0: conv = 0.
            tritData.append([float(data[1]),conv])
    tritData = np.array(tritData)

    # NOTE: check sandbox/th1.py for examples of manually filling TH1D's and verifying wl.GetHisto and wl.npTH1D.

    # ROOT output
    h1 = TH1D("h1","photoelectric",nB,eLo,eHi)         # [cm^2 / kg]
    h2 = TH1D("h2","axioelectric",nB,eLo,eHi)          # [cm^2 / kg]
    h3 = TH1D("h3","axion flux, gae=1",nB,eLo,eHi)     # [cts / (keV cm^2 d)]
    h4 = TH1D("h4","convolved flux",nB,eLo,eHi)        # [cts / (keV d kg)]
    h5 = TH1D("h5","tritium",nB,eLo,eHi)               # [cts] (normalized to 1)

    # manually fill ROOT histos (don't normalize yet)
    for iB in range(nB+1):
        ctr = (iB + 0.5)*epb + eLo
        bLo, bHi = ctr - epb/2, ctr + epb/2
        with warnings.catch_warnings():
            warnings.simplefilter("ignore",category=RuntimeWarning)

            # if ma>0, we ignore entries with E <= m.

            # photoelectric x-section [cm^2 / kg]
            idx = np.where((phoData[:,0] >= bLo) & (phoData[:,0] < bHi))
            pho = np.mean(phoData[idx][:,1]) * 1000
            if np.isnan(pho) or len(phoData[idx][:,1]) == 0: pho = 0.
            if phoData[idx][:,1].any() <= ma: pho = 0.
            h1.SetBinContent(iB+1,pho)

            # axioelectric x-section [cm^2 / kg]
            if ctr > ma: axio = pho * sig_ae(ctr, ma)
            else: axio=0.
            h2.SetBinContent(iB+1,axio)

            # axion flux [flux / (cm^2 d keV)]
            idx = np.where((axData[:,0] >= bLo) & (axData[:,0] < bHi))
            flux = np.mean(axData[idx][:,1]) * redondoScale
            if np.isnan(flux): flux = 0.
            h3.SetBinContent(iB+1, flux)
            # YES, adding 1 here. keeps the 6.6 keV line in the proper place for all binnings.
            # it must have to do w/ the way i'm reading in the data from the text files ...

            # axion flux PDF [flux / (keV d kg)]
            axConv = axio * flux
            h4.SetBinContent(iB+1, axConv)

            # tritium
            idx = np.where((tritData[:,0] >= bLo) & (tritData[:,0] <= bHi))
            trit = np.mean(tritData[idx][:,1])
            if np.isnan(trit): trit = 0.
            h5.SetBinContent(iB+1, trit)

    # Pb210 (from separate file)
    tf2 = TFile("../data/Pb210PDFs.root")
    h6 = tf2.Get("hPb210TDL") # with TDL
    h7 = tf2.Get("hPb210") # without TDL
    h6.SetName("h6")
    h7.SetName("h7")

    if makePlots:

        # === 1. verify the numpy histogram and ROOT histogram give the same output. OK

        x, h210, xpb = wl.npTH1D(h7)
        iE = np.where((x > 45) & (x < 48))
        plt.plot(x[iE], h210[iE], ls='steps', lw=3, c='b')
        plt.xlabel("Energy (keV)", ha='right', x=1)
        plt.tight_layout()
        plt.savefig("./plots/sf-pk210.pdf")

        from ROOT import TCanvas
        c = TCanvas()
        h7.GetXaxis().SetTitle("Energy (keV)")
        h7.GetXaxis().SetRangeUser(45, 48)
        h7.Draw('hist')
        c.Print('./plots/sf-pb210th1d.pdf')

        # === 2. print ROOT histos to match w/ numpy histos

        c.Clear(); h1.Draw("hist"); c.Print("./plots/root-sigGe.pdf")
        c.Clear(); h2.Draw("hist"); c.Print("./plots/root-sigAe.pdf")
        c.Clear(); h3.Draw("hist"); c.Print("./plots/root-axFlux.pdf")
        c.Clear(); h4.Draw("hist"); c.Print("./plots/root-axPDF.pdf")
        c.Clear(); h5.Draw("hist"); c.Print("./plots/root-trit.pdf")
        c.Clear(); h6.Draw("hist"); c.Print("./plots/root-pb210TDL.pdf")
        c.Clear(); h7.Draw("hist"); c.Print("./plots/root-pb210.pdf")

    gROOT.cd(td.GetPath())
    h1.Write()
    h2.Write()
    h3.Write()
    h4.Write()
    h5.Write()
    h6.Write()
    h7.Write()
    tf.Close()


def runFit():
    from ROOT import TFile, TH1D

    tf = TFile("../data/latDS%s.root" % ''.join([str(d) for d in dsList]))
    tt = tf.Get("skimTree")
    tCut = "isEnr==1" if enr is True else "isEnr==0"
    hitE = ROOT.RooRealVar("trapENFCal", "Energy", eLo, eHi, "keV")
    hEnr = ROOT.RooRealVar("isEnr", "isEnr", 0, 1, "")
    # hitW = ROOT.RooRealVar("weight", "weight", 1, 1000, "")
    fData = ROOT.RooDataSet("data", "data", tt, ROOT.RooArgSet(hitE, hEnr), tCut)
    # # fData = ROOT.RooDataSet("data", "data", tt, ROOT.RooArgSet(hitE, hEnr, hitW), "", "weight")
    nData = fData.numEntries()

    tf2 = TFile("../data/specPDFs.root")

    trVal = 1000

    trTH1D = tf2.Get("h5")
    trNum = ROOT.RooRealVar("amp-trit", "amp-trit", trVal)
    trDataHist = ROOT.RooDataHist("tr", "tr", ROOT.RooArgList(hitE), RF.Import(trTH1D))
    hitE.setRange(eLo, eHi)
    trPdf = ROOT.RooHistPdf("trPdf", "trPdf", ROOT.RooArgSet(hitE), trDataHist, 0)
    trExt = ROOT.RooExtendPdf("ext-trit", "ext-trit", trPdf, trNum)

    # rooplot before fitting
    fSpec = hitE.frame(RF.Range(eLo,eHi), RF.Bins(int((eHi-eLo)/epb)))

    # wouter's note: DON'T DELETE
    # "the default behavior is when you plot a p.d.f. on an empty frame it is
    # plotted with unit normalization. When you plot it on a frame with data in
    # it, it will normalize to the number of events in that dataset."
    fData.plotOn(fSpec)
    # trExt.plotOn(fSpec) # default normalization: 1 if no data, nData if data
    # trExt.plotOn(fSpec, RF.Normalization(1, ROOT.RooAbsReal.Relative)) # this is relative to total counts, don't use it
    # trExt.plotOn(fSpec, RF.LineColor(ROOT.kBlue), RF.Normalization(trVal/epb, ROOT.RooAbsReal.NumEvent)) # equivalent to the Raw way

    # use this one (you have to divide by epb in numpy when you plot, but not when you integrate)
    trExt.plotOn(fSpec, RF.LineColor(ROOT.kMagenta), RF.Normalization(trVal, ROOT.RooAbsReal.Raw))

    from ROOT import TCanvas
    c = TCanvas("c","c", 1400, 1000)
    fSpec.SetTitle("")
    fSpec.Draw()
    c.Print("../plots/spectrum-before.pdf")


    # replicate the rooplot with numpy (no weights)

    tCut = "isEnr" if enr else "!isEnr"
    tCut += " && trapENFCal >= %.1f && trapENFCal <= %.1f" % (eLo, eHi)
    n = tt.Draw("trapENFCal", tCut, "goff")
    trapE = tt.GetV1()
    trapE = [trapE[i] for i in range(n)]
    x, hData = wl.GetHisto(trapE, eLo, eHi, epb)
    # plt.plot(x, hData, ls='steps', c='b') # normal histo
    hErr = np.asarray([np.sqrt(h) for h in hData]) # statistical error
    plt.errorbar(x, hData, yerr=hErr, c='k', ms=10, linewidth=0.8, fmt='.', capsize=2, zorder=1) # pretty convincing rooplot fake

    x, y, xpb = wl.npTH1D(trTH1D)
    # yn = np.divide(y, np.sum(y))
    # plt.step(x, np.cumsum(yn))

    # yn = np.divide(y, np.sum(y[np.where((x >= eLo) & (x <= eHi))] * xpb)) < -- this makes the norm != 1
    # yn = np.divide(y, np.sum(y[np.where((x >= eLo) & (x <= eHi))])) # < -- this one works w/ no data
    # yn = np.multiply(yn, nData) # <-- this one works w/ data
    # plt.plot(x, yn, ls='steps', c='b', lw=2, label="tritium, nData %d  sum %d  max %.1f" % (nData, int(np.sum(yn)), np.amax(yn)))

    # NOTE: to plot to a set number of counts, divide by epb.  to integrate, don't.
    plt.plot(x1, y1 * trVal / epb, ls='steps', c='m', lw=2, label="trit init value: %d  int: %d" % (trVal, np.sum(y1*trVal)))

    plt.xlabel("Energy (keV)", ha='right', x=1)
    plt.ylabel("Counts / %.1f keV" % epb, ha='right', y=1)
    plt.legend(loc=1)
    plt.tight_layout()
    # plt.show()
    plt.savefig("../plots/sf3-mplplot.pdf")



if __name__=="__main__":
    main(sys.argv[1:])