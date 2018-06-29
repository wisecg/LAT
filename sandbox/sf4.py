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

eLo, eHi, epb = 5, 50, 0.3
pLo, pHi, ppb = 0, 30, 0.05
nB = int((eHi-eLo)/epb)
nBP = int((pHi-pLo)/ppb)

# dsList = [0,1,2,3,4,"5A","5B","5C"]
dsList = [1,2,3,4,"5B","5C"]

enr = True
eff = True

def main(argv):

    initialize()
    # loadDataMJD()
    # getUnscaledPDFs(makePlots=True)
    plotModel()
    plotFit()


def initialize():
    """ For this dsList, load the efficiency curves and exposure into globals """
    global effLim, effMax, xEff, detEff

    # load efficiency correction
    f1 = np.load('%s/data/lat-expo-efficiency-all-e95.npz' % dsi.latSWDir)
    xEff = f1['arr_0']
    totEnrEff, totNatEff = f1['arr_1'].item(), f1['arr_2'].item()
    detEff = np.zeros(len(xEff))
    for ds in dsList:
        if enr: detEff += totEnrEff[ds]
        else: detEff += totNatEff[ds]

    # load exposure
    f2 = np.load("%s/data/expo-totals-e95.npz"  % dsi.latSWDir)
    dsExpo, detExpo = f2['arr_0'].item(), f2['arr_1'].item()
    detExp = 0
    for d in dsExpo:
        if d in dsList:
            if enr: detExp += dsExpo[d][0]
            else:   detExp += dsExpo[d][1]

    # normalize the efficiency
    detEff = np.divide(detEff, detExp)
    effLim, effMax = xEff[-1], detEff[-1]

    # diagnostic plot
    # plt.axhline(np.amax(detEff), c='orange', label="%.2f" % np.amax(detEff))
    # plt.plot(xEff, detEff)
    # plt.legend(loc=4)
    # plt.xlabel("Energy (keV)", ha='right', x=1)
    # plt.tight_layout()
    # plt.show()
    # exit()


def loadDataMJD():
    """ Load MJD data based on the global variable 'dsList'.
    RooFit can't handle the vector<double> format for energies.
    So save a few select branches into a new file.
    """
    from array import array
    from ROOT import TChain, TFile, TTree

    # load the data
    tt = TChain("skimTree")
    for ds in dsList:
        tt.Add("%s/final95t/final95t_DS%s.root" % (dsi.cutDir, ds))

    # declare output
    fName = "%s/data/latDS%s.root" % (dsi.latSWDir, ''.join([str(d) for d in dsList]))
    fOut = TFile(fName,"RECREATE")
    tOut = TTree("skimTree", "skimTree")
    run = array('i',[0])
    iEvt = array('i',[0])
    iHit = array('i',[0])
    chan = array('i',[0])
    hitE = array('d',[0.])
    isEnr = array('i',[0])
    weight = array('d',[0.])
    tOut.Branch("run", run, "run/I")
    tOut.Branch("iEvent", iEvt, "iEvent/I")
    tOut.Branch("iHit", iHit, "iHit/I")
    tOut.Branch("channel", chan, "channel/I")
    tOut.Branch("trapENFCal", hitE, "trapENFCal/D")
    tOut.Branch("isEnr", isEnr, "isEnr/I")
    tOut.Branch("weight", weight, "weight/D")

    for iE in range(tt.GetEntries()):
        tt.GetEntry(iE)
        run[0] = tt.run
        iEvt[0] = tt.iEvent
        for iH in range(tt.channel.size()):
            iHit[0] = tt.iHit.at(iH)
            chan[0] = tt.channel.at(iH)
            hitE[0] = tt.trapENFCal.at(iH)

            # calculate weight based on 1/efficiency
            if hitE[0] > effLim:
                weight[0] = 1/effMax
            else:
                idx = (np.abs(xEff-hitE[0])).argmin()
                weight[0] = 1/np.interp(hitE[0], xEff[idx:idx+1], detEff[idx:idx+1])
            # if hitE[0] < effLim:
                # print("%.2f  %.2f " % (hitE[0], weight[0]))

            if "P" in tt.detName.at(iH): isEnr[0] = 1
            elif "B" in tt.detName.at(iH): isEnr[0] = 0
            else:
                print("WTF, error")
                exit(0)
            tOut.Fill()

    tOut.Write()
    fOut.Close()

    # verify
    f2 = TFile(fName)
    t2 = f2.Get("skimTree")
    t2.Scan("run:channel:isEnr:trapENFCal:weight")


def getUnscaledPDFs(ma=0, makePlots=False):
    """ Generate a set of TH1D's to be turned into RooDataHist objects.
    Be careful they have the same axis limits and binning as the RooDataSet.
    Takes axion mass (in keV) as a parameter.
    """
    from ROOT import TFile, TH1D, gROOT

    # output files
    rOut = "%s/data/specPDFs.root" % dsi.latSWDir
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
    with open("%s/data/redondoFlux.txt" % dsi.latSWDir) as f1: # 23577 entries
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
    with open("%s/data/ge76peXS.txt" % dsi.latSWDir) as f2: # 2499 entries, 0.01 kev intervals
        lines = f2.readlines()
        for line in lines:
            data = line.split()
            phoData.append([float(data[0]),float(data[1])])
    phoData = np.array(phoData)

    # === 3. tritium
    tritData = []
    with open("%s/data/TritiumSpectrum.txt" % dsi.latSWDir) as f3: # 20000 entries
        lines = f3.readlines()[1:]
        for line in lines:
            data = line.split()
            conv = float(data[2]) # raw spectrum convolved w/ ge cross section
            if conv < 0: conv = 0.
            tritData.append([float(data[1]),conv])
    tritData = np.array(tritData)

    # NOTE: check sandbox/th1.py for examples of manually filling TH1D's and verifying wl.GetHisto and wl.npTH1D.

    # ROOT output
    h1 = TH1D("h1","photoelectric",nBP,pLo,pHi)         # [cm^2 / kg]
    h2 = TH1D("h2","axioelectric",nBP,pLo,pHi)          # [cm^2 / kg]
    h3 = TH1D("h3","axion flux, gae=1",nBP,pLo,pHi)     # [cts / (keV cm^2 d)]
    h4 = TH1D("h4","convolved flux",nBP,pLo,pHi)        # [cts / (keV d kg)]
    h5 = TH1D("h5","tritium",nBP,pLo,pHi)               # [cts] (normalized to 1)

    # manually fill ROOT histos (don't normalize yet)
    for iB in range(nBP+1):
        ctr = (iB + 0.5)*ppb + pLo
        bLo, bHi = ctr - ppb/2, ctr + ppb/2
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
    tf2 = TFile("%s/data/Pb210PDFs.root" % dsi.latSWDir)
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
        plt.savefig("%s/plots/sf-pk210.pdf" % dsi.latSWDir)

        from ROOT import TCanvas
        c = TCanvas()
        h7.GetXaxis().SetTitle("Energy (keV)")
        h7.GetXaxis().SetRangeUser(45, 48)
        h7.Draw('hist')
        c.Print('%s/plots/sf-pb210th1d.pdf' % dsi.latSWDir)

        # === 2. print ROOT histos to match w/ numpy histos

        c.Clear(); h1.Draw("hist"); c.Print("%s/plots/root-sigGe.pdf" % dsi.latSWDir)
        c.Clear(); h2.Draw("hist"); c.Print("%s/plots/root-sigAe.pdf" % dsi.latSWDir)
        c.Clear(); h3.Draw("hist"); c.Print("%s/plots/root-axFlux.pdf" % dsi.latSWDir)
        c.Clear(); h4.Draw("hist"); c.Print("%s/plots/root-axPDF.pdf" % dsi.latSWDir)
        c.Clear(); h5.Draw("hist"); c.Print("%s/plots/root-trit.pdf" % dsi.latSWDir)
        c.Clear(); h6.Draw("hist"); c.Print("%s/plots/root-pb210TDL.pdf" % dsi.latSWDir)
        c.Clear(); h7.Draw("hist"); c.Print("%s/plots/root-pb210.pdf" % dsi.latSWDir)

    gROOT.cd(td.GetPath())
    h1.Write()
    h2.Write()
    h3.Write()
    h4.Write()
    h5.Write()
    h6.Write()
    h7.Write()
    tf.Close()


def getBkgPDF(eff=False):
    from ROOT import TH1D
    bkg = TH1D("bkg","flat BG",nB,eLo,eHi)
    for iB in range(nBP+1):
        bkg.SetBinContent(iB, 1) # the initial amplitude doesn't matter b/c we normalize to 1 (no width option)
    bkg.Scale(1 / bkg.Integral(bkg.FindBin(eLo), bkg.FindBin(eHi)))

    if eff:
        # efficiency correct doesn't renormalize to 1
        bkg = getEffCorrTH1D(bkg, eLo, eHi, nB)

    return bkg


def peakPDF(pkE, sig, name, eff=False):
    """ Make a TH1D of the peaks so we can apply the efficiency correction
    the same way we do for the continuous pdfs.
    """
    from ROOT import TH1D

    def gaus(x, a, x0, sigma):
        return a * np.exp(-(x - x0)**2 / (2 * sigma**2))

    xLo = pkE-2 if pkE-2 > eLo else eLo
    xHi = pkE+2 if pkE+2 < eHi else eHi

    # normalize based on global eLo, eHi, epb (not necessarily 1)
    x = np.arange(xLo, xHi, ppb)
    hP = gaus(x, 1, pkE, sig)
    hP = np.divide(hP, np.sum(hP[np.where((x >= eLo) & (x <= eHi))]) * epb)
    # print("peakgen sum", np.sum(hP))

    # make the TH1D
    nB = int((xHi-xLo)/ppb)
    hPR = TH1D(name,name,nB,xLo,xHi)
    for iB in range(nB):
        hPR.SetBinContent(iB, hP[iB])

    # efficiency correct (doesn't normalize)
    if eff:
        hPR = getEffCorrTH1D(hPR, xLo, xHi, nB)

    return hPR


def normPDF(x, y, xLo, xHi):
    """ Normalize a numpy pdf to 1 in the global
    energy range (not the range the pdf was generated) """
    idx = np.where((x>=xLo)&(x<=xHi))
    return x[idx], np.divide(y, np.sum(y[idx]))[idx]


def getEffCorr(x, h, inv=False):
    """ Returns numpy arrays, uses global xEff and detEff functions.
    If inv is False, "corrects" a PDF by weighting it DOWN
    If inv is True, "uncorrects" a PDF by weighting it UP.
    """
    hc = []
    for i in range(len(x)):
        idx = (np.abs(xEff-x[i])).argmin()
        if not inv:
            hc.append(h[i] * np.interp(x[i], xEff[idx:idx+1], detEff[idx:idx+1]))
        else:
            hc.append(h[i] / np.interp(x[i], xEff[idx:idx+1], detEff[idx:idx+1]))
    hc = np.asarray(hc)
    # plt.close()
    # plt.step(x, h)
    # plt.step(x, hc)
    # plt.show()
    # exit()
    return hc


def getEffCorrTH1D(h, xLo, xHi, xNB):
    """ Returns a copy of a ROOT TH1D, uses global xEff and detEff functions.
    Doesn't normalize.
    """
    from ROOT import TH1D
    nB = h.GetNbinsX()
    hEff = TH1D(h.GetName(), h.GetTitle(),xNB, xLo, xHi)
    for i in range(nB+1):
        binE = h.GetXaxis().GetBinCenter(i)
        hBin = h.GetBinContent(i)
        idx = (np.abs(xEff-binE)).argmin()
        hBinC = hBin * np.interp(binE, xEff[idx:idx+1], detEff[idx:idx+1])
        hEff.SetBinContent(i, hBinC)
    return hEff


def getTotalModel(pdfs, eLo, eHi, epb, smooth=True, amp=True):

    xT = np.arange(eLo, eHi, epb)
    yT = np.zeros(len(xT))

    for iB in range(len(xT)-1):

        for xp, yp, pb, v in pdfs:
            idx = np.where((xp>=xT[iB]) & (xp<xT[iB+1]))
            if len(idx[0]) > 0:
                if amp:
                    yT[iB] += v * np.average(yp[idx]) * (epb/pb)
                else:
                    yT[iB] += np.average(yp[idx]) * (epb/pb)

            if iB == len(xT)-2:  # fill the last bin w/ the 2nd to last bin
                yT[iB+1] = yT[iB]

    if smooth:
        xTS = np.arange(eLo, eHi, 0.01)
        yTS = np.interp(xTS, xT, yT)
        return xTS, yTS

    return xT, yT


def getSigma(E, ds=None):
    """ Get the MJ energy resolution.
    This is just used to set an intial guess on the resolution
    because we let sigma float for the peaks.
    """
    if ds==0:   p = [0.147, 0.0173, 0.0003]
    elif ds==1: p = [0.136, 0.0174, 0.00028]
    elif ds==2: p = [0.143, 0.0172, 0.000284]
    elif ds==3: p = [0.162, 0.0172, 0.000297]
    elif ds==4: p = [0.218, 0.015, 0.00035]
    elif ds=='5A': p = [0.2121, 0.01838, 0.00031137]
    elif ds=='5B': p = [0.18148, 0.01690, 0.00031873]
    else:
        # use DS5B numbers for any other DS
        p = [0.18148, 0.01690, 0.00031873]

    return np.sqrt(p[0]**2 + p[1]**2 * E + p[2]**2 * E**2)


def plotModel():
    from ROOT import TFile, TH1D

    tf = TFile("%s/data/latDS%s.root" % (dsi.latSWDir,''.join([str(d) for d in dsList])))
    tt = tf.Get("skimTree")
    tCut = "isEnr==1" if enr is True else "isEnr==0"
    hitE = ROOT.RooRealVar("trapENFCal", "Energy", eLo, eHi, "keV")
    hEnr = ROOT.RooRealVar("isEnr", "isEnr", 0, 1, "")
    # hitW = ROOT.RooRealVar("weight", "weight", 1, 1000, "")
    fData = ROOT.RooDataSet("data", "data", tt, ROOT.RooArgSet(hitE, hEnr), tCut)
    # fData = ROOT.RooDataSet("data", "data", tt, ROOT.RooArgSet(hitE, hEnr, hitW), "", "weight")
    nData = fData.numEntries()
    fitWorkspace = ROOT.RooWorkspace("fitWorkspace","Fit Workspace")
    getattr(fitWorkspace,'import')(hitE)
    getattr(fitWorkspace,'import')(fData)
    # getattr(fitWorkspace,'import')(fWeight)

    tf2 = TFile("%s/data/specPDFs.root" % dsi.latSWDir)
    pdfList = ROOT.RooArgList("shapes")

    # tritium
    nTr = 1000
    hTr = tf2.Get("h5")
    if eff: hTr = getEffCorrTH1D(hTr, pLo, pHi, nBP)
    trNum = ROOT.RooRealVar("amp-trit", "amp-trit", nTr, 0., 50000.)
    trDH = ROOT.RooDataHist("tr", "tr", ROOT.RooArgList(hitE), RF.Import(hTr))
    hitE.setRange(eLo, eHi)
    trPdf = ROOT.RooHistPdf("trPdf", "trPdf", ROOT.RooArgSet(hitE), trDH, 2)
    trExt = ROOT.RooExtendPdf("ext-trit", "ext-trit", trPdf, trNum)
    pdfList.add(trExt)

    # flat bg
    nBk = 1000
    hBkg = getBkgPDF(eff)
    bkgNum = ROOT.RooRealVar("amp-bkg", "amp-bkg", nBk, 0., 10000.)
    bkgDH = ROOT.RooDataHist("bkg", "bkg", ROOT.RooArgList(hitE), RF.Import(hBkg))
    hitE.setRange(eLo, eHi)
    bkgPdf = ROOT.RooHistPdf("bkgPdf", "bkgPdf", ROOT.RooArgSet(hitE), bkgDH, 2)
    bkgExt = ROOT.RooExtendPdf("ext-bkg", "ext-bkg", bkgPdf, bkgNum)
    pdfList.add(bkgExt)

    # 68ge peak
    nPk = 100
    hPk = peakPDF(10.37, getSigma(10.37), "68GeK", eff)
    pkDH = ROOT.RooDataHist("pk", "pk", ROOT.RooArgList(hitE), RF.Import(hPk))
    hitE.setRange(eLo, eHi)
    pkPdf = ROOT.RooHistPdf("pkPdf", "pkPdf", ROOT.RooArgSet(hitE), pkDH, 2)
    pkNum = ROOT.RooRealVar("amp-68GeK", "amp-68GeK", nPk, 0.0, 1000.)
    pkExt = ROOT.RooExtendPdf("ext-68GeK", "ext-68GeK", pkPdf, pkNum)
    pdfList.add(pkExt)

    model = ROOT.RooAddPdf("model","total pdf",pdfList)

    # rooplot before fitting
    fSpec = hitE.frame(RF.Range(eLo,eHi), RF.Bins(int((eHi-eLo)/epb)))

    # wouter's note: DON'T DELETE
    # "the default behavior is when you plot a p.d.f. on an empty frame it is
    # plotted with unit normalization. When you plot it on a frame with data in
    # it, it will normalize to the number of events in that dataset."
    # (then after you do a fit, the pdf normalization changes again ...)
    fData.plotOn(fSpec)
    # bkgExt.plotOn(fSpec)
    # pkExt.plotOn(fSpec)

    # 1 -- individual components at their initial fit values
    # use this one for one component (you have to divide by (bin width of orig pdf) in numpy when you plot, but not when you integrate)
    trExt.plotOn(fSpec, RF.LineColor(ROOT.kMagenta), RF.Normalization(nTr, ROOT.RooAbsReal.Raw))
    bkgExt.plotOn(fSpec, RF.LineColor(ROOT.kGreen), RF.Normalization(nBk, ROOT.RooAbsReal.Raw))
    pkExt.plotOn(fSpec, RF.LineColor(ROOT.kBlue), RF.Normalization(nPk, ROOT.RooAbsReal.Raw))

    # 2 -- the model, normalized according to the total number of counts in a really fking stupid way
    # model.plotOn(fSpec, RF.LineColor(ROOT.kRed))
    # model.plotOn(fSpec, RF.Components("ext-trit"), RF.LineColor(ROOT.kMagenta), RF.Name("ext-trit"))
    # model.plotOn(fSpec, RF.Components("ext-bkg"), RF.LineColor(ROOT.kGreen), RF.Name("ext-bkg"))
    # model.plotOn(fSpec, RF.Components("ext-68GeK"), RF.LineColor(ROOT.kBlue), RF.Name("ext-68GeK"))


    from ROOT import TCanvas
    c = TCanvas("c","c", 1400, 1000)
    fSpec.SetTitle("")
    fSpec.Draw()
    c.Print("%s/plots/spectrum-before.pdf" % dsi.latSWDir)
    c.Clear()

    # === replicate the rooplot with numpy (no weights) ===

    tCut = "isEnr" if enr else "!isEnr"
    tCut += " && trapENFCal >= %.1f && trapENFCal <= %.1f" % (eLo, eHi)
    n = tt.Draw("trapENFCal", tCut, "goff")
    trapE = tt.GetV1()
    trapE = [trapE[i] for i in range(n)]
    x, hData = wl.GetHisto(trapE, eLo, eHi, epb)
    # plt.plot(x, hData, ls='steps', c='b') # normal histo
    hErr = np.asarray([np.sqrt(h) for h in hData]) # statistical error
    plt.errorbar(x, hData, yerr=hErr, c='k', ms=5, linewidth=0.5, fmt='.', capsize=1, zorder=1) # pretty convincing rooplot fake

    # get (eff-corrected) histos and normalize them to the global energy range
    x1, y1, _ = wl.npTH1D(hTr)
    x2, y2, _ = wl.npTH1D(hBkg)
    x3, y3, _ = wl.npTH1D(hPk)
    x1, y1 = normPDF(x1, y1, eLo, eHi)
    x2, y2 = normPDF(x2, y2, eLo, eHi)
    x3, y3 = normPDF(x3, y3, eLo, eHi)

    # === 1. plot individual components of the model
    # *** NOTE: to plot, divide by (bin width when generated).  to integrate, don't. ***
    plt.plot(x1, y1 * nTr / ppb, ls='steps', c='m', lw=2, label="trit init: %d int: %d" % (nTr, np.sum(y1 * nTr)))
    plt.plot(x2, y2 * nBk / epb, ls='steps', c='g', lw=2, label="bkg init: %d int: %d" % (nBk, np.sum(y2 * nBk)))
    plt.plot(x3, y3 * nPk / ppb, ls='steps', c='b', lw=2, label="68GeK init: %d int: %d" % (nPk, np.sum(y3 * nPk)))

    # === 2. replicate the stupid way a rooplot normalizes multiple pdf's based on the number of data counts (before a fit)
    # nModel = 3
    # yTot = np.add(y1, y2, y3)
    # yTot *= nData/nModel
    # plt.plot(x1, yTot, ls='steps', c='r', lw=2, label="total, nData %d  sum %d  max %.1f" % (nData, int(np.sum(yTot)), np.amax(yTot)))
    # plt.plot(x1, y1 * nData/2, ls='steps', c='b', lw=2, label="tritium")
    # plt.plot(x2, y2 * nData/2, ls='steps', c='g', lw=2, label="bkg")

    # === 3. check a peak which was generated with a different binning than the global binning
    # print(np.sum(y3)) # this is 1
    # print(np.sum(y3 * (epb/xpb3))) # this is bigger than 1, but matches the way rooplot normalizes it when plotted by itself.  fk rooplot
    # print(np.sum(y3 * nPk))
    # plt.plot(x3, y3 * nPk / xpb3, ls='steps', label="pk init value: %d  int: %d" % (nPk, np.sum(y3 * nPk)))

    plt.xlabel("Energy (keV)", ha='right', x=1)
    plt.ylabel("Counts / %.1f keV" % epb, ha='right', y=1)
    plt.legend(loc=1)
    plt.xlim(eLo, eHi)
    plt.ylim(ymin=0)
    plt.tight_layout()
    # plt.show()
    plt.savefig("%s/plots/sf4-mplplot.pdf" % dsi.latSWDir)


    # === alright, now run the fit and check the plot again
    minimizer = ROOT.RooMinimizer( model.createNLL(fData, RF.NumCPU(2,0), RF.Extended(True)) )
    minimizer.setPrintLevel(-1)
    minimizer.setStrategy(2)
    minimizer.migrad()
    fitResult = minimizer.save()

    # according to the internet, covQual==3 is a good indicator that it converged
    print("Fitter is done. Fit Cov Qual:", fitResult.covQual())

    # save workspace to a TFile
    getattr(fitWorkspace,'import')(fitResult)
    getattr(fitWorkspace,'import')(model)
    tf3 = TFile("%s/data/fitWorkspace.root" % dsi.latSWDir,"RECREATE")
    fitWorkspace.Write()
    tf3.Close()


def plotFit():
    from ROOT import TFile, TCanvas, TH1D

    f = TFile("%s/data/fitWorkspace.root" % dsi.latSWDir)
    fitWorkspace = f.Get("fitWorkspace")
    fData = fitWorkspace.allData().front()
    fitResult = fitWorkspace.allGenericObjects().front()
    nPars = fitResult.floatParsFinal().getSize()
    hitE = fitWorkspace.var("trapENFCal")
    model = fitWorkspace.pdf("model")
    # fitWorkspace.Print()

    # plot data
    fSpec = hitE.frame(RF.Range(eLo,eHi), RF.Bins(nB))
    fData.plotOn(fSpec)

    # plot model and components
    model.plotOn(fSpec, RF.LineColor(ROOT.kRed), RF.Name("FullModel"))

    c = TCanvas("c","c", 1400, 1000)
    fSpec.SetTitle("")
    fSpec.Draw()
    c.Print("%s/plots/spectrum-after.pdf" % dsi.latSWDir)
    c.Clear()

    # get fit results
    fitVals = {}
    for i in range(nPars):
        fp = fitResult.floatParsFinal()
        name = fp.at(i).GetName()
        fitVal = fp.at(i).getValV()
        fitErr = fp.at(i).getError()
        print("%s  fitVal %.2f  error %.2f" % (name, fitVal, fitErr))
        if name == "amp-68GeK": nPk = fitVal
        if name == "amp-bkg":   nBk = fitVal
        if name == "amp-trit":  nTr = fitVal

    # === duplicate the rooplot in matplotlib ===
    plt.close()

    tf = TFile("%s/data/latDS%s.root" % (dsi.latSWDir,''.join([str(d) for d in dsList])))
    tt = tf.Get("skimTree")
    tCut = "isEnr==1" if enr is True else "isEnr==0"
    tCut = "isEnr" if enr else "!isEnr"
    tCut += " && trapENFCal >= %.1f && trapENFCal <= %.1f" % (eLo, eHi)
    n = tt.Draw("trapENFCal", tCut, "goff")
    trapE = tt.GetV1()
    trapE = [trapE[i] for i in range(n)]
    x, hData = wl.GetHisto(trapE, eLo, eHi, epb)
    # plt.plot(x, hData, ls='steps', c='b') # normal histo
    hErr = np.asarray([np.sqrt(h) for h in hData]) # statistical error
    plt.errorbar(x, hData, yerr=hErr, c='k', ms=5, linewidth=0.5, fmt='.', capsize=1, zorder=1) # pretty convincing rooplot fake

    # plot the model components and total
    tf2 = TFile("%s/data/specPDFs.root" % dsi.latSWDir)

    # get (eff-corrected) histos and normalize to 1 in the global energy range
    hTr = tf2.Get("h5")
    if eff: hTr = getEffCorrTH1D(hTr, pLo, pHi, nBP)
    hBkg = getBkgPDF(eff)
    hPk = peakPDF(10.37, getSigma(10.37), "68GeK", eff)
    x1, y1, xpb1 = wl.npTH1D(hTr)
    x2, y2, xpb2 = wl.npTH1D(hBkg)
    x3, y3, xpb3 = wl.npTH1D(hPk)
    x1, y1 = normPDF(x1, y1, eLo, eHi)
    x2, y2 = normPDF(x2, y2, eLo, eHi)
    x3, y3 = normPDF(x3, y3, eLo, eHi)
    nTot = nTr + nBk + nPk
    if abs(nTr - np.sum(y1 * nTr)) > 3: print("Error in trit:  nTr %d  cts in curve %d" % (nTr, np.sum(y1*nTr)))
    if abs(nTr - np.sum(y2 * nTr)) > 3: print("Error in bkg:  nTr %d  cts in curve %d" % (nTr, np.sum(y2*nTr)))
    if abs(nTr - np.sum(y3 * nTr)) > 3: print("Error in peak:  nTr %d  cts in curve %d" % (nTr, np.sum(y3*nTr)))

    # === reverse the efficiency correction to get the "true" number of counts
    y1c = nTr * getEffCorr(x1, y1, inv=True)
    y2c = nBk * getEffCorr(x2, y2, inv=True)
    y3c = nPk * getEffCorr(x3, y3, inv=True)
    nTotC = np.sum(y1c) + np.sum(y2c) + np.sum(y3c)

    # === plot total model
    pdfs = [[x1,y1,xpb1,nTr],[x2,y2,xpb2,nBk],[x3,y3,xpb3,nPk]]
    xT, yT = getTotalModel(pdfs, eLo, eHi, epb, smooth=True)
    plt.step(xT, yT, c='b', lw=2, label="Raw (no eff. corr): %d cts" % nTot)

    # === plot components of the (uncorrected) model
    # *** NOTE: to plot after the fit, multiply by (global bin width / bin width when generated).  to integrate, don't. ***
    plt.step(x1, y1 * nTr * (epb/ppb), c='m', lw=2, alpha=0.7, label="Tritium: %d cts" % nTr)
    plt.step(x2, y2 * nBk * (epb/epb), c='g', lw=2, alpha=0.7, label="Bkg: %d cts" % nBk)
    plt.step(x3, y3 * nPk * (epb/ppb), c='c', lw=2, alpha=0.7, label="68GeK %d cts" % nPk)

    # === plot efficiency corrected final model
    pdfs = [[x1,y1c,xpb1,nTr],[x2,y2c,xpb2,nBk],[x3,y3c,xpb3,nPk]]
    xTc, yTc = getTotalModel(pdfs, eLo, eHi, epb, smooth=True, amp=False)
    plt.step(xTc, yTc, c='r', lw=3, label="Efficiency corrected: %d cts" % nTotC)

    # === plot components of the corrected model
    # plt.step(x1, y1c * (epb/ppb), c='orange', lw=2, alpha=0.7, label="trit fit: %d corr: %d" % (nTr, np.sum(y1c)))
    # plt.step(x2, y2c * (epb/epb), c='orange', lw=2, alpha=0.7, label="bkg fit: %d corr: %d" % (nBk, np.sum(y2c)))
    # plt.step(x3, y3c * (epb/ppb), c='orange', lw=2, alpha=0.7, label="peak fit: %d corr: %d" % (nPk, np.sum(y3c)))

    plt.xlabel("Energy (keV)", ha='right', x=1)
    plt.ylabel("Counts / %.1f keV" % epb, ha='right', y=1)
    plt.legend(loc=1, fontsize=12)
    plt.xlim(eLo, eHi)
    plt.ylim(ymin=0)
    plt.tight_layout()
    # plt.show()
    plt.savefig("%s/plots/sf4-mplafter.pdf" % dsi.latSWDir)


if __name__=="__main__":
    main(sys.argv[1:])