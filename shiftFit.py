#!/usr/bin/env python3
import sys, warnings, time, dsi
import numpy as np
from scipy.optimize import curve_fit
from scipy.interpolate import spline
import matplotlib
import matplotlib.pyplot as plt
plt.style.use('%s/pltReports.mplstyle' % dsi.latSWDir)
sys.argv.append("-b")
import waveLibs as wl
import ROOT
from ROOT import RooFit as RF
from ROOT import gROOT
gROOT.ProcessLine("gErrorIgnoreLevel = 3001;")
gROOT.ProcessLine("RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);")
ROOT.Math.MinimizerOptions.SetDefaultMinimizer("Minuit") # the PLC doesn't work w/ minuit2

def main():
    initialize(makePlots=False)

    # getUnscaledPDFs(makePlots=True)
    # plotPDF()
    # fitPeaks()
    # getPeakFluxRF() # do this 2 ways, since the roofit model has to be super constrained to look ok
    # getPeakFluxPY(makePlots=False) # does a sideband analysis and gaussian peak fitting.  this gives the final flux results
    # loadShiftedData() # this leaves something wierd in memory that causes other functions to segfault
    # plotShiftedData()

    # getNLLOffset()
    combineProfiles()
    exit()

    # lower than 2.0 introduces big outliers, higher than 3.5 introduces a peak
    eLo, eHi, epb = 2.0, 3.5, 0.05
    global nPks, pkModel, axPeaks

    nPks = 4
    pkModel = ["axSi_a","axSi_b","axS_a","axS_b"]
    axPeaks = [sigVals[name][0] for name in pkModel]
    pk0 = 4-nPks
    axPeaks = axPeaks[pk0:]

    fitShiftModel(eLo, eHi, epb, makePlots=False)
    plotShiftModel(eLo, eHi, epb)
    getShiftProfile(eLo, eHi, epb)
    plotShiftProfile(eLo, eHi, epb, makePlots=True)


def initialize(makePlots=False):
    global dsList, enr, opt, useWeight, floatWidth, zeroPk
    global pkModel, dMu, sigVals, simEffCorr, nPks, dsList, axPeaks
    global effLim, effMax, xEff, detEff, dsExpo, detExp, bkgModelPeaks

    dsList = [1,2,3,4,"5A","5B","5C"]
    enr = True
    useWeight = True
    simEffCorr = False
    floatWidth = True # only set true for systematic uncertainty cross check

    zeroPk = False # if true, fix the amplitude of the peak to 0
    if zeroPk:
        print("Warning, fixing amplitude to 0")

    nPks = 4 # keep only this many of the axion peaks (from 1 to 4)

    dMu = 0.05 # how much we let the peak energies float
    sigVals = {
        # expo bkg:  amp, lo, hi,    tau, lo, hi
        "exp1":     [1000, 0, 100000,   -0.8, -2, 0.1],
        "exp2":     [100, 0, 2000,   -2, -3, 0],

        # flat bkg:  amp, lo, hi
        "pol":      [200, 1, 10000],

        # peaks:       mu,  amp,  lo,   hi     sig,   lo,  hi      # peak values are confirmed by the fit to redondo's data
        "axSi_a":   [1.86,  0.1,   0, 1.,     0.03, 0.02,  0.05],  # tight-ish constraints, forcing the curve to match
        "axSi_b":   [2.00,  0.004, 0, 0.08,   0.01, 0.005, 0.02],  # the _b peaks just always go the maximum allowed
        "axS_a":    [2.45,  0.1,   0, 1.,     0.03, 0.02, 0.05],   # which is why we go with the python sideband analysis.
        "axS_b":    [2.62,  0.004, 0, 0.08,   0.01, 0.005, 0.02],

        # shifted peak
        "sPk":      [2.62,  10,  -0.1, 100.] # was 100, -0.1, 300
        }

    pkModel = ["axSi_a","axSi_b","axS_a","axS_b"]
    axPeaks = [sigVals[name][0] for name in pkModel]
    pk0 = 4-nPks
    axPeaks = axPeaks[pk0:]

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

    # special -- load slowness fraction
    fS = np.load('%s/data/efficiency-corr250.npz' % dsi.latSWDir)
    hTotSim, hSurfSim, xTotSim = fS['arr_0'], fS['arr_1'], fS['arr_2']

    # suppress a dumb divide by 0 warning for a bin I don't care about
    # print(np.geterr())
    np.seterr(divide='ignore', invalid='ignore')
    hFracSim = np.divide(hSurfSim, hTotSim, dtype=float)
    np.seterr(divide='warn', invalid='warn') # default

    idx = np.where((xTotSim >0) & (xTotSim <= 50.01))
    x, h = xTotSim[idx], hFracSim[idx]

    xS = np.arange(0, 50, 0.01)
    hS = spline(x, h, xS)
    idx = np.where((hS>0) & (detEff>0))
    effCorr = detEff[idx] / (1 - hS[idx])

    if makePlots:

        fig, p1 = plt.subplots(1, 1)

        p1.plot(xEff[idx], detEff[idx], c='b', label="Measured Efficiency")
        p1.plot(xEff[idx], effCorr, c='r', label="Sim-Corrected Efficiency")
        p1.axvline(1.5, c='g', lw=1, label="1.5 keV")
        p1.plot(np.nan, np.nan, '-m', label="Sim. Slow Pulse Fraction")

        p1a = p1.twinx()
        p1a.plot(xS[idx], hS[idx], 'm', lw=3)
        p1a.set_ylabel('Slow Fraction', color='m', ha='right', y=1)
        p1a.set_yticks(np.arange(0, 1.1, 0.2))
        p1a.tick_params('y', colors='m')

        p1.set_xlabel("Energy (keV)", ha='right', x=1)
        p1.set_ylabel("Efficiency", ha='right', y=1)
        p1.set_ylim(0, 1)
        p1.legend(loc=1, bbox_to_anchor=(0., 0.7, 0.97, 0.2))
        plt.tight_layout()
        # plt.show()
        plt.savefig("%s/plots/sf-sim-eff-corr.pdf" % dsi.latSWDir)
        plt.close()

    # *** replaces the measured efficiency w/ the simulated slow-pulse-corrected efficiency ***
    if simEffCorr:
        print("WARNING: using sim-corrected efficiency")
        detEff = effCorr
        xEff = xEff[idx]


def getUnscaledPDFs(makePlots=False):
    """ Generate a set of TH1D's to be turned into RooDataHist objects.
    Be careful they have the same axis limits and binning as the RooDataSet.
    """
    from ROOT import TFile, TH1D, gROOT

    pLo, pHi, ppb = 0, 30, 0.03 # requires ppb=0.03, the fit parameters are optimized for it
    nB = int((pHi-pLo)/ppb)

    # output file
    rOut = "%s/data/specPDFs-sf7.root" % dsi.latSWDir
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
    h1 = TH1D("h1","photoelectric",nB,pLo,pHi)         # [cm^2 / kg]
    h2 = TH1D("h2","axioelectric",nB,pLo,pHi)          # [cm^2 / kg]
    h3 = TH1D("h3","axion flux, gae=1",nB,pLo,pHi)     # [cts / (keV cm^2 d)]
    h4 = TH1D("h4","convolved flux",nB,pLo,pHi)        # [cts / (keV d kg)]
    h5 = TH1D("h5","tritium",nB,pLo,pHi)               # [cts] (normalized to 1)

    # manually fill ROOT histos (don't normalize yet)
    for iB in range(nB+1):
        ctr = (iB + 0.5)*ppb + pLo
        bLo, bHi = ctr - ppb/2, ctr + ppb/2
        with warnings.catch_warnings():
            warnings.simplefilter("ignore",category=RuntimeWarning)

            # if ma>0, we ignore entries with E <= m.
            ma=0 # this used to be a parameter but it's deprecated.

            # photoelectric x-section [cm^2 / kg]
            idx = np.where((phoData[:,0] >= bLo) & (phoData[:,0] < bHi))
            pho = np.mean(phoData[idx][:,1]) * 1000
            if np.isnan(pho) or len(phoData[idx][:,1]) == 0: pho = 0.
            if phoData[idx][:,1].any() <= ma: pho = 0.
            h1.SetBinContent(iB+1,pho)

            # axioelectric x-section [cm^2 / kg]
            if ctr > ma: axio = pho * wl.sig_ae(ctr, ma)
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


def plotPDF():

    from ROOT import TFile
    tf = TFile("%s/data/specPDFs-sf7.root" % dsi.latSWDir)

    hA = tf.Get("h3")
    xR, hR, xpb = wl.npTH1D(hA)
    plt.step(xR, hR)

    aLo, aHi = 6.4, 6.8
    nexp = hA.Integral(hA.FindBin(aLo), hA.FindBin(aHi), "width")
    # nexpy = np.sum(hR[ np.where((xR >= aLo) & (xR < aHi+xpb/2)) ]) * xpb # have to be careful about bin endpoints
    nexpy = np.sum(hR[ np.where((xR > aLo) & (xR <= aHi+xpb/2)) ]) * xpb
    print("%.1f-%.1f  %.2e [cts / cm^2 d], py: %.2e" % (aLo, aHi, nexp, nexpy))
    plt.xlim(aLo, aHi)
    plt.show()

    # === axion PDF ===
    plt.close()
    xR1, hR1, xpb = wl.npTH1D(tf.Get("h4"))
    plt.plot(xR1, hR1, ls='steps', c='b', lw=3, label=r"$\Phi_a$, %.2f keV/bin, $\mathregular{g_{ae}=1}$" % xpb)

    plt.axvline(1.85, c='g', lw=2, alpha=0.5, label=r"1.85 keV, Si ($\mathregular{K_{\alpha 1,\alpha 2}}$)")
    plt.axvline(2.00, c='m', lw=2, alpha=0.5, label=r"2.00 keV, Si ($\mathregular{K_{\beta}}$)")
    plt.axvline(2.45, c='r', lw=2, alpha=0.5, label=r"2.45 keV, S ($\mathregular{K_{\alpha 1,\alpha 2}}$)")
    plt.axvline(2.62, c='k', lw=2, alpha=0.5, label=r"2.62 keV, S ($\mathregular{K_{\beta}}$)")
    plt.axvline(6.67, c='orange', lw=2, alpha=0.8, label=r"6.67 keV, S ($\mathregular{K_{\beta}}$)")

    plt.xlabel("Energy (keV)", ha='right', x=1)
    plt.ylabel("Flux / (keV d kg)", ha='right', y=1)
    plt.xlim(0,10)
    plt.legend()
    plt.tight_layout()
    # plt.show()
    plt.savefig("%s/plots/sf7-axPDF.pdf" % dsi.latSWDir)


def fitPeaks():
    from ROOT import TFile, TH1D, TCanvas, TLegend, gStyle

    bkModel = ["exp1","exp2"]
    sigModel = pkModel + bkModel
    eLo, eHi, epb = 1.5, 3.0, 0.05
    nB = int((eHi-eLo)/epb)

    # treat the axion PDF as data
    tf = TFile("%s/data/specPDFs-sf7.root" % dsi.latSWDir)
    hAx = tf.Get("h4")
    hAx.Scale(1/1e42) # can't handle the huge scale factor
    hitE = ROOT.RooRealVar("hitE","Energy (keV)",eLo,eHi)
    hAxDH = ROOT.RooDataHist("hAx", "hAx", ROOT.RooArgList(hitE), hAx)

    # === background model ===

    pkVars = []
    for name in pkModel:
        mu, sig, amp = sigVals[name][0], 0.01, sigVals[name][1]
        if not eLo < mu < eHi: continue
        pN = ROOT.RooRealVar("amp-"+name, "amp-"+name, amp, sigVals[name][2], sigVals[name][3])
        pM = ROOT.RooRealVar("mu-"+name, "mu-"+name, mu, mu - dMu, mu + dMu)
        pS = ROOT.RooRealVar("sig-"+name, "sig-"+name, sigVals[name][4], sigVals[name][5], sigVals[name][6])
        pG = ROOT.RooGaussian("gaus-"+name, "gaus-"+name, hitE, pM, pS)
        pE = ROOT.RooExtendPdf("ext-"+name, "ext-"+name, pG, pN)
        pkVars.append([pE, name, mu, sig, amp, pN, pM, pS, pG])

    bkVars = []
    for name in bkModel:
        bkN = ROOT.RooRealVar("amp-"+name,"amp-"+name, sigVals[name][0], sigVals[name][1], sigVals[name][2])
        bkT = ROOT.RooRealVar("tau-"+name,"tau-"+name, sigVals[name][3], sigVals[name][4], sigVals[name][5])
        bkE = ROOT.RooExponential("expo-"+name,"expo-"+name, hitE, bkT)
        bkP = ROOT.RooExtendPdf("ext-"+name,"ext-"+name, bkE, bkN)
        bkVars.append([bkP,name,bkN,bkT,bkE])

    sigVars = bkVars + pkVars

    # this is separate b/c all the RooVars have to remain in memory
    pdfList = ROOT.RooArgList("shapes")
    for bkg in sigVars:
        pdfList.add(bkg[0])
    model = ROOT.RooAddPdf("model", "total PDF", pdfList)


    # === make a rooplot of the initial guess ===

    c = TCanvas("c","c",800,600)
    leg = TLegend(0.83,0.5,0.97,0.9)
    gStyle.SetPalette(ROOT.kRainBow)
    nCol = float(gStyle.GetNumberOfColors())

    fSpec = hitE.frame(RF.Range(eLo, eHi), RF.Bins(nB))
    ROOT.RooAbsData.plotOn( hAxDH, fSpec ) # can't just use plotOn with a RooDataHist b/c roofit sucks

    nTot = 0
    for i, ext in enumerate(sigVars):
        extPDF, name = ext[0], ext[1]
        col = gStyle.GetColorPalette(int(nCol/len(sigModel) * i))
        extPDF.plotOn(fSpec, RF.LineColor(col), RF.Normalization(sigVals[name][1], ROOT.RooAbsReal.Raw), RF.Name(name))
        leg.AddEntry(fSpec.findObject(name), name, "l")
        nTot += sigVals[name][1]

    model.plotOn(fSpec, RF.LineColor(ROOT.kRed), RF.Name("fmodel"), RF.Normalization(nTot, ROOT.RooAbsReal.Raw))

    fSpec.SetTitle("")
    fSpec.Draw()
    leg.Draw("same")
    c.Print("%s/plots/sf7-before.pdf" % dsi.latSWDir)

    minimizer = ROOT.RooMinimizer( model.createNLL(hAxDH, RF.NumCPU(2,0), RF.Extended(True)) )
    minimizer.setPrintLevel(-1)
    minimizer.setStrategy(2)
    minimizer.migrad()
    fitRes = minimizer.save()

    # according to the internet, covQual==3 is a good indicator that it converged
    print("Fitter is done. Fit Cov Qual:", fitRes.covQual())

    # save workspace to a TFile
    fitWS = ROOT.RooWorkspace("fitWS","Fit Workspace")
    getattr(fitWS,'import')(hitE)
    getattr(fitWS,'import')(hAxDH)
    getattr(fitWS,'import')(fitRes)
    getattr(fitWS,'import')(model)
    tf3 = TFile("%s/data/fitWS-axPks.root" % dsi.latSWDir,"RECREATE")
    fitWS.Write()
    tf3.Close()


def getPeakFluxRF():

    from ROOT import TFile, TCanvas, TH1D, TLegend, gStyle

    bkModel = ["exp1","exp2"]
    sigModel = pkModel + bkModel
    eLo, eHi, epb = 1.5, 3.0, 0.05
    nB = int((eHi-eLo)/epb)

    f = TFile("%s/data/fitWS-axPks.root" % dsi.latSWDir)
    fitWS = f.Get("fitWS")
    hAxDH = fitWS.allData().front()
    fitRes = fitWS.allGenericObjects().front()
    nPars = fitRes.floatParsFinal().getSize()
    hitE = fitWS.var("hitE")
    model = fitWS.pdf("model")

    # === get fit results: {name : [nCts, err]} ===
    fitVals = {}
    print("fit vals:")
    for i in range(nPars):
        fp = fitRes.floatParsFinal()
        name = fp.at(i).GetName()
        fitVal, fitErr = fp.at(i).getValV(), fp.at(i).getError()
        fitVals[name] = [fitVal, fitErr]
        print("%-10s" % name, wl.niceList(fitVals[name], "%.3f"))

    for name in fitVals:
        if "amp-" in name:
            print("%-10s  %.3f ± %-5.3f  Flux: %.3e ± %.3e" % (name, fitVals[name][0], fitVals[name][1],  fitVals[name][0]*1e42, fitVals[name][1]*1e42))

    # === make a rooplot of the fit ===

    leg = TLegend(0.83,0.5,0.97,0.9)
    gStyle.SetPalette(ROOT.kRainBow)
    nCol = float(gStyle.GetNumberOfColors())

    fSpec = hitE.frame(RF.Range(eLo,eHi), RF.Bins(nB))

    hAxDH.plotOn(fSpec)

    for i, name in enumerate(sigModel):
        pdfName = "ext-"+name
        col = gStyle.GetColorPalette(int(nCol/len(sigModel) * i))
        model.plotOn(fSpec, RF.Components(pdfName), RF.LineColor(col), RF.LineStyle(ROOT.kDashed), RF.Name(name))
        leg.AddEntry(fSpec.findObject(name), name, "l")

    chiSquare = fSpec.chiSquare(nPars)
    model.plotOn(fSpec, RF.LineColor(ROOT.kRed), RF.Name("fmodel"))
    # leg.AddEntry(fSpec.findObject("fmodel"),"Full Model, #chi^{2}/NDF = %.3f" % chiSquare, "l")
    leg.AddEntry(fSpec.findObject("fmodel"),"Full Model", "l")

    c = TCanvas("c","c", 1400, 1000)
    fSpec.SetTitle("")
    fSpec.Draw()
    leg.Draw("same")
    c.Print("%s/plots/sf7-after.pdf" % dsi.latSWDir)
    c.Clear()


def getPeakFluxPY(makePlots=False):
    """
    Flux results (h3), eLo, eHi, epb = 1.5, 3, 0.05
    1.86  5.534e+38 ± 2.243e+37
    2.00  1.340e+38 ± 1.803e+37
    2.45  4.671e+38 ± 1.990e+37
    2.62  8.858e+37 ± 1.828e+37
    bkg-mu   17.170 ± 0.901
    bkg-tau  -4.603 ± 0.464
    bkg-b    -5.776 ± 1.064
    Total peak flux: 1.24e+39
    Integral 1.5-3.0  7.52e+39 [cts / cm^2 d], py: 7.52e+39
    Pct in peaks: 16.523%

    PDF results (h4), eLo, eHi, epb = 1.5, 3, 0.05
    1.86  1.438e+41
    2.00  3.244e+40
    2.45  1.061e+41
    2.62  2.104e+40
    bkg-mu   7.272 ± 0.108
    bkg-tau  -1.173 ± 0.021
    bkg-b    0.037 ± 0.020
    Total peak flux: 3.03e+41
    Integral 1.5-3.0  1.84e+42 [cts / cm^2 d], py: 1.84e+42
    Pct in peaks: 16.463%
    """
    from ROOT import TFile

    eLo, eHi, epb = 1.5, 3.0, 0.05
    nB = int((eHi-eLo)/epb)

    tf = TFile("%s/data/specPDFs-sf7.root" % dsi.latSWDir)

    # ---- two options for this plot ----
    rHist, scale = tf.Get("h3"), 1e39  # axion flux (reproduce frank's numbers)
    # rHist, scale = tf.Get("h4"), 1e42  # axion PDF
    # -----------------------------------

    x0, h0, xpb = wl.npTH1D(rHist, opt="")
    idx = np.where((x0>=eLo) & (x0<=eHi))
    xA, hA = x0[idx], h0[idx] / scale # scale the histo

    # sideband analysis to get the background
    apb = 0.03
    pkExc = [1.86 - 4*apb, 2.00 + 3*apb, 2.45-3*apb, 2.62+3*apb]
    idxSB = np.where( (xA < pkExc[0]) | ((xA > pkExc[1]) & (xA < pkExc[2]) | (xA > pkExc[3])) )

    # fit to one exp (a fit to twoExp has super huge errors, the two are probably degenerate)
    init = (1, -2, -1) # mu, tau, b
    poptBk, pcovBk = curve_fit(wl.expFunc, xA[idxSB], hA[idxSB], p0=init)
    perrBk = np.sqrt(np.diag(pcovBk))

    xpbF = 0.005
    xF = np.arange(eLo, eHi, xpbF)
    hS = hA - wl.expFunc(xA, *poptBk)
    hF = wl.expFunc(xF, *poptBk)

    # set up the initial guesses and fit pars roughly the same as in roofit
    pars = []
    pars.extend([1.86, 0.025, 0.72]) # mu, sig, amp
    pars.extend([2.00, 0.005, 0.01])
    pars.extend([2.45, 0.023, 0.58])
    pars.extend([2.62, 0.005, 0.01])
    bLo = [
        1.86-dMu, 0.02, 0,
        2.00-dMu, 0.005, 0,
        2.45-dMu, 0.02, 0,
        2.62-dMu, 0.004, 0
        ]
    bHi = [
        1.86+dMu, 0.05, 1,
        2.00+dMu, 0.02, 0.03,
        2.45+dMu, 0.05, 1,
        2.62+dMu, 0.02, 0.03,
        ]
    bnds = (tuple(bLo),tuple(bHi))
    po, pcov = curve_fit(wl.nGaus, xA, hS, p0=pars, bounds=bnds)
    pe = np.sqrt(np.diag(pcov))

    fitVals = {}
    fitVals["axSi_a"] = [po[0],pe[0], po[1],pe[1], po[2],pe[2]]  # mu, muE, sig, sigE, amp, ampE
    fitVals["axSi_b"] = [po[3],pe[3], po[4],pe[4], po[5],pe[5]]
    fitVals["axS_a"]  = [po[6],pe[6], po[7],pe[7], po[8],pe[8]]
    fitVals["axS_b"]  = [po[9],pe[9], po[10],pe[10], po[11],pe[11]]

    # make peak gaussians
    yP1 = wl.gaus(xF, po[0], po[1], po[2])
    yP2 = wl.gaus(xF, po[3], po[4], po[5])
    yP3 = wl.gaus(xF, po[6], po[7], po[8])
    yP4 = wl.gaus(xF, po[9], po[10], po[11])

    # figure out the error
    yS1 = wl.gaus(xF, po[0], po[1]-pe[1], po[2]-pe[2]) # mu doesn't make a difference.  max error is - -.
    # print(np.sum(yS1), np.sum(yP1), np.sum(yP1) / (1 - pe[2]/po[2]))
    yE1 = pe[2]/po[2] # take the error to be dominated by the amplitude
    yE2 = pe[5]/po[5]
    yE3 = pe[8]/po[8]
    yE4 = pe[11]/po[11]

    # print results

    # NOTE: when plotting in mpl, to match roofit's rooplot, do:  y *= (xpbF/xpb) * scale,
    #          but when getting the number of counts, do:         y *= xpb * scale.
    #
    # (this is not quite the same rule as for when the PDF's are from TH1D's, as in the continuum fit.)

    nCts1 = np.sum(yP1) * xpb * scale
    nCts2 = np.sum(yP2) * xpb * scale
    nCts3 = np.sum(yP3) * xpb * scale
    nCts4 = np.sum(yP4) * xpb * scale
    nTot = nCts1 + nCts2 + nCts3 + nCts4

    nE1 = nCts1 * yE1
    nE2 = nCts2 * yE2
    nE3 = nCts3 * yE3
    nE4 = nCts4 * yE4

    print("%.2f  %.3e ± %.3e" % (po[0], nCts1, nE1) )
    print("%.2f  %.3e ± %.3e" % (po[3], nCts2, nE2) )
    print("%.2f  %.3e ± %.3e" % (po[6], nCts3, nE3) )
    print("%.2f  %.3e ± %.3e" % (po[9], nCts4, nE4) )
    print("bkg-mu   %.3f ± %.3f" % (poptBk[0], perrBk[0]))
    print("bkg-tau  %.3f ± %.3f" % (poptBk[1], perrBk[1]))
    print("bkg-b    %.3f ± %.3f" % (poptBk[2], perrBk[2]))
    print("Total peak flux: %.2e" % (nTot))

    # save flux output for g_ae calculation
    axFluxes = {"%.2f"%po[0]:nCts1, "%.2f"%po[3]:nCts2, "%.2f"%po[6]:nCts3, "%.2f"%po[9]:nCts4}
    np.savez("%s/data/sf7-pkFluxes.npz" % dsi.latSWDir, axFluxes)

    # check results are consistent w/ axion flux
    nexp = rHist.Integral(rHist.FindBin(eLo), rHist.FindBin(eHi), "width")
    nexpy = np.sum(h0[ np.where((x0 > eLo) & (x0 <= eHi+xpb/2)) ]) * xpb # have to include the last bin
    print("Integral %.1f-%.1f  %.2e [cts / cm^2 d], py: %.2e" % (eLo, eHi, nexp, nexpy))
    print("Pct in peaks: %.3f%%" % (100*nTot/nexpy))

    if makePlots:

        plt.plot(xA, hA, '.k', ms=6, label=r'Axion PDF, $\mathregular{g_{ae}}$=1')
        # plt.plot(xA[idxSB], hA[idxSB], ".", c='b', lw=3, label='Sideband')
        plt.plot(xA, hS, ls='steps-mid', c='k', alpha=0.7, label='Bkg. Subt.')
        plt.plot(xF, hF + (yP1 + yP2 + yP3 + yP4), c='r', alpha=0.7, label="Total Model")

        i1 = np.where(yP1>0.01)
        i2 = np.where(yP2>0.01)
        i3 = np.where(yP3>0.01)
        i4 = np.where(yP4>0.01)

        sigLabels = {
            "axSi_a": r"Si ($\mathregular{K_{\alpha 1,\alpha 2}}$)",
            "axSi_b": r"Si ($\mathregular{K_{\beta}}$)",
            "axS_a": r"S ($\mathregular{K_{\alpha 1,\alpha 2}}$)",
            "axS_b": r"S ($\mathregular{K_{\beta}}$)"
            }

        plt.plot(xF[i1], yP1[i1], c='g', lw=3, label=r"%.2f keV, $\phi$ = %.2e, %s" % (po[0], nCts1, sigLabels["axSi_a"]))
        plt.plot(xF[i2], yP2[i2], c='orange', lw=3, label=r"%.2f keV, $\phi$ = %.2e, %s" % (po[3], nCts2, sigLabels["axSi_b"]))
        plt.plot(xF[i3], yP3[i3], c='m', lw=3, label=r"%.2f keV, $\phi$ = %.2e, %s" % (po[6], nCts3, sigLabels["axS_a"]))
        plt.plot(xF[i4], yP4[i4], c='b', lw=3, label=r"%.2f keV, $\phi$ = %.2e, %s" % (po[9], nCts4, sigLabels["axS_b"]))
        # plt.plot(xF, expFunc(xF, *poptBk))

        plt.xlabel("Energy (keV)", ha='right', x=1)

        if rHist.GetName()=="h3":
            plt.ylabel(r"Flux / %.0e (keV $\mathregular{cm^2}$ d)" % scale, ha='right', y=1)
            plt.ylim(-0.1, 7)
            plt.legend(loc=1, fontsize=12, bbox_to_anchor=(0.5, 0.75))
            plt.tight_layout()
            plt.savefig("%s/plots/sf7-axPeakFluxFit.pdf" % (dsi.latSWDir))

        elif rHist.GetName()=="h4":
            plt.ylabel("Flux / %.0e (keV d kg)" % scale, ha='right', y=1)
            plt.legend(loc=1, fontsize=12)
            plt.ylim(ymax=2.5)
            plt.tight_layout()
            plt.savefig("%s/plots/sf7-axPeakPDFFit.pdf" % (dsi.latSWDir))


def loadShiftedData():
    from ROOT import TFile, TTree, TList
    from array import array

    print("Getting shifted data for %d peaks:" % nPks, axPeaks)

    fName = "%s/data/latDS%s.root" % (dsi.latSWDir, ''.join([str(d) for d in dsList]))
    tf1 = TFile(fName)
    treeIn = tf1.Get("skimTree")

    run, iEvent, iHit = array('i',[0]), array('i',[0]), array('i',[0]),
    channel, hitE = array('i',[0]), array('d',[0])
    isEnr, weight = array('i',[0]), array('d',[0])
    treeIn.SetBranchAddress("run",run)
    treeIn.SetBranchAddress("iEvent",iEvent)
    treeIn.SetBranchAddress("iHit",iHit)
    treeIn.SetBranchAddress("channel",channel)
    treeIn.SetBranchAddress("trapENFCal",hitE)
    treeIn.SetBranchAddress("isEnr",isEnr)
    treeIn.SetBranchAddress("weight",weight)

    # create ths shifted tree
    tfName = "%s/data/latDS%s_shifted_%dpks.root" % (dsi.latSWDir, ''.join([str(d) for d in dsList]), nPks)
    tf2 = TFile(tfName, "RECREATE")
    tList = TList()
    tVec = [0,0,0,0]

    # get the axion peak list from globals
    for i in range(len(axPeaks)):
        shift = axPeaks[-1] - axPeaks[i] # line up everything with the last peak

        tTitle = "pkE: %.2f keV, shift %.2f keV" % (axPeaks[i], shift)
        tVec[i] = TTree("t%d" % i, tTitle)
        ene, wt = array('d',[0.]), array('d',[0.])
        enr1 = array('i',[0])
        tVec[i].Branch("trapENFCal",ene,"trapENFCal/D")
        tVec[i].Branch("weight",wt,"weight/D")
        tVec[i].Branch("isEnr",enr1,"isEnr/I")
        for j in range(treeIn.GetEntries()):
            treeIn.GetEntry(j)
            ene[0], enr1[0] = hitE[0], isEnr[0]

            # calculate weight based on 1/efficiency
            if ene[0] > effLim:
                wt[0] = 1 / effMax
            else:
                idx = (np.abs(xEff-hitE[0])).argmin()
                wt[0] = 1/np.interp(hitE[0], xEff[idx:idx+1], detEff[idx:idx+1])
            # if hitE[0] < effLim:
                # print("%.2f  %.2f " % (hitE[0], wt[0]))

            ene[0] += shift
            tVec[i].Fill()

        tList.Add(tVec[i])
        tVec[i].Write()
        print("Tree %d - %d entries. %s" % (i, tVec[i].GetEntries(), tTitle))

    tShift = ROOT.TTree.MergeTrees(tList)
    tShift.SetName("mergeTree")
    tShift.SetTitle("mergeTree")
    print("Trees merged, with %d entries total." % tShift.GetEntries())
    print("Writing file:", tfName)
    tShift.Write()
    tf2.Close()
    tf1.Close()


def plotShiftedData():

    from ROOT import TFile, TTree, TList

    xLo, xHi, xpb = 1.8, 15, 0.05

    tfName = "%s/data/latDS%s_shifted_%dpks.root" % (dsi.latSWDir, ''.join([str(d) for d in dsList]), nPks)
    tf = TFile(tfName)
    tCut = "isEnr" if enr is True else "!isEnr"
    # print("Plotting shifted data for %d peaks" % nPks, axPeaks)

    # merged tree
    ttM = tf.Get("mergeTree")
    n = ttM.Draw("trapENFCal:weight",tCut,"goff")
    hE, hW = ttM.GetV1(), ttM.GetV2()
    hE = [hE[i] for i in range(n)]
    hW = [hW[i] for i in range(n)] if useWeight else None
    xM, yM = wl.GetHisto(hE, xLo, xHi, xpb, wts=hW)

    # create shifted & merged spectrum
    lab = "Spectrum, Merged & Weighted" if useWeight else "Merged Spectrum"
    plt.step(xM, yM, c='k', lw=2, label=lab)

    # loop over individual trees (the last one will be unshifted)
    tKeys = [key.GetName() for key in tf.GetListOfKeys() if key.GetName() != "mergeTree"]
    cmap = plt.cm.get_cmap('jet',len(tKeys))
    for i, tKey in enumerate(tKeys):
        tt = tf.Get(tKey)
        n = tt.Draw("trapENFCal:weight",tCut,"goff")
        hE, hW = tt.GetV1(), tt.GetV2()
        hE = [hE[i] for i in range(n)]
        hW = [hW[i] for i in range(n)] if useWeight else None
        x, y = wl.GetHisto(hE, xLo, xHi, xpb, wts=hW)
        plt.step(x, y, lw=2, c=cmap(i), label=tt.GetTitle())

    plt.axvline(axPeaks[-1], c='r', lw=4, alpha=0.7, label="Shift Peak: %.2f keV" % axPeaks[-1])
    plt.legend(fontsize=12)
    plt.xlabel("Shifted Energy (keV)", ha='right', x=1)
    plt.ylabel("Counts / %.2f keV" % xpb, ha='right', y=1)
    plt.xlim(xLo, xHi)
    plt.ylim(ymin=0)
    plt.tight_layout()
    # plt.show()
    # plt.savefig("%s/plots/sf7-shiftSpec-%dpks.pdf" % (dsi.latSWDir, nPks))


    # === try fitting the data to some different backgrounds ===
    # since the axion fit is suuper sensitive to choice of bkg and energy range
    plt.close()

    xLo, xHi = 1.9, 3.5 # lower than 1.9 introduces big outliers, higher than 3.5 introduces a peak

    xP = np.arange(xLo, xHi, xpb) + xpb/2
    idx = np.where((yM>0) & (xM>=xLo) & (xM<=xHi))
    xM, yM = xM[idx], yM[idx]
    yE = np.asarray([np.sqrt(y) for y in yM]) # statistical error

    # one expo - wenqin says this is much too strong an assumption
    # it gives the best result on g_ae, but probably b/c there aren't enough free parameters
    # to fit the background properly
    # p0 = (100, -0.5)
    # popt,_ = curve_fit(wl.oneExp, xM, yM, p0=p0)
    # yP = wl.oneExp(xP, *popt)
    # chi2R = np.sum(((wl.oneExp(xM, *popt)-yM)/yE)**2) / (len(xM)-len(popt)) # chi2/ndf
    # print("%d pks  oneExp  %.1f-%.1f kev chi2/ndf %.4f" % (nPks, xLo, xHi, chi2R), wl.niceList(list(popt)))
    # plt.plot(xP, yP, 'r', label=r"one exp, $\chi^2$/ndf = %.2f" % (chi2R))

    # two expos - lower chi2 than pol3 for this energy range and 4 peaks
    p0 = (100, -0.5, -0.8)
    popt,_ = curve_fit(wl.twoExp, xM, yM, p0=p0)
    yP = wl.twoExp(xP, *popt)
    chi2R = np.sum(((wl.twoExp(xM, *popt)-yM)/yE)**2) / (len(xM)-len(popt)) # chi2/ndf
    print("%d pks  twoExp  %.1f-%.1f kev chi2/ndf %.4f" % (nPks, xLo, xHi, chi2R), wl.niceList(list(popt)))
    plt.plot(xP, yP, 'g', label=r"two exp, $\chi^2$/ndf = %.2f" % (chi2R))

    # pol3 - lowest chi2R for this energy range and 4 peaks
    p0 = (200, -100, 8, 1)
    popt,_ = curve_fit(wl.nPol, xM, yM, p0=p0)
    yP = wl.nPol(xP, *popt)
    chi2R = np.sum(((wl.nPol(xM, *popt)-yM)/yE)**2) / (len(xM)-len(popt)) # chi2/ndf
    print("%d pks  pol%d  %.1f-%.1f kev chi2/ndf %.4f" % (nPks, len(p0)-1, xLo, xHi, chi2R), wl.niceList(list(popt)))
    plt.plot(xP, yP, 'b', label=r"pol3, $\chi^2$/ndf = %.2f" % (chi2R))

    plt.errorbar(xM, yM, yerr=yE, c='k', ms=5, linewidth=0.5, fmt='.', capsize=1, zorder=1)
    plt.xlim(xLo, xHi)
    plt.legend()
    plt.xlabel("Shifted Energy (keV)", ha='right', x=1)
    plt.ylabel("Counts / %.2f keV" % xpb, ha='right', y=1)
    plt.tight_layout()
    # plt.show()
    plt.savefig("%s/plots/sf7-polDataFit-%dpks.pdf" % (dsi.latSWDir, nPks))


def getSigma(E, opt=""):
    """ Get the MJ energy resolution, IN SIGMA (source tables are in terms of FWHM)
    If multiple DS are selected, weight the curve by DS exposure.
    Uses the global variable 'dsList'.
    """
    # HG resolutions IN FWHM, from the energy unidoc.
    eRes = {
        0 :    {"nat": [1.260e-1, 1.790e-2, 2.370e-4], "enr": [1.500e-1, 1.750e-2, 2.820e-4], "both": [1.470e-1, 1.730e-2, 3.000e-4]},
        1 :    {"nat": [1.470e-1, 1.770e-2, 2.010e-4], "enr": [1.340e-1, 1.750e-2, 2.820e-4], "both": [1.360e-1, 1.740e-2, 2.800e-4]},
        2 :    {"nat": [1.410e-1, 1.800e-2, 1.680e-4], "enr": [1.420e-1, 1.720e-2, 2.860e-4], "both": [1.430e-1, 1.720e-2, 2.840e-4]},
        3 :    {"nat": [1.800e-1, 1.820e-2, 2.090e-4], "enr": [1.580e-1, 1.710e-2, 3.090e-4], "both": [1.620e-1, 1.720e-2, 2.970e-4]},
        4 :    {"nat": [2.140e-1, 1.540e-2, 3.970e-4], "enr": [2.170e-1, 1.490e-2, 3.190e-4], "both": [2.180e-1, 1.500e-2, 3.500e-4]},
        "5A" : {"nat": [2.248e-1, 1.894e-2, 2.794e-4], "enr": [2.660e-1, 2.215e-2, 2.868e-4], "both": [2.592e-1, 2.057e-2, 3.086e-4]},
        "5B" : {"nat": [1.650e-1, 1.760e-2, 2.828e-4], "enr": [1.815e-1, 1.705e-2, 3.153e-4], "both": [1.815e-1, 1.690e-2, 3.187e-4]},
        "5C" : {"nat": [1.565e-1, 1.810e-2, 2.201e-4], "enr": [1.361e-1, 1.740e-2, 2.829e-4], "both": [1.519e-1, 1.718e-2, 2.762e-4]}
    }

    if len(dsList)==1:
        p = eRes[dsList[0]][opt]
        return np.sqrt(p[0]**2 + p[1]**2 * E + p[2]**2 * E**2)/2.355
    else:
        # weight the curve by exposure
        sig, expTot = 0, 0
        for ds in dsList:
            if opt=="enr": exp = dsExpo[ds][0]
            if opt=="nat": exp = dsExpo[ds][1]
            if opt=="both": exp = dsExpo[ds][0] + dsExpo[ds][1]
            p = eRes[ds][opt]
            sig += np.sqrt(p[0]**2 + p[1]**2 * E + p[2]**2 * E**2) * exp
            expTot += exp
        sig /= expTot
        return sig/2.355


def fitShiftModel(eLo, eHi, epb, makePlots=True):

    from ROOT import TFile, TH1D, TCanvas, TLegend, gStyle

    bkModel = ["exp1", "pol"]
    sigModel = ["sPk"] + bkModel
    nB = int((eHi-eLo)/epb)

    # === load data into workspace ===

    tfName = "%s/data/latDS%s_shifted_%dpks.root" % (dsi.latSWDir, ''.join([str(d) for d in dsList]), nPks)
    print("Fitting shifted data for %d peaks" % nPks, axPeaks)

    tf = TFile(tfName)
    tt = tf.Get("mergeTree")
    tCut = "isEnr" if enr is True else "!isEnr"
    hitE = ROOT.RooRealVar("trapENFCal", "Energy", eLo, eHi, "keV")
    hEnr = ROOT.RooRealVar("isEnr", "isEnr", 0, 1, "")
    hitW = ROOT.RooRealVar("weight", "weight", 0, 10000, "")

    if useWeight:
        fData = ROOT.RooDataSet("data", "data", tt, ROOT.RooArgSet(hitE, hEnr, hitW), tCut, "weight")
    else:
        fData = ROOT.RooDataSet("data", "data", tt, ROOT.RooArgSet(hitE, hEnr), tCut)


    # === signal model: 1 peak, 2 exponentials ===

    # NOTE: since we overlap everything at the highest-E peak, we are limited to the MJD resolution at that value.
    # we also don't allow the peak mean to float (since trapENFCal seems OK to 0.01 kev in this region)
    # or sigma, i guess

    name = "sPk"
    pkVars = []
    opt = "enr" if enr else "nat"
    mu, sig, amp = sigVals[name][0], getSigma(sigVals[name][0], opt), sigVals[name][1]
    # print("e-region: sigma: %.2f  lo %.2f  mean %.2f  hi %.2f" % (sig, mu-3*sig, mu, mu+3*sig))

    if zeroPk:
        pN = ROOT.RooRealVar("amp-"+name, "amp-"+name, 0)
    else:
        pN = ROOT.RooRealVar("amp-"+name, "amp-"+name, amp, sigVals[name][2], sigVals[name][3])

    pM = ROOT.RooRealVar("mu-"+name, "mu-"+name, mu)

    if floatWidth:
        pS = ROOT.RooRealVar("sig-"+name, "sig-"+name, sig, sig - 0.3*sig, sig + 0.3*sig) # << systematic check, floating width
        print("Warning, using floating width. sig %.3f, 0.3*sig: %.3f" % (sig, 0.3*sig))

    else:
        pS = ROOT.RooRealVar("sig-"+name, "sig-"+name, sig) # << fixed value, used in main fit
        print("Fixed mean and sigma:",mu, sig)

    pG = ROOT.RooGaussian("gaus-"+name, "gaus-"+name, hitE, pM, pS)
    pE = ROOT.RooExtendPdf("ext-"+name, "ext-"+name, pG, pN)
    pkVars.append([pE, name, mu, sig, amp, pN, pM, pS, pG])

    bkVars = []
    for name in bkModel:
        if "exp" in name:
            bkN = ROOT.RooRealVar("amp-"+name,"amp-"+name, sigVals[name][0], sigVals[name][1], sigVals[name][2])
            bkT = ROOT.RooRealVar("tau-"+name,"tau-"+name, sigVals[name][3], sigVals[name][4], sigVals[name][5])
            bkE = ROOT.RooExponential("expo-"+name,"expo-"+name, hitE, bkT)
            bkP = ROOT.RooExtendPdf("ext-"+name,"ext-"+name, bkE, bkN)
            bkVars.append([bkP,name,bkN,bkT,bkE])
        if "pol" in name:
            # this sucks, it keeps giving me negative values
            name="pol"
            bkN = ROOT.RooRealVar("amp-"+name,"amp-"+name, sigVals[name][0], sigVals[name][1], sigVals[name][2])
            # bkC1 = ROOT.RooRealVar("c1","c1", -1, -1000, 1000)
            # bkC2 = ROOT.RooRealVar("c2","c2", -10)
            # bkC3 = ROOT.RooRealVar("c3","c3", 0)
            bkPo = ROOT.RooPolynomial("pol-"+name,"pol-"+name, hitE, ROOT.RooArgList()) # y=1+c1x+c2x^2+c3x^3
            bkEP = ROOT.RooExtendPdf("ext-"+name,"ext-"+name,bkPo, bkN)
            bkVars.append([bkEP,name,bkPo,bkN])

    sigVars = bkVars + pkVars

    # this is separate b/c all the RooVars have to remain in memory
    pdfList = ROOT.RooArgList("shapes")
    for bkg in sigVars:
        pdfList.add(bkg[0])
    model = ROOT.RooAddPdf("model", "total PDF", pdfList)

    if makePlots:

        # === make a rooplot of the initial guess ===

        c = TCanvas("c","c",800,600)
        leg = TLegend(0.83,0.5,0.97,0.9)
        gStyle.SetPalette(ROOT.kRainBow)
        nCol = float(gStyle.GetNumberOfColors())

        fSpec = hitE.frame(RF.Range(eLo, eHi), RF.Bins(nB))
        fData.plotOn(fSpec)

        nTot = 0
        for i, ext in enumerate(sigVars):
            extPDF, name = ext[0], ext[1]
            col = gStyle.GetColorPalette(int(nCol/len(sigModel) * i))
            extPDF.plotOn(fSpec, RF.LineColor(col), RF.Normalization(sigVals[name][1], ROOT.RooAbsReal.Raw), RF.Name(name))
            leg.AddEntry(fSpec.findObject(name), name, "l")
            nTot += sigVals[name][1]

        model.plotOn(fSpec, RF.LineColor(ROOT.kRed), RF.Name("fmodel"), RF.Normalization(nTot, ROOT.RooAbsReal.Raw))

        fSpec.SetTitle("")
        fSpec.Draw()
        leg.Draw("same")
        c.Print("%s/plots/sf7-shift-before-%dpks.pdf" % (dsi.latSWDir, nPks))

    # ==== ok, now run the fit ===

    minimizer = ROOT.RooMinimizer( model.createNLL(fData, RF.NumCPU(2,0), RF.Extended(True)) )
    minimizer.setPrintLevel(-1)
    minimizer.setStrategy(2)
    minimizer.migrad()
    fitRes = minimizer.save()

    # according to the internet, covQual==3 is a good indicator that it converged
    print("Fitter is done. Fit Cov Qual:", fitRes.covQual())

    # save workspace to a TFile
    fitWS = ROOT.RooWorkspace("fitWS","Fit Workspace")
    getattr(fitWS,'import')(hitE)
    getattr(fitWS,'import')(fData)
    getattr(fitWS,'import')(hitW)
    getattr(fitWS,'import')(fitRes)
    getattr(fitWS,'import')(model)
    tf3 = TFile("%s/data/fitWS-axShift-%dpks.root" % (dsi.latSWDir, nPks),"RECREATE")
    fitWS.Write()
    tf3.Close()


def plotShiftModel(eLo, eHi, epb, plotProfileResults=True):
    from ROOT import TFile, TCanvas, TH1D, TLegend, gStyle

    print("Plotting shifted fit for %d peaks" % nPks, axPeaks)

    bkModel = ["exp1","pol"]
    sigModel = ["sPk"] + bkModel
    nB = int((eHi-eLo)/epb)

    f = TFile("%s/data/fitWS-axShift-%dpks.root" % (dsi.latSWDir, nPks))
    fitWS = f.Get("fitWS")
    fData = fitWS.allData().front()
    fitRes = fitWS.allGenericObjects().front()
    nPars = fitRes.floatParsFinal().getSize()
    hitE = fitWS.var("trapENFCal")
    model = fitWS.pdf("model")
    minNLL = fitRes.minNll()
    print("minNLL:",minNLL)

    # if zeroPk:
    # exit()

    # === get fit results: {name : [nCts, err]} ===
    fitVals = {}
    print("fit vals:")
    for i in range(nPars):
        fp = fitRes.floatParsFinal()
        name = fp.at(i).GetName()
        fitVal, fitErr = fp.at(i).getValV(), fp.at(i).getError()
        fitVals[name] = [fitVal, fitErr]
        # print("%-10s" % name, wl.niceList(fitVals[name], "%.3f"))

    profileVars = ["sPk"]
    if plotProfileResults and not zeroPk:
        from ROOT import RooStats as RS
        for pName in profileVars:
            fitVar = "amp-"+pName
            fitVal = fitVals[fitVar][0]
            thisVar = fitWS.var(fitVar)
            pCL = 0.9
            plc = RS.ProfileLikelihoodCalculator(fData, model, ROOT.RooArgSet(thisVar))
            plc.SetConfidenceLevel(0.90)
            interval = plc.GetInterval()
            lower = interval.LowerLimit(thisVar)
            upper = interval.UpperLimit(thisVar)
            print("%.0f%% CL Upper Limit, %s: %.2f" % (100*pCL, fitVar, upper))
            fitVals[fitVar][0] = upper
            thisVar.setVal(upper)

    for name in fitVals:
        # if "amp-" in name:
        print("%-10s  %.3f ± %-5.3f" % (name, fitVals[name][0], fitVals[name][1]))

    # === make a rooplot of the fit ===

    leg = TLegend(0.83,0.5,0.97,0.9)
    gStyle.SetPalette(ROOT.kRainBow)
    nCol = float(gStyle.GetNumberOfColors())

    fSpec = hitE.frame(RF.Range(eLo,eHi), RF.Bins(nB))

    fData.plotOn(fSpec)

    for i, name in enumerate(sigModel):
        pdfName = "ext-"+name
        col = gStyle.GetColorPalette(int(nCol/len(sigModel) * i))
        model.plotOn(fSpec, RF.Components(pdfName), RF.LineColor(col), RF.LineStyle(ROOT.kDashed), RF.Name(name))
        leg.AddEntry(fSpec.findObject(name), name, "l")

    chiSquare = fSpec.chiSquare(nPars)
    model.plotOn(fSpec, RF.LineColor(ROOT.kRed), RF.Name("fmodel"))
    # leg.AddEntry(fSpec.findObject("fmodel"),"Full Model, #chi^{2}/NDF = %.3f" % chiSquare, "l")
    leg.AddEntry(fSpec.findObject("fmodel"),"Full Model", "l")

    c = TCanvas("c","c", 1400, 1000)
    fSpec.SetTitle("")
    fSpec.Draw()
    leg.Draw("same")
    c.Print("%s/plots/sf7-shift-after-%dpks.pdf" % (dsi.latSWDir, nPks))
    c.Clear()

    # === make a pyplot of the fit.  be careful, the data might be weighted ===

    plt.close()

    tfName = "%s/data/latDS%s_shifted_%dpks.root" % (dsi.latSWDir, ''.join([str(d) for d in dsList]), nPks)
    tf = TFile(tfName)
    tt = tf.Get("mergeTree")
    print("mergeTree has",tt.GetEntries())
    tCut = "isEnr" if enr else "!isEnr"
    # tCut += " && trapENFCal >= %.1f && trapENFCal <= %.1f" % (eLo, eHi)
    tCut = "trapENFCal > 1 && trapENFCal < 15 && isEnr"
    n = tt.Draw("trapENFCal:weight", tCut, "goff")
    trapE, weight = tt.GetV1(), tt.GetV2()
    trapE = [trapE[i] for i in range(n)]
    weight = [weight[i] for i in range(n)]

    # plot data
    x, hData = wl.GetHisto(trapE, eLo, eHi, epb, wts=weight)
    hErr = np.asarray([np.sqrt(h) for h in hData]) # statistical error
    print("Avg. error: ± %.2f cts" % (np.mean(hErr)*2))

    plt.errorbar(x, hData, yerr=hErr, c='k', ms=5, linewidth=0.5, fmt='.', capsize=1, zorder=1)

    lab = "Shifted & Weighted Data" if useWeight else "Shifted Data"
    lab += ", DS%s-%s, %.3f kg-y" % (str(dsList[0]), str(dsList[-1]), detExp/365.25)
    plt.plot(np.nan, np.nan, ".k", ms=5, label=lab)

    # plot the fit components
    pdfs, pdfsCorr = [], []
    nTot = 0

    xpbF = 0.005
    xF = np.arange(eLo, eHi, xpbF)
    yF = np.zeros(len(xF))

    # sum peak
    opt = "enr" if enr else "nat"
    mu, sig = sigVals["sPk"][0], getSigma(sigVals["sPk"][0], opt)

    if zeroPk:
        pCts = 0
    else:
        pCts = fitVals["amp-sPk"][0]
    yP = wl.gaus(xF, mu, sig, pCts) # note: to match pCts, do np.sum(yP * xpbF)
    # plt.plot(xF, yP * epb, c='b', lw=3, label="%.2f keV sum peak, %.2f cts (%.0f%% CL)" % (axPeaks[-1], pCts, 100*pCL))
    nTot += pCts
    yF = np.add(yF, yP * epb)

    # two exponential bkg's, added together
    # bTot, bTotCts = np.zeros(len(xF)), 0
    # for i, name in enumerate(bkModel):
    #     nCts, c = fitVals["amp-"+name][0], fitVals["tau-"+name][0]
    #     yB = wl.oneExp(xF, nCts, c)
    #     yB /= np.sum(yB)
    #     bTot = np.add(bTot, yB * nCts * (epb/xpbF))
    #     bTotCts += nCts
    #     nTot += nCts

    # one expo and one flat bkg
    bTot, bTotCts = np.zeros(len(xF)), 0

    nCts, c = fitVals["amp-exp1"][0], fitVals["tau-exp1"][0]
    yB = wl.oneExp(xF, nCts, c)
    yB /= np.sum(yB)
    bTot = np.add(bTot, yB * nCts * (epb/xpbF))
    bTotCts += nCts
    nTot += nCts

    nCts = fitVals["amp-pol"][0]
    yB = np.ones(len(xF))
    yB /= np.sum(yB)
    bTot = np.add(bTot, yB * nCts * (epb/xpbF))
    bTotCts += nCts
    nTot += nCts

    plt.plot(xF, bTot, c='m', lw=3, label="Expo Bkg, %.2f cts" % bTotCts)
    yF = np.add(yF, bTot)
    print("Total bkg cts:",bTotCts)

    # total model
    plt.plot(xF, yF, c='r', lw=3, label="Total Model. Peak Cts: %.2f (90%% CL)" % pCts)

    plt.xlabel("Energy (keV)", ha='right', x=1)
    plt.ylabel("Counts / %.2f keV" % epb, ha='right', y=1)
    plt.legend(loc=4, fontsize=12)
    plt.minorticks_on()
    plt.xlim(eLo, eHi)
    # plt.ylim(ymin=0)
    plt.tight_layout()
    # plt.show()
    plt.savefig("%s/plots/sf7-after-mpl-%dpks.pdf" % (dsi.latSWDir, nPks))


def getShiftProfile(eLo, eHi, epb):

    from ROOT import TFile, TCanvas
    from ROOT import RooStats as RS

    f = TFile("%s/data/fitWS-axShift-%dpks.root" % (dsi.latSWDir, nPks))
    fitWS = f.Get("fitWS")
    fData = fitWS.allData().front()
    hitE = fitWS.var("trapENFCal")
    model = fitWS.pdf("model")
    fitRes = fitWS.allGenericObjects().front()
    fPars = fitRes.floatParsFinal()
    nPars = fPars.getSize()
    minNLL = fitRes.minNll()

    # === get fit results: {name : [nCts, err]} ===
    fitVals = {}
    for i in range(nPars):
        fp = fitRes.floatParsFinal()
        name = fp.at(i).GetName()
        fitVal, fitErr = fp.at(i).getValV(), fp.at(i).getError()
        if "amp" in name:
            fitVals[name] = [fitVal, fitErr, name.split('-')[1]]
    for f in fitVals:
        print(f, fitVals[f])


    if floatWidth:
        tOut = TFile("%s/data/rs-plc-shift-%dpks-eLo%.1f-float.root" % (dsi.latSWDir, nPks, eLo), "RECREATE")
        print("Warning, saving the profile for a floating signal width")
    else:
        tOut = TFile("%s/data/rs-plc-shift-%dpks-eLo%.1f.root" % (dsi.latSWDir, nPks, eLo), "RECREATE")


    start = time.clock()

    name = "amp-sPk"
    fitVal = fitVals[name][0]
    thisVar = fitWS.var(name)

    pCL = 0.9
    plc = RS.ProfileLikelihoodCalculator(fData, model, ROOT.RooArgSet(thisVar))
    plc.SetConfidenceLevel(0.90)
    interval = plc.GetInterval()
    lower = interval.LowerLimit(thisVar)
    upper = interval.UpperLimit(thisVar)
    plot = RS.LikelihoodIntervalPlot(interval)
    plot.SetNPoints(100)
    plot.Draw("tf1")

    # from ROOT import TCanvas
    # c = TCanvas("c","c",800,600)
    # plot.Draw("tf1")
    # c.Print("%s/plots/profile-sPk-test.pdf" % dsi.latSWDir)

    pName = "hP"
    hProfile = plot.GetPlottedObject()
    hProfile.SetName(pName)
    hProfile.SetTitle("PL %.2f  %s  lo %.3f  mid %.3f  hi %.3f" % (pCL, name, lower, fitVal, upper))
    hProfile.Write()
    print(hProfile.GetTitle())

    tOut.Close()

    # print("elapsed:",time.clock()-start)


def plotShiftProfile(eLo, eHi, epb, makePlots=False):
    from ROOT import TFile

    # use the results from GetPeakFluxPY
    f = np.load("%s/data/sf7-pkFluxes.npz" % dsi.latSWDir)
    axFluxes = f['arr_0'].item()
    nFlux = np.sum([axFluxes[str("%.2f"%pk)] for pk in axPeaks])
    Nexp_P = detExp * nFlux

    tf = TFile("%s/data/rs-plc-shift-%dpks-eLo%.1f.root" % (dsi.latSWDir, nPks, eLo))
    hP = tf.Get("hP")
    hT = hP.GetTitle().split()
    pars = []
    for v in hT:
        try: pars.append(float(v))
        except: pass
    pcl, intLo, bestFit, intHi = pars

    xP, yP, xpb = wl.npTH1D(hP)
    xP, yP = xP[1:] - xpb/2, yP[1:]

    # sanity check 2 - get the maximum of the confint, same as in RooStats:
    # Double_t Yat_Xmax = 0.5*ROOT::Math::chisquared_quantile(fInterval->ConfidenceLevel(),1);
    from scipy.stats import chi2
    df = 1 # this is an assumption
    cx = np.linspace(chi2.ppf(0.85, df), chi2.ppf(0.92, df), 100)
    cy = chi2.cdf(cx, df)
    chi2max = cx[np.where(cy>=0.9)][0] * 0.5
    pyCts = xP[np.where(yP>chi2max)][0]

    # Nobs = intHi
    Nobs = pyCts
    gae = np.power(Nobs/Nexp_P, 1/4)
    print("Expo (kg-d) %.2f  Nobs: %.2f  Nexp_P %.2e  eLo %.1f  eHi %.1f  g_ae U.L. %.4e" % (detExp, Nobs, Nexp_P, eLo, eHi, gae))

    if makePlots:
        plt.close()
        plt.axhline(chi2max, c='m', lw=2, label=r"$\chi^2\mathregular{/2\ (90\%\ C.L.)}}$")

        plt.plot(xP, yP, '-b', lw=4, label=r"w/ Measured Eff. $\mathregular{g_{ae} \leq}$ %.2e" % gae)
        # plt.plot([intHi,intHi],[0,chi2max], '-b', lw=2, alpha=0.5) # this isn't always at the 90% value, the PLC sucks
        plt.plot([intLo,intLo],[0,chi2max], '-b', lw=2, alpha=0.5)
        plt.plot([pyCts,pyCts],[0,chi2max], '-b', lw=2, alpha=0.5)

        plt.xlabel(r"$\mathregular{N_{obs}}$", ha='right', x=1)
        plt.ylabel(r"-log $\mathregular{\lambda(\mu_{axion})}$", ha='right', y=1)
        plt.legend(loc=2)
        plt.ylim(0, chi2max*3)
        # plt.xlim(500, 1000)
        plt.tight_layout()
        # plt.show()
        plt.savefig("%s/plots/sf-axion-profile-shift-%dpks.pdf" % (dsi.latSWDir, nPks))


def getNLLOffset():

    n1 = 2 * np.log(-3539.630380/-3539.4261122)
    n2 = 2 * np.log(-8437.22376865 / -8437.2228011)
    n3 = 2 * np.log(-15069.657386 / -15069.6567999)
    n4 = 2 * np.log(-22784.546954 / -22784.5460849)

    print("n1:",n1)
    print("n2:",n2)
    print("n3:",n3)
    print("n4:",n4)

def combineProfiles():

    from ROOT import TFile
    from scipy.stats import chi2

    print("i;m here")

    # expected axion counts, gae=1
    f = np.load("%s/data/sf7-pkFluxes.npz" % dsi.latSWDir)
    axFluxes = f['arr_0'].item()

    # 1. profiles for each number of peaks, with g_ae in legend, and 2.0 eLo

    # === 1. profiles for each number of peaks, with g_ae in legend, and 2.0 eLo ===

    eLo = 2.0
    cols = ['k','g','b','r']
    for i in range(1,5):

        tf = TFile("%s/data/rs-plc-shift-%dpks-eLo%.1f.root" % (dsi.latSWDir, i, eLo))
        hP = tf.Get("hP")
        hT = hP.GetTitle().split()
        pars = []
        for v in hT:
            try: pars.append(float(v))
            except: pass
        pcl, intLo, bestFit, intHi = pars
        print(i,"---",pars)

        xP, yP, xpb = wl.npTH1D(hP)
        xP, yP = xP[1:] - xpb/2, yP[1:]
        tf.Close()

        # observed counts
        df = 1
        cx = np.linspace(chi2.ppf(0.85, df), chi2.ppf(0.92, df), 100)
        cy = chi2.cdf(cx, df)
        chi2max = cx[np.where(cy>=0.9)][0] * 0.5
        pyCts = xP[np.where(yP>chi2max)][0]-xpb/2 # sometimes the interval reported by the PLC isn't the true 90% value in the curve??
        # Nobs = intHi
        Nobs = pyCts

        # expected counts
        pk0 = 4-i
        axPeaks = [sigVals[name][0] for name in pkModel]
        axPeaks = axPeaks[pk0:]
        nFlux = np.sum([axFluxes[str("%.2f"%pk)] for pk in axPeaks])
        Nexp_P = detExp * nFlux

        gae = np.power(Nobs/Nexp_P, 1/4)
        # print("Expo (kg-d) %.2f  Nobs: %.2f  Nexp_P %.2e  eLo %.1f  eHi %.1f  g_ae U.L. %.4e" % (detExp, Nobs, Nexp_P, eLo, eHi, gae))

        plt.axhline(chi2max, c='m', lw=2, label=r"$\chi^2\mathregular{/2\ (90\%\ C.L.)}}$")

        plt.plot(xP, yP, '-', c=cols[i-1], lw=2, label=r"nPks=%d: $\mathregular{g_{ae} \leq}$ %.2e" % (i,gae))
        plt.plot([pyCts,pyCts],[0,chi2max], '-', c=cols[i-1], lw=2, alpha=0.5)

    plt.xlabel(r"$\mathregular{N_{obs}}$", ha='right', x=1)
    plt.ylabel(r"-log $\mathregular{\lambda(\mu_{axion})}$", ha='right', y=1)
    plt.legend(loc=1, fontsize=14)
    plt.ylim(0, chi2max*3)
    # plt.xlim(500, 1000)
    plt.tight_layout()
    # plt.show()
    print("i saved")
    plt.savefig("%s/plots/sf-axion-profile-shift-allPks.pdf" % dsi.latSWDir)


    # === 2. profile for the nPks=4 curve, with a 1.9 and 2.0 eLo ===
    # ===    AND, profile for the n=4 curve, with a fixed and floating signal peak width

    plt.close()

    # df = 1
    # cx = np.linspace(chi2.ppf(0.85, df), chi2.ppf(0.92, df), 100)
    # cy = chi2.cdf(cx, df)
    # chi2max = cx[np.where(cy>=0.9)][0] * 0.5
    chi2max = 1.355
    plt.axhline(chi2max, c='m', lw=2, label=r"$\chi^2\mathregular{/2\ (90\%\ C.L.)}}$")

    nPks = 4
    pk0 = 4-nPks
    axPeaks = [sigVals[name][0] for name in pkModel]
    axPeaks = axPeaks[pk0:]
    nFlux = np.sum([axFluxes[str("%.2f"%pk)] for pk in axPeaks])
    Nexp_P = detExp * nFlux

    eLo = 2.0
    tf1 = TFile("%s/data/rs-plc-shift-%dpks-eLo%.1f.root" % (dsi.latSWDir, nPks, eLo))
    hP1 = tf1.Get("hP")
    hT1 = hP1.GetTitle().split()
    xP1, yP1, xpb1 = wl.npTH1D(hP1)
    xP1, yP1 = xP1[1:] - xpb1/2, yP1[1:]
    pyCts1 = xP1[np.where(yP1>chi2max)][0]-xpb1/2
    gae = np.power(pyCts1/Nexp_P, 1/4)
    tf1.Close()
    plt.plot(xP1, yP1, '-b', lw=2, label=r"nPks: %d, eLo: %.1f keV, nCts: %.1f, $\mathregular{g_{ae} \leq}$ %.2e" % (nPks,eLo,pyCts1,gae))
    plt.plot([pyCts1,pyCts1],[0,chi2max], '-b', lw=2, alpha=0.5)

    eLo = 1.9
    tf2 = TFile("%s/data/rs-plc-shift-%dpks-eLo%.1f.root" % (dsi.latSWDir, nPks, eLo))
    hP2 = tf2.Get("hP")
    hT12 = hP2.GetTitle().split()
    xP2, yP2, xpb2 = wl.npTH1D(hP2)
    xP2, yP2 = xP2[1:] - xpb2/2, yP2[1:]
    pyCts2 = xP2[np.where(yP2>chi2max)][0]-xpb2/2
    gae = np.power(pyCts2/Nexp_P, 1/4)
    tf2.Close()
    plt.plot(xP2, yP2, '-r', lw=2, label=r"nPks: %d, eLo: %.1f keV, nCts: %.1f, $\mathregular{g_{ae} \leq}$ %.2e" % (nPks,eLo,pyCts2,gae))
    plt.plot([pyCts2,pyCts2],[0,chi2max], '-r', lw=2, alpha=0.5)

    eLo = 2.0
    tf3 = TFile("%s/data/rs-plc-shift-%dpks-eLo%.1f-float.root" % (dsi.latSWDir, nPks, eLo))
    hP3 = tf3.Get("hP")
    hT3 = hP3.GetTitle().split()
    xP3, yP3, xpb3 = wl.npTH1D(hP3)
    xP3, yP3 = xP3[1:] - xpb3/2, yP3[1:]
    pyCts3 = xP3[np.where(yP3>chi2max)][0]-xpb3/2
    gae = np.power(pyCts3/Nexp_P, 1/4)
    tf1.Close()
    plt.plot(xP3, yP3, '-g', lw=2, label=r"Floating $\sigma$, nPks: %d, eLo: %.1f keV, nCts: %.1f, $\mathregular{g_{ae} \leq}$ %.2e" % (nPks,eLo,pyCts3,gae))
    plt.plot([pyCts3,pyCts3],[0,chi2max], '-g', lw=2, alpha=0.5)
    plt.plot([0,0],[0,chi2max], '-g', lw=2, alpha=0.5)

    plt.xlabel(r"$\mathregular{N_{obs}}$", ha='right', x=1)
    plt.ylabel(r"-log $\mathregular{\lambda(\mu_{axion})}$", ha='right', y=1)
    plt.legend(loc=1, fontsize=13)
    plt.ylim(0, chi2max*3)
    # plt.xlim(500, 1000)
    plt.tight_layout()
    # plt.show()
    plt.savefig("%s/plots/sf-axion-profile-shift-systematics.pdf" % dsi.latSWDir)


    # plt.close()
    # plt.axhline(chi2max, c='m', lw=2, label=r"$\chi^2\mathregular{/2\ (90\%\ C.L.)}}$")
    #
    # plt.xlabel(r"$\mathregular{N_{obs}}$", ha='right', x=1)
    # plt.ylabel(r"-log $\mathregular{\lambda(\mu_{axion})}$", ha='right', y=1)
    # plt.legend(loc=1, fontsize=14)
    # plt.ylim(0, chi2max*3)
    # # plt.xlim(500, 1000)
    # plt.tight_layout()
    #
    # # plt.show()
    # plt.savefig("%s/plots/sf-axion-profile-shift-widthSystematic.pdf" % dsi.latSWDir)




if __name__=="__main__":
    main()