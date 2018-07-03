#!/usr/bin/env python3
import sys, warnings, time
import numpy as np
from scipy.interpolate import spline
import waveLibs as wl
import dsi
import matplotlib.pyplot as plt
plt.style.use('../pltReports.mplstyle')
sys.argv.append("-b")
import ROOT
from ROOT import RooFit as RF

eLo, eHi, epb = 1.5, 3, 0.03
pLo, pHi, ppb = 0, 30, 0.03 # fitPeaks doesn't work for ppb!=0.03, the fit parameters are optimized for it
dsList = [1,2,3,4,"5A","5B","5C"]
enr = True
opt = "enr"

nB = int((eHi-eLo)/epb)
nBP = int((pHi-pLo)/ppb)

pkModel = ["axSi_a","axSi_b","axS_a","axS_b"]
bkModel = ["bkg1","bkg2"]
sigModel = pkModel + bkModel

dMu = 0.05 # was using 0.005
sigVals = {
    # expo bkg:      amp, lo, hi,  tau, lo, hi
    "bkg1":     [-1, 0.9, 0, 5,    -0.8, -5, 0],
    "bkg2":     [-1, 0.9, 0, 5,    -0.9, -5, 0],
    # peaks:       mu,  amp,  lo,   hi     sig,   lo,  hi
    "axSi_a":   [1.86,  0.1,   0, 0.15,   0.02, 0.01,  0.03],  # super tight constraints, forcing the curve to match
    "axSi_b":   [2.00,  0.004, 0, 0.015,   0.01, 0.005, 0.02],
    "axS_a":    [2.45,  0.1,   0, 0.13,   0.025, 0.015, 0.025],
    "axS_b":    [2.62,  0.004, 0, 0.015,   0.01, 0.005, 0.02],
    # "axSi_a":   [1.88,  0.1,   0, 1.,   0.02, 0.01,  0.1],  # loose constraints, just to see what roofit wants to do
    # "axSi_b":   [2.00,  0.004, 0, 0.015,   0.01, 0.005, 0.02],  # these totally suck
    # "axS_a":    [2.45,  0.1,   0, 1.,   0.025, 0.02, 0.1],
    # "axS_b":    [2.62,  0.004, 0, 0.015,   0.01, 0.005, 0.02],  # the _b peaks just always go the maximum allowed
    }

from ROOT import gROOT
gROOT.ProcessLine("gErrorIgnoreLevel = 3001;")
gROOT.ProcessLine("RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);")

def main():

    # initialize()
    getUnscaledPDFs(makePlots=False)
    # plotPDF()
    # fitPeaks()
    # getPeakFluxRF() # do this 2 ways, since the roofit model has to be super constrained to look ok,
    getPeakFluxPY() # does a sideband analysis and gaussian peak fitting


def initialize(makePlots=False):
    """ For this dsList, load the efficiency curves and exposure into globals.
    Also tweak the bkgModelPeaks list to only include peaks in the energy range.
    """
    global effLim, effMax, xEff, detEff, dsExpo, detExp, bkgModelPeaks

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


def getUnscaledPDFs(ma=0, makePlots=False):
    """ Generate a set of TH1D's to be turned into RooDataHist objects.
    Be careful they have the same axis limits and binning as the RooDataSet.
    Takes axion mass (in keV) as a parameter.
    """
    from ROOT import TFile, TH1D, gROOT

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

    plt.close()
    xR1, hR1, _ = wl.npTH1D(tf.Get("h4"))
    plt.plot(xR1, hR1, ls='steps', c='b', lw=3, label=r"$\Phi_a$, %.2f keV/bin, $\mathregular{g_{ae}=1}$" % ppb)

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


def getSigma(E, opt=""):
    """ Get the MJ energy resolution.
    If multiple DS are selected, weight the curve by DS exposure.
    Uses the global variable 'dsList'.
    """

    # HG resolutions, from the energy unidoc.
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
        return np.sqrt(p[0]**2 + p[1]**2 * E + p[2]**2 * E**2)
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
        return sig


def fitPeaks():
    from ROOT import TFile, TH1D, TCanvas, TLegend, gStyle

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
        bkN = ROOT.RooRealVar("amp-"+name,"amp-"+name, sigVals[name][1], sigVals[name][2], sigVals[name][3])
        bkT = ROOT.RooRealVar("tau-"+name,"tau-"+name, sigVals[name][4], sigVals[name][5], sigVals[name][6])
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

    fSpec = hitE.frame(RF.Range(eLo, eHi)) #, RF.Bins(nB))
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


def getPeakFluxPY():

    from ROOT import TFile
    from scipy.optimize import curve_fit

    # tf = TFile("%s/data/specPDFs-sf7.root" % )



    from ROOT import TFile
    from scipy.optimize import curve_fit
    opt = "enr" if enr else "nat"

    tf = TFile("%s/data/specPDFs-sf7.root" % dsi.latSWDir)
    xR, hR, _ = wl.npTH1D(tf.Get("h4"))
    # plt.plot(xR, hR, ls='steps', c='b', lw=3, label=r"$\Phi_a$, %.2f keV/bin, $\mathregular{g_{ae}=1}$" % ppb)

    aLo, aHi, apb = 1.5, 3.0, 0.03

    idx = np.where((xR>=aLo)&(xR<=aHi))
    # plt.plot(xR[idx], hR[idx] / 1e42, alpha=0.3, ls='steps-mid', c='k', lw=3, label=r"$\Phi_a$, %.2f keV/bin, $\mathregular{g_{ae}=1}$" % ppb)

    pkE = [1.85, 2.00, 2.45, 2.62]
    pkX = [1.85 - 3*apb, 2.00 + 3*apb, 2.45-3*apb, 2.62+3*apb] # x for eXclude

    for i in range(len(pkE)):
        # plt.axvline(pkE[i], c='g')
        plt.axvline(pkX[i], c='r')

    # get sidebands
    xR, hR = xR[idx], hR[idx]/1e42

    idxSB = np.where( (xR<pkX[0]) | ((xR > pkX[1])&(xR < pkX[2]) | (xR>pkX[3])))
    # plt.plot(xR[idxSB], hR[idxSB], ".", c='b', lw=3) # this is ugly, don't show it, but it worked.

    xF = np.arange(aLo, aHi, 0.01)

    def expFunc(x, a, mu, tau, b):
        return a * np.exp((x-mu)/tau) + b

    def doubleExpFunc(x, a1, tau1, a2, tau2):
        return a1 * np.exp(x/tau1) + a2 * np.exp(x/tau2)

    def pol1(x, a, b):
        return a*x + b

    def pol2(x, a, b, c):
        return a*x**2 + b*x + c

    init=(1, 2, -0.8, 0.2)
    popt, pcov = curve_fit(expFunc, xR[idxSB], hR[idxSB], p0=init)
    # plt.plot(xF, expFunc(xF, *popt), c='g')


    hSubt = hR - expFunc(xR, *popt)

    def gauss_function(x, a, x0, sigma):
        return a * np.exp(-(x - x0)**2 / (2 * sigma**2))

    def twoGauss(x, a1, mu1, sig1, a2, mu2, sig2):
        return gauss_function(x, a1, mu1, sig1) + gauss_function(x, a2, mu2, sig2)
        # return a1 * np.exp(-(x - mu1)**2 / (2 * sig1**2)) + a2 * np.exp(-(x - mu2)**2 / (2 * sig2**2))

    plt.axvline(1.85, c='g', lw=2, alpha=0.5, label=r"1.85 keV, Si ($\mathregular{K_{\alpha 1,\alpha 2}}$)")
    # plt.axvline(2.00, c='m', lw=2, alpha=0.5, label=r"2.00 keV, Si ($\mathregular{K_{\beta}}$)")
    # plt.axvline(2.45, c='r', lw=2, alpha=0.5, label=r"2.45 keV, S ($\mathregular{K_{\alpha 1,\alpha 2}}$)")
    # plt.axvline(2.62, c='k', lw=2, alpha=0.5, label=r"2.62 keV, S ($\mathregular{K_{\beta}}$)")

    # plt.step(xR, hSubt, c='b')

    # init = (0.5, 1.85, 0.01)
    # bnds = ((0,1.8499,0),(np.inf,1.8501,np.inf))
    # plt.plot(xF, gauss_function(xF, *init), c='g')
    # popt, pcov = curve_fit(gauss_function, xR, hSubt, p0=init, bounds=bnds)
    # plt.plot(xF, gauss_function(xF, *popt), c='g')

    init = (0.5, 1.85, 0.01, 0.3, 2.00, 0.01)
    bnds = ((0,1.8499,0, 0,1.9999,0),(np.inf,1.8501,0.03, np.inf,2.0001,0.03))
    popt, pcov = curve_fit(twoGauss, xR, hSubt, p0=init, bounds=bnds)
    # plt.plot(xF, twoGauss(xF, *popt), c='g')

    plt.step(xR, hSubt * 1e42, c='b')
    plt.plot(xF, twoGauss(xF, *popt) * 1e42, c='g')






    # plt.ylim(min(hR/1e42), max(hR/1e42))


    plt.show()


def loadShiftedDataMJD():
    """ Sometimes this segfaults.  No idea why.  Just re-run and it's fine. """
    from ROOT import TFile, TTree, TList
    from array import array

    fName = "%s/data/latDS%s.root" % (dsi.latSWDir, ''.join([str(d) for d in dsList]))
    f1 = TFile(fName)
    tree0 = f1.Get("skimTree")

    # axPeaks = [1.739, 1.836, 2.307, 2.464]
    axPeaks = [1.85, 2.00, 2.45, 2.62] # shifted to match Redondo paper

    ene0, wt0, shift = array('d',[0.]), array('d',[0.]), array('d',[0.])

    tree0.SetBranchAddress("trapENFCal",ene0)
    tree0.SetBranchAddress("weight",wt0)

    tList = TList()
    tVec = [0.,0.,0.,0.]
    tmp = TFile("%s/data/latDS%s_shifted.root" % (dsi.latSWDir, ''.join([str(d) for d in dsList])),"RECREATE")

    for i in range(len(axPeaks)):
        ene, wt = array('d',[0.]), array('d',[0.])
        shift = axPeaks[len(axPeaks)-1] - axPeaks[i]
        tVec[i] = TTree("t%d" % i, "t%d" %i)
        tVec[i].Branch("energy_keV",ene,"energy_keV/D")
        tVec[i].Branch("weight",wt,"weight/D")
        for j in range(tree0.GetEntries()):
            tree0.GetEntry(j)
            ene[0], wt[0] = ene0[0], wt0[0]
            ene[0] += shift
            tVec[i].Fill()

        print("axPeaks: %d  tree %d - %d entries  shift %.3f keV" % (len(axPeaks),len(tVec),tVec[i].GetEntries(),shift))
        tList.Add(tVec[i])
        tVec[i].Write()

    tree1 = ROOT.TTree.MergeTrees(tList)
    tree1.SetName("mergeTree")
    tree1.SetTitle("mergeTree")
    print("Trees merged, with %d entries total." % tree1.GetEntries())
    tree1.Write()


def plotShiftedData():

    # TODO: load the trees here

    # create diagnostic shifted & merged spectrum
    xLo, xHi, xpb = 1, 5, 0.1
    nBX = int((xHi-xLo)/xpb)
    h0 = wl.H1D(tree0, nB, xLo, xHi, "trapENFCal","")
    x, y, xpb = wl.npTH1D(h0)
    # plt.step(x, y)
    # plt.show()

    # # -- create diagnostic shifted & merged spectrum --
    # malbekExposure = 89.5
    # binSize = 0.04
    # eLo, eHi = 0.8, 5.
    # nBins = int((eHi-eLo)/binSize + 0.5)
    #
    # h0 = H1D(tree0,nBins,eLo,eHi,"energy_keV","weight")
    # h0.SetLineColor(ROOT.kBlack)
    # h0.SetLineWidth(2)
    # h0.Scale(1./malbekExposure)
    #
    # h1 = H1D(tVec[0],nBins,eLo,eHi,"energy_keV","weight")
    # h1.SetLineColorAlpha(ROOT.kBlue,0.5)
    # h1.Scale(1./malbekExposure)
    #
    # h2 = H1D(tVec[1],nBins,eLo,eHi,"energy_keV","weight")
    # h2.SetLineColorAlpha(ROOT.kGreen,0.5)
    # h2.Scale(1./malbekExposure)
    #
    # h3 = H1D(tVec[2],nBins,eLo,eHi,"energy_keV","weight")
    # h3.SetLineColorAlpha(ROOT.kMagenta,0.5)
    # h3.Scale(1./malbekExposure)
    #
    # h4 = H1D(tVec[3],nBins,eLo,eHi,"energy_keV","weight")
    # h4.SetLineColor(ROOT.kCyan)
    # h4.Scale(1./malbekExposure)
    #
    # hAdd = TH1D("hAdd","",nBins,eLo,eHi)
    # hAdd.Add(h1)
    # hAdd.Add(h2)
    # hAdd.Add(h3)
    # hAdd.Add(h4)
    # hAdd.SetLineColor(ROOT.kBlack)
    # hAdd.SetLineWidth(2)
    # hAdd.SetMinimum(0)
    # hAdd.SetMaximum(2)
    # ymax = hAdd.GetMaximum()
    # hAdd.GetXaxis().SetTitle("Energy (keV)")
    # hAdd.GetYaxis().SetTitle("Counts / kg-d")
    #
    # c = TCanvas("c","Bob Ross's Canvas",1100,800)
    # hAdd.Draw("hist")
    # h0.Draw("hist same")
    # h1.Draw("hist same")
    # h2.Draw("hist same")
    # h3.Draw("hist same")
    # h4.Draw("hist same")
    #
    # # Draw lines around the gaussian fit ROI, using 2 sigma of the energy resolution
    # fitWin = 3 * getSigma(axPeaks[3])
    # print("Malbek resolution at %.2f keV is %.2f" % (axPeaks[3], getSigma(axPeaks[3])))
    # print("Set 3-sigma fit region to %.2f - %.2f" % (axPeaks[3] - fitWin, axPeaks[3] + fitWin))
    #
    # l1 = ROOT.TLine(axPeaks[3],0.,axPeaks[3],ymax)
    # l1.SetLineColor(ROOT.kRed)
    # l1.SetLineWidth(2)
    # l1.Draw("same")
    #
    # l2 = ROOT.TLine(axPeaks[3]+fitWin,0.,axPeaks[3]+fitWin,ymax)
    # l2.SetLineColorAlpha(ROOT.kRed, 0.5)
    # l2.SetLineWidth(2)
    # l2.Draw("same")
    #
    # l3 = ROOT.TLine(axPeaks[3]-fitWin,0.,axPeaks[3]-fitWin,ymax)
    # l3.SetLineColorAlpha(ROOT.kRed, 0.5)
    # l3.SetLineWidth(2)
    # l3.Draw("same")
    #
    # leg = ROOT.TLegend(0.6,0.6,0.85,0.85)
    # leg.AddEntry(hAdd,"Shifted+Summed","l")
    # leg.AddEntry(h4,"peak-%.3f (unshifted)" % (axPeaks[3]),"l")
    # leg.AddEntry(h3,"peak-%.3f" % (axPeaks[2]),"l")
    # leg.AddEntry(h2,"peak-%.3f" % (axPeaks[1]),"l")
    # leg.AddEntry(h1,"peak-%.3f" % (axPeaks[0]),"l")
    # leg.AddEntry(l1,"3-sigma fit region","l")
    # leg.Draw("same")
    #
    # c.Print("./plots/shiftSpec.pdf")


if __name__=="__main__":
    main()