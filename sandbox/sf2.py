#!/usr/bin/env python3
import sys
import warnings

import matplotlib.pyplot as plt
plt.style.use('./pltReports.mplstyle')

import numpy as np
import ROOT
from ROOT import RooFit as RF
from ROOT import RooStats as RS

import waveLibs as wl
import dsi

# from ROOT import gROOT
# gROOT.ProcessLine("gErrorIgnoreLevel = 3001;") # suppress ROOT messages
# gROOT.ProcessLine("RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);")
# gROOT.ProcessLine("RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);")
# ROOT.Math.MinimizerOptions.SetDefaultMinimizer("Minuit") # the PLC doesn't work w/ minuit2

dsList = [0,1,2,3,4,"5A","5B","5C"]
eLo, eHi, epb = 1, 50, 0.2
enr = False

# must be specified here to be included in the fit
bkgModel = ["68GeK","55Fe","65ZnK","trit","flat"]#,"Pb210"]#,"axion"]

# list of all available peaks (could add fit constraints here)
pkList = {
    "68GeK": 10.37,
    "68Ga": 9.66,
    "65ZnK": 8.98,
    "55Fe": 6.54,
    "54Mn": 5.99,
    "49V": 4.97,
    "68GeL":1.29,
    "65ZnL":1.10
    }

# make efficiency and exposure global
xEff, detEff, effLim, effMax = -1, -1, -1, -1
detExp, detExpo = -1, -1 # total exposure, and by-detector exposure

def main(argv):

    # global eLo, eHi, epb, dsList, enr
    # eLo, eHi, epb, enr = ... # adjust the global values here
    # expoDict = {} todo, make exposures accessible
    initialize()

    # === routines ===
    # loadDataMJD()
    # getUnscaledPDFs(makePlots=False)
    # scalePDFs(eff=True, makePlots=True)
    # plotPDFs()
    # axionPeaks()
    # peakPDF(pkList["68GeK"], getSigma(pkList["68GeK"]))
    # runFit()
    # compareData()
    plotFit()
    # plotFitRF()


def initialize():
    """ For this dsList, load the efficiency curves and exposure into globals """
    global effLim, effMax, xEff, detEff

    # load efficiency correction
    f1 = np.load('./data/lat-expo-efficiency-all-e95.npz')
    xEff = f1['arr_0']
    totEnrEff, totNatEff = f1['arr_1'].item(), f1['arr_2'].item()
    detEff = np.zeros(len(xEff))
    for ds in dsList:
        if enr: detEff += totEnrEff[ds]
        else: detEff += totNatEff[ds]

    # load exposure
    f2 = np.load("./data/expo-totals-e95.npz")
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


def getUnscaledPDFs(ma=0, makePlots=False):
    """
    Generate a set of TH1D's to be turned into RooDataHist objects.
    Takes axion mass (in keV) as a parameter.
    - Binning is in 0.05 keV intervals - hopefully that's fine enough.
    TODO:
    - BDM/ALP PDF. Get from Kris's analysis in GAT
    """
    from ROOT import TFile, TH1D, gROOT

    # output files
    rOut = "./data/specPDFs.root"
    tf = TFile(rOut,"RECREATE")
    td = gROOT.CurrentDirectory()

    # energy limits
    xLo, xHi, xpb = 0, 30, 0.05
    nB = int((xHi-xLo)/xpb)

    print("Generating unscaled PDFs, xLo %.1f  xHi %.1f  xpb %.2f: %s" % (xLo, xHi, xpb, rOut))

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
    with open("./data/redondoFlux.txt") as f1: # 23577 entries
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
    with open("./data/ge76peXS.txt") as f2: # 2499 entries, 0.01 kev intervals
        lines = f2.readlines()
        for line in lines:
            data = line.split()
            phoData.append([float(data[0]),float(data[1])])
    phoData = np.array(phoData)

    # === 3. tritium
    tritData = []
    with open("./data/TritiumSpectrum.txt") as f3: # 20000 entries
        lines = f3.readlines()[1:]
        for line in lines:
            data = line.split()
            conv = float(data[2]) # raw spectrum convolved w/ ge cross section
            if conv < 0: conv = 0.
            tritData.append([float(data[1]),conv])
    tritData = np.array(tritData)

    # NOTE: check sandbox/th1.py for examples of manually filling TH1D's and verifying wl.GetHisto and wl.npTH1D.

    # ROOT output
    h1 = TH1D("h1","photoelectric",nB,xLo,xHi)         # [cm^2 / kg]
    h2 = TH1D("h2","axioelectric",nB,xLo,xHi)          # [cm^2 / kg]
    h3 = TH1D("h3","axion flux, gae=1",nB,xLo,xHi)     # [cts / (keV cm^2 d)]
    h4 = TH1D("h4","convolved flux",nB,xLo,xHi)        # [cts / (keV d kg)]
    h5 = TH1D("h5","tritium",nB,xLo,xHi)               # [cts] (normalized to 1)

    # manually fill ROOT histos (don't normalize yet)
    for iB in range(nB+1):
        ctr = (iB + 0.5)*xpb + xLo
        bLo, bHi = ctr - xpb/2, ctr + xpb/2
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
    tf2 = TFile("./data/Pb210PDFs.root")
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
    fName = "./data/latDS%s.root" % ''.join([str(d) for d in dsList])
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
    return hc


def getEffCorrTH1D(h, xLo, xHi):
    """ Returns a copy of a ROOT TH1D, uses global xEff and detEff functions. """
    from ROOT import TH1D
    nB = h.GetNbinsX()
    hEff = TH1D(h.GetName()+"-e", h.GetTitle()+"-eff",nB, xLo, xHi)
    for i in range(nB+1):
        binE = h.GetXaxis().GetBinCenter(i)
        hBin = h.GetBinContent(i)
        idx = (np.abs(xEff-binE)).argmin()
        hBinC = hBin * np.interp(binE, xEff[idx:idx+1], detEff[idx:idx+1])
        hEff.SetBinContent(i, hBinC)
    return hEff


def scalePDFs(eff=False, makePlots=False):
    """ Create scaled (normalized to 1) PDFs (both ROOT and numpy output).
    RooFit does that automatically, but I also need the numpy histos to be correct.
    Optionally, apply the DS-specific efficiency correction, using the global variable 'dsList'.
    Be very careful when you apply efficiency -- RooFit will automatically normalize those too.
    """
    from ROOT import TFile, TH1D, gROOT

    rOut = "./data/scaledPDFs.root"
    tfOut = TFile(rOut,"RECREATE")
    td = gROOT.CurrentDirectory()

    npOut = "./data/scaledPDFs.npz"

    print("Generating scaled PDFs: %s" % (rOut))

    # load unscaled PDFs
    tf = TFile("./data/specPDFs.root")
    h4 = tf.Get('h4') # axion pdf
    h5 = tf.Get('h5') # tritium
    h6 = tf.Get('h6') # Pb210TDL
    h7 = tf.Get('h7') # Pb210

    x4, np4, xpb4 = wl.npTH1D(h4)
    x5, np5, xpb5 = wl.npTH1D(h5)
    x6, np6, xpb6 = wl.npTH1D(h6)
    x7, np7, xpb7 = wl.npTH1D(h7)

    # create a flat BG PDF over the current energy range
    nB = int((eHi-eLo)/0.05)
    h8 = TH1D("h8","flat BG",nB,eLo,eHi)
    for iB in range(nB+1):
        h8.SetBinContent(iB, 1) # the initial amplitude doesn't matter b/c we normalize to 1
    x8, np8, xpb8 = wl.npTH1D(h8)

    # normalize PDFs, in the energy range we're using for the fit (eLo, eHi)

    # root
    h4.Scale(1/h4.Integral(h4.FindBin(eLo), h4.FindBin(eHi), 'width')) # axion
    h5.Scale(1/h5.Integral(h5.FindBin(eLo), h5.FindBin(eHi), 'width')) # tritium
    h6.Scale(1/h6.Integral(h6.FindBin(eLo), h6.FindBin(eHi), 'width')) # pb210-tdl
    h7.Scale(1/h7.Integral(h7.FindBin(eLo), h7.FindBin(eHi), 'width')) # pb210
    h8.Scale(1/h8.Integral(h8.FindBin(eLo), h8.FindBin(eHi), 'width')) # flat bg

    # numpy
    np4n = np.divide(np4, np.sum(np4[np.where((x4 >= eLo) & (x4 <= eHi))] * xpb4))
    np5n = np.divide(np5, np.sum(np5[np.where((x5 >= eLo) & (x5 <= eHi))] * xpb5))
    np6n = np.divide(np6, np.sum(np6[np.where((x6 >= eLo) & (x6 <= eHi))] * xpb6))
    np7n = np.divide(np7, np.sum(np7[np.where((x7 >= eLo) & (x7 <= eHi))] * xpb7))
    np8n = np.divide(np8, np.sum(np8[np.where((x8 >= eLo) & (x8 <= eHi))] * xpb8))

    if makePlots:

        # === 1. make sure I normalized correctly
        xR, yR, xpbR = wl.npTH1D(h4)
        plt.plot(x4, np4n, ls='steps', lw=3, label="numpy normalized")
        plt.plot(xR, yR, ls='steps', lw=2, label="root normalized")
        plt.xlabel("Energy (keV)", ha='right', x=1)
        plt.ylabel("Counts (norm)", ha='right', y=1)
        plt.legend()
        plt.tight_layout()
        plt.savefig("./plots/sf-axionPDF-norm.pdf")
        # return

    # apply efficiency correction (to the scaled TH1D)
    if eff:
        h4e = getEffCorrTH1D(h4, h4.GetXaxis().GetXmin(), h4.GetXaxis().GetXmax()) # axion
        h5e = getEffCorrTH1D(h5, h5.GetXaxis().GetXmin(), h5.GetXaxis().GetXmax()) # tritium
        h6e = getEffCorrTH1D(h6, h6.GetXaxis().GetXmin(), h6.GetXaxis().GetXmax()) # pb210-tdl
        h7e = getEffCorrTH1D(h7, h7.GetXaxis().GetXmin(), h7.GetXaxis().GetXmax()) # pb210
        h8e = getEffCorrTH1D(h8, h8.GetXaxis().GetXmin(), h8.GetXaxis().GetXmax()) # flat

        np4ne = getEffCorr(x4, np4n)
        np5ne = getEffCorr(x5, np5n)
        np6ne = getEffCorr(x6, np6n)
        np7ne = getEffCorr(x7, np7n)
        np8ne = getEffCorr(x8, np8n)

        if makePlots:

            # === 2. make sure the efficiency correction is right
            plt.close()
            xR, yR, xpbR = wl.npTH1D(h4e)
            plt.plot(xR, yR, ls='steps', label="root eff-corr")
            plt.plot(x4, np4ne, ls='steps', label="np eff-corr")
            plt.legend()
            # plt.show()
            plt.savefig("./plots/sf-effCorr-axion.pdf")

        # finally, normalize and scale the efficiency-corrected histos

        # this is what RooFit will essentially use (automatically)
        h4e.Scale(1/h4e.Integral(h4e.FindBin(eLo), h4e.FindBin(eHi), 'width')) # axion
        h5e.Scale(1/h5e.Integral(h5e.FindBin(eLo), h5e.FindBin(eHi), 'width')) # tritium
        h6e.Scale(1/h6e.Integral(h6e.FindBin(eLo), h6e.FindBin(eHi), 'width')) # pb210-tdl
        h7e.Scale(1/h7e.Integral(h7e.FindBin(eLo), h7e.FindBin(eHi), 'width')) # pb210
        h8e.Scale(1/h8e.Integral(h8e.FindBin(eLo), h8e.FindBin(eHi), 'width')) # flat

        np4nen = np.divide(np4ne, np.sum(np4ne[np.where((x4 >= eLo) & (x4 <= eHi))] * xpb4))
        np5nen = np.divide(np5ne, np.sum(np5ne[np.where((x5 >= eLo) & (x5 <= eHi))] * xpb5))
        np6nen = np.divide(np6ne, np.sum(np6ne[np.where((x6 >= eLo) & (x6 <= eHi))] * xpb6))
        np7nen = np.divide(np7ne, np.sum(np7ne[np.where((x7 >= eLo) & (x7 <= eHi))] * xpb7))
        np8nen = np.divide(np8ne, np.sum(np8ne[np.where((x8 >= eLo) & (x8 <= eHi))] * xpb8))

        if makePlots:

            # === 3. make sure I normalized the efficiency-corrected histos correctly
            plt.close()
            xR, yR, xpbR = wl.npTH1D(h4e)
            plt.plot(xR, yR, ls='steps', lw=3, label="root normalized")
            plt.plot(x4, np4nen, ls='steps', lw=2, label="numpy normalized")
            plt.xlabel("Energy (keV)", ha='right', x=1)
            plt.ylabel("Counts (norm)", ha='right', y=1)
            plt.legend()
            plt.tight_layout()
            # plt.show()
            plt.savefig("./plots/sf-effCorr-axion-norm.pdf")

            # === 4. plot the different corrections together
            plt.close()
            plt.plot(x4, np4n, ls='steps', lw=3, label="normalized")
            plt.plot(x4, np4ne, ls='steps', lw=3, label="norm, eff.corr")
            plt.plot(x4, np4nen, ls='steps', lw=3, label="norm, eff.corr, norm'd again")
            plt.xlabel("Energy (keV)", ha='right', x=1)
            plt.ylabel("Counts (norm)", ha='right', y=1)
            plt.xlim(0, 12)
            plt.legend()
            plt.tight_layout()
            # plt.show()
            plt.savefig("./plots/sf-effCorr-axion-norm-2.pdf")

            # === 5. plot the efficiency-corrected flat bg
            plt.close()
            plt.plot(x8, np8n, ls='steps', lw=3, label="Flat BG, normalized")
            plt.plot(x8, np8nen, ls='steps', lw=3, label="Flat BG, efficiency corrected, normalized")

            # reverse the normalization
            np8nenc = getEffCorr(x8, np8nen, inv=True)
            plt.plot(x8, np8nenc, ls='steps', lw=3, label="Flat BG, reverted efficiency correction")

            plt.xlabel("Energy (keV)", ha='right', x=1)
            plt.ylabel("Counts (norm)", ha='right', y=1)
            plt.ylim(0, 0.05)
            plt.legend()
            plt.tight_layout()
            # plt.show()
            plt.savefig("./plots/sf-effCorr-flat.pdf")


    # save ROOT output (used by the fitter)
    gROOT.cd(td.GetPath())
    if eff:
        h4e.Write()
        h5e.Write()
        h6e.Write()
        h7e.Write()
        h8e.Write()
    else:
        h4.Write()
        h5.Write()
        h6.Write()
        h7.Write()
        h8.Write()
    tfOut.Close()

    # save numpy output
    pdfRaw = {
        # this is the unscaled pdf
        "axion" : [x4, np4, xpb4],
        "trit"  : [x5, np5, xpb5],
        "Pb210" : [x6, np6, xpb6],
        "Pb210-TDL" : [x7, np7, xpb7],
        "flat" : [x8, np8, xpb8]
    }
    pdfNorm = {
        # this is the normalized pdf
        "axion" : [x4, np4n, xpb4],
        "trit"  : [x5, np5n, xpb5],
        "Pb210" : [x6, np6n, xpb6],
        "Pb210-TDL" : [x7, np7n, xpb7],
        "flat" : [x8, np8n, xpb8]
    }
    pdfEff = {
        # this is the efficiency corrected normalized pdf
        "axion" : [x4, np4ne, xpb4],
        "trit"  : [x5, np5ne, xpb5],
        "Pb210" : [x6, np6ne, xpb6],
        "Pb210-TDL" : [x7, np7ne, xpb7],
        "flat" : [x8, np8ne, xpb8]
    }
    pdfEffN = {
        # this is the efficiency corrected normalized pdf, normalized again
        # b/c that's what roofit does automatically in a fit
        "axion" : [x4, np4nen, xpb4],
        "trit"  : [x5, np5nen, xpb5],
        "Pb210" : [x6, np6nen, xpb6],
        "Pb210-TDL" : [x7, np7nen, xpb7],
        "flat" : [x8, np8nen, xpb8]
    }
    np.savez(npOut, pdfRaw, pdfNorm, pdfEff, pdfEffN)


def plotPDFs():
    """ Final consistency check on PDFs before we use them in the fitter. """
    from ROOT import TFile

    tfU = TFile("./data/specPDFs.root")
    tfS = TFile("./data/scaledPDFs.root")
    f = np.load("./data/scaledPDFs.npz")
    pdfRaw, pdfNorm, pdfEff, pdfEffN = f['arr_0'].item(), f['arr_1'].item(), f['arr_2'].item(), f['arr_3'].item()

    # === 1. Ge photoelectric XS
    plt.close()
    xR, hR, _ = wl.npTH1D(tfU.Get("h1"))
    plt.semilogy(xR, hR, ls='steps', c='b', lw=3, label=r"$\sigma_{ge}$")
    plt.xlabel("Energy (keV)", ha='right', x=1)
    plt.ylabel(r"$\mathregular{cm^2/kg}$", ha='right', y=1)
    plt.legend()
    plt.tight_layout()
    # plt.show()
    plt.savefig("./plots/sf-gexs.pdf")

    # === 2. axioelectric XS
    plt.close()
    xR, hR, _ = wl.npTH1D(tfU.Get("h2"))
    plt.plot(xR, hR, ls='steps', c='b', lw=3, label=r"$\sigma_{ae}$")
    plt.xlabel("Energy (keV)", ha='right', x=1)
    plt.ylabel(r"$\mathregular{cm^2/kg}$", ha='right', y=1)
    plt.legend()
    plt.tight_layout()
    # plt.gca().yaxis.set_label_coords(-0.06, 1)
    # plt.show()
    plt.savefig("./plots/sf-axs.pdf")

    # === 3. axion flux, gae=1
    plt.close()
    xR, hR, _ = wl.npTH1D(tfU.Get("h3"))
    plt.plot(xR, hR, ls='steps', c='b', lw=3, label=r"$\Phi_a$")
    plt.xlabel("Energy (keV)", ha='right', x=1)
    plt.ylabel(r"Flux / (keV cm${}^2$ d)]", ha='right', y=1)
    plt.xlim(0,12)
    plt.legend()
    plt.tight_layout()
    # plt.show()
    plt.savefig("./plots/sf-axFlux.pdf")

    # === 4. solar axion PDF, gae=1 -- this is what we integrate to get N_exp
    plt.close()
    xR1, hR1, _ = wl.npTH1D(tfU.Get("h4"))
    xr, hr, xbr = pdfRaw["axion"]
    plt.plot(xR1, hR1, ls='steps', c='b', lw=3, label=r"$\Phi_a$")
    plt.plot(xr, hr, ls='steps', c='r', lw=1)
    plt.xlabel("Energy (keV)", ha='right', x=1)
    plt.ylabel("Flux / (keV d kg)", ha='right', y=1)
    plt.xlim(0,10)
    plt.legend()
    plt.tight_layout()
    # plt.show()
    plt.savefig("./plots/sf-axPDF.pdf")

    # === 5. normalized, efficiency corrected solar axion PDF, gae=1
    plt.close()
    xR2, hR2, _ = wl.npTH1D(tfS.Get("h4-e")) # this is normalized w/ efficiency - not what we want for N_exp
    xn, hn, xbn = pdfNorm["axion"]
    xe, he, xbe = pdfEff["axion"]
    xen, hen, xben = pdfEffN["axion"]
    # plt.plot(xR2, hR2, ls='steps', c='r', lw=3, label="root, eff. corrected, norm'd")
    plt.plot(xn, hn, ls='steps', c='m', lw=2, label="Normalized PDF")
    # plt.plot(xe, he, ls='steps', c='g', lw=2, label="numpy, eff. corrected")
    plt.plot(xen, hen, ls='steps', c='b', lw=2, label="Normalized Eff.Corr. PDF")

    henc = getEffCorr(xen, hen, inv=True)
    plt.plot(xen, henc, ls='steps', c='g', lw=2, label="Final")

    plt.xlabel("Energy (keV)", ha='right', x=1)
    plt.ylabel("Counts (norm)", ha='right', y=1)
    plt.xlim(0,10)
    plt.legend()
    plt.tight_layout()
    # plt.show()
    plt.savefig("./plots/sf-axPDF-effnorm.pdf")
    # return

    # # === 5. tritium
    plt.close()
    xn, hn, xbn = pdfNorm["trit"]
    xen, hen, xben = pdfEffN["trit"]
    henc = getEffCorr(xen, hen, inv=True)

    plt.plot(xn, hn, ls='steps', c='r', lw=3, label="Tritium, norm. PDF")
    plt.plot(xen, hen, ls='steps', c='b', lw=3, label="w/ normalized efficiency correction")
    plt.plot(xen, henc, ls='steps', c='g', lw=2, label="Final")

    plt.xlabel("Energy (keV)", ha='right', x=1)
    plt.ylabel("Counts (norm)", ha='right', y=1)
    plt.legend()
    plt.xlim(0,20)
    plt.tight_layout()
    # plt.show()
    plt.savefig('./plots/sf-trit.pdf')

    # === 6. Pb210-TDL PDF
    plt.close()

    xn, hn, xbn = pdfNorm["Pb210-TDL"]
    xen, hen, xben = pdfEffN["Pb210-TDL"]
    henc = getEffCorr(xen, hen, inv=True)

    plt.plot(xn, hn, ls='steps', c='r', lw=2, label="Pb210TDL, norm. PDF")
    plt.plot(xen, hen, ls='steps', c='b', lw=3, label="w/ normalized efficiency correction")
    plt.plot(xen, henc, ls='steps', c='g', lw=2, label="Final")

    plt.xlabel("Energy (keV)", ha='right', x=1)
    plt.ylabel("Counts (norm)", ha='right', y=1)
    plt.legend(loc=2)
    plt.tight_layout()
    # plt.show()
    plt.savefig("./plots/sf-pb210.pdf")


def axionPeaks():
    """ placeholder: get the axion peak energies here. """
    from ROOT import TH1D

    gae = 1
    gRat = (gae / 5.11e-11)
    redondoScale = 1e19 * gRat**2 # convert table to [flux / (keV cm^2 d)]

    xLo, xHi, xpb = 0, 10, 0.1
    nB = int((xHi-xLo)/xpb)

    axData = []
    with open("./data/redondoFlux.txt") as f1: # 23577 entries
        lines = f1.readlines()[11:]
        for line in lines:
            data = line.split()
            axData.append([float(data[0]),float(data[1])])
    axData = np.array(axData)

    hAx = TH1D("hAx","axion flux, gae=1",nB,xLo,xHi)     # [cts / (keV cm^2 d)]
    for iB in range(nB+1):
        ctr = (iB + 0.5)*xpb + xLo
        bLo, bHi = ctr - xpb/2, ctr + xpb/2
        # print(iB, bLo, bHi)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore",category=RuntimeWarning)

            # axion flux [flux / (cm^2 d keV)]
            idx = np.where((axData[:,0] >= bLo) & (axData[:,0] < bHi))
            flux = np.mean(axData[idx][:,1]) * redondoScale
            if np.isnan(flux): flux = 0.
            hAx.SetBinContent(iB+1, flux)

    xA, hA, xpb2 = wl.npTH1D(hAx)
    if xpb!=xpb2:
        print("WTF, error")
        exit(1)

    # make sure the saved PDF matches
    f = np.load("./data/specPDFs.npz")
    pdfs = f['arr_0'].item()
    xP, aP, xpb3 = pdfs["axFlux"]

    plt.plot(axData[:,0], axData[:,1] * redondoScale, '.k', ms=1, alpha=0.3, label='Raw Flux Data')
    plt.plot(xA, hA, ls='steps', c='r', lw=1, label='Flux, %.1f keV bins' % xpb)
    plt.plot(xP, aP, ls='steps', c='m', lw=1)
    plt.axvline(1, c='m', lw=1, label="1.0 keV")

    # try to locate the shifted peaks
    msList = []
    msThresh = 10
    maxtab,_ = wl.peakdet(hA, msThresh)
    for iMax in range(len(maxtab)):
        idx = int(maxtab[iMax][0])
        if not 1 <= xA[idx] <= 7 : continue
        val = maxtab[iMax][1]
        xPk = xA[idx]-xpb/2
        msList.append(xPk)
        print("%d  idx %d  ene %.2f  val %.2e  thresh %.2f" % (iMax, idx, xPk, val, msThresh))
        plt.axvline(xPk, lw=1, c='g', label="%.2f" % xPk)

    plt.legend()
    plt.xlabel("Energy (keV)", ha='right', x=1)
    plt.tight_layout()
    # plt.show()
    plt.savefig('./plots/sf-axFlux-binned.pdf')


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


def gaus(x, a, x0, sigma):
    return a * np.exp(-(x - x0)**2 / (2 * sigma**2))


def peakPDF(pkE, sig):
    """ To apply the efficiency correction the same way as the continuous
    PDFs, let's make the peaks TH1D's that the fitter can quickly regenerate.
    Floats: mu, sig.
    """
    from ROOT import TH1D
    # pkE = 6.54
    # pkE = 3
    # sig = getSigma(pkE)

    xLo = pkE-2 if pkE-2 > eLo else eLo
    xHi = pkE+2 if pkE+2 < eHi else eHi
    xpb = 0.05

    # normalize based on global eLo, eHi
    x = np.arange(xLo, xHi, xpb)
    hP = gaus(x, 1, pkE, sig)
    hP = np.divide(hP, np.sum(hP[np.where((x >= eLo) & (x <= eHi))] * xpb))

    hSum = np.sum(hP[np.where((x >= eLo) & (x <= eHi))]) * xpb
    print(hSum)

    # efficiency correct and normalize
    hPe = getEffCorr(x, hP)
    hPe = np.divide(hPe, np.sum(hPe[np.where((x >= eLo) & (x <= eHi))] * xpb))

    # efficiency inflate
    hPf = getEffCorr(x, hPe, inv=True)

    # plt.plot(x, hP, ls='steps', lw=3, c='r', label="norm PDF")
    # plt.plot(x, hPe, ls='steps', lw=2, c='g', label="normalized eff. corr")
    # plt.plot(x, hPf, ls='steps', lw=2, c='m', label="Final")
    # plt.legend(loc=1)
    # plt.show()

    # make a TH1D
    nB = int((xHi-xLo)/xpb)
    hPR = TH1D("68GeK","68GeK",nB,xLo,xHi)
    for iB in range(nB):
        hPR.SetBinContent(iB, hPe[iB])

    # xR, hR, xpbR = wl.npTH1D(hPR)
    # plt.plot(xR, hR, ls='steps', lw=2, label='ROOT normalized eff corr')
    # plt.plot(x, hPe, ls='steps', lw=2, label='numpy normalized eff corr')
    # plt.legend(loc=1)
    # plt.show()

    return hPR


class pkModel:

    def __init__(self,fEnergy,name,ene,amp=0):
        """ TODO: add in constraints to the likelihood function instead of hard bounds
        https://root.cern.ch/root/html/tutorials/roofit/rf604_constraints.C.html
        """
        sig = getSigma(ene)
        muLo, muHi = ene * 0.99, ene * 1.01
        sgLo, sgHi = sig * 0.66, sig * 1.33
        ampLo, ampHi = -0.1, 300    # letting it go slightly negative helps the plc calculator
        if eHi > 10.: ampHi = 500. # 68GeK is around 1700, all the rest are under 500

        # specifying a min/max value allows the var to float
        if amp!=0:
            self.pkNum = ROOT.RooRealVar("amp-"+name, "amp-"+name, amp)
            self.pkMu = ROOT.RooRealVar("mu-"+name, "mu-"+name, ene)
            self.pkSig = ROOT.RooRealVar("sig-"+name, "sig-"+name, sig)
        else:
            self.pkNum = ROOT.RooRealVar("amp-"+name, "amp-"+name, ampLo, ampHi)
            self.pkMu = ROOT.RooRealVar("mu-"+name, "mu-"+name, ene, muLo, muHi)
            self.pkSig = ROOT.RooRealVar("sig-"+name, "sig-"+name, sig, sgLo, sgHi)

        # pdf and extended pdf
        self.pkGaus = ROOT.RooGaussian("gaus-"+name, "gaus-"+name, fEnergy, self.pkMu, self.pkSig)
        self.pkExt = ROOT.RooExtendPdf("ext-"+name, "ext-"+name, self.pkGaus, self.pkNum)

    def GetPkExt(self): return self.pkExt
    def GetAll(self): return self.pkNum, self.pkMu, self.pkSig, self.pkGaus, self.pkExt


def runFit():
    from ROOT import TFile

    # rescale the pdf's for this ds and this efficiency
    scalePDFs(eff=True, makePlots=False)

    # load data and create a workspace
    f1 = TFile("./data/latDS%s.root" % ''.join([str(d) for d in dsList]))
    t = f1.Get("skimTree")
    fEnergy = ROOT.RooRealVar("trapENFCal","Energy",eLo,eHi,"keV")
    fEnr = ROOT.RooRealVar("isEnr","isEnr",0,1,"")
    fWeight = ROOT.RooRealVar("weight","weight",1,10,"")

    tCut = "isEnr==1" if enr is True else "isEnr==0"
    fData = ROOT.RooDataSet("data", "data", t, ROOT.RooArgSet(fEnergy,fEnr), tCut)
    # fData = ROOT.RooDataSet("data", "data", t, ROOT.RooArgSet(fEnergy,fEnr,fWeight),"","weight")

    fitWorkspace = ROOT.RooWorkspace("fitWorkspace","Fit Workspace")
    getattr(fitWorkspace,'import')(fEnergy)
    getattr(fitWorkspace,'import')(fData)
    # getattr(fitWorkspace,'import')(fWeight)

    # ==== background model ====
    pdfList = ROOT.RooArgList("shapes")

    # load special PDFs
    f2 = TFile("./data/scaledPDFs.root")

    # ge68-k - with delta-E shift (broken)
    # pkE = pkList["68GeK"]
    # pkTH1D = peakPDF(pkE, getSigma(pkE))
    # dE = ROOT.RooRealVar("dE","dE",0., -0.2, 0.2)
    # xEf = ROOT.RooFormulaVar("xEf","@0-@1",ROOT.RooArgList(fEnergy, dE))
    # xE = fData.addColumn(xEf)
    # pk = ROOT.RooDataHist("dx","dx", ROOT.RooArgList(xE), RF.Import(pkTH1D))
    # fEnergy.setRange(eLo, eHi)
    # pkPdf = ROOT.RooHistPdf("pkPdf", "pkPdf", ROOT.RooArgSet(xE), pk, 2)
    # pkNum = ROOT.RooRealVar("amp-68GeK", "amp-68GeK", 10, -0.1, 1000)
    # pkExt = ROOT.RooExtendPdf("ext-68GeK", "ext-68GeK", pkPdf, pkNum)
    # if "68GeK" in bkgModel:
    #     pdfList.add(pkExt)

    # flat BG
    bkgTH1D = f2.Get("h8-e")
    bkgNum = ROOT.RooRealVar("amp-flat", "amp-flat", 10., -1, 20000.)
    bkgDataHist = ROOT.RooDataHist("bkg", "bkg", ROOT.RooArgList(fEnergy), RF.Import(bkgTH1D))
    fEnergy.setRange(eLo, eHi)
    bkgPdf = ROOT.RooHistPdf("bkgPdf", "bkgPdf", ROOT.RooArgSet(fEnergy), bkgDataHist, 2)
    bkgExt = ROOT.RooExtendPdf("ext-flat", "ext-flat", bkgPdf, bkgNum)
    if "flat" in bkgModel:
        pdfList.add(bkgExt)

    # tritium
    trTH1D = f2.Get("h5-e")
    trNum = ROOT.RooRealVar("amp-trit", "amp-trit", 1000., -1, 20000.)
    trDataHist = ROOT.RooDataHist("tr", "tr", ROOT.RooArgList(fEnergy), RF.Import(trTH1D))
    fEnergy.setRange(eLo, eHi)
    trPdf = ROOT.RooHistPdf("trPdf", "trPdf", ROOT.RooArgSet(fEnergy), trDataHist, 2)
    trExt = ROOT.RooExtendPdf("ext-trit", "ext-trit",trPdf,trNum)
    if "trit" in bkgModel:
        pdfList.add(trExt)

    # ge68-k - no delta-E shift
    pkE = pkList["68GeK"]
    pkTH1D = peakPDF(pkE, getSigma(pkE))
    pkDataHist = ROOT.RooDataHist("pk", "pk", ROOT.RooArgList(fEnergy), RF.Import(pkTH1D))
    fEnergy.setRange(eLo, eHi)
    pkPdf = ROOT.RooHistPdf("pkPdf", "pkPdf", ROOT.RooArgSet(fEnergy), pkDataHist, 2)
    pkNum = ROOT.RooRealVar("amp-68GeK", "amp-68GeK", 10, -0.1, 1000)
    pkExt = ROOT.RooExtendPdf("ext-68GeK", "ext-68GeK",pkPdf,pkNum)
    if "68GeK" in bkgModel:
        pdfList.add(pkExt)

    # gaussian peak list
    # (the trick seems to be that you can't overwrite the pkModel object in the loop.)
    # pks = []
    # for pk in pkList:
    #     if pk not in bkgModel: continue
    #     pks.append( pkModel(fEnergy, pk, pkList[pk]) )
    # for pk in pks:
    #     pdfList.add(pk.GetPkExt())

    # axion continuum
    # axTH1D = f2.Get("h4-e")
    # axNum = ROOT.RooRealVar("amp-axion", "amp-axion", 10., 0., 2000.)
    # # axNum = ROOT.RooRealVar("amp-axion","amp-axion",16.904) # can hardcode the profile upper limit
    # axDataHist = ROOT.RooDataHist("ax", "ax", ROOT.RooArgList(fEnergy), RF.Import(axTH1D))
    # fEnergy.setRange(eLo, eHi) # have to reset after loading a histo w/ different bounds
    # axPdf = ROOT.RooHistPdf("axPdf", "axPdf", ROOT.RooArgSet(fEnergy), axDataHist, 2)
    # axExt = ROOT.RooExtendPdf("ext-axion", "ext-axion", axPdf, axNum)
    # if "axion" in bkgModel:
    #     pdfList.add(axExt)

    # Pb210 continuum (w/ TDL)
    # pbTH1D = f2.Get("h6")
    # pbNum = ROOT.RooRealVar("amp-Pb210", "amp-Pb210", 10., 0., 1000.)
    # pbDataHist = ROOT.RooDataHist("pb", "pb", ROOT.RooArgList(fEnergy), RF.Import(pbTH1D))
    # fEnergy.setRange(eLo, eHi)
    # pbPdf = ROOT.RooHistPdf("pbPdf", "pbPdf", ROOT.RooArgSet(fEnergy), pbDataHist, 2)
    # pbExt = ROOT.RooExtendPdf("ext-Pb210", "ext-Pb210",pbPdf,pbNum)
    # if "Pb210" in bkgModel:
    #     pdfList.add(pbExt)

    # create total model pdf
    model = ROOT.RooAddPdf("model","total pdf",pdfList)

    # === efficiency ===
    # effFile = TFile("./data/lat-expo-efficiency.root")
    # effHist = effFile.Get("hDS5B_Norm_Enr")
    # x, y, xpb = wl.npTH1D(effHist)
    # plt.plot(x, y, ls='steps')
    # plt.show()
    # effRooHist = ROOT.RooDataHist("eff","Efficiency", ROOT.RooArgList(fEnergy), RF.Import(effHist))
    # fEnergy.setRange(eLo, eHi)
    # effPdf = ROOT.RooHistPdf("effPdf","effPdf", ROOT.RooArgSet(fEnergy), effRooHist, 0)
    # modelEff = ROOT.RooProdPdf("modelEff","model with efficiency", model, effPdf)

    # run fitter
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
    # getattr(fitWorkspace,'import')(modelEff)
    f2 = TFile("./data/fitWorkspace.root","RECREATE")
    fitWorkspace.Write()
    f2.Close()


def getRooArgDict(arglist):
    """ Convert a RooArgList into a python dict.
    NOTE: I'm doing it this way because this method doesn't work:
    pkValF = fitResult.floatParsFinal().find("pk_gaus") # can't cast to RooRealVar
    """
    pkVals = {}
    for i in range(arglist.getSize()):
        pkVals[ arglist.at(i).GetName() ] = arglist.at(i).getValV()
    return pkVals


def compareData():
    """ Make sure I can match the data points of a RooPlot """
    from ROOT import TFile, TCanvas

    # === matplotlib style
    tf2 = TFile("./data/latDS%s.root" % ''.join([str(d) for d in dsList]))
    tt = tf2.Get("skimTree")
    tCut = "isEnr" if enr else "!isEnr"
    tCut += " && trapENFCal >= %.1f && trapENFCal <= %.1f" % (eLo, eHi)
    n = tt.Draw("trapENFCal", tCut, "goff")
    hitE = tt.GetV1()
    hitE = [hitE[i] for i in range(n)]
    x, hData = wl.GetHisto(hitE, eLo, eHi, epb)
    hErr = np.asarray([np.sqrt(h) for h in hData]) # statistical error
    plt.figure(figsize=(7,5))
    # plt.plot(x, hData, ls='steps', c='b') # normal histo
    plt.errorbar(x, hData, yerr=hErr, c='k', ms=10, linewidth=0.8, fmt='.', capsize=2) # pretty convincing rooplot fake
    plt.xlabel("Energy (keV)", ha='right', x=1)
    plt.ylabel("Events / %.1f keV" % epb, ha='right', y=1)
    plt.xlim(eLo, eHi)
    plt.tight_layout()
    # plt.show()
    plt.savefig("./plots/sf-data-np.pdf")

    # === roofit style
    f = TFile("./data/fitWorkspace.root")
    fitWorkspace = f.Get("fitWorkspace")
    fData = fitWorkspace.allData().front()
    fitResult = fitWorkspace.allGenericObjects().front()
    nPars = fitResult.floatParsFinal().getSize()
    fEnergy = fitWorkspace.var("trapENFCal")
    modelPDF = fitWorkspace.pdf("model")
    nBins = int((eHi-eLo)/epb + 0.5)
    fSpec = fEnergy.frame(RF.Range(eLo,eHi), RF.Bins(nBins))
    fData.plotOn(fSpec)
    c = TCanvas()
    fSpec.Draw()
    c.Print("./plots/sf-data-rp.pdf")


def plotFit():

    from ROOT import TFile

    tf1 = TFile("./data/fitWorkspace.root")
    fitWorkspace = tf1.Get("fitWorkspace")
    fitResult = fitWorkspace.allGenericObjects().front()

    # -- print fit results --
    print("-- spec-fit RESULTS -- ")
    # print("%-10s = %.3f" % ("chiSq",chiSquare))
    fitValsFinal = getRooArgDict( fitResult.floatParsFinal() )

    bgr = {b:{} for b in bkgModel} # results

    for name in sorted(fitValsFinal):

        # find the part of the model this corresponds to
        for b in bkgModel:
            if b in name:
                break

        fitVal = fitValsFinal[name]
        error = fitWorkspace.var(name).getError()

        if "amp" in name:
            print("%-10s = best %-7.3f  error %.3f (w/o profile)" % (name, fitVal, error))
            bgr[b]["amp"] = [fitVal, error]
        elif "mu" in name:
            # compare the energy offset
            pkName = name[3:]
            pct = 100*(1 - fitVal/pkList[pkName])
            print("%-10s : fit %-6.3f  lit %-6.3f  (%.3f%%)" % (name, fitVal, pkList[pkName], pct))
            bgr[b]["mu"] = [fitVal, error]
        elif "sig" in name:
            # compare the sigma difference
            pkName = name[4:]
            pct = 100*(1 - fitVal/getSigma(pkList[pkName]))
            print("%-10s : fit %-6.3f  func %-6.3f  (%.3f%%)" % (name, fitVal, getSigma(pkList[pkName]), pct))
            bgr[b]["sig"] = [fitVal, error]
        elif "bkg" in name:
            print("%s = %.4f" % (name, fitVal))
            bgr[b]["bkg"] = [fitVal, error]

    for b in bgr:
        print(b, bgr[b])

    tf2 = TFile("./data/latDS%s.root" % ''.join([str(d) for d in dsList]))
    tt = tf2.Get("skimTree")
    tCut = "isEnr" if enr else "!isEnr"
    tCut += " && trapENFCal >= %.1f && trapENFCal <= %.1f" % (eLo, eHi)
    n = tt.Draw("trapENFCal", tCut, "goff")
    hitE = tt.GetV1()
    hitE = [hitE[i] for i in range(n)]
    x, hData = wl.GetHisto(hitE, eLo, eHi, epb)

    # plt.plot(x, hData, ls='steps', c='b') # normal histo
    hErr = np.asarray([np.sqrt(h) for h in hData]) # statistical error
    plt.errorbar(x, hData, yerr=hErr, c='k', ms=10, linewidth=0.8, fmt='.', capsize=2) # pretty convincing rooplot fake

    # plot tritium
    f = np.load("./data/scaledPDFs.npz")
    pdfRaw, pdfNorm, pdfEff, pdfEffN = f['arr_0'].item(), f['arr_1'].item(), f['arr_2'].item(), f['arr_3'].item()

    # xn, hn, xbn = pdfNorm["trit"]
    xen, hen, xben = pdfEffN["trit"]
    henc = getEffCorr(xen, hen, inv=True)

    # plt.plot(xn, hn, ls='steps', c='r', lw=3, label="Tritium, norm. PDF")
    # plt.plot(xen, hen * bgr["trit"]["amp"][0] * 4, ls='steps', c='b', lw=3, label="w/ normalized efficiency correction")
    plt.plot(xen, henc * bgr["trit"]["amp"][0], ls='steps', c='b', lw=3, label="w/ normalized efficiency correction")
    # plt.plot(xen, henc, ls='steps', c='g', lw=2, label="Final")

    plt.xlabel("Energy (keV)", ha='right', x=1)
    plt.ylabel("Events / %.1f keV" % epb, ha='right', y=1)
    plt.xlim(eLo, eHi)
    plt.tight_layout()
    plt.show()


def plotFitRF():
    from ROOT import TFile, gStyle, TLegend, TCanvas, TPad

    # load workspace
    f = TFile("./data/fitWorkspace.root")
    fitWorkspace = f.Get("fitWorkspace")
    fData = fitWorkspace.allData().front()
    fitResult = fitWorkspace.allGenericObjects().front()
    nPars = fitResult.floatParsFinal().getSize()
    fEnergy = fitWorkspace.var("trapENFCal")
    # fEnergy = fitWorkspace.var("xE")
    modelPDF = fitWorkspace.pdf("model")
    fitWorkspace.Print()
    # return

    # get a list of pdf names
    pdfNames = []
    pdfList = fitWorkspace.allPdfs() # goddamn RooArgSet
    itr = pdfList.createIterator()
    var = itr.Next()
    while var :
        name = var.GetName()
        if "ext" in name:
            pdfNames.append(name)
        var = itr.Next()
    pdfNames = sorted(pdfNames)

    # -- create spectrum rooplot --
    nCol = float(gStyle.GetNumberOfColors())
    binSize = 0.2
    nBins = int((eHi-eLo)/binSize + 0.5)
    fSpec = fEnergy.frame(RF.Range(eLo,eHi), RF.Bins(nBins))
    fData.plotOn(fSpec)

    modelPDF.plotOn(fSpec, RF.LineColor(ROOT.kRed), RF.Name("FullModel"))
    chiSquare = fSpec.chiSquare(nPars)

    # draw components
    leg = TLegend(0.83,0.1,0.97,0.9)
    leg.AddEntry(fSpec.findObject("FullModel"),"model #chi^{2}=%.3f" % chiSquare,"l")
    # for idx in range(len(pdfNames)):
    #     name = pdfNames[idx]
    #     lineCol = gStyle.GetColorPalette(int(nCol / len(pdfNames) * idx))
    #     modelPDF.plotOn(fSpec, RF.Components(name), RF.LineColor(lineCol), RF.Name(name))
    #     plotName = name
    #     if "ext-" in name:
    #         plotName = plotName[4:]
    #     leg.AddEntry(fSpec.findObject(name), plotName, "l")

    # create normalized residual ("pull") rooplot
    # (draw full model again so residuals are calculated against it)
    modelPDF.plotOn(fSpec, RF.LineColor(ROOT.kRed), RF.Name("FullModel"))
    res = fSpec.pullHist()
    fRes = fEnergy.frame(RF.Title(" "))
    fRes.addPlotable(res,"P")

    # -- print fit results --
    print("-- spec-fit RESULTS -- ")
    print("%-10s = %.3f" % ("chiSq",chiSquare))
    fitValsFinal = getRooArgDict( fitResult.floatParsFinal() )

    for name in sorted(fitValsFinal):
        fitVal = fitValsFinal[name]

        print("%s  fitVal %.2e" % (name, fitVal))

        # if "amp-" in name:
        #     error = fitWorkspace.var(name).getError()
        #     print("%-10s = best %-7.3f  error %.3f (w/o profile)" % (name, fitVal, error))
        # elif "mu-" in name:
        #     # compare the energy offset
        #     pkName = name[3:]
        #     pct = 100*(1 - fitVal/pkList[pkName])
        #     print("%-10s : fit %-6.3f  lit %-6.3f  (%.3f%%)" % (name, fitVal, pkList[pkName], pct))
        # elif "sig-" in name:
        #     # compare the sigma difference
        #     pkName = name[4:]
        #     pct = 100*(1 - fitVal/getSigma(pkList[pkName]))
        #     print("%-10s : fit %-6.3f  func %-6.3f  (%.3f%%)" % (name, fitVal, getSigma(pkList[pkName]), pct))
        # else:
        #     print("%s = %.4f" % (name, fitVal))
        #     continue

    # -- make spectrum plot w/ residual, w/ all the formatting crap --
    gStyle.SetOptStat(0);
    gStyle.SetPalette(ROOT.kRainBow) # https://root.cern.ch/doc/master/classTColor.html#C06
    c = TCanvas("c","Bob Ross's Canvas", 1400, 1000)
    p1 = TPad("p1","spectrum",0.,0.3,1.,1.)
    p2 = TPad("p2","residual",0.,0.,1.,0.3)
    p1.Draw()
    p2.Draw()
    p1.SetRightMargin(0.2)
    p1.SetBottomMargin(0.)
    p2.SetRightMargin(0.2)
    p2.SetTopMargin(0.)
    p2.SetBottomMargin(0.4)

    p1.cd()
    # p1.SetGrid()
    fSpec.SetTitle("")
    fSpec.GetXaxis().SetLabelSize(0.)
    fSpec.GetYaxis().SetLabelSize(0.05)
    fSpec.GetYaxis().SetTitleSize(0.06)
    fSpec.GetYaxis().SetTitleOffset(0.6)
    fSpec.GetYaxis().SetTitle("Events / (%.1f keV)" % binSize)

    fSpec.Draw()
    leg.Draw("same")

    p2.cd()
    fRes.GetXaxis().SetLabelSize(0.1)
    fRes.GetYaxis().SetLabelSize(0.1)
    fRes.GetXaxis().SetTitleSize(0.15)
    fRes.GetYaxis().SetTitle("Normalized Resid.")
    fRes.GetYaxis().SetTitleOffset(0.3)
    fRes.GetYaxis().SetTitleSize(0.1)
    fRes.Draw()
    zeroLine = ROOT.TLine(eLo,0.,eHi,0.)
    zeroLine.SetLineColor(ROOT.kBlack)
    zeroLine.Draw("same")
    c.Print("./plots/spectrum.pdf")

    # -- correlation matrix --
    # c2 = TCanvas("c2","Bob Ross's Canvas", 1100, 800)
    # c2.SetLeftMargin(0.15)
    # c2.SetRightMargin(0.12)
    # gStyle.SetPalette(ROOT.kDarkBodyRadiator)
    # corrMat = fitResult.correlationHist("Correlation Matrix")
    # corrMat.SetTitle("")
    # corrMat.GetXaxis().SetLabelSize(0.02)
    # corrMat.Draw("colz")
    # c2.Print("./plots/corrMatrix.pdf")


if __name__=="__main__":
    main(sys.argv[1:])