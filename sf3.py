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

def main(argv):

    # global eLo, eHi, epb, dsList, enr
    # eLo, eHi, epb, enr = ... # adjust the global values here
    # expoDict = {} todo, make exposures accessible
    initialize()

    # === routines ===
    # loadDataMJD()
    # getUnscaledPDFs(makePlots=False)
    runFit()
    plotFit()


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
    xLo, xHi, xpb = 0, 20, 0.05
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


def runFit():
    from ROOT import TFile, TH1D

    # load data and create a workspace
    f1 = TFile("./data/latDS%s.root" % ''.join([str(d) for d in dsList]))
    t = f1.Get("skimTree")
    fEnergy = ROOT.RooRealVar("trapENFCal","Energy",eLo,eHi,"keV")
    fEnr = ROOT.RooRealVar("isEnr","isEnr",0,1,"")
    fWeight = ROOT.RooRealVar("weight","weight",1,10,"")

    tCut = "isEnr==1" if enr is True else "isEnr==0"
    fData = ROOT.RooDataSet("data", "data", t, ROOT.RooArgSet(fEnergy,fEnr), tCut)
    # fData = ROOT.RooDataSet("data", "data", t, ROOT.RooArgSet(fEnergy,fEnr,fWeight),tCut,"weight")

    fitWorkspace = ROOT.RooWorkspace("fitWorkspace","Fit Workspace")
    getattr(fitWorkspace,'import')(fEnergy)
    getattr(fitWorkspace,'import')(fData)
    # getattr(fitWorkspace,'import')(fWeight)

    # ==== background model ====
    pdfList = ROOT.RooArgList("shapes")

    # load special PDFs
    f2 = TFile("./data/specPDFs.root")

    # tritium
    trTH1D = f2.Get("h5")
    trNum = ROOT.RooRealVar("amp-trit", "amp-trit", 100.0, 0.0, 5000.)
    trDataHist = ROOT.RooDataHist("tr", "tr", ROOT.RooArgList(fEnergy), RF.Import(trTH1D))
    fEnergy.setRange(eLo, eHi)
    trPdf = ROOT.RooHistPdf("trPdf", "trPdf", ROOT.RooArgSet(fEnergy), trDataHist, 2)
    trExt = ROOT.RooExtendPdf("ext-trit", "ext-trit",trPdf,trNum)
    pdfList.add(trExt)

    # flat bg
    # nB = int((eHi-eLo)/0.05)
    # bkgTH1D = TH1D("h8","flat BG",nB,eLo,eHi)
    # for iB in range(nB+1):
    #     bkgTH1D.SetBinContent(iB, 100) # the initial amplitude doesn't matter b/c we normalize to 1
    # bkgNum = ROOT.RooRealVar("amp-flat", "amp-flat", 50.0, 0.0, 5000.)
    # bkgDataHist = ROOT.RooDataHist("bkg", "bkg", ROOT.RooArgList(fEnergy), RF.Import(bkgTH1D))
    # fEnergy.setRange(eLo, eHi)
    # bkgPdf = ROOT.RooHistPdf("bkgPdf", "bkgPdf", ROOT.RooArgSet(fEnergy), bkgDataHist, 2)
    # bkgExt = ROOT.RooExtendPdf("ext-flat", "ext-flat", bkgPdf, bkgNum)
    # pdfList.add(bkgExt)

    # linear background function
    bkgNum = ROOT.RooRealVar("bkgNum","bkgNum",1.,5000.)
    bkgPol = ROOT.RooChebychev("bkgPol","bkgPol",fEnergy,ROOT.RooArgList()) # more stable than roopolynomial
    bkgExt = ROOT.RooExtendPdf("bkg-ext", "bkg-ext", bkgPol, bkgNum);
    pdfList.add(bkgExt);


    # nBkg = ROOT.RooRealVar("Bkg", "Bkg", 50.0, 0.0, 5000.)
    # bkgP = ROOT.RooPolynomial("Bkg", "Flat background function", fEnergy, ROOT.RooArgList())
    # bkgE = ROOT.RooExtendPdf("bkg-ext","Extended bkg function", bkgP, nBkg)
    # pdfList.add(bkgE)

    # create total model pdf
    model = ROOT.RooAddPdf("model","total pdf",pdfList)

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
    f2 = TFile("./data/fitWorkspace.root","RECREATE")
    fitWorkspace.Write()
    f2.Close()


def plotFit():
    from ROOT import TFile, gStyle, TLegend, TCanvas, TPad

    # load workspace
    f = TFile("./data/fitWorkspace.root")
    fitWorkspace = f.Get("fitWorkspace")
    fData = fitWorkspace.allData().front()
    fitResult = fitWorkspace.allGenericObjects().front()
    nPars = fitResult.floatParsFinal().getSize()
    fEnergy = fitWorkspace.var("trapENFCal")
    modelPDF = fitWorkspace.pdf("model")
    fitWorkspace.Print()

    binSize = 0.2
    nBins = int((eHi-eLo)/binSize + 0.5)
    fSpec = fEnergy.frame(RF.Range(eLo,eHi), RF.Bins(nBins))
    fData.plotOn(fSpec)
    modelPDF.plotOn(fSpec, RF.LineColor(ROOT.kRed), RF.Name("FullModel"))

    # plot the components
    pdfNames = []
    pdfList = fitWorkspace.allPdfs() # goddamn RooArgSet
    itr = pdfList.createIterator()
    var = itr.Next()
    while var :
        name = var.GetName()
        print("name:",name)
        if "ext" in name:
            pdfNames.append(name)
        var = itr.Next()
    pdfNames = sorted(pdfNames)

    gStyle.SetPalette(ROOT.kRainBow)
    nCol = float(gStyle.GetNumberOfColors())
    for idx in range(len(pdfNames)):
        name = pdfNames[idx]
        lineCol = gStyle.GetColorPalette(int(nCol / len(pdfNames) * idx))
        modelPDF.plotOn(fSpec, RF.Components(name), RF.LineColor(lineCol), RF.Name(name))

    fitValsFinal = getRooArgDict( fitResult.floatParsFinal() )
    for name in sorted(fitValsFinal):
        fitVal = fitValsFinal[name]
        print("%s  fitVal %.2e" % (name, fitVal))

    # draw the final rooplot
    c = TCanvas("c","c", 1400, 1000)
    fSpec.SetTitle("")
    fSpec.Draw()
    c.Print("./plots/spectrum.pdf")
    return


    # == now try to reproduce the plot in numpy

    # plot the data
    tf2 = TFile("./data/latDS%s.root" % ''.join([str(d) for d in dsList]))
    tt = tf2.Get("skimTree")
    tCut = "isEnr" if enr else "!isEnr"
    tCut += " && trapENFCal >= %.1f && trapENFCal <= %.1f" % (eLo, eHi)
    n = tt.Draw("trapENFCal", tCut, "goff")
    hitE = tt.GetV1()
    hitE = [hitE[i] for i in range(n)]
    x, hData = wl.GetHisto(hitE, eLo, eHi, epb)
    hErr = np.asarray([np.sqrt(h) for h in hData])
    # plt.errorbar(x, hData, yerr=hErr, c='k', ms=10, linewidth=0.8, fmt='.', capsize=2)


    # get the fit values



    # plot the scaled pdfs

    # f = np.load("./data/scaledPDFs.npz")
    # pdfRaw, pdfNorm, pdfEff, pdfEffN = f['arr_0'].item(), f['arr_1'].item(), f['arr_2'].item(), f['arr_3'].item()
    # x, h, _ = pdfNorm["trit"]
    # scale = 3000 # this is like, pretty close
    # plt.plot(x, h * scale, ls='steps', c='b', lw=3)

    # load unscaled PDFs
    tf = TFile("./data/specPDFs.root")
    h5 = tf.Get('h5') # tritium
    x5, np5, xpb5 = wl.npTH1D(h5)
    np5n = np.divide(np5, np.sum(np5[np.where((x5 >= eLo) & (x5 <= eHi))] * xpb5))

    plt.step(x5, np5n/np.amax(np5n)) # pin max to 1
    plt.show()
    return

    # plt.plot(x5, np5n * 100*30, ls='steps') # why is 3000 so close?
    plt.plot(x5, np5n * 10*250, ls='steps')

    plt.xlabel("Energy (keV)", ha='right', x=1)
    plt.ylabel("Events / %.1f keV" % epb, ha='right', y=1)
    plt.xlim(eLo, eHi)
    plt.tight_layout()
    plt.show()

    # plt.savefig("./plots/sf3-spectrum.pdf")


def getRooArgDict(arglist):
    """ Convert a RooArgList into a python dict.
    NOTE: I'm doing it this way because this method doesn't work:
    pkValF = fitResult.floatParsFinal().find("pk_gaus") # can't cast to RooRealVar
    """
    pkVals = {}
    for i in range(arglist.getSize()):
        pkVals[ arglist.at(i).GetName() ] = arglist.at(i).getValV()
    return pkVals


if __name__=="__main__":
    main(sys.argv[1:])