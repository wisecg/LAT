#!/usr/bin/env python3
import sys
import warnings

import matplotlib.pyplot as plt
plt.style.use('./pltReports.mplstyle')

import numpy as np
import ROOT
from ROOT import gROOT
from ROOT import RooFit as RF
from ROOT import RooStats as RS

import waveLibs as wl
import dsi

# gStyle.SetOptStat(0)
gROOT.ProcessLine("gErrorIgnoreLevel = 3001;") # suppress ROOT messages
gROOT.ProcessLine("RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);")
# # gROOT.ProcessLine("RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);")
# ROOT.Math.MinimizerOptions.SetDefaultMinimizer("Minuit") # the PLC doesn't work w/ minuit2

dsList = [0,1,2,3,4,"5A","5B","5C"]
eLo, eHi, epb = 1, 50, 0.2
enr = True

# must be specified here to be included in the fit
bkgModel = ["68GeK","55Fe","65ZnK","trit","bkg","Pb210","axion"]

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
    # "51Cr": 5.46,
    # "5678Co":7.11,  "56Ni": 7.71,
    # "734As":11.10,
    # "41Ca": 3.31,
    # "36Cl": 2.307,
    # "ax_Si_Ka1a2":1.739, "ax_Si_Kb1":1.836,
    # "ax_S_Ka1a2":2.307,
    # "ax_S_Kb1":2.464
    }


def main(argv):

    # global eLo, eHi, epb, dsList, enr
    # eLo, eHi, epb, enr = ... # adjust the global values here
    # expoDict = {} todo, make exposures accessible

    # === routines ===
    # loadDataMJD()
    # generatePDFs(makePlots=False)
    # axionPeaks()
    # plotPDFs()
    runFit()
    # compareData()
    # plotFit()
    plotFitRF()


def loadDataMJD():
    """ RooFit can't handle the vector<double> format for energies.
        So save a few select branches only (can add branches if necessary) into a new file.
    """
    from array import array
    from ROOT import TChain, TFile, TTree

    tt = TChain("skimTree")
    for ds in dsList:
        tt.Add("%s/final95t/final95t_DS%s.root" % (dsi.cutDir, ds))

    fName = "./data/latDS%s.root" % ''.join([str(d) for d in dsList])
    fOut = TFile(fName,"RECREATE")
    tOut = TTree("skimTree", "skimTree")
    run = array('i',[0])
    iEvt = array('i',[0])
    iHit = array('i',[0])
    chan = array('i',[0])
    hitE = array('d',[0.])
    # weight = array('d',[0.])
    isEnr = array('i',[0])
    tOut.Branch("run", run, "run/I")
    tOut.Branch("iEvent", iEvt, "iEvent/I")
    tOut.Branch("iHit", iHit, "iHit/I")
    tOut.Branch("channel", chan, "channel/I")
    tOut.Branch("trapENFCal", hitE, "trapENFCal/D")
    # tOut.Branch("weight", weight, "weight/D")
    tOut.Branch("isEnr", isEnr, "isEnr/I")

    for iE in range(tt.GetEntries()):
        tt.GetEntry(iE)
        run[0] = tt.run
        iEvt[0] = tt.iEvent
        for iH in range(tt.channel.size()):
            iHit[0] = tt.iHit.at(iH)
            chan[0] = tt.channel.at(iH)
            hitE[0] = tt.trapENFCal.at(iH)
            # weight[0] = 1.
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
    t2.Scan("run:channel:isEnr:trapENFCal")


def sig_ae(E,m):
    """ E, m are in units of keV.  must multiply result by sig_pe """
    beta = (1 - m**2./E**2.)**(1./2)
    return (1 - (1./3.)*beta**(2./3.)) * (3. * E**2.) / (16. * np.pi * (1./137.) * 511.**2. * beta)


def generatePDFs(ma=0, makePlots=False):
    """ Generate a set of TH1D's to be turned into RooDataHist objects.
    Takes axion mass (in keV) as a parameter.
    Full suite of diganostic plots is in ./sandbox/specPlots.py
    Make the binning 0.02 keV intervals - hopefully that's fine enough.
    TODO:
        - BDM/ALP PDF. Get from Kris's analysis in GAT
        - Pb210 pdf (46 kev line + continuum)
    """
    from ROOT import TFile, TH1D, gROOT

    # output files
    tf = TFile("./data/specPDFs.root","RECREATE")
    td = gROOT.CurrentDirectory()
    nf = "./data/specPDFs.npz"

    # energy limits
    xLo, xHi, xpb = 0, 30, 0.05
    nB = int((xHi-xLo)/xpb)

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

    # === 2. ge photoelectric xs
    phoData = []
    with open("./data/ge76peXS.txt") as f2: # 2499 entries, 0.01 kev intervals
        lines = f2.readlines()
        for line in lines:
            data = line.split()
            phoData.append([float(data[0]),float(data[1])])
    phoData = np.array(phoData)

    # plt.plot(tritData[:,0], tritData[:,1])
    # plt.show()
    # exit()

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
    tritData[:,1] *= 1/np.max(tritData[:,1]) # normalize the max to 1

    # === 4. Pb210
    tf2 = TFile("./data/Pb210PDFs.root")
    hPb210 = tf2.Get("hPb210")       # without TDL
    hPb210TDL = tf2.Get("hPb210TDL") # with TDL

    # NOTE: check sandbox/th1.py for examples of manually filling TH1D's and verifying wl.GetHisto and wl.npTH1D.

    # ROOT output
    h1 = TH1D("h1","photoelectric",nB,xLo,xHi)         # [cm^2 / kg]
    h2 = TH1D("h2","axioelectric",nB,xLo,xHi)          # [cm^2 / kg]
    h3 = TH1D("h3","axion flux, gae=1",nB,xLo,xHi)     # [cts / (keV cm^2 d)]
    h4 = TH1D("h4","convolved flux",nB,xLo,xHi)        # [cts / (keV d kg)]
    h5 = TH1D("h5","tritium",nB,xLo,xHi)               # [cts] (normalized to 1)

    # manually fill ROOT histos output
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

            # axion flux convolved [flux / (keV d kg)]
            axConv = axio * flux
            h4.SetBinContent(iB+1, axConv)

            # tritium
            idx = np.where((tritData[:,0] >= bLo) & (tritData[:,0] <= bHi))
            trit = np.mean(tritData[idx][:,1])
            if np.isnan(trit): trit = 0.
            h5.SetBinContent(iB+1, trit)


    if makePlots:
        from ROOT import TCanvas

        # === 1. verify the numpy histogram and ROOT histogram give the same output. OK

        # x, h210, xpb = wl.npTH1D(hPb210)
        # iE = np.where((x > 45) & (x < 48))
        # plt.plot(x[iE], h210[iE], ls='steps', lw=3, c='b')
        # plt.xlabel("Energy (keV)", ha='right', x=1)
        # plt.tight_layout()
        # plt.savefig("./plots/sf-pk210.pdf")
        #
        c = TCanvas()
        # hPb210.GetXaxis().SetTitle("Energy (keV)")
        # hPb210.GetXaxis().SetRangeUser(45, 48)
        # hPb210.Draw('hist')
        # c.Print('./plots/sf-pb210th1d.pdf')

        # === 2. print ROOT histos to match w/ numpy histos

        c.Clear(); h1.Draw("hist"); c.Print("./plots/root-sigGe.pdf")
        c.Clear(); h2.Draw("hist"); c.Print("./plots/root-sigAe.pdf")
        c.Clear(); h3.Draw("hist"); c.Print("./plots/root-axFlux.pdf")
        c.Clear(); h4.Draw("hist"); c.Print("./plots/root-axPDF.pdf")
        c.Clear(); h5.Draw("hist"); c.Print("./plots/root-trit.pdf")

    # save numpy output
    pdfs = {
        "sig_ge":   wl.npTH1D(h1),        # x, h, xpb. [cm^2 / kg]
        "sig_ae":   wl.npTH1D(h2),        # [cm^2 / kg]
        "axFlux":   wl.npTH1D(h3),        # [cts / (keV cm^2 d)]
        "axPDF":    wl.npTH1D(h4),        # [cts / (keV d kg)]
        "trit":     wl.npTH1D(h5),        # [cts] (normalized)
        "Pb210":    wl.npTH1D(hPb210),    # [cts] (normalized)
        "Pb210TDL": wl.npTH1D(hPb210TDL)  # [cts] (normalized)
    }
    np.savez(nf, pdfs)

    gROOT.cd(td.GetPath())
    h1.Write()
    h2.Write()
    h3.Write()
    h4.Write()
    h5.Write()
    # hPb210.Write()
    hPb210TDL.Write()
    tf.Close()


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


def plotPDFs():

    f = np.load("./data/specPDFs.npz")
    pdfs = f['arr_0'].item()

    # === 1. Ge photoelectric XS
    plt.close()
    x, n1, xpb = pdfs["sig_ge"]
    plt.semilogy(x, n1, ls='steps', c='b', lw=3, label=r"$\sigma_{ge}$")
    plt.xlabel("Energy (keV)", ha='right', x=1)
    plt.ylabel(r"$\mathregular{cm^2/kg}$", ha='right', y=1)
    plt.legend()
    plt.tight_layout()
    # plt.show()
    plt.savefig("./plots/sf-gexs.pdf")

    # === 2. axioelectric XS
    plt.close()
    x, n1, xpb = pdfs["sig_ae"]
    plt.plot(x, n1, ls='steps', c='b', lw=3, label=r"$\sigma_{ae}$")
    plt.xlabel("Energy (keV)", ha='right', x=1)
    plt.ylabel(r"$\mathregular{cm^2/kg}$", ha='right', y=1)
    plt.legend()
    plt.tight_layout()
    # plt.gca().yaxis.set_label_coords(-0.06, 1)
    # plt.show()
    plt.savefig("./plots/sf-axs.pdf")

    # === 3. axion flux, gae=1
    plt.close()
    x, n1, xpb = pdfs["axFlux"]
    plt.plot(x, n1, ls='steps', c='b', lw=3, label=r"$\Phi_a$")
    plt.xlabel("Energy (keV)", ha='right', x=1)
    plt.ylabel(r"Flux / (keV cm${}^2$ d)]", ha='right', y=1)
    plt.xlim(0,10)
    plt.legend()
    plt.tight_layout()
    # plt.show()
    plt.savefig("./plots/sf-axFlux.pdf")

    # === 4. solar axion PDF, gae=1
    plt.close()
    x, n1, xpb = pdfs["axPDF"]
    plt.plot(x, n1, ls='steps', c='b', lw=3, label=r"$\Phi_a$")
    plt.xlabel("Energy (keV)", ha='right', x=1)
    plt.ylabel("Flux / (keV d kg)", ha='right', y=1)
    plt.xlim(0,10)
    plt.legend()
    plt.tight_layout()
    # plt.show()
    plt.savefig("./plots/sf-axPDF.pdf")

    # === 5. tritium
    plt.close()
    x, n1, xpb = pdfs["trit"]
    plt.plot(x, n1, ls='steps', c='b', lw=3, label="Tritium")
    plt.xlabel("Energy (keV)", ha='right', x=1)
    plt.ylabel("Counts (norm)", ha='right', y=1)
    plt.legend()
    plt.xlim(0,20)
    plt.tight_layout()
    plt.savefig('./plots/sf-trit.pdf')

    # === 6. Pb210 PDF (with and without TDL)
    plt.close()
    x, n1, xpb = pdfs["Pb210"]
    x, n2, xpb = pdfs["Pb210TDL"]
    iE = np.where(x > 1)
    x, n1, n2 = x[iE], n1[iE], n2[iE]
    plt.plot(x, n1, ls='steps', c='b', lw=3, alpha=0.7, label="Pb210, no TDL")
    plt.plot(x, n2, ls='steps', c='r', lw=2, alpha=0.7, label="Pb210, with TDL")
    plt.axvline(1., c='g', lw=1, label="1.0 keV")
    plt.xlabel("Energy (keV)", ha='right', x=1)
    plt.ylabel("Counts (norm)", ha='right', y=1)
    plt.legend()
    plt.tight_layout()
    plt.savefig("./plots/sf-pb210.pdf")


def getSigma(E):
    # pg. 92 of graham's thesis (he ended up letting sigma float)
    sig_e, F, expval = 0.0698, 0.21, 0.00296
    return np.sqrt(sig_e * sig_e + expval * F * E)


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

    # load data and create a workspace
    f1 = TFile("./data/latDS%s.root" % ''.join([str(d) for d in dsList]))
    t = f1.Get("skimTree")
    fEnergy = ROOT.RooRealVar("trapENFCal","Energy",eLo,eHi,"keV")
    fEnr = ROOT.RooRealVar("isEnr","isEnr",0,1,"")
    # fWeight = ROOT.RooRealVar("weight","weight",1,10,"")

    tCut = "isEnr==1" if enr is True else "isEnr==0"
    fData = ROOT.RooDataSet("data", "data", t, ROOT.RooArgSet(fEnergy,fEnr), tCut)

    fitWorkspace = ROOT.RooWorkspace("fitWorkspace","Fit Workspace")
    getattr(fitWorkspace,'import')(fEnergy)
    getattr(fitWorkspace,'import')(fData)
    # getattr(fitWorkspace,'import')(fWeight)

    # ==== background model ====
    pdfList = ROOT.RooArgList("shapes")

    # linear background function (basically flat)
    bkgNum = ROOT.RooRealVar("bkgNum","bkgNum",1.,5000.)
    bkgSlo = ROOT.RooRealVar("bkgSlo","bkgSlo",-1.,1.)
    bkgPol = ROOT.RooChebychev("bkgPol","bkgPol",fEnergy,ROOT.RooArgList(bkgSlo)) # more stable than roopolynomial
    bkgExt = ROOT.RooExtendPdf("bkgExt", "bkgExt", bkgPol, bkgNum)
    if "bkg" in bkgModel:
        pdfList.add(bkgExt)

    # testing
    # pk1 = pkModel(fEnergy,"68GeK",pkList["68GeK"])
    # pkNum, pkMu, pkSig, pkGaus, pkExt = pk1.GetAll()
    # pdfList.add(pkExt)

    # gaussian peak list
    # (the trick seems to be that you can't overwrite the pkModel object in the loop.)
    pks = []
    for pk in pkList:
        if pk not in bkgModel: continue
        pks.append( pkModel(fEnergy, pk, pkList[pk]) )
    for pk in pks:
        pdfList.add(pk.GetPkExt())

    # load special PDFs
    f2 = TFile("./data/specPDFs.root")

    # axion continuum
    axTH1D = f2.Get("h4")
    axNum = ROOT.RooRealVar("amp-axion", "amp-axion", 5000., 0., 20000.)
    # axNum = ROOT.RooRealVar("amp-axion","amp-axion",16.904) # hardcode the profile upper limit
    axDataHist = ROOT.RooDataHist("ax", "ax", ROOT.RooArgList(fEnergy), RF.Import(axTH1D))
    fEnergy.setRange(eLo, eHi) # have to reset after loading a histo w/ different bounds
    axPdf = ROOT.RooHistPdf("axPdf", "axPdf", ROOT.RooArgSet(fEnergy), axDataHist, 2)
    axExt = ROOT.RooExtendPdf("ext-axion", "ext-axion", axPdf, axNum)
    if "axion" in bkgModel:
        pdfList.add(axExt)

    # tritium
    trTH1D = f2.Get("h5")
    trNum = ROOT.RooRealVar("amp-trit", "amp-trit", 100., 0., 5000.)
    trDataHist = ROOT.RooDataHist("tr", "tr", ROOT.RooArgList(fEnergy), RF.Import(trTH1D))
    fEnergy.setRange(eLo, eHi)
    trPdf = ROOT.RooHistPdf("trPdf", "trPdf", ROOT.RooArgSet(fEnergy), trDataHist, 2)
    trExt = ROOT.RooExtendPdf("ext-trit", "ext-trit",trPdf,trNum)
    if "trit" in bkgModel:
        pdfList.add(trExt)

    # Pb210 continuum (w/ TDL)
    pbTH1D = f2.Get("hPb210TDL")
    pbNum = ROOT.RooRealVar("amp-Pb210", "amp-Pb210", 100., 0., 1000.)
    pbDataHist = ROOT.RooDataHist("pb", "pb", ROOT.RooArgList(fEnergy), RF.Import(pbTH1D))
    fEnergy.setRange(eLo, eHi)
    pbPdf = ROOT.RooHistPdf("pbPdf", "pbPdf", ROOT.RooArgSet(fEnergy), pbDataHist, 2)
    pbExt = ROOT.RooExtendPdf("ext-Pb210", "ext-Pb210",pbPdf,pbNum)
    if "Pb210" in bkgModel:
        pdfList.add(pbExt)

    # create total model pdf
    model = ROOT.RooAddPdf("model","total pdf",pdfList)

    # === efficiency ===
    effFile = TFile("./data/lat-expo-efficiency.root")
    effHist = effFile.Get("hDS5B_Norm_Enr")
    # x, y, xpb = wl.npTH1D(effHist)
    # plt.plot(x, y, ls='steps')
    # plt.show()

    effRooHist = ROOT.RooDataHist("eff","Efficiency", ROOT.RooArgList(fEnergy), RF.Import(effHist))
    fEnergy.setRange(eLo, eHi)
    effPdf = ROOT.RooHistPdf("effPdf","effPdf", ROOT.RooArgSet(fEnergy), effRooHist, 0)
    modelEff = ROOT.RooProdPdf("modelEff","model with efficiency", model, effPdf)

    # modelEff = ROOT.RooProdPdf("modelEff","model with efficiency", ROOT.RooArgList(model, effPdf))

    # run fitter
    # minimizer = ROOT.RooMinimizer( model.createNLL(fData, RF.NumCPU(2,0), RF.Extended(True)) )
    minimizer = ROOT.RooMinimizer( modelEff.createNLL(fData, RF.NumCPU(2,0), RF.Extended(False)) )
    minimizer.setPrintLevel(-1)
    minimizer.setStrategy(2)
    minimizer.migrad()
    fitResult = minimizer.save()

    # according to the internet, covQual==3 is a good indicator that it converged
    print("Fitter is done. Fit Cov Qual:", fitResult.covQual())

    # save workspace to a TFile
    getattr(fitWorkspace,'import')(fitResult)
    # getattr(fitWorkspace,'import')(model)
    getattr(fitWorkspace,'import')(modelEff)
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


def gaus(x, a, x0, sigma):
    return a * np.exp(-(x - x0)**2 / (2 * sigma**2))


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

    # plot peaks
    # xF = np.arange(eLo, eHi, 0.001)
    # for b in bkgModel:
    #     if b in pkList.keys():
    #         amp, mu, sig = bgr[b]["amp"][0], bgr[b]["mu"][0], bgr[b]["sig"][0]
    #         plt.plot(xF, gaus(xF, amp, mu, sig))

    # plot tritium
    f = np.load("./data/specPDFs.npz")
    pdfs = f['arr_0'].item()
    xP, aP, xpb3 = pdfs["trit"]

    tNorm = np.sum(aP)

    plt.plot(xP, aP)

    plt.plot(xP, aP * bgr["trit"]["amp"][0]/ tNorm )
    # plt.plot(xP, aP * bgr["trit"]["amp"][0])

    plt.xlabel("Energy (keV)", ha='right', x=1)
    plt.ylabel("Events / %.1f keV" % epb, ha='right', y=1)
    plt.xlim(eLo, eHi)
    plt.tight_layout()
    plt.show()

    # yF = gaus(xF, 5.385, 10.412, 0.115)
    # plt.plot(xF, yF, '-r')

    # plt.xlabel("Energy (keV)", ha='right', x=1)
    # plt.ylabel("Counts", ha='right', y=1)
    # plt.tight_layout()
    # plt.show()


def plotFitRF():
    from ROOT import TFile, gStyle, TLegend, TCanvas, TPad

    # load workspace
    f = TFile("./data/fitWorkspace.root")
    fitWorkspace = f.Get("fitWorkspace")
    fData = fitWorkspace.allData().front()
    fitResult = fitWorkspace.allGenericObjects().front()
    nPars = fitResult.floatParsFinal().getSize()
    fEnergy = fitWorkspace.var("trapENFCal")
    # modelPDF = fitWorkspace.pdf("model")
    # modelPDF = fitWorkspace.pdf("modelEff")
    fitWorkspace.Print()
    return

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
        if "amp-" in name:
            error = fitWorkspace.var(name).getError()
            print("%-10s = best %-7.3f  error %.3f (w/o profile)" % (name, fitVal, error))
        elif "mu-" in name:
            # compare the energy offset
            pkName = name[3:]
            pct = 100*(1 - fitVal/pkList[pkName])
            print("%-10s : fit %-6.3f  lit %-6.3f  (%.3f%%)" % (name, fitVal, pkList[pkName], pct))
        elif "sig-" in name:
            # compare the sigma difference
            pkName = name[4:]
            pct = 100*(1 - fitVal/getSigma(pkList[pkName]))
            print("%-10s : fit %-6.3f  func %-6.3f  (%.3f%%)" % (name, fitVal, getSigma(pkList[pkName]), pct))
        else:
            print("%s = %.4f" % (name, fitVal))
            continue

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