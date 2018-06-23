#!/usr/bin/env python3
import sys
import warnings

import numpy as np
import ROOT

import waveLibs as wl
import dsi


# import sys, time, warnings
# sys.argv.append("-b") # kill all interactive crap
# import ROOT
# from ROOT import TFile, TTree, TChain, TCanvas, TH1D, TLegend, TPad, TLine, TGraph
# from ROOT import gStyle, gPad, gROOT, std
from ROOT import RooFit as RF
# from ROOT import RooStats as RS
# import numpy as np
# from array import array
#
# gStyle.SetOptStat(0)
# gROOT.ProcessLine("gErrorIgnoreLevel = 3001;") # suppress ROOT messages
# gROOT.ProcessLine("RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);")
# # gROOT.ProcessLine("RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);")
# ROOT.Math.MinimizerOptions.SetDefaultMinimizer("Minuit") # the PLC doesn't work w/ minuit2

def main(argv):

    # keep params we know we want to adjust up top
    dsList = [0,1,2,3,4,"5A","5B","5C"]
    # expoDict = {} todo, make exposures accessible
    global enr
    enr = False

    global eLo, eHi, epb
    eLo, eHi, epb = 1, 50, 0.2
    global pkList
    pkList = {
        "55Fe":6.54, "68Ga":9.66, "68GeK":10.37
    }
    # loadDataMJD(dsList)
    # generatePDFs()
    # runFit(dsList)
    plotSpectrum()


def loadDataMJD(dsList):
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


def generatePDFs(ma=0):
    """ Generate a set of TH1D's to be turned into RooDataHist objects.
        Takes axion mass (in keV) as a parameter.
        Full suite of diganostic plots is in ./sandbox/specPlots.py
        Make the binning 0.02 keV intervals - hopefully that's fine enough.
        TODO:
            - BDM/ALP PDF. Get from Kris's analysis in GAT
            - Pb210 pdf (46 kev line + continuum)
    """
    from ROOT import TFile, TH1D

    tf = TFile("./data/inputHists.root","RECREATE")

    # axion flux scale.
    # NOTE: to do the fit and set a new limit, we set g_ae=1.
    # To plot an expected flux, we would use a real value.
    # Redondo note: I calculated the flux using gae = 0.511*10^-10
    # for other values of gae use: FLUX = Table*[gae/(0.511*10^-10)]^2
    gae = 1
    gRat = (gae / 5.11e-11)
    redondoScale = 1e19 * gRat**2 # convert table to [flux / (keV cm^2 d)]

    axData, phoData, tritData = [], [], []
    with open("./data/redondoFlux.txt") as f1: # 23577 entries
        lines = f1.readlines()[11:]
        for line in lines:
            data = line.split()
            axData.append([float(data[0]),float(data[1])])
    with open("./data/ge76peXS.txt") as f2: # 2499 entries, 0.01 kev intervals
        lines = f2.readlines()
        for line in lines:
            data = line.split()
            phoData.append([float(data[0]),float(data[1])])
    with open("./data/TritiumSpectrum.txt") as f3: # 20000 entries
        lines = f3.readlines()[1:]
        for line in lines:
            data = line.split()
            conv = float(data[2]) # raw spectrum convolved w/ ge cross section
            if conv < 0: conv = 0.
            tritData.append([float(data[1]),conv])
    axData, phoData, tritData = np.array(axData), np.array(phoData), np.array(tritData)

    def sig_ae(E,m):
        """ E, m are in units of keV.  must multiply result by sig_pe """
        beta = (1 - m**2./E**2.)**(1./2)
        return (1 - (1./3.)*beta**(2./3.)) * (3. * E**2.) / (16. * np.pi * (1./137.) * 511.**2. * beta)

    eLo, eHi, epb = 0.1, 30, 0.02

    # root output
    nBins = int((eHi-eLo)/epb)
    h1 = TH1D("h1","photoelectric",nBins,eLo,eHi)         # [cm^2 / kg]
    h2 = TH1D("h2","axioelectric",nBins,eLo,eHi)          # [cm^2 / kg]
    h3 = TH1D("h3","axion flux, gae=1",nBins,eLo,eHi)     # [cts / (keV cm^2 d)]
    h4 = TH1D("h4","convolved flux",nBins,eLo,eHi)        # [cts / (keV d kg)]
    h5 = TH1D("h5","tritium",nBins,eLo,eHi)               # [cts] (normalized to 1)

    # TODO: numpy output (compare wl.GetHisto with wl.npTH1D and use it to find the axion peaks)
    # xn, h1n = wl.GetHisto()

    for i in range(nBins):
        ene = i * epb + eLo
        eneLo, eneHi = ene - epb/2., ene + epb/2.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore",category=RuntimeWarning)

            # if ma>0, we ignore entries with E <= m.

            # photoelectric x-section [cm^2 / kg]
            idx = np.where((phoData[:,0] >= eneLo) & (phoData[:,0] <= eneHi))
            pho = np.mean(phoData[idx][:,1]) * 1000
            if np.isnan(pho) or len(phoData[idx][:,1]) == 0: pho = 0.
            if phoData[idx][:,1] <= ma: pho = 0.
            h1.SetBinContent(i,pho)

            # axioelectric x-section [cm^2 / kg]
            if ene > ma: axio = pho * sig_ae(ene, ma)
            else: axio=0.
            h2.SetBinContent(i,axio)

            # axion flux [cts / (cm^2 d keV)]
            idx = np.where((axData[:,0] >= eneLo) & (axData[:,0] <= eneHi))
            flux = np.mean(axData[idx][:,1]) * redondoScale
            if np.isnan(flux): flux = 0.
            h3.SetBinContent(i,flux)

            # axion flux convolved [cts / (keV d kg)]
            axConv = axio * flux
            h4.SetBinContent(i, axConv)

            # tritium
            idx = np.where((tritData[:,0] >= eneLo) & (tritData[:,0] <= eneHi))
            trit = np.mean(tritData[idx][:,1])
            if np.isnan(trit): trit = 0.
            h5.SetBinContent(i, trit)

    h1.Write()
    h2.Write()
    h3.Write()
    h4.Write()
    h5.Write()
    tf.Close()


def getSigma(E):
    # pg. 92 of graham's thesis (he ended up letting sigma float)
    sig_e, F, expval = 0.0698, 0.21, 0.00296
    return np.sqrt(sig_e * sig_e + expval * F * E)


class pkModel:

    def __init__(self,fEnergy,name,ene,amp=0):
        """ TODO: add in constraints to the likelihood function instead of hard bounds
        https://root.cern.ch/root/html/tutorials/roofit/rf604_constraints.C.html
        """
        import ROOT
        sig = getSigma(ene)
        muLo, muHi = ene * 0.99, ene * 1.01
        sgLo, sgHi = sig * 0.66, sig * 1.33
        ampLo, ampHi = -0.1, 500    # letting it go slightly negative helps the plc calculator
        if eHi > 10.: ampHi = 2000. # 68GeK is around 1700, all the rest are under 500

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


def runFit(dsList):
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

    # background model
    pdfList = ROOT.RooArgList("shapes")

    # linear background function
    bkgNum = ROOT.RooRealVar("bkgNum","bkgNum",1.,5000.)
    bkgSlo = ROOT.RooRealVar("bkgSlo","bkgSlo",-1.,1.)
    bkgPol = ROOT.RooChebychev("bkgPol","bkgPol",fEnergy,ROOT.RooArgList(bkgSlo)) # more stable than roopolynomial
    bkgExt = ROOT.RooExtendPdf("bkgExt", "bkgExt", bkgPol, bkgNum);
    pdfList.add(bkgExt);

    # cosmogenic (and axion) peaks (declared in main function)
    peaks = {}
    for key in pkList:
        ene = pkList[key]
        if eLo <= ene <= eHi:
            if key=="68GeL":
                thisPk = pkModel(fEnergy, key, ene) # reminder: can set peak-specific "amp" here
            else:
                thisPk = pkModel(fEnergy, key, ene)
            peaks[key] = thisPk
            pdfList.add(thisPk.GetPkExt())

    # axion continuum
    axNum = ROOT.RooRealVar("amp-axion", "amp-axion", 100., 0., 10000.)
    # axNum = ROOT.RooRealVar("amp-axion","amp-axion",16.904) # hardcode the profile upper limit
    f2 = TFile("./data/inputHists.root")
    axTH1D = f2.Get("h4")
    axDataHist = ROOT.RooDataHist("ax", "ax", ROOT.RooArgList(fEnergy), RF.Import(axTH1D))
    fEnergy.setRange(eLo, eHi) # have to reset after loading a histo w/ different bounds
    axPdf = ROOT.RooHistPdf("axPdf", "axPdf", ROOT.RooArgSet(fEnergy), axDataHist, 2)
    axExt = ROOT.RooExtendPdf("ext-axion", "ext-axion", axPdf, axNum)
    pdfList.add(axExt)

    # tritium
    trNum = ROOT.RooRealVar("amp-trit", "amp-trit", 100., 0., 10000.)
    trTH1D = f2.Get("h5")
    trDataHist = ROOT.RooDataHist("tr", "tr", ROOT.RooArgList(fEnergy), RF.Import(trTH1D))
    fEnergy.setRange(eLo, eHi)
    trPdf = ROOT.RooHistPdf("trPdf", "trPdf", ROOT.RooArgSet(fEnergy), trDataHist, 2)
    trExt = ROOT.RooExtendPdf("ext-trit", "ext-trit",trPdf,trNum)
    pdfList.add(trExt)

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


def getRooArgDict(arglist):
    """ Convert a RooArgList into a python dict.
    NOTE: I'm doing it this way because this method doesn't work:
    pkValF = fitResult.floatParsFinal().find("pk_gaus") # can't cast to RooRealVar
    """
    pkVals = {}
    for i in range(arglist.getSize()):
        pkVals[ arglist.at(i).GetName() ] = arglist.at(i).getValV()
    return pkVals


def plotSpectrum():
    from ROOT import TFile, gStyle, TLegend, TCanvas, TPad

    # load workspace
    f = TFile("./data/fitWorkspace.root")
    fitWorkspace = f.Get("fitWorkspace")
    fData = fitWorkspace.allData().front()
    fitResult = fitWorkspace.allGenericObjects().front()
    nPars = fitResult.floatParsFinal().getSize()
    fEnergy = fitWorkspace.var("trapENFCal")
    modelPDF = fitWorkspace.pdf("model")
    # fitWorkspace.Print()

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
    # print(type(fData.plotOn(fSpec)))
    # return

    modelPDF.plotOn(fSpec, RF.LineColor(ROOT.kRed), RF.Name("FullModel"))
    chiSquare = fSpec.chiSquare(nPars)

    # draw components
    leg = TLegend(0.83,0.1,0.97,0.9)
    leg.AddEntry(fSpec.findObject("FullModel"),"model #chi^{2}=%.3f" % chiSquare,"l")
    for idx in range(len(pdfNames)):
        name = pdfNames[idx]
        lineCol = gStyle.GetColorPalette(int(nCol / len(pdfNames) * idx))
        modelPDF.plotOn(fSpec, RF.Components(name), RF.LineColor(lineCol), RF.Name(name))
        plotName = name
        if "ext-" in name:
            plotName = plotName[4:]
        leg.AddEntry(fSpec.findObject(name), plotName, "l")

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
    c2 = TCanvas("c2","Bob Ross's Canvas", 1100, 800)
    c2.SetLeftMargin(0.15)
    c2.SetRightMargin(0.12)
    gStyle.SetPalette(ROOT.kDarkBodyRadiator)
    corrMat = fitResult.correlationHist("Correlation Matrix")
    corrMat.SetTitle("")
    corrMat.GetXaxis().SetLabelSize(0.02)
    corrMat.Draw("colz")
    c2.Print("./plots/corrMatrix.pdf")


if __name__=="__main__":
    main(sys.argv[1:])