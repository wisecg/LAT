#!/usr/bin/env python
"""
specFit.py
Likelihood analysis of MALBEK and MJD low energy spectra.
Adapted from GPXFitter by Brian Zhu, LANL.  Much thanks to Brian, Graham, and Lukas.
Clint Wiseman, USC

v1. 5 Oct. 2017 (ported GPXFitter, performed MALBEK study)

Inputs:
  - MALBEK data from Graham, with rise time cuts and efficiency correction applied
  - MJD data from LAT
Background Model:
  - Polynomial (linear) background function
  - Cosmogenic peaks (CoGeNT (Clint's thesis proposal, slide 22), Graham's thesis, pg 78)
  - It seems 41Ca and 36Cl also MAY be contributing to MALBEK [2] [3]:
Axion Spectrum (Redondo)
  - Can either look for discrete peaks [3]
  - Or fit the entire spectrum [4].
  - We use Mucal.c to convolve the axion flux with the photoelectric cross section [5]

Ref 1: http://www.nndc.bnl.gov/chart/decaysearchdirect.jsp?nuc=41CA&unc=nds
Ref 2: http://www.kayelaby.npl.co.uk/atomic_and_nuclear_physics/4_2/4_2_1.html
Ref 3: Solar x-ray table (but no axion fluxes): http://xdb.lbl.gov/Section1/Table_1-2.pdf
Ref 4: Redondo's solar axion flux: http://wwwth.mpp.mpg.de/members/redondo/material.html
Ref 5: Mucal on the web: http://www.csrri.iit.edu/mucal.html
"""
import sys, time, warnings
sys.argv.append("-b") # kill all interactive crap
import ROOT
from ROOT import TFile, TTree, TCanvas, TH1D, TLegend, TPad, TLine
from ROOT import gStyle, gPad, gROOT, std
from ROOT import RooFit as RF
from ROOT import RooStats as RS
import numpy as np

gStyle.SetOptStat(0)
gROOT.ProcessLine("gErrorIgnoreLevel = 3001;") # suppress ROOT messages
# gROOT.ProcessLine("RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);")
gROOT.ProcessLine("RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);")
ROOT.Math.MinimizerOptions.SetDefaultMinimizer("Minuit") # the PLC doesn't work w/ minuit2

# ===================================================================================================
def main(argv):

    global eLo, eHi
    eLo, eHi = 1.5, 8. # graham's axion region
    # eLo, eHi = 8., 12.
    # eLo, eHi = 0.8, 5.

    global pkList
    pkList = {
        "65ZnL":1.10, "68GeL":1.29, "49V":   4.97,  # "51Cr": 5.46,
        "54Mn": 5.99, "55Fe": 6.54, # "5678Co":7.11,  "56Ni": 7.71,
        "65ZnK":8.98, "68Ga": 9.66, "68GeK": 10.37, # "734As":11.10,
        "41Ca": 3.31,
        "36Cl": 2.307,
        # "ax_Si_Ka1a2":1.739, "ax_Si_Kb1":1.836,
        # "ax_S_Ka1a2":2.307,
        # "ax_S_Kb1":2.464
        }

    getSpecPDFs()
    # runFit()
    # plotSpectrum()
    # plotProfiles()
    # calculate_g_ae()
    # plotContours()
    # plotMCSpectrum()
    # runMCStudy()

# ===================================================================================================
def getSigma(E):
    # pg. 92 of graham's thesis (he ended up letting sigma float)
    sig_e, F, expval = 0.0698, 0.21, 0.00296
    return np.sqrt(sig_e * sig_e + expval * F * E)


def getRooArgDict(arglist):
    """ Convert a RooArgList into a python dict.
    NOTE: I'm doing it this way because this method doesn't work:
    pkValF = fitResult.floatParsFinal().find("pk_gaus") # can't cast to RooRealVar
    """
    pkVals = {}
    for i in range(arglist.getSize()):
        pkVals[ arglist.at(i).GetName() ] = arglist.at(i).getValV()
    return pkVals


def getSpecPDFs():
    """ Return a set of TH1D's to be turned into RooDataHist objects.
        Make the binning 0.01 keV intervals - hopefully that's fine enough.
        TODO: WIMP pdf?  Get PyWIMP from Graham
    """
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

    # -- open an output file --
    f4 = TFile("./data/inputHists.root","RECREATE")

    # set the energy range from 0 to 50 (hardcoded)
    kevPerBin = 0.02
    nBins = int((50.)/kevPerBin + 0.5)
    hPhoto = TH1D("hPhoto","hPhoto",nBins, 0., 50.)
    hAxio = TH1D("hAxio","hAxio",nBins, 0., 50.)
    hTritium = TH1D("hTritium","hTritium",nBins, 0., 50.)
    hAxionFlux = TH1D("hAxionFlux","hAxionFlux",nBins, 0., 50.)
    hAxionConv = TH1D("hAxionConv","hAxionConv",nBins, 0., 50.)
    # hBremFlux = TH1D("hBremFlux","hBremFlux",nBins, 0., 50.)
    # hBremConv = TH1D("hBremConv","hBremConv",nBins ,0., 50.)

    for i in range(nBins):

        # make this energy the midpoint of the bin
        ene = float(i)/nBins/kevPerBin
        eneLo, eneHi = ene - kevPerBin/2., ene + kevPerBin/2.

        # get the average value of column 1 for a range of energies (column 0)
        idx = np.where((tritData[:,0] >= eneLo) & (tritData[:,0] <= eneHi))
        trit = np.mean(tritData[idx][:,1])
        if np.isnan(trit): trit = 0.
        hTritium.SetBinContent(i, trit)

        # this gives a "mean of empty slice" for the first entry, just suppress it
        idx = np.where((phoData[:,0] >= eneLo) & (phoData[:,0] <= eneHi))
        pho = 0.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore",category=RuntimeWarning)
            pho = np.mean(phoData[idx][:,1])
        if np.isnan(pho) or len(phoData[idx][:,1]) == 0: pho = 0.
        pho = pho * 1000. # kg
        hPhoto.SetBinContent(i, pho)

        idx = np.where((axData[:,0] >= eneLo) & (axData[:,0] <= eneHi))
        axFlux = np.mean(axData[idx][:,1])
        if np.isnan(axFlux): axFlux = 0.
        hAxionFlux.SetBinContent(i, axFlux)

        axo = pho * np.power(ene,2.) * 2.088e-5
        hAxio.SetBinContent(i, axo)

        axConv = axFlux * axo # scaling happens later
        hAxionConv.SetBinContent(i, axConv)

        # make brems curve: x^0.89 Exp[-0.7 x - 1.26 Sqrt[x]]
        # bremFlux = np.power(ene,0.89) * np.exp(-0.8 * ene - 1.26 * np.sqrt(ene))
        # bremConv = bremFlux * pho
        # hBremFlux.SetBinContent(i, bremFlux)
        # hBremConv.SetBinContent(i, bremConv)

        # print ene, pho, trit, axFlux, axConv, bremFlux, bremConv


    # -- diagnostic plots --

    # reproduce graham's thesis plot, pg. 105
    c = TCanvas("c","Bob Ross' Canvas",800,600)
    c.SetLogy()
    c.SetGrid()
    hAxionFlux.SetTitle(" ")

    hAxionFlux.Scale(1e19/(np.power(0.511e-10,2.) * 86400)) # to match graham's thesis plot
    # hAxionFlux.SetMinimum(2e32)
    # hAxionFlux.SetMaximum(1.2e35)
    hAxionFlux.GetXaxis().SetRangeUser(0.1,13.)
    hAxionFlux.GetXaxis().SetTitle("Energy (keV)")
    hAxionFlux.GetYaxis().SetTitle("Flux (cts / (cm^{2} s^{1} keV) )")
    hAxionFlux.SetLineColor(ROOT.kBlue)
    hAxionFlux.Draw("hist")
    c.Print("./plots/axionFlux.pdf")
    c.Clear()

    # check the units of the photoelectric cross section are in cm^2/g
    hPhoto.SetTitle(" ")
    hPhoto.GetXaxis().SetTitle("Energy (keV)")
    hPhoto.GetYaxis().SetTitle("#sigma_{pe} (cm^{2}/g)")
    hPhoto.Scale(1/1000.) # undo the conversion to kg
    # phoTest = hPhoto.GetBinContent(hPhoto.GetXaxis().FindBin(1.74)) # should be about 4412
    hPhoto.Draw("hist")
    c.Print("./plots/photoElectric.pdf")

    # axioelectric cross section
    hAxio.GetXaxis().SetRangeUser(0.5, 12.)
    hAxio.SetTitle(" ")
    hAxio.Scale(120.5 / 1000.) # convert to barns/atom (mucal)
    hAxio.GetXaxis().SetTitle("Energy (keV)")
    hAxio.GetYaxis().SetTitle("#sigma_{ae} (barns/atom)")
    hAxio.Draw("hist")
    c.Print("./plots/axioElectric.pdf")

    # final axion spectrum, with integral for expected counts
    hAxionConv.GetXaxis().SetRangeUser(0.5, 12.)
    hAxionConv.Scale((89.5 * 1e19 * 1000) / (np.power(0.511e-10, 2.)))

    hAxionConv.SetTitle(" ")
    hAxionConv.GetXaxis().SetTitle("Energy (keV)")
    hAxionConv.GetYaxis().SetTitle("Counts / keV")
    hAxionConv.Draw("hist")

    # ymin, ymax = 2e44, 5e47
    # hAxionConv.SetMinimum(ymin)
    # hAxionConv.SetMaximum(ymax)
    # eLoLine = TLine(eLo, ymin, eLo, ymax)
    # eLoLine.SetLineColor(ROOT.kRed)
    # eLoLine.SetLineWidth(2)
    # eLoLine.Draw("same")
    # eHiLine = TLine(eHi, ymin, eHi, ymax)
    # eHiLine.SetLineColor(ROOT.kRed)
    # eHiLine.SetLineWidth(2)
    # eHiLine.Draw("same")

    c.Print("./plots/axionConv.pdf")

    xax = hAxionConv.GetXaxis()
    N_expected = hAxionConv.Integral(xax.FindBin(eLo), xax.FindBin(eHi))
    # print "N_expected:",N_expected

    # reproduce the numbers in frank's tables
    eTest = 1.74
    phoTest = hPhoto.GetBinContent(hPhoto.GetXaxis().FindBin(eTest)) # should be about 4412 YES
    phoConv = phoTest * 72.64 / 6.022e23 # should be about 5.3e-19 YES
    axoTest = hAxio.GetBinContent(hAxio.GetXaxis().FindBin(eTest))
    axoTest = axoTest * (1000/120.5)

    (1000/120.5) / (6.022e23 * 1000)
    print axoTest


    # write to file.
    hPhoto.Write()
    hTritium.Write()
    hAxionFlux.Write()
    hAxionConv.Write()
    # hBremFlux.Write()
    # hBremConv.Write()
    f4.Close()

    return N_expected


class pkModel:
    def __init__(self,ene,name,fEnergy):
        """ TODO: add in constraints to the likelihood function instead of hard bounds
        https://root.cern.ch/root/html/tutorials/roofit/rf604_constraints.C.html
        """
        sig = getSigma(ene)
        muLo, muHi = ene * 0.99, ene * 1.01
        sgLo, sgHi = sig * 0.66, sig * 1.33
        ampLo, ampHi = -0.1, 500   # letting it go slightly negative helps the plc calculator
        if eHi > 10.: ampHi = 2000. # 68GeK is around 1700, all the rest are under 500

        # specifying a min/max value allows the var to float
        self.pkNum = ROOT.RooRealVar("amp-"+name, "amp-"+name, -0.1, 500.)
        self.pkMu = ROOT.RooRealVar("mu-"+name, "mu-"+name, ene, muLo, muHi)
        self.pkSig = ROOT.RooRealVar("sig-"+name, "sig-"+name, sig, sgLo, sgHi)

        # pdf and extended pdf
        self.pkGaus = ROOT.RooGaussian("gaus-"+name, "gaus-"+name, fEnergy, self.pkMu, self.pkSig)
        self.pkExt = ROOT.RooExtendPdf("ext-"+name, "ext-"+name, self.pkGaus, self.pkNum)

    def GetPkExt(self): return self.pkExt


def runFit():

    # load malbek data and create a workspace
    f1 = TFile("./data/malbek_data.root")
    t = f1.Get("malbek_wrt")
    fEnergy = ROOT.RooRealVar("energy_keV","Energy",eLo,eHi,"keV")
    fWeight = ROOT.RooRealVar("weight","weight",1,10,"")
    fData = ROOT.RooDataSet("data","data", t, ROOT.RooArgSet(fEnergy,fWeight),"","weight")
    fitWorkspace = ROOT.RooWorkspace("fitWorkspace","Fit Workspace")
    getattr(fitWorkspace,'import')(fEnergy)
    getattr(fitWorkspace,'import')(fData)
    getattr(fitWorkspace,'import')(fWeight)

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
            thisPk = pkModel(ene, key, fEnergy)
            peaks[key] = thisPk
            pdfList.add(thisPk.GetPkExt())

    # axion continuum
    axNum = ROOT.RooRealVar("amp-axion","amp-axion",100., 0., 10000.)
    # axNum = ROOT.RooRealVar("amp-axion","amp-axion",16.904) # hardcode the profile upper limit
    f2 = TFile("./data/inputHists.root")
    axTH1D = f2.Get("hAxionConv")
    axDataHist = ROOT.RooDataHist("ax", "ax", ROOT.RooArgList(fEnergy), RF.Import(axTH1D))
    fEnergy.setRange(eLo, eHi) # have to reset after loading a histo w/ different bounds
    axPdf = ROOT.RooHistPdf("axPdf","axPdf", ROOT.RooArgSet(fEnergy), axDataHist, 2)
    axExt = ROOT.RooExtendPdf("ext-axion","ext-axion", axPdf, axNum)
    pdfList.add(axExt)

    # create total model pdf
    model = ROOT.RooAddPdf("model","total pdf",pdfList)

    # run fitter
    minimizer = ROOT.RooMinimizer( model.createNLL(fData, RF.NumCPU(2,0), RF.Extended(True)) )
    minimizer.setPrintLevel(-1)
    minimizer.setStrategy(2)
    minimizer.migrad()
    fitResult = minimizer.save()

    # according to the internet, covQual==3 is a good indicator that it converged
    print "Fitter is done. Fit Cov Qual:", fitResult.covQual()

    # save workspace to a TFile
    getattr(fitWorkspace,'import')(fitResult)
    getattr(fitWorkspace,'import')(model)
    f2 = TFile("./data/fitWorkspace.root","RECREATE")
    fitWorkspace.Write()
    f2.Close()


def plotSpectrum():

    # load workspace
    f = TFile("./data/fitWorkspace.root")
    fitWorkspace = f.Get("fitWorkspace")
    fData = fitWorkspace.allData().front()
    fitResult = fitWorkspace.allGenericObjects().front()
    nPars = fitResult.floatParsFinal().getSize()
    fEnergy = fitWorkspace.var("energy_keV")
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
    print "-- SPECFIT RESULTS -- "
    print "%-10s = %.3f" % ("chiSq",chiSquare)
    fitValsFinal = getRooArgDict( fitResult.floatParsFinal() )

    for name in sorted(fitValsFinal):
        fitVal = fitValsFinal[name]
        if "amp-" in name:
            error = fitWorkspace.var(name).getError()
            print "%-10s = best %-7.3f  error %.3f (w/o profile)" % (name, fitVal, error)
        elif "mu-" in name:
            # compare the energy offset
            pkName = name[3:]
            pct = 100*(1 - fitVal/pkList[pkName])
            print "%-10s : fit %-6.3f  lit %-6.3f  (%.3f%%)" % (name, fitVal, pkList[pkName], pct)
        elif "sig-" in name:
            # compare the sigma difference
            pkName = name[4:]
            pct = 100*(1 - fitVal/getSigma(pkList[pkName]))
            print "%-10s : fit %-6.3f  func %-6.3f  (%.3f%%)" % (name, fitVal, getSigma(pkList[pkName]), pct)
        else:
            print "%s = %.4f" % (name, fitVal)
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


def plotProfiles():

    f = TFile("./data/fitWorkspace.root")
    fitWorkspace = f.Get("fitWorkspace")
    fData = fitWorkspace.allData().front()
    fEnergy = fitWorkspace.var("energy_keV")
    model = fitWorkspace.pdf("model")
    fitResult = fitWorkspace.allGenericObjects().front()
    fitValsFinal = getRooArgDict( fitResult.floatParsFinal() )
    print "Fit Cov Qual:", fitResult.covQual()

    nameList = ["amp-axion"]
    # nameList = []
    # for key in sorted(fitValsFinal):
    #     if "amp-" in key:
    #         nameList.append(key)

    print "Generating profiles ..."
    c = TCanvas("c","c",800,600)
    for name in nameList:

        fitVal = fitValsFinal[name]
        thisVar = fitWorkspace.var(name)

        # Brian & Clint method
        plc = RS.ProfileLikelihoodCalculator(fData, model, ROOT.RooArgSet(thisVar))
        plc.SetConfidenceLevel(0.90)
        interval = plc.GetInterval()
        lower = interval.LowerLimit(thisVar)
        upper = interval.UpperLimit(thisVar)
        print "%-10s = lo %-7.3f  best %-7.3f  hi %.3f" % (name, lower, fitVal, upper)
        plot = RS.LikelihoodIntervalPlot(interval)
        # p1, pltHi = fitVal - 1.5*(fitVal - lower), fitVal + 1.5*(upper - fitVal)
        # plot.SetRange(pltLo,pltHi)
        plot.SetTitle("%s: lo %.3f  fit %.3f  hi %.3f" % (name,lower,fitVal,upper))
        plot.SetNPoints(50)
        plot.Draw("tf1")

        """
        # Lukas method
        # note: lukas uses ModelConfig, to explicitly set observables, constraints, and nuisance parameters
        # https://root-forum.cern.ch/t/access-toy-datasets-generated-by-roostat-frequentistcalculator/24465/8

        mc = RS.ModelConfig('mc', fitWorkspace)
        mc.SetPdf( model )
        mc.SetParametersOfInterest( ROOT.RooArgSet(thisVar) )
        mc.SetObservables( ROOT.RooArgSet(fEnergy) )

        # lukas example
        # mc.SetConstraintParameters( ROOT.RooArgSet(mean, sigma) )
        # mc.SetNuisanceParameters( ROOT.RooArgSet(mean, sigma, n_bkg) )
        # mc.SetGlobalObservables( ROOT.RooArgSet(mean_obs, sigma_obs) )

        # need to make some RooArgSets from RooRealVars
        # constraints, nuisances, globalobs = ROOT.RooArgSet(), ROOT.RooArgSet(), ROOT.RooArgSet()
        # for parName in sorted(fitValsFinal):
        #     rrv = fitWorkspace.var(parName)
        #     if parName != name:
        #         constraints.add(rrv)
        #         # .. etc.

        # pl = RS.ProfileLikelihoodCalculator(fData, mc)
        # pl.SetConfidenceLevel(0.683) # lukas used 0.90
        # interval = pl.GetInterval()
        # plot = RS.LikelihoodIntervalPlot(interval)
        # plot.SetNPoints(50)
        # plot.Draw("")
        """
        c.Print("./plots/profile_%s.pdf" % name)


def calculate_g_ae():

    N_expected = getSpecPDFs()
    N_observed = 16.904 # profile U.L.
    # N_observed = 40.8
    print "g_ae U.L.:", np.power(N_observed/N_expected, 1./4.)


def plotContours():

    par1, par2 = "amp-41Ca", "amp-36Cl"

    f = TFile("./data/fitWorkspace.root")
    fitWorkspace = f.Get("fitWorkspace")
    fData = fitWorkspace.allData().front()
    model = fitWorkspace.pdf("model")

    # redeclare the minimizer since it can't be saved into the workspace for some reason
    minimizer = ROOT.RooMinimizer( model.createNLL(fData, RF.NumCPU(2,0), RF.Extended(True)) )
    # minimizer.setMinimizerType("Minuit2")
    minimizer.setPrintLevel(-1)
    minimizer.setStrategy(2)
    minimizer.migrad()
    fitResult = minimizer.save()

    c = TCanvas("c","c",1100,800)
    fContour = minimizer.contour(fitWorkspace.var(par1), fitWorkspace.var(par2), 1,2,3)
    fContour.SetTitle("Contour of %s vs %s" % (par1, par2))

    # set plot range
    mux = fitWorkspace.var(par1).getValV()
    upx = fitWorkspace.var(par1).getErrorHi()
    lox = fitWorkspace.var(par1).getErrorLo()
    muy = fitWorkspace.var(par2).getValV()
    upy = fitWorkspace.var(par2).getErrorHi()
    loy = fitWorkspace.var(par2).getErrorLo()
    fContour.GetXaxis().SetRangeUser(mux + 4*lox, mux + 4*upx)
    fContour.GetYaxis().SetRangeUser(muy + 4*loy, muy + 4*upy)
    fContour.Draw()

    c.Print("./plots/contour-%s-vs-%s.pdf" % (par1, par2))


def plotMCSpectrum():

    # load workspace
    f = TFile("./data/fitWorkspace.root","UPDATE")
    fitWorkspace = f.Get("fitWorkspace")
    fData = fitWorkspace.allData().front()
    fitResult = fitWorkspace.allGenericObjects().front()
    fEnergy = fitWorkspace.var("energy_keV")
    model = fitWorkspace.pdf("model")
    pdfList = fitWorkspace.allPdfs() # RooArgSet

    # create a toy dataset with the same number of expected events as the best-fit model
    # by sampling from the best-fit model PDF.
    # makes it look like the real data points got slightly randomized.
    print model.expectedEvents(pdfList)
    fMCData = model.generate(ROOT.RooArgSet(fEnergy), RF.Name("Toy MC Data"), RF.Extended())
    # getattr(fitWorkspace,'import')(fMCData)
    # fitWorkspace.Write("",TObject.kOverwrite)
    # f.Close()

    # -- plot the data vs. the toy data --
    c = TCanvas("c","c",1100,800)
    binSize = 0.2
    nBins = int((eHi-eLo)/binSize + 0.5)
    fSpec = fEnergy.frame(RF.Range(eLo,eHi), RF.Bins(nBins))
    fData.plotOn(fSpec)
    fMCData.plotOn(fSpec, RF.MarkerColor(ROOT.kBlue))
    model.plotOn(fSpec, RF.LineColor(ROOT.kRed), RF.Name("FullModel"))
    fSpec.SetTitle(" ")
    fSpec.Draw()
    c.Print("./plots/toyMC.pdf")


def runMCStudy():

    # load workspace
    f = TFile("./data/fitWorkspace.root")
    fitWorkspace = f.Get("fitWorkspace")
    fEnergy = fitWorkspace.var("energy_keV")
    model = fitWorkspace.pdf("model")
    fitResult = fitWorkspace.allGenericObjects().front()
    fitValsFinal = getRooArgDict( fitResult.floatParsFinal() )
    pdfList = fitWorkspace.allPdfs() # RooArgSet
    nExpected = int(model.expectedEvents(pdfList))

    # thank you, lindsey gray. http://www.hep.wisc.edu/~lgray/toyMC.py
    nSamples = 10
    args = [
        # RF.Binned(True),
        RF.Silence(),
        RF.FitOptions( RF.PrintEvalErrors(0), RF.NumCPU(2) ) # not sure numCPU is doing anything
        ]
    mcStudy = ROOT.RooMCStudy(model, ROOT.RooArgSet(fEnergy), *args)
    start = time.time()
    mcStudy.generateAndFit(nSamples, nExpected) # nSamples (nToys), nEvtPerSample (nEvents)
    print "MC Study Complete. %i samples, %.2f seconds." % (nSamples, time.time()-start)


    parName = "amp-36Cl"
    parVal = fitResult.floatParsFinal().find(parName).getValV()
    parErr = fitResult.floatParsFinal().find(parName).getError()
    parVar = fitWorkspace.var(parName)

    # make 4 rooplot frames
    fParam = mcStudy.plotParam( parVar, RF.Bins(50))
    fError = mcStudy.plotError( parVar, RF.FrameRange(parErr-0.5*parErr, parErr+0.5*parErr), RF.Bins(50) )
    fPull = mcStudy.plotPull( parVar, RF.FrameRange(-5, 5), RF.Bins(50) )
    fNLL = mcStudy.plotNLL( RF.Bins(50) )

    line = TLine()
    line.SetLineColor(ROOT.kBlue)
    line.SetLineWidth(2)
    line.SetLineStyle(3)

    c = TCanvas("c","c",1100,800)
    c.Divide(2,2)
    c.cd(1)
    fParam.SetTitle(parName)
    fParam.Draw()
    line.DrawLine(parVal, fParam.GetMinimum(), parVal, fParam.GetMaximum()) # best fit

    c.cd(2)
    fError.SetTitle(parName + " Error")
    fError.Draw()
    line.DrawLine(parErr, fError.GetMinimum(), parErr, fError.GetMaximum())

    c.cd(3)
    fPull.SetTitle(parName + " Pull")
    fPull.Draw()

    c.cd(4)
    fNLL.SetTitle("-log(likelihood)")
    fNLL.Draw()
    line.DrawLine( fitResult.minNll(), fNLL.GetMinimum(), fitResult.minNll(), fNLL.GetMaximum() )

    c.Print("plots/mcStudy.pdf")


    """
    // Study done for Jason to convince him I know what I'm doing
    // double parVal0 = dynamic_cast<RooRealVar*>(fFitResult->floatParsFinal().find(Form("%s", argN.c_str())))->getValV();
    // double parErrHi0 = dynamic_cast<RooRealVar*>(fFitResult->floatParsFinal().find(Form("%s", argN.c_str())))->getErrorHi();
    // double parErrLo0 = dynamic_cast<RooRealVar*>(fFitResult->floatParsFinal().find(Form("%s", argN.c_str())))->getErrorLo();

    // TH1D *hMean = new TH1D("hMean", "Tritium Mean", 200, parVal0+2.5*parErrLo0, parVal0+2.5*parErrHi0);
    // TH1D *hErrLo = new TH1D("hErrLo", "Tritium Error Low", 200, parErrLo0+parErrLo0/4, parErrLo0+parErrHi0/4);
    // TH1D *hErrHi = new TH1D("hErrHi", "Tritium Error High", 200, parErrHi0+parErrLo0/4, parErrHi0+parErrHi0/4);

    // for(int i = 0; i < nMC; i++)
    // {
    //     const RooFitResult *fFitMCResult = fMCStudy->fitResult(i);
    //     double parVal2 = dynamic_cast<RooRealVar*>(fFitMCResult->floatParsFinal().find(Form("%s", argN.c_str())))->getValV();
    //     double parErrHi2 = dynamic_cast<RooRealVar*>(fFitMCResult->floatParsFinal().find(Form("%s", argN.c_str())))->getErrorHi();
    //     double parErrLo2 = dynamic_cast<RooRealVar*>(fFitMCResult->floatParsFinal().find(Form("%s", argN.c_str())))->getErrorLo();

    //     hMean->Fill(parVal2);
    //     hErrLo->Fill(parErrLo2);
    //     hErrHi->Fill(parErrHi2);
    // }

    // TCanvas *cM = new TCanvas("cM", "cM", 800, 600);
    // hMean->Draw();
    // l1.DrawLine(parVal0, 0, parVal0, hMean->GetBinContent(hMean->GetMaximumBin()) );
    // cM->SaveAs("./plots/MeanTest.pdf");

    // TCanvas *cLo = new TCanvas("cLo", "cLo", 800, 600);
    // hErrLo->Draw();
    // l1.DrawLine(parErrLo0, 0, parErrLo0, hErrLo->GetBinContent(hErrLo->GetMaximumBin()));
    // cLo->SaveAs("./plots/MeanTest_Low.pdf");

    // TCanvas *cHi = new TCanvas("cHi", "cHi", 800, 600);
    // hErrHi->Draw();
    // l1.DrawLine(parErrHi0, 0, parErrHi0, hErrHi->GetBinContent(hErrHi->GetMaximumBin()));
    // cHi->SaveAs("./plots/MeanTest_High.pdf");

    // double NLLmean = dynamic_cast<RooRealVar*>(fFitResult->floatParsFinal().find("NLL"))->getValV();
    // double NLLmean = dynamic_cast<RooRealVar*>(fFitResult->floatParsFinal().find("NLL"))->getError();

    /*
    // Extract nLL variables from MCStudy -- maybe useful later
    RooDataSet MCFitData = fMCStudy->fitParDataSet();
    // RooDataSet *NLL = static_cast<RooDataSet*>(MCFitData.reduce(RooArgSet()));
    MCFitData.Print("v");
    const RooArgSet* row = MCFitData.get();
    row->Print("v");
    shared_ptr<TCanvas> cMCNLL( make_shared<TCanvas>("cMCNLL", "cMCNLL", 1100, 800) );
    RooRealVar *NLL = static_cast<RooRealVar*>(row->find("NLL"));
    RooFormulaVar sNLL("sNLL", "Shifted NLL", Form("NLL - %f", fFitResult->minNll()), *NLL);
    RooRealVar *SNLL = static_cast<RooRealVar*>(MCFitData.addColumn(sNLL));

    RooPlot* frameMC = SNLL->frame(Range((0.2*fFitResult->minNll()), -(0.2*fFitResult->minNll()) ));
    MCFitData.plotOn(frameMC);
    frameMC->Draw();
    cMCNLL->SaveAs("Test.pdf");
    */
    }
    """


if __name__=="__main__":
    main(sys.argv[1:])