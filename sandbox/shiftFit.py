#!/usr/bin/env python
"""

"""
import sys, time, random
sys.argv.append("-b") # kill all interactive crap
import ROOT
from ROOT import TFile, TTree, TList, TCanvas, TH1D, TLegend, gStyle, gROOT, TFeldmanCousins
from ROOT import RooFit as RF
from ROOT import RooStats as RS
import numpy as np
from array import array

def main():

    gStyle.SetOptStat(0)
    gROOT.ProcessLine("gErrorIgnoreLevel = 3001;") # suppress ROOT messages
    # gROOT.ProcessLine("RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);")
    gROOT.ProcessLine("RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);")
    ROOT.Math.MinimizerOptions.SetDefaultMinimizer("Minuit") # the PLC doesn't work w/ minuit2

    global eLo, eHi
    # eLo, eHi = 2.22, 2.70 # results from splitData()
    eLo, eHi = 2.22, 3. # help constrain the flat BG

    # compareData()
    # splitData()
    # runFit()
    # plotFit()
    calculate_g_ae()
    # plotProfiles()



# ==========================================================================
def compareData():
    """ Compare Frank's table with Graham's data.
    NOTES:
    fs idx 10 = 1.74 kev
    fs idx 12 = 1.84 kev
    fs idx 24 = 2.31 kev
    fs idx 28 = 2.46 kev
    linear fit gives 0.0396 ~ 0.04 bins/kev.
    that means that fLo = 1.84 - n*0.04, n = 12, *** fLo = 1.36 ***
    and fHi = 1.36 + len(fs)*0.04, len(fs)=39,   *** fHi = 2.92 ***
    """
    malbekExposure = 89.5 # kg-d

    # -- frank's table --
    col1 = [2.40, 3.07, 2.23, 4.34, 1.83, 2.48, 2.79, 2.23, 2.23, 1.96, 2.51, 1.95, 1.96, 2.51, 2.79, 1.68, 2.79, 2.51, 4.75, 1.40, 2.23]
    col2 = [5.86, 1.92, 2.40, 3.07, 2.23, 4.34, 1.84, 2.48, 2.79, 2.23, 2.23, 1.95, 2.51, 1.96, 1.96, 2.51, 2.79, 1.68, 2.79, 2.51, 4.75]
    col3 = [1.96, 2.51, 2.79, 1.68, 2.79, 2.51, 4.75, 1.40, 2.23, 3.07, 4.47, 3.63, 1.96, 3.07, 2.23, 2.79, 1.12, 3.91, 1.40, 1.68, 3.07]
    col4 = [2.79, 2.51, 4.75, 1.40, 2.23, 3.07, 4.47, 3.63, 1.95, 3.07, 2.23, 2.79, 1.12, 3.91, 1.40, 1.68, 3.07, 1.12, 2.23, 2.23, 2.79]
    col5 = [13.01, 10.01, 12.17, 10.49, 9.08, 12.40, 13.85, 9.74, 9.20, 10.33, 11.44, 10.32, 7.55, 11.45, 8.38, 8.66, 9.77, 9.22, 11.17,
            7.82, 12.84] # verified the sum is correct.
    axionCol = 10 # 0-indexed

    # frank's spectrum (from comparing when #'s in columns 1-4 repeat)
    fs = [5.86, 1.92, 2.40, 3.07, 2.23, 4.34, 1.84, 2.48, 2.79, 2.23, 2.23, 1.95, 2.51, 1.96, 1.96, 2.51, 2.79, 1.68, 2.79, 2.51, 4.75,
          1.40, 2.23, 3.07, 4.47, 3.63, 1.96, 3.07, 2.23, 2.79, 1.12, 3.91, 1.40, 1.68, 3.07, 1.12, 2.23, 2.23, 2.79]

    # frank didn't adjust for kev/bin of 0.04 keV bins times exposure (0.04 * 89.5 = 3.58)
    fs = [val * 3.58 for val in fs]

    # set binning and endpoints, and fill a histogram
    fBins = len(fs) # 39
    binSize = 0.04
    fLo = 1.84 - 12 * binSize # 1.3648
    fHi = fLo + fBins * binSize # 2.9092
    hF = TH1D("hF","",fBins,fLo,fHi)
    for i in range(fBins):
        hF.SetBinContent(i,fs[i])
    hF.SetLineColor(ROOT.kRed)
    hF.Scale(1./malbekExposure)
    print "Frank's data: bins %d  eLo %.3f  eHi %.3f" % (fBins, fLo, fHi)


    # -- graham's data --
    f1 = TFile("./data/malbek_data.root")
    t1 = f1.Get("malbek_wrt")

    # duplicate the frank region only
    hG1 = TH1D("hG1","",fBins,fLo,fHi)
    t1.Project("hG1","energy_keV","weight")
    hG1.Scale(1./malbekExposure)
    hG1.SetLineColor(ROOT.kGreen)

    # get histo for a larger energy region
    eLo, eHi = 1.0, 3.5
    nBins = int((eHi-eLo)/binSize + 0.5)
    hG2 = TH1D("hG2","",nBins, eLo, eHi)
    t1.Project("hG2","energy_keV","weight")
    hG2.SetLineColor(ROOT.kBlue)
    hG2.GetXaxis().SetTitle("Energy (keV)")
    hG2.Scale(1./malbekExposure)

    # -- plot all three --
    c = TCanvas("c","c",800,600)
    gStyle.SetOptStat(0)

    hG2.Draw("hist")
    hG1.Draw("hist same")
    hF.Draw("hist same")

    leg = TLegend(0.5,0.7,0.85,0.85)
    leg.AddEntry(hF,"frank, %.2f-%.2f, %.3f b/kev" % (fLo, fHi, binSize),"l")
    leg.AddEntry(hG1,"graham-1, %.2f-%.2f" % (fLo, fHi),"l")
    leg.AddEntry(hG2,"graham-2, %.2f-%.2f" % (eLo, eHi),"l")
    leg.Draw("same")

    c.Print("./plots/testSpec.pdf")


def splitData():
    """ Sometimes this segfaults.  No idea why.  Just re-run and it's fine. """

    f1 = TFile("./data/malbek_data.root")
    tree0 = f1.Get("malbek_wrt")

    axPeaks = [1.739, 1.836, 2.307, 2.464]

    # declare as arrays s/t ROOT branches work.  access values with e.g. ene[0] = 1.
    ene0, wt0, shift = array('d',[0.]), array('d',[0.]), array('d',[0.])

    tree0.SetBranchAddress("energy_keV",ene0)
    tree0.SetBranchAddress("weight",wt0)

    tList = TList()
    tVec = [0.,0.,0.,0.]
    tmp = TFile("./data/malbek_splitData.root","RECREATE")

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

        print "axPeaks: %d  tree %d - %d entries  shift %.3f keV" % (len(axPeaks),len(tVec),tVec[i].GetEntries(),shift)
        tList.Add(tVec[i])
        tVec[i].Write()

    tree1 = ROOT.TTree.MergeTrees(tList)
    tree1.SetName("mergeTree")
    tree1.SetTitle("mergeTree")
    print "Trees merged, with %d entries total." % tree1.GetEntries()
    tree1.Write()

    # -- create diagnostic shifted & merged spectrum --
    malbekExposure = 89.5
    binSize = 0.04
    eLo, eHi = 0.8, 5.
    nBins = int((eHi-eLo)/binSize + 0.5)

    h0 = H1D(tree0,nBins,eLo,eHi,"energy_keV","weight")
    h0.SetLineColor(ROOT.kBlack)
    h0.SetLineWidth(2)
    h0.Scale(1./malbekExposure)

    h1 = H1D(tVec[0],nBins,eLo,eHi,"energy_keV","weight")
    h1.SetLineColorAlpha(ROOT.kBlue,0.5)
    h1.Scale(1./malbekExposure)

    h2 = H1D(tVec[1],nBins,eLo,eHi,"energy_keV","weight")
    h2.SetLineColorAlpha(ROOT.kGreen,0.5)
    h2.Scale(1./malbekExposure)

    h3 = H1D(tVec[2],nBins,eLo,eHi,"energy_keV","weight")
    h3.SetLineColorAlpha(ROOT.kMagenta,0.5)
    h3.Scale(1./malbekExposure)

    h4 = H1D(tVec[3],nBins,eLo,eHi,"energy_keV","weight")
    h4.SetLineColor(ROOT.kCyan)
    h4.Scale(1./malbekExposure)

    hAdd = TH1D("hAdd","",nBins,eLo,eHi)
    hAdd.Add(h1)
    hAdd.Add(h2)
    hAdd.Add(h3)
    hAdd.Add(h4)
    hAdd.SetLineColor(ROOT.kBlack)
    hAdd.SetLineWidth(2)
    hAdd.SetMinimum(0)
    hAdd.SetMaximum(2)
    ymax = hAdd.GetMaximum()
    hAdd.GetXaxis().SetTitle("Energy (keV)")
    hAdd.GetYaxis().SetTitle("Counts / kg-d")

    c = TCanvas("c","Bob Ross's Canvas",1100,800)
    hAdd.Draw("hist")
    h0.Draw("hist same")
    h1.Draw("hist same")
    h2.Draw("hist same")
    h3.Draw("hist same")
    h4.Draw("hist same")

    # Draw lines around the gaussian fit ROI, using 2 sigma of the energy resolution
    fitWin = 3 * getSigma(axPeaks[3])
    print "Malbek resolution at %.2f keV is %.2f" % (axPeaks[3], getSigma(axPeaks[3]))
    print "Set 3-sigma fit region to %.2f - %.2f" % (axPeaks[3] - fitWin, axPeaks[3] + fitWin)

    l1 = ROOT.TLine(axPeaks[3],0.,axPeaks[3],ymax)
    l1.SetLineColor(ROOT.kRed)
    l1.SetLineWidth(2)
    l1.Draw("same")

    l2 = ROOT.TLine(axPeaks[3]+fitWin,0.,axPeaks[3]+fitWin,ymax)
    l2.SetLineColorAlpha(ROOT.kRed, 0.5)
    l2.SetLineWidth(2)
    l2.Draw("same")

    l3 = ROOT.TLine(axPeaks[3]-fitWin,0.,axPeaks[3]-fitWin,ymax)
    l3.SetLineColorAlpha(ROOT.kRed, 0.5)
    l3.SetLineWidth(2)
    l3.Draw("same")

    leg = ROOT.TLegend(0.6,0.6,0.85,0.85)
    leg.AddEntry(hAdd,"Shifted+Summed","l")
    leg.AddEntry(h4,"peak-%.3f (unshifted)" % (axPeaks[3]),"l")
    leg.AddEntry(h3,"peak-%.3f" % (axPeaks[2]),"l")
    leg.AddEntry(h2,"peak-%.3f" % (axPeaks[1]),"l")
    leg.AddEntry(h1,"peak-%.3f" % (axPeaks[0]),"l")
    leg.AddEntry(l1,"3-sigma fit region","l")
    leg.Draw("same")

    c.Print("./plots/shiftSpec.pdf")


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
        # self.pkNum = ROOT.RooRealVar("amp-"+name, "amp-"+name, -0.01, 500.)
        self.pkNum = ROOT.RooRealVar("amp-"+name, "amp-"+name, 23.45)
        self.pkMu = ROOT.RooRealVar("mu-"+name, "mu-"+name, ene, muLo, muHi)
        self.pkSig = ROOT.RooRealVar("sig-"+name, "sig-"+name, sig, sgLo, sgHi)

        # pdf and extended pdf
        self.pkGaus = ROOT.RooGaussian("gaus-"+name, "gaus-"+name, fEnergy, self.pkMu, self.pkSig)
        self.pkExt = ROOT.RooExtendPdf("ext-"+name, "ext-"+name, self.pkGaus, self.pkNum)

    def GetPkExt(self): return self.pkExt


def runFit():
    """ model: flat BG + gaussian peak at 2.464 keV. """

    axPeak = 2.464

    # load data
    f1 = TFile("./data/malbek_splitData.root")
    t = f1.Get("mergeTree")
    fEnergy = ROOT.RooRealVar("energy_keV","Energy",eLo,eHi,"keV")
    fWeight = ROOT.RooRealVar("weight","weight",1,10,"")
    fData = ROOT.RooDataSet("data","data", t, ROOT.RooArgSet(fEnergy,fWeight),"","weight")
    fitWorkspace = ROOT.RooWorkspace("fitWorkspace","Fit Workspace")
    getattr(fitWorkspace,'import')(fEnergy)
    getattr(fitWorkspace,'import')(fData)
    getattr(fitWorkspace,'import')(fWeight)

    # background model
    pdfList = ROOT.RooArgList("shapes")

    # flat background function
    bkgNum = ROOT.RooRealVar("bkgNum","bkgNum",1.,5000.)
    bkgLin = ROOT.RooPolynomial("bkgLin", "bkgLin", fEnergy, ROOT.RooArgList() )
    bkgExt = ROOT.RooExtendPdf("bkgExt", "bkgExt", bkgLin, bkgNum);
    pdfList.add(bkgExt);

    axionSumPeak = pkModel(axPeak, "axion", fEnergy)
    pdfList.add(axionSumPeak.GetPkExt())

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
    f2 = TFile("./data/splitWorkspace.root","RECREATE")
    fitWorkspace.Write()
    f2.Close()


def plotFit():

    axPeak = 2.464

    # load workspace
    f = TFile("./data/splitWorkspace.root")
    fitWorkspace = f.Get("fitWorkspace")
    fData = fitWorkspace.allData().front()
    fitResult = fitWorkspace.allGenericObjects().front()
    nPars = fitResult.floatParsFinal().getSize()
    fEnergy = fitWorkspace.var("energy_keV")
    modelPDF = fitWorkspace.pdf("model")
    # fitWorkspace.Print()

    pdfNames = ["ext-axion"]

    # -- create spectrum rooplot --
    nCol = float(gStyle.GetNumberOfColors())
    binSize = 0.04
    nBins = int((eHi-eLo)/binSize + 0.5)
    fSpec = fEnergy.frame(RF.Range(eLo,eHi), RF.Bins(nBins))
    fData.plotOn(fSpec)
    modelPDF.plotOn(fSpec, RF.LineColor(ROOT.kRed), RF.Name("FullModel"))
    chiSquare = fSpec.chiSquare(nPars)
    modelPDF.plotOn(fSpec, RF.Components(pdfNames[0]), RF.LineColor(ROOT.kBlue), RF.Name(pdfNames[0]))

    # -- make a fake Gaussian --
    # gaus = ROOT.TF1("g","gaus",-3, eHi)
    # gaus.SetParameters(10., axPeak, getSigma(axPeak))
    # gaus.SetParameter(0,10.)
    # gaus.SetParameter(1, axPeak)
    # gaus.SetParameter(2, getSigma(axPeak))
    # erf
    # rrvGaus = ROOT.RooRealVar("ax","ax",-3.,3.)
    # rarGaus = RF.bindFunction("gaus", ROOT.TMath.Erf, rrvGaus)
    # rarGaus.Print()
    # frame2 = rrvGaus.frame(RF.Title("mygaus"))
    # rarGaus.plotOn(frame2, RF.LineColor(ROOT.kGreen), RF.LineStyle(ROOT.kDashed), RF.Name("axGaus"))
    # c0 = TCanvas("c0","",800,600)
    # frame2.Draw()
    # c0.Print("./plots/testGaus.pdf")


    # -- print fit results --
    print "-- SHIFTFIT RESULTS -- "
    print "%-10s = %.3f" % ("chiSq",chiSquare)
    fitValsFinal = getRooArgDict( fitResult.floatParsFinal() )

    bkgVal = 0.
    for name in sorted(fitValsFinal):
        fitVal = fitValsFinal[name]
        if "amp-" in name:
            error = fitWorkspace.var(name).getError()
            print "%-10s = best %-7.3f  error %.3f (w/o profile)" % (name, fitVal, error)
        elif "mu-" in name:
            # compare the energy offset
            pkName = name[3:]
            pct = 100*(1 - fitVal/axPeak)
            print "%-10s : fit %-6.3f  lit %-6.3f  (%.3f%%)" % (name, fitVal, axPeak, pct)
        elif "sig-" in name:
            # compare the sigma difference
            pkName = name[4:]
            pct = 100*(1 - fitVal/getSigma(axPeak))
            print "%-10s : fit %-6.3f  func %-6.3f  (%.3f%%)" % (name, fitVal, getSigma(axPeak), pct)
        else:
            print "%s = %.4f (%.4f)" % (name, fitVal, fitVal/nBins)
            bkgVal = fitVal/nBins
            continue

    # -- make spectrum plot --
    c = TCanvas("c","Bob Ross's Canvas", 1400, 1000)
    c.SetRightMargin(0.2)
    fSpec.SetTitle(" ")
    fSpec.Draw()

    ymax = fSpec.GetMaximum()
    l1 = ROOT.TLine(axPeak,0.,axPeak,ymax)
    l1.SetLineColor(ROOT.kBlue)
    l1.SetLineWidth(2)
    l1.Draw("same")

    leg = TLegend(0.83,0.6,0.97,0.9)
    leg.AddEntry(fSpec.findObject("FullModel"),"model #chi^{2}=%.3f" % chiSquare,"l")
    leg.AddEntry(fSpec.findObject("FullModel"),"cts/bin=%.2f" % bkgVal, "")
    leg.AddEntry(fSpec.findObject(pdfNames[0]),"axion gaussian", "l")
    leg.AddEntry(l1,"axion-%.2f" % axPeak, "l")
    leg.Draw("same")

    c.Print("./plots/shiftFit.pdf")

    # -- get FC Limit --
    # Calculate the confidence interval for the axion peak, which is too low to use profile likelihood.
    # TFeldmanCousins version:
    # https://root.cern.ch/root/html/tutorials/math/FeldmanCousins.C.html
    # RooStats version:
    # https://root.cern.ch/root/html/tutorials/roostats/StandardFeldmanCousinsDemo.C.html
    # FC Paper: https://arxiv.org/pdf/physics/9711021.pdf
    f = TFeldmanCousins(0.95)
    N_obs = 1.
    N_bkg = bkgVal/3.
    ul = f.CalculateUpperLimit(N_obs, N_bkg)
    ll = f.GetLowerLimit()
    print "For %.2f events observed, and %.2f background events," % (N_obs, N_bkg)
    print "F-C method gives UL: %.2f and LL %.2f (90%% CL)" % (ul, ll)


def calculate_g_ae():
    """ Repeat Frank's MALBEK calculation w/ his data and binning. """

    # axion inputs
    axPeaks = [1.739, 1.836, 2.307, 2.464]             # keV
    axFlux = [4.95e+38, 4.95e+38, 3.94e+38, 2.27e+38]  # cm^2/day
    gePhoto = [5.32e-19, 4.60e-19, 2.47e-19, 2.06e-19] # cm^2/atom
    def sigAe(E,idx): return 2.088e-5 * np.power(E,2.) * gePhoto[idx]

    # exposure and expected counts
    N_avogadro = 6.0221409e+23
    ge_molar_mass = 5.5624 # from 404 pm 15 grams Ge
    malbek_livetime = 221.5 # days
    exposure = N_avogadro * ge_molar_mass * malbek_livetime # atom-days
    N_expected = exposure * sum(axFlux[i] * sigAe(axPeaks[i],i) for i in range(4))

    # trick to make list printing prettier
    class prettyfloat(float):
        def __repr__(self): return "%.2e" % self

    # check frank's tables
    print "Tab3-Col1", map(prettyfloat, [np.power(ax, 2.) * 2.088e-5 for ax in axPeaks])
    print "Tab3-Col3", map(prettyfloat, [sigAe(axPeaks[i],i) for i in range(4)])
    rates = [axFlux[i] * sigAe(axPeaks[i],i) for i in range(4)]
    print "Tab4-Col4:", map(prettyfloat, rates)
    print "Total 4-pk rate: %.2e" % sum(rates)
    print "Expected axion counts: %.2e" % N_expected

    # what we saw in the data, to a 95% confidence interval
    # N_observed = 53.76 # his method
    # N_observed = 10. # a more reasonable guess
    # N_observed = 1.26 # feldman-cousins fit result
    N_observed = 23.45 # profile result
    print "Observed counts: ", N_observed

    # g_ae must be less than this value
    g_ae_upper = np.power(N_observed / N_expected, 1./4.)
    print "Frank's upper bound g_ae: %.2e" % g_ae_upper

    print "MALBEK resolution at %.2f keV is %.2f keV" % (axPeaks[3], getSigma(axPeaks[3]))


def plotProfiles():

    f = TFile("./data/splitWorkspace.root")
    fitWorkspace = f.Get("fitWorkspace")
    fData = fitWorkspace.allData().front()
    fEnergy = fitWorkspace.var("energy_keV")
    model = fitWorkspace.pdf("model")
    fitResult = fitWorkspace.allGenericObjects().front()
    fitValsFinal = getRooArgDict( fitResult.floatParsFinal() )
    print "Fit Cov Qual:", fitResult.covQual()

    name = "amp-axion"

    print "Generating profile ..."
    c = TCanvas("c","c",800,600)
    fitVal = fitValsFinal[name]
    thisVar = fitWorkspace.var(name)

    # Brian & Clint method
    plc = RS.ProfileLikelihoodCalculator(fData, model, ROOT.RooArgSet(thisVar))
    plc.SetConfidenceLevel(0.683)
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
    c.Print("./plots/shiftProfile_%s.pdf" % name)


# ==========================================================================
def H1D(tree,bins,xlo,xhi,drawStr,cutStr,xTitle="",yTitle=""):
    nameStr = str(random.uniform(1.,2.))
    h1 = TH1D(nameStr,nameStr,bins,xlo,xhi)
    tree.Project(nameStr,drawStr,cutStr)
    h1.SetTitle("")
    if xTitle!="": h1.GetXaxis().SetTitle(xTitle)
    if yTitle!="": h1.GetYaxis().SetTitle(yTitle)
    return h1


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


if __name__=="__main__":
    main()
