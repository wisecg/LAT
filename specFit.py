#!/usr/bin/env python3
import sys, warnings, time
import numpy as np
from scipy.stats import chi2
from scipy.interpolate import spline
import waveLibs as wl
import dsi
import matplotlib.pyplot as plt
plt.style.use('%s/pltReports.mplstyle' % dsi.latSWDir)
sys.argv.append("-b")
import ROOT
from ROOT import TFile
from ROOT import RooFit as RF
from ROOT import gROOT
gROOT.ProcessLine("gErrorIgnoreLevel = 3001;")
gROOT.ProcessLine("RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);")

def main(argv):
    initialize(makePlots=False)

    # loadDataMJD()
    # getUnscaledPDFs(makePlots=True)
    # plotPDFs()
    # testFunc()
    # getHistList()

    # fitModel(makePlots=True)
    # plotFit(plotRate=False, plotProfileResults=True)
    # getProfile()
    # plotProfile(makePlots=True)
    # combineProfileLoHi()

    # getProfileM1()
    # plotProfileM1()
    # combineProfileM1()

    # hadronicCurve()

    gaeProj()


def initialize(makePlots=False):

    global dsList, enr, eff, res, simEffCorr, use90, eLo, eHi, epb, pLo, pHi, ppb, nB, nBP
    global bkgModelHists, bkgModelPeaks, profileVars, bkgVals, sigLabels
    global effLim, effMax, xEff, detEff, dsExpo, detExp, bkgModelPeaks
    global useLoEff, useHiEff, constrainPeak, constrainMu

    # dsList = [0,1,2,3,4,"5A","5B","5C"]
    # dsList = [1,2,3,4,"5A","5B","5C"]
    dsList = [1,2,3,4,"5B","5C"]
    # dsList = [0]

    enr = True

    # special switches
    eff = True # should almost always be true
    res = True # should almost always be true
    simEffCorr = False
    use90 = False
    useLoEff = False # this is the one that really makes us eat shit
    useHiEff = False
    constrainPeak = True # if True, will constrain width of all peaks in "bkgModelPeaks" to a gaussian
    constrainMu = True # use only for 14.4 study

    print(dsList)
    print("Enriched" if enr else "Natural")

    # full spectrum axion fit
    eLo, eHi, epb = 1.5, 20, 0.2
    bkgModelHists = ["trit","flat","49V","54Mn","55Fe","57Co","65Zn","68Ga","68Ge","axion"]
    bkgModelPeaks = []
    profileVars = ["axion"]

    # reference fit
    # eLo, eHi, epb = 1.5, 20, 0.2
    # bkgModelHists = ["trit","flat","49V","54Mn","55Fe","57Co","65Zn","68Ga","68Ge"]
    # bkgModelPeaks = []
    # profileVars = []

    # # enr/nat resolution study
    # eLo, eHi, epb = 5, 20, 0.2
    # bkgModelHists = ["trit","flat"]
    # # bkgModelPeaks = ["54Mn","55Fe","57Co","65Zn","68Ga","68Ge"]
    # bkgModelPeaks = ["55Fe","65Zn","68Ga","68Ge"]
    # profileVars = []

    # # 14.4 keV axion
    # eLo, eHi, epb = 5, 30, 0.2
    # # bkgModelHists = ["trit","flat","54Mn","55Fe","57Co","65Zn","68Ga","68Ge","axFe_M1"] # fixed peak
    # # bkgModelPeaks = []
    # bkgModelHists = ["trit","flat","54Mn","55Fe","57Co","65Zn","68Ga","68Ge"] # floating width peak
    # bkgModelPeaks = ["axFe_M1"]
    # profileVars = ["axFe_M1"]

    # binning of continuum histograms
    pLo, pHi, ppb = 0, 30, 0.05
    nB = int((eHi-eLo)/epb)
    nBP = int((pHi-pLo)/ppb)

    bkgVals = {
        # key,      [muE, init guess amp, ampLo, ampHi]
        "flat":     [-1, 1000, 0, 10000],
        "trit":     [-1, 1000, 0, 50000],
        "trit_s":   [-1, 1, 0, 1000],
        "210Pb_c":  [-1, 50, 0, 10000],
        "axion":    [-1, 1, 0, 10000],
        # CDMS paper, Table 3. https://arxiv.org/pdf/1806.07043.pdf
        "49V":      [4.97, 3, 0, 1000],
        "51Cr":     [5.46, 0, 0, 1000], # not in cdms
        "54Mn":     [5.99, 3, 0, 1000],
        "55Fe":     [6.54, 50, 0, 1000],
        "57Co":     [7.11, 3, 0, 1000],
        "65Zn":     [8.98, 10, 0, 1000],
        "65Zn_L":   [1.10, 0, 0, 1000],
        "68Ga":     [9.66, 3, 0, 1000],
        "68Ge":     [10.37, 50, 0, 1000],
        "68Ge_L":   [1.29, 1, 0, 1000],
        "73As":     [11.10, 0, 0, 1000], # not in cdms
        "210Pb_46": [46.54, 50, 0, 1000],
        "210Pb_10": [10.8, 1, 0, 1000],
        "axSi_a":   [1.86, 0, 0, 1000], # k_a1,a2 (lab 1.739, diff 0.12)
        "axSi_b":   [2.00, 0, 0, 1000], # k_b     (lab 1.836, diff 0.16)
        "axS_a":    [2.45, 0, 0, 1000], # k_a1,a2 (lab 2.307, diff 0.14)
        "axS_b":    [2.62, 0, 0, 1000], # k_b     (lab 2.464, diff 0.16)
        "axFe_a":   [6.68, 0, 0, 1000], # k_a   (combination 6.668, 6.701, heliumlike (k_a))
        "axFe_b":   [6.96, 0, 0, 1000], # k_b   (combination 6.952, 6.973, hydrogenlike (k_b))
        "axFe_M1":  [14.4, 1, 0, 1000]
        }

    sigLabels = {
        "axSi_a": r"Si ($\mathregular{K_{\alpha 1,\alpha 2}}$)",
        "axSi_b": r"Si ($\mathregular{K_{\beta}}$)",
        "axS_a": r"S ($\mathregular{K_{\alpha 1,\alpha 2}}$)",
        "axS_b": r"S ($\mathregular{K_{\beta}}$)",
        "axFe_M1": "57Fe (M1)",
        "axion": "Axion Cont."
        }

    # load efficiency correction(s)
    if use90:
        f1 = np.load("%s/data/lat-expo-efficiency-all.npz" % dsi.latSWDir) # << just for debugging
        print("Warning, I'm using the 90\% value")
    else:
        f1 = np.load('%s/data/lat-expo-efficiency-all-e95.npz' % dsi.latSWDir) # << use this one

    xEff = f1['arr_0']
    totEnrEff, totNatEff = f1['arr_1'].item(), f1['arr_2'].item()
    if not use90:
        totEnrEffLo, totEnrEffHi = f1['arr_9'].item(), f1['arr_10'].item()  # this doesn't work with use90
        totNatEffLo, totNatEffHi = f1['arr_11'].item(), f1['arr_12'].item()

    detEff = np.zeros(len(xEff))
    detEffLo = np.zeros(len(xEff))
    detEffHi = np.zeros(len(xEff))
    for ds in dsList:
        if enr:
            detEff += totEnrEff[ds]
            if not use90:
                detEffLo += totEnrEffLo[ds]
                detEffHi += totEnrEffHi[ds]
        else:
            detEff += totNatEff[ds]
            if not use90:
                detEffLo += totNatEffLo[ds]
                detEffHi += totNatEffHi[ds]

    # load exposure
    if use90:
        f2 = np.load("%s/data/expo-totals-e90.npz" % dsi.latSWDir)
        print("Warning, I'm using the 90\% value")
    else:
        f2 = np.load("%s/data/expo-totals-e95.npz"  % dsi.latSWDir) # << use this one

    if use90:
        dsExpo = f2['arr_0'].item()
    else:
        dsExpo, dsUnc = f2['arr_0'].item(), f2['arr_1'].item()
    detExp = 0
    detUncVals = []
    for d in dsExpo:
        if d in dsList:
            if enr:
                detExp += dsExpo[d][0]
                if not use90:
                    detUncVals.append(dsUnc[d][0])
            else:
                detExp += dsExpo[d][1]
                if not use90:
                    detUncVals.append(dsUnc[d][1])

    if not use90:
        detUnc = np.sqrt(np.sum(detUncVals))

    # normalize the efficiency
    detEff = np.divide(detEff, detExp)
    if not use90:
        detEffLo = np.divide(detEffLo, detExp)
        detEffHi = np.divide(detEffHi, detExp)

    effLim, effMax = xEff[-1], detEff[-1] # used in interpolating for weighted data

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

        # === 1. low, bestfit, and high efficiencies
        plt.close()

        plt.plot(xEff[idx], detEff[idx], c='b', label="Measured Efficiency (%.2f%% at %d keV)" % (100*detEff[-1], xEff[-1]))
        plt.plot(xEff[idx], detEffLo[idx], c='g', label="Lower Limit (%.2f%% at %d keV)" % (100*detEffLo[-1], xEff[-1]))
        plt.plot(xEff[idx], detEffHi[idx], c='r', label="Upper Limit (%.2f%% at %d keV)" % (100*detEffHi[-1], xEff[-1]))
        plt.legend(loc=4)
        plt.xlabel("Energy (keV)", ha='right', x=1)
        plt.ylabel("Efficiency", ha='right', y=1)
        plt.tight_layout()
        # plt.show()
        plt.savefig("%s/plots/sf-effBounds.pdf" % dsi.latSWDir)

        # === 2. sim-corrected efficiency plot
        plt.close()
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

    # *** uses the lower or upper limit of the efficiency ***
    if useLoEff:
        print("WARNING: using lower-limit efficiency")
        detEff = detEffLo

    if useHiEff:
        print("WARNING: using upper-limit efficiency")
        detEff = detEffHi

    # fix the bkg peaks list
    bkgPeaksInRange = []
    for pk in bkgModelPeaks:
        opt = "enr" if enr else "nat"
        mu, sig = bkgVals[pk][0], getSigma(bkgVals[pk][0],opt)
        if mu+5*sig <= eLo or mu-5*sig >= eHi: continue
        bkgPeaksInRange.append(pk)
    bkgModelPeaks = bkgPeaksInRange
    if len(bkgModelPeaks)>0:
        print("Floating peaks in %.1f -- %.1f keV range:" % (eLo, eHi), bkgModelPeaks)


def loadDataMJD():
    """ Load MJD data based on the global variable 'dsList'.
    RooFit can't handle the vector<double> format for energies.
    So save a few select branches into a new file.
    """
    from array import array
    from ROOT import TChain, TTree

    # load the data
    tt = TChain("skimTree")
    for ds in dsList:
        if use90:
            tt.Add("%s/final90/final90_DS%s.root" % (dsi.cutDir, ds)) # this is just for debugging
        else:
            tt.Add("%s/final95t/final95t_DS%s.root" % (dsi.cutDir, ds)) # << use this one
    if use90:
        print("Warning, I'm using the 90\% value")

    # declare output
    if use90:
        fName = "%s/data/latDS%s-90.root" % (dsi.latSWDir, ''.join([str(d) for d in dsList]))
    else:
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
    # f2 = TFile(fName)
    # t2 = f2.Get("skimTree")
    # t2.Scan("run:channel:isEnr:trapENFCal:weight")


def getUnscaledPDFs(ma=0, makePlots=False):
    """ Generate a set of TH1D's to be turned into RooDataHist objects.
    Be careful they have the same axis limits and binning as the RooDataSet.
    Takes axion mass (in keV) as a parameter.
    """
    from ROOT import TH1D

    # output file
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


def plotPDFs():
    """ Final consistency check on PDFs before we use them in the fitter. """

    tf = TFile("%s/data/specPDFs.root" % dsi.latSWDir)

    # === 1. Ge photoelectric XS
    plt.close()
    xP, yP, _ = wl.npTH1D(tf.Get("h1"))
    plt.semilogy(xP, yP, ls='steps', c='b', lw=3, label=r"$\sigma_{ge}$")
    plt.xlabel("Energy (keV)", ha='right', x=1)
    plt.ylabel(r"$\mathregular{cm^2 / kg}$", ha='right', y=1)
    plt.legend()
    plt.xlim(0, 15)
    plt.ylim(ymin=1e4)
    plt.tight_layout()
    # plt.show()
    plt.savefig("%s/plots/sf-gexs.pdf" % dsi.latSWDir)

    # === 2. axioelectric XS (mA = 0 and mA = 5)
    plt.close()
    xR, hR, _ = wl.npTH1D(tf.Get("h2"))
    plt.plot(xR, hR, ls='steps', c='b', lw=3, label=r"$\sigma_{ae}, \mathregular{m_a \sim 0}$")

    idx = np.where(xP >= 5.1)
    hR2 = np.asarray([wl.sig_ae(xP[idx][i], 5) * yP[idx][i] for i in range(len(xP[idx]))])
    plt.plot(xP[idx], hR2, ls='steps', c='r', lw=3, label=r"$\sigma_{ae}, \mathregular{m_a = 5\ keV}$")

    plt.xlabel("Energy (keV)", ha='right', x=1)
    plt.ylabel(r"$\mathregular{\sigma_{ae}\ (cm^2\ /\ kg)}$", ha='right', y=1)
    plt.legend()
    plt.xlim(0, 12)
    plt.tight_layout()
    # plt.show()
    plt.savefig("%s/plots/sf-axs.pdf" % dsi.latSWDir)

    # === 3. axion flux, gae=1
    plt.close()
    xR, hR, _ = wl.npTH1D(tf.Get("h3"))
    plt.plot(xR, hR, ls='steps', c='b', lw=3, label=r"$\mathregular{\Phi_a}$, %.2f keV/bin" % ppb)

    plt.axvline(1.85, c='g', lw=2, alpha=0.5, label=r"1.85 keV, Si ($\mathregular{K_{\alpha 1,\alpha 2}}$)")
    plt.axvline(2.00, c='m', lw=2, alpha=0.5, label=r"2.00 keV, Si ($\mathregular{K_{\beta}}$)")
    plt.axvline(2.45, c='r', lw=2, alpha=0.5, label=r"2.45 keV, S ($\mathregular{K_{\alpha 1,\alpha 2}}$)")
    plt.axvline(2.62, c='k', lw=2, alpha=0.5, label=r"2.62 keV, S ($\mathregular{K_{\beta}}$)")
    plt.axvline(6.67, c='orange', lw=2, alpha=0.8, label=r"6.67 keV, S ($\mathregular{K_{\beta}}$)")

    plt.xlabel("Energy (keV)", ha='right', x=1)
    plt.ylabel(r"Flux / (keV cm${}^2$ d)]", ha='right', y=1)
    plt.xlim(0,12)
    plt.legend()
    plt.tight_layout()
    # plt.show()
    plt.savefig("%s/plots/sf-axFlux.pdf" % dsi.latSWDir)

    # === 4. solar axion PDF, gae=1 -- this is what we integrate to get N_exp

    plt.close()
    hB = tf.Get("h4")
    xR1, hR1, _ = wl.npTH1D(hB)
    xR2, hR2, _ = wl.npTH1D(getResConv(hB,"axion"))
    plt.plot(xR1, hR1, ls='steps', c='b', lw=3, label=r"$\Phi_a$, %.2f keV/bin, $\mathregular{g_{ae}=1}$" % ppb)
    plt.plot(xR2, hR2, ls='steps', c='r', lw=3, label=r"Convolved with $\sigma_{MJ}$")
    plt.xlabel("Energy (keV)", ha='right', x=1)
    plt.ylabel("Flux / (keV d kg)", ha='right', y=1)
    plt.xlim(0,10)
    plt.legend()
    plt.tight_layout()
    # plt.show()
    # exit()
    plt.savefig("%s/plots/sf-axPDF.pdf" % dsi.latSWDir)

    # === 5. tritium (show example of efficiency correction)
    plt.close()
    hT = tf.Get("h5")

    xT, yT, xpb = wl.npTH1D(hT)

    xT, yT = normPDF(xT, yT, 0.1, 20)

    xT2, yT2, _ = wl.npTH1D(getResConv(hT,"trit"))
    xT3, yT3 = normPDF(xT2, yT2, 0.1, 20)

    hC = getEffCorrTH1D(hT, pLo, pHi, nBP)
    xC, yC, xpb = wl.npTH1D(hC)

    xC, yC = normPDF(xC, yC, 0.1, 20)

    plt.plot(xT, yT, "-", c='b', label="Raw")
    plt.step(xT3, yT3, "-", c='g', label=r"Convolved with $\sigma_{MJ}$")
    opt = "Enr" if enr else "Nat"
    plt.plot(xC, yC, "-", c='r', label="Eff.Corr, DS%d-%s, %s, Norm'd" % (dsList[0], dsList[-1], opt))
    plt.axvline(1, c='g', lw=1, label="1.0 keV")

    plt.xlabel("Energy (keV)", ha='right', x=1)
    plt.ylabel("Counts (norm)", ha='right', y=1)
    plt.legend()
    plt.xlim(0,20)
    plt.tight_layout()
    # plt.show()
    # exit()
    plt.savefig('%s/plots/sf-trit.pdf' % dsi.latSWDir)

    # === 6. Pb210-TDL PDF
    plt.close()
    # hT = tf.Get("h6")

    # Pb210 (from separate file)
    tf2 = TFile("%s/data/Pb210PDFs.root" % dsi.latSWDir)
    hT = tf2.Get("hPb210TDL") # with TDL
    # h7 = tf2.Get("hPb210") # without TDL

    xT, yT, xpb = wl.npTH1D(hT)
    xT, yT = normPDF(xT, yT, 1, 50)

    plt.step(xT, yT, c='b', lw=2, label=r"$\mathregular{{}^{210}Pb}$, prelim. simulation")
    plt.axvline(1.0, c='g', lw=1, label="1.0 keV")

    plt.xlabel("Energy (keV)", ha='right', x=1)
    plt.ylabel("Counts (norm)", ha='right', y=1)
    plt.legend(loc=2)
    plt.xlim(0.5,50)
    plt.tight_layout()
    # plt.show()
    plt.savefig("%s/plots/sf-pb210.pdf" % dsi.latSWDir)


def getBkgPDF(eff=False):
    from ROOT import TH1D
    bkg = TH1D("bkg","flat BG",nBP,eLo,eHi)
    for iB in range(nBP+1):
        bkg.SetBinContent(iB, 1) # the initial amplitude doesn't matter b/c we normalize to 1 (no width option)
    bkg.SetBinContent(nBP+1, 1)
    bkg.Scale(1 / bkg.Integral(bkg.FindBin(eLo), bkg.FindBin(eHi)))

    if eff:
        # efficiency correct doesn't renormalize to 1
        bkg = getEffCorrTH1D(bkg, eLo, eHi, nBP)

    return bkg


def peakPDF(pkE, sig, name, eff=False):
    """ Make a TH1D of the peaks so we can apply the efficiency correction
    the same way we do for the continuous pdfs.
    """
    from ROOT import TH1D

    xLo = pkE-5*sig
    xHi = pkE+5*sig
    if xHi < eLo or xLo > eHi:
        return None

    # make a simple gaussian
    x = np.arange(eLo, eHi, ppb)
    hP = wl.gauss_function(x, 1, pkE, sig)

    # make the TH1D
    nB2 = int((eHi-eLo)/ppb)
    hPk = TH1D(name,name,nB2,eLo,eHi)
    for iB in range(nB2):
        hPk.SetBinContent(iB, hP[iB])
    hPk.Scale(1 / hPk.Integral(hPk.FindBin(eLo), hPk.FindBin(eHi)))

    # efficiency correct (doesn't normalize)
    if eff:
        hPk = getEffCorrTH1D(hPk, eLo, eHi, nB2)

    return hPk


def normPDF(x, y, xLo, xHi):
    """ Normalize a numpy pdf to 1 in the global
    energy range (not the range the pdf was generated) """
    idx = np.where((x >= xLo)&(x <= xHi + epb/2))
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
    hEff = TH1D(h.GetName()+"_e", h.GetTitle()+"_e", xNB, xLo, xHi)
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

    for iB in range(1, len(xT)):
        for xp, yp, pb, v in pdfs:
            idx = np.where((xp >= xT[iB-1]) & (xp < xT[iB]))
            if len(idx[0]) > 0:
                if amp:
                    yT[iB] += v * np.average(yp[idx]) * (epb/pb)
                else:
                    yT[iB] += np.average(yp[idx]) * (epb/pb)
            if iB == 0:
                yT[iB] = yT[iB+1]

    if smooth:

        xTS = np.linspace(eLo, eHi, 10000)
        yTS = np.zeros(len(xTS))
        for xp, yp, pb, v in pdfs:
            if amp: ypS = spline(xp - pb/2, yp * v * (epb/pb), xTS)
            else:   ypS = spline(xp - pb/2, yp * (epb/pb), xTS)
            yTS += ypS
        # plt.step(xT, yT, 'b')
        # plt.plot(xTS, yTS, 'r')
        # plt.show()
        # exit()
        return xTS, yTS

    return xT, yT


def getSigma(E, opt=""):
    """ Get the MJ energy resolution, in sigma.
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


def getResConv(h, hName):
    """ Convolve continuum PDF's (axion, tritium) with the MJ energy resolution """
    from ROOT import TH1D

    np.random.seed(0)
    def jitter(vals,mu,sig):
        return [np.random.normal(mu, sig) for v in vals]
    opt = "enr" if enr else "nat"

    if hName == "axFlux": hScale = 1e35
    elif hName == "axion": hScale = 1e38
    elif hName == "trit": hScale = 1e-9
    else:
        print("unknown name!")
        exit()

    x1, y1, xpb = wl.npTH1D(h)
    y1 /= hScale

    xj = []
    for i in range(len(x1)):
        tmp = np.full(int(y1[i]), x1[i]) # get n copies of value x
        jit = jitter(tmp, x1[i], getSigma(x1[i], opt))
        xj.extend(jit)
    x2, y2 = wl.GetHisto(xj, pLo, pHi, ppb)

    # correct the scale of y2 (taking integer copies above messes it up)
    int1 = np.sum(y1[np.where((x1>=5)&(x1<=10))]) * ppb
    int2 = np.sum(y2[np.where((x2>=5)&(x2<=10))]) * ppb
    ySca = int1/int2

    # # make sure you're getting the resolution you expect
    # from scipy.optimize import curve_fit
    # def axBkg(x, mu, sig, amp, B, C):
    #     return amp * np.exp(-(x - mu)**2 / (2 * sig**2)) / (sig * np.sqrt(2*np.pi)) + B*x + C
    # idx = np.where((x2>=6)&(x2<=7.5))
    # popt, pcov = curve_fit(axBkg, x2[idx], y2[idx], p0=(6.6, getSigma(6.6, opt), 1000, 10, 500))
    # xF = np.arange(6, 7.5, 0.005)
    # plt.plot(xF, axBkg(xF, *popt), c='g')
    # print("init sigma",getSigma(6.6, opt), "fit sigma",popt[1])
    # # print(popt)
    # plt.step(x1, y1, c='b', lw=2)
    # plt.step(x2, y2, c='r', lw=2)
    # plt.ylim(0, 2500)
    # plt.xlim(6, 7.5)
    # plt.show()

    y1 = np.multiply(y1, hScale)
    y2 = np.multiply(y2, ySca * hScale)
    int1 = np.sum(y1[np.where((x1>=5)&(x1<=10))]) * ppb
    int2 = np.sum(y2[np.where((x2>=5)&(x2<=10))]) * ppb
    # print("int1 %.4e  int2 %.4e" % (int1, int2))

    # now create a TH1D and return it
    hOut = TH1D(hName,hName,nBP,pLo,pHi)
    for iB in range(nBP+1):
        ctr = (iB + 0.5)*ppb + pLo
        bLo, bHi = ctr - ppb/2, ctr + ppb/2
        hOut.SetBinContent(iB, y2[iB])

    # check it
    # x3, y3, xpb3 = wl.npTH1D(hOut,"") # don't shift up
    # plt.step(x2, y2, c='r')
    # plt.step(x3, y3, c='g')
    # plt.show()
    # exit()

    return hOut


def getHistList(hName=None):
    """ Loads TH1D's for the background model. """
    from ROOT import TH1D
    print("Getting hists ...")

    tf2 = TFile("%s/data/specPDFs.root" % dsi.latSWDir)

    if hName is not None:
        if hName == "axion":
            hB = tf2.Get("h4")
            if res: hB = getResConv(hB, "axion")
            if eff: hB = getEffCorrTH1D(hB, pLo, pHi, nBP)
            hB.SetDirectory(0)
            return hB

    hList = []
    for name in [n for n in bkgVals if n in bkgModelHists]:

        pars = bkgVals[name]
        if name=="flat":
            hB = getBkgPDF(eff)
        elif name=="trit":
            hB = tf2.Get("h5")
            # if res: hB = getResConv(hB, name) # trit does not need to be res-corrected, marino phd pg.185
            if eff: hB = getEffCorrTH1D(hB, pLo, pHi, nBP)
        elif name=="axion":
            hB = tf2.Get("h4")
            if res: hB = getResConv(hB, name)
            if eff: hB = getEffCorrTH1D(hB, pLo, pHi, nBP)
        elif name=="210Pb_c":
            hB = tf2.Get("h6") # with preliminary TDL
            if eff: hB = getEffCorrTH1D(hB, hB.GetXaxis().GetXmin(), hB.GetXaxis().GetXmax(), hB.GetNbinsX())
        else:
            opt = "enr" if enr else "nat"
            hB = peakPDF(pars[0], getSigma(pars[0], opt), name, eff)

        if hB is None:
            continue

        hB.SetDirectory(0) # detach "hB" from the "tf2" file

        hList.append([hB, name])

    return hList


def testFunc():
    """ Low-priority to-do: it would be neat to be able to efficiency correct the RooGaussian
    for the peaks.  I was looking at RooGenericPdf but got sidetracked. """

    from ROOT import TCanvas

    tf = TFile("%s/data/latDS%s.root" % (dsi.latSWDir,''.join([str(d) for d in dsList])))
    tt = tf.Get("skimTree")
    tCut = "isEnr==1" if enr is True else "isEnr==0"
    hitE = ROOT.RooRealVar("trapENFCal", "Energy", eLo, eHi, "keV")
    hEnr = ROOT.RooRealVar("isEnr", "isEnr", 0, 1, "")
    fData = ROOT.RooDataSet("data", "data", tt, ROOT.RooArgSet(hitE, hEnr), tCut)

    name = "68Ge"
    opt = "enr" if enr else "nat"
    mu, sig, amp = bkgVals[name][0], getSigma(bkgVals[name][0], opt), bkgVals[name][1]
    pN = ROOT.RooRealVar("amp-"+name, "amp-"+name, amp)
    pM = ROOT.RooRealVar("mu-"+name, "mu-"+name, mu)
    pS = ROOT.RooRealVar("sig-"+name, "sig-"+name, sig)
    pG = ROOT.RooGaussian("gaus-"+name, "gaus-"+name, hitE, pM, pS)
    pE = ROOT.RooExtendPdf("ext-"+name, "ext-"+name, pG, pN)

    fSpec = hitE.frame(RF.Range(eLo,eHi), RF.Bins(nB))
    fData.plotOn(fSpec)

    pE.plotOn(fSpec, RF.Normalization(amp, ROOT.RooAbsReal.Raw))

    c = TCanvas("c","c", 1400, 1000)
    fSpec.SetTitle("")
    fSpec.Draw()
    c.Print("%s/plots/pk-before.pdf" % dsi.latSWDir)
    c.Clear()

    xP = np.arange(eLo, eHi, ppb)
    hP = wl.gauss_function(xP, 20, mu, sig)
    hP /= np.sum(hP)
    print(np.sum(hP)*amp)
    plt.plot(xP, hP * amp / ppb, c='b', lw=2, label="name %d" % amp)
    plt.show()


def fitModel(makePlots=False):
    from ROOT import TH1D, TCanvas, TLegend, gStyle

    # === load data into workspace ===
    if use90:
        tf = TFile("%s/data/latDS%s-90.root" % (dsi.latSWDir,''.join([str(d) for d in dsList]))) # this is for debugging
        print("Warning, I'm using the 90\% value")
    else:
        tf = TFile("%s/data/latDS%s.root" % (dsi.latSWDir,''.join([str(d) for d in dsList]))) # << use this one

    tt = tf.Get("skimTree")
    tCut = "isEnr==1" if enr is True else "isEnr==0"
    hitE = ROOT.RooRealVar("trapENFCal", "Energy", eLo, eHi, "keV")
    hEnr = ROOT.RooRealVar("isEnr", "isEnr", 0, 1, "")
    # hitW = ROOT.RooRealVar("weight", "weight", 1, 1000, "")
    fData = ROOT.RooDataSet("data", "data", tt, ROOT.RooArgSet(hitE, hEnr), tCut)
    # fData = ROOT.RooDataSet("data", "data", tt, ROOT.RooArgSet(hitE, hEnr, hitW), "", "weight")
    fitWS = ROOT.RooWorkspace("fitWS","Fit Workspace")
    getattr(fitWS,'import')(hitE)
    getattr(fitWS,'import')(fData)
    # getattr(fitWS,'import')(hitW)

    tf2 = TFile("%s/data/specPDFs.root" % dsi.latSWDir)
    pdfList = ROOT.RooArgList("shapes")

    # === background model ===
    histModel = []
    for h in getHistList():
        hB, name = h[0], h[1]
        pars = bkgVals[name]
        bkN = ROOT.RooRealVar("amp-"+name, "amp-"+name, pars[1], pars[2], pars[3])
        bkDH = ROOT.RooDataHist("dh-"+name, "dh-"+name, ROOT.RooArgList(hitE), RF.Import(hB))
        bkPDF = ROOT.RooHistPdf("pdf-"+name, "pdf-"+name, ROOT.RooArgSet(hitE), bkDH, 2)
        bkExt = ROOT.RooExtendPdf("ext-"+name, "ext-"+name, bkPDF, bkN)
        hitE.setRange(eLo, eHi)
        histModel.append([bkExt, name, bkN, bkDH, bkPDF, hB])

    peakModel = []
    for name in bkgModelPeaks:
        opt = "enr" if enr else "nat"
        mu, sig, amp = bkgVals[name][0], getSigma(bkgVals[name][0], opt), bkgVals[name][1]

        if constrainPeak:
            pN = ROOT.RooRealVar("amp-"+name, "amp-"+name, amp, bkgVals[name][2], bkgVals[name][3])

            if constrainMu:
                pM = ROOT.RooRealVar("mu-"+name, "mu-"+name, mu) # NOTE: hardcoded mean
            else:
                pM = ROOT.RooRealVar("mu-"+name, "mu-"+name, mu, mu-0.2, mu+0.2)

            # this is what i used for the 14.4 axion study where only the axion peak was floating
            pS = ROOT.RooRealVar("sig-"+name, "sig-"+name, sig, 0.01, 0.5)
            pSC = ROOT.RooGaussian("scon-"+name,"scon-"+name, pS, RF.RooConst(sig), RF.RooConst(0.1)) # "floating width"
            # pSC = ROOT.RooGaussian("scon-"+name,"scon-"+name, pS, RF.RooConst(sig), RF.RooConst(0.001)) # "fixed width"

            # this is what i used for the cosmogenic peak study
            # pS = ROOT.RooRealVar("sig-"+name, "sig-"+name, sig, 0.01, 0.5)
            # pSC = ROOT.RooGaussian("scon-"+name,"scon-"+name, pS, RF.RooConst(sig), RF.RooConst(0.1)) # gaussian constraint pdf mu, sig

            pG1 = ROOT.RooGaussian("g-"+name, "g-"+name, hitE, pM, pS)
            pG = ROOT.RooProdPdf("gaus-"+name, "gaus-"+name, ROOT.RooArgList(pG1, pSC))
            pE = ROOT.RooExtendPdf("ext-"+name, "ext-"+name, pG, pN)
            peakModel.append([pE, name, mu, sig, amp, pN, pM, pS, pG, pSC, pG1])
        else:
            mu, sig, amp = bkgVals[name][0], getSigma(bkgVals[name][0], opt), bkgVals[name][1]
            pN = ROOT.RooRealVar("amp-"+name, "amp-"+name, amp, bkgVals[name][2], bkgVals[name][3])
            pM = ROOT.RooRealVar("mu-"+name, "mu-"+name, mu, mu - 0.2, mu + 0.2)
            # pS = ROOT.RooRealVar("sig-"+name, "sig-"+name, sig, sig - 0.1, sig + 0.1)
            pS = ROOT.RooRealVar("sig-"+name, "sig-"+name, sig, sig-0.001, sig+0.001) # hardcoded width, use only for 14.4 comparison study
            pG = ROOT.RooGaussian("gaus-"+name, "gaus-"+name, hitE, pM, pS)
            pE = ROOT.RooExtendPdf("ext-"+name, "ext-"+name, pG, pN)
            peakModel.append([pE, name, mu, sig, amp, pN, pM, pS, pG])

    bkgModel = histModel + peakModel

    # exit()

    # this is separate b/c all the RooVars have to remain in memory
    for bkg in bkgModel:
        pdfList.add(bkg[0])

    model = ROOT.RooAddPdf("model", "total PDF", pdfList)

    if makePlots:

        # === make a rooplot of the initial guess parameters, don't let roofit normalize automatically

        leg = TLegend(0.83,0.5,0.97,0.9)
        gStyle.SetPalette(ROOT.kRainBow)
        nCol = float(gStyle.GetNumberOfColors())

        fSpec = hitE.frame(RF.Range(eLo,eHi), RF.Bins(nB))

        fData.plotOn(fSpec)

        # wouter's note: DON'T DELETE
        # "the default behavior is when you plot a p.d.f. on an empty frame it is
        # plotted with unit normalization. When you plot it on a frame with data in
        # it, it will normalize to the number of events in that dataset."
        # (then after you do a fit, the pdf normalization changes again ...)
        nData = fData.numEntries()

        nTot = 0
        bkgModel = histModel + peakModel
        for i, ext in enumerate(bkgModel):
            extPDF, name = ext[0], ext[1]
            col = gStyle.GetColorPalette(int(nCol/len(bkgModel) * i))
            extPDF.plotOn(fSpec, RF.LineColor(col), RF.Normalization(bkgVals[name][1], ROOT.RooAbsReal.Raw), RF.Name(name))
            leg.AddEntry(fSpec.findObject(name), name, "l")
            nTot += bkgVals[name][1]

        model.plotOn(fSpec, RF.LineColor(ROOT.kRed), RF.Name("fmodel"), RF.Normalization(nTot, ROOT.RooAbsReal.Raw))
        leg.AddEntry(fSpec.findObject("fmodel"),"Full Model","l")

        c = TCanvas("c","c", 1400, 1000)
        fSpec.SetTitle("")
        fSpec.Draw()
        leg.Draw("same")
        c.Print("%s/plots/sf-before.pdf" % dsi.latSWDir)
        c.Clear()

        # === make a pyplot of the same thing

        tCut = "isEnr" if enr else "!isEnr"
        tCut += " && trapENFCal >= %.1f && trapENFCal <= %.1f" % (eLo, eHi)
        n = tt.Draw("trapENFCal", tCut, "goff")
        trapE = tt.GetV1()
        trapE = [trapE[i] for i in range(n)]
        x, hData = wl.GetHisto(trapE, eLo, eHi, epb)
        hErr = np.asarray([np.sqrt(h) for h in hData])

        # comment this out to show just the bkg model
        # plt.errorbar(x, hData, yerr=hErr, c='k', ms=5, linewidth=0.5, fmt='.', capsize=1, zorder=1) # pretty convincing rooplot fake

        cmap = plt.cm.get_cmap('jet',len(bkgModel))
        pdfs = []
        for i, bkg in enumerate(bkgModel):
            name = bkg[1]

            if name in bkgModelHists:
                x, y, xpb = wl.npTH1D(bkg[5])
                x, y = normPDF(x, y, eLo, eHi)
                nCts = bkgVals[bkg[1]][1]

            elif name in bkgModelPeaks:
                mu, sig, nCts = bkg[2], bkg[3], bkg[4]
                x = np.arange(eLo, eHi, ppb)
                xpb = ppb
                y = wl.gauss_function(x, 20, mu, sig)
                y /= np.sum(y)

            xS = np.arange(eLo, eHi, 0.01)
            yS = spline(x - xpb/2, y, xS)
            plt.plot(xS, yS * nCts / xpb, c=cmap(i), lw=2, label="%s init cts: %d" % (name, nCts))

            pdfs.append([x, y, xpb, nCts])

        xT, yT = getTotalModel(pdfs, eLo, eHi, epb, smooth=True)
        plt.step(xT, yT / epb, c='r', lw=2, label="Total Model")

        plt.xlabel("Energy (keV)", ha='right', x=1)
        plt.ylabel("Counts / %.1f keV" % epb, ha='right', y=1)
        plt.legend(loc=1, fontsize=11)
        # plt.minorticks_on()
        plt.xlim(eLo, eHi)
        plt.ylim(ymin=0)
        plt.tight_layout()
        # plt.show()
        plt.savefig("%s/plots/sf-before-mpl.pdf" % dsi.latSWDir)


    # === alright, now run the fit and output to the workspace
    minimizer = ROOT.RooMinimizer( model.createNLL(fData, RF.NumCPU(2,0), RF.Extended(True)) )
    minimizer.setPrintLevel(-1)
    minimizer.setStrategy(2)
    minimizer.migrad()
    fitRes = minimizer.save()

    # according to the internet, covQual==3 is a good indicator that it converged
    print("Fitter is done. Fit Cov Qual:", fitRes.covQual())

    # save workspace to a TFile
    getattr(fitWS,'import')(fitRes)
    getattr(fitWS,'import')(model)
    tf3 = TFile("%s/data/fitWS.root" % dsi.latSWDir,"RECREATE")
    fitWS.Write()
    tf3.Close()

    fitVals = {}
    nPars = fitRes.floatParsFinal().getSize()
    for i in range(nPars):
        fp = fitRes.floatParsFinal()
        name = fp.at(i).GetName()
        fitVal, fitErr = fp.at(i).getValV(), fp.at(i).getError()
        fitVals[name] = [fitVal, fitErr]
    np.savez("%s/data/sf6-results.npz" % dsi.latSWDir, fitVals)


def plotFit(plotRate=False, plotProfileResults=False):

    from ROOT import TCanvas, TH1D, TLegend, gStyle

    f = TFile("%s/data/fitWS.root" % dsi.latSWDir)
    fitWS = f.Get("fitWS")
    fData = fitWS.allData().front()
    fitRes = fitWS.allGenericObjects().front()
    nPars = fitRes.floatParsFinal().getSize()
    hitE = fitWS.var("trapENFCal")
    model = fitWS.pdf("model")
    # fitWS.Print()

    # === get fit results: {name : [nCts, err]} ===
    fitVals = {}
    for i in range(nPars):
        fp = fitRes.floatParsFinal()
        name = fp.at(i).GetName()
        fitVal, fitErr = fp.at(i).getValV(), fp.at(i).getError()
        fitVals[name] = [fitVal, fitErr]
        print(name, wl.niceList(fitVals[name]))

    # compare the diffs w/ the lit values if you're floating any peaks
    for name in bkgModelPeaks:
        litE = bkgVals[name][0]

        # mu, sig, amp = fitVals["mu-"+name], fitVals["sig-"+name], fitVals["amp-"+name] # floating mean
        # mu = mu[0]
        # print("using floating mean")

        sig, amp = fitVals["sig-"+name], fitVals["amp-"+name] # fixed mean
        mu = litE
        print("using fixed mean")

        dMu = litE - mu
        opt = "enr" if enr else "nat"
        dSig = getSigma(litE,opt) - sig[0]
        print(name)
        print("E   --  lit %.4f  fit %.4f  diff %.4f  (%.2f%%)" % (litE, mu, dMu, 100*dMu/litE))
        print("Sig -- func %.4f  fit %.4f  diff %.4f  (%.2f%%)" % (getSigma(litE, opt), sig[0], dSig, 100*dSig/sig[0]))

    if plotProfileResults:
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

    # === make a rooplot of the fit ===

    bkgModel = bkgModelHists + bkgModelPeaks

    leg = TLegend(0.83,0.5,0.97,0.9)
    gStyle.SetPalette(ROOT.kRainBow)
    nCol = float(gStyle.GetNumberOfColors())

    fSpec = hitE.frame(RF.Range(eLo,eHi), RF.Bins(nB))

    fData.plotOn(fSpec)

    for i, name in enumerate(bkgModel):
        pdfName = "ext-"+name
        col = gStyle.GetColorPalette(int(nCol/len(bkgModel) * i))
        model.plotOn(fSpec, RF.Components(pdfName), RF.LineColor(col), RF.LineStyle(ROOT.kDashed), RF.Name(name))
        leg.AddEntry(fSpec.findObject(name), name, "l")

    chiSquare = fSpec.chiSquare(nPars)
    model.plotOn(fSpec, RF.LineColor(ROOT.kRed), RF.Name("fmodel"))
    leg.AddEntry(fSpec.findObject("fmodel"),"Full Model, #chi^{2}/NDF = %.3f" % chiSquare, "l")

    c = TCanvas("c","c", 1400, 1000)
    fSpec.SetTitle("")
    fSpec.Draw()
    leg.Draw("same")
    c.Print("%s/plots/sf-after.pdf" % dsi.latSWDir)
    c.Clear()

    # === duplicate the rooplot in matplotlib ===
    plt.close()

    if use90:
        tf = TFile("%s/data/latDS%s-90.root" % (dsi.latSWDir,''.join([str(d) for d in dsList]))) # this is for debugging
        print("Warning, I'm using the 90\% value")
    else:
        tf = TFile("%s/data/latDS%s.root" % (dsi.latSWDir,''.join([str(d) for d in dsList]))) # << use this one

    tt = tf.Get("skimTree")
    tCut = "isEnr==1" if enr is True else "isEnr==0"
    tCut = "isEnr" if enr else "!isEnr"
    tCut += " && trapENFCal >= %.1f && trapENFCal <= %.1f" % (eLo, eHi)
    n = tt.Draw("trapENFCal", tCut, "goff")
    trapE = tt.GetV1()
    trapE = [trapE[i] for i in range(n)]
    x, hData = wl.GetHisto(trapE, eLo, eHi, epb)

    if plotRate:
        hErr = np.asarray([np.sqrt(h)/detExp for h in hData]) # statistical error / exposure
        plt.errorbar(x, hData/detExp, yerr=hErr, c='k', ms=5, linewidth=0.5, fmt='.', capsize=1, zorder=1)
    else:
        hErr = np.asarray([np.sqrt(h) for h in hData]) # statistical error
        dCts = np.sum( hData[np.where((x>=eLo)&(x<=eHi))] )
        lab = "DS%s-%s, %.2f kg-y, %d cts" % (str(dsList[0]), str(dsList[-1]), detExp/365.25, dCts)
        # lab = "DS1-4,5B,5C, %.2f kg-y, %d cts" % (detExp/365.25, dCts)
        plt.errorbar(x, hData, yerr=hErr, c='k', ms=5, linewidth=0.5, fmt='.', capsize=1, zorder=1, label=lab)

    # plot the bkg components
    cmap = plt.cm.get_cmap('jet',len(bkgModel))
    pdfs, pdfsCorr = [], []
    nTot, nTotC = 0, 0
    hList = getHistList()

    for i, name in enumerate(bkgModel):

        if name in bkgModelHists:
            for h in hList:
                if h[1]==name: hB = h[0]
            x, y, xpb = wl.npTH1D(hB)
            x, y = normPDF(x, y, eLo, eHi)
            nCts, nErr = fitVals["amp-"+name][0], fitVals["amp-"+name][1]

            yc = nCts * getEffCorr(x, y, inv=True)
            ec = np.sum(yc)/nCts
            print(name,"eff corr: %.2f ± %.2f (factor %.2f)" %(np.sum(yc),nErr*ec,ec))

            # print the efficiency-corrected counts
            lab = "%s: %.1f ± %.1f cts" % (name, np.sum(yc), nErr*ec)
            if name in profileVars:
                lab = "%s %.1f cts (90%%CL)" % (sigLabels[name], np.sum(yc))

        elif name in bkgModelPeaks:
            sig, nCts = fitVals["sig-"+name][0], fitVals["amp-"+name][0]
            mu = bkgVals[name][0]
            x = np.arange(eLo, eHi, ppb)
            xpb = ppb
            y = wl.gauss_function(x, 20, mu, sig)
            y /= np.sum(y)
            lab = r"%s (%.2f keV): %.2f keV, $\mathregular{\sigma=}$%.2f" % (name, bkgVals[name][0], mu, sig)
            yc = nCts * getEffCorr(x, y, inv=True)
            ec = np.sum(yc)/nCts
            print(name,"eff corr: %.2f ± %.2f (factor %.2f)" %(np.sum(yc),nErr*ec,ec))

        # plot a smoothed version, like roofit does
        xS = np.arange(eLo, eHi, 0.001)
        yS = spline(x - xpb/2, y, xS)

        if plotRate:
            plt.plot(xS, yS * nCts * (epb/xpb) / detExp, c=cmap(i), lw=2, label="%s %.3f ± %.3f" % (name, nCts/detExp, nErr/detExp))
        else:
            plt.plot(xS, yS * nCts * (epb/xpb), c=cmap(i), lw=2, label=lab)

        pdfs.append([x, y, xpb, nCts])
        pdfsCorr.append([x, yc, xpb, nCts])
        nTot += nCts
        nTotC += np.sum(yc) # reverse the efficiency correction to get the "true" number of counts

    # get the fit model, and the efficiency-corrected final model
    xT, yT = getTotalModel(pdfs, eLo, eHi, epb, smooth=True)
    xTc, yTc = getTotalModel(pdfsCorr, eLo, eHi, epb, smooth=True, amp=False)

    if plotRate:
        plt.plot(xT, yT / detExp, 'r', lw=2, alpha=0.7, label="Rate: %.2f c/kev/kg-d" % (nTot/((eHi-eLo)*detExp)) )
        plt.plot(xTc, yTc / detExp, c='m', lw=3, alpha=0.7, label="Eff.Corr: %.2f " % (nTotC/((eHi-eLo)*detExp))  )
        plt.ylabel("Counts / keV / kg-d", ha='right', y=1)

    else:
        plt.plot(xT, yT, 'r', lw=2, alpha=0.7, label="Full Model PDF")
        plt.plot(xTc, yTc, c='m', lw=3, alpha=0.7, label="Model w/ Eff: %d cts" % nTotC)
        plt.ylabel("Counts / %.1f keV" % epb, ha='right', y=1)



    # get an efficiency-corrected number of counts for this energy region
    # to help figure out the background index
    n = tt.Draw("trapENFCal:weight", tCut, "goff")
    trapE, weight = tt.GetV1(), tt.GetV2()
    trapE = [trapE[i] for i in range(n)]
    weight = [weight[i] for i in range(n)]
    x, hDataW = wl.GetHisto(trapE, eLo, eHi, epb, wts=weight)

    nReg = np.sum(hData)
    nWtd = np.sum(hDataW)

    print("Data counts, %.2f-%.2f: %d, eff-corrected: %.2f, expo: %.3f kg-y" % (eLo, eHi, nReg, nWtd, detExp/365.25))

    print("Efficiency corrected rate, %.1f-%.1f keV:  %.3f c/(kev-kg-d)" % (eLo, eHi, nWtd/detExp/(eHi-eLo)))

    idx = np.where((x>=1.5) & (x<=8))
    nWtd2 = np.sum(hDataW[idx])
    print("Efficiency corrected rate, %.1f-%.1f keV:  %.3f c/(kev-kg-d)" % (1.5, 8, nWtd2/detExp/(8-1.5)))

    # plt.plot(np.nan, np.nan, ".w", label="Data w/ Eff: %d cts" % nWtd)

    plt.xlabel("Energy (keV)", ha='right', x=1)
    plt.legend(loc=1, fontsize=11)
    plt.minorticks_on()
    plt.xlim(eLo, eHi)
    plt.ylim(ymin=0)
    plt.tight_layout()
    # plt.show()
    plt.savefig("%s/plots/sf-after-mpl.pdf" % dsi.latSWDir)


def getProfile(idx=None, update=False):
    from ROOT import TCanvas
    from ROOT import RooStats as RS

    f = TFile("%s/data/fitWS.root" % dsi.latSWDir)
    fitWS = f.Get("fitWS")
    fData = fitWS.allData().front()
    hitE = fitWS.var("trapENFCal")
    model = fitWS.pdf("model")
    fitRes = fitWS.allGenericObjects().front()
    fPars = fitRes.floatParsFinal()
    nPars = fPars.getSize()

    # === get fit results: {name : [nCts, err]} ===
    fitVals = {}
    for i in range(nPars):
        fp = fitRes.floatParsFinal()
        name = fp.at(i).GetName()
        fitVal, fitErr = fp.at(i).getValV(), fp.at(i).getError()
        if "amp" in name:
            fitVals[name] = [fitVal, fitErr, name.split('-')[1]]

    # for f in fitVals:
        # print(f, fitVals[f])

    # === get "true" counts (reverse efficiency correction based on fit value) ===
    hAx = getHistList("axion")
    x, y, xpb = wl.npTH1D(hAx)
    x, y = normPDF(x, y, eLo, eHi)
    nCts, nErr, _ = fitVals["amp-axion"] # fit result
    yc = nCts * getEffCorr(x, y, inv=True)
    nCorr = np.sum(yc)
    effCorr = nCorr/nCts
    print("nCts %d  nCorr %.2f  effCorr %.2f" % (nCts, nCorr, effCorr))

    # thesis plot, don't delete
    plt.close()
    plt.step(x, y * nCts * (epb/xpb), c='r', lw=2, label="PDF w/ Efficiency: %.1f" % (nCts))
    plt.step(x, yc * (epb/xpb), c='b', lw=2, label="True counts, DS%s-%s: %.1f" % (dsList[0], dsList[-1], nCorr))
    plt.axvline(eLo, c='g', lw=1, label="%.2f keV" % eLo)
    plt.xlabel("Energy (keV)", ha='right', x=1)
    plt.ylabel("Counts / keV", ha='right', y=1)
    plt.legend(loc=1)
    plt.tight_layout()
    # plt.show()
    plt.xlim(xmax=10)
    plt.savefig("%s/plots/lat-axionPDFCorr.pdf" % dsi.latSWDir)

    idx = 1 if simEffCorr else 0
    tMode = "UPDATE" if update else "RECREATE"

    if use90:
        tOut = TFile("%s/data/rs-plc-%d-fs90.root" % (dsi.latSWDir, idx), tMode)
        print("Warning, saving the profile curve for eff=90...")
    elif useLoEff:
        tOut = TFile("%s/data/rs-plc-%d-effLo.root" % (dsi.latSWDir, idx), tMode)
        print("Warning, saving the profile curve for eff=lower limit...")
    elif useHiEff:
        tOut = TFile("%s/data/rs-plc-%d-effHi.root" % (dsi.latSWDir, idx), tMode)
        print("Warning, saving the profile curve for eff=upper limit...")
    else:
        tOut = TFile("%s/data/rs-plc-%d.root" % (dsi.latSWDir, idx), tMode) # << use this one by default

    name = "amp-axion"
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

    from ROOT import TCanvas
    c = TCanvas("c","c",800,600)
    plot.Draw("tf1")
    c.Print("%s/plots/profile-axion-test.pdf" % dsi.latSWDir)

    pName = "hP_%d" % idx if idx is not None else "hP"
    hProfile = plot.GetPlottedObject()
    hProfile.SetName(pName)
    hProfile.SetTitle("PL %.2f  %s  lo %.3f  mid %.3f  hi %.3f  eC %.3f" % (pCL, name, lower, fitVal, upper, effCorr))
    hProfile.Write()
    print(hProfile.GetTitle())

    tOut.Close()


def plotProfile(makePlots=False):

    getSimEff = True
    get90 = True

    tf2 = TFile("%s/data/specPDFs.root" % dsi.latSWDir)
    hA = tf2.Get("h4")
    if res: hA = getResConv(hA, "axion")
    Nexp = hA.Integral(hA.FindBin(eLo), hA.FindBin(eHi), "width") * detExp

    tf0 = TFile("%s/data/rs-plc-0.root" % dsi.latSWDir) # no correction
    hP0 = tf0.Get("hP_0")
    hT0 = hP0.GetTitle().split()
    pars = []
    for v in hT0:
        try: pars.append(float(v))
        except: pass
    pcl0, intLo0, bestFit0, intHi0, effCorr0 = pars
    print(pars)

    # apply the efficiency correction
    Nobs = intHi0 * effCorr0
    gae0 = np.power(Nobs/Nexp, 1/4)
    print("No-corr: Expo (kg-d) %.2f  Nobs: %.2f  Nexp %.2e  eLo %.1f  eHi %.1f  g_ae U.L. %.4e" % (detExp, Nobs, Nexp, eLo, eHi, gae0))

    if getSimEff:
        tf1 = TFile("%s/data/rs-plc-1.root" % dsi.latSWDir) # sim eff correction
        hP1 = tf1.Get("hP_1")
        hT1 = hP1.GetTitle().split()
        pars = []
        for v in hT1:
            try: pars.append(float(v))
            except: pass
        pcl1, intLo1, bestFit1, intHi1, effCorr1 = pars
        Nobs = intHi1 * effCorr1
        gae1 = np.power(Nobs/Nexp, 1/4)
        print("With corr: Expo (kg-d) %.2f  Nobs: %.2f  Nexp %.2e  eLo %.1f  eHi %.1f  g_ae U.L. %.4e" % (detExp, Nobs, Nexp, eLo, eHi, gae1))

    if get90:
        tf2 = TFile("%s/data/rs-plc-0-fs90.root" % dsi.latSWDir) # 90% fitslo cut
        hP2 = tf2.Get("hP_0")
        hT2 = hP2.GetTitle().split()
        pars = []
        for v in hT2:
            try: pars.append(float(v))
            except: pass
        pcl2, intLo2, bestFit2, intHi2, effCorr2 = pars
        Nobs = intHi2 * effCorr2
        gae2 = np.power(Nobs/Nexp, 1/4)
        print("With corr: Expo (kg-d) %.2f  Nobs: %.2f  Nexp %.2e  eLo %.1f  eHi %.1f  g_ae U.L. %.4e" % (detExp, Nobs, Nexp, eLo, eHi, gae2))

    # sanity check - graham's number
    # malbExp  = 89.5 # kg-d
    # Nobs = 106  # fig 5.14
    # Nexp = hA.Integral(hA.FindBin(1.5), hA.FindBin(8), "width") * malbExp # [cts / (keV d kg)] * [kg d]
    # gae = np.power(Nobs/Nexp, 1/4)
    # print("Malbek: Expo (kg-d) %.2f  Nobs: %.2f  Nexp %.2e" % (malbExp, Nobs, Nexp))
    # print("g_ae U.L.:  %.4e" % gae)

    # sanity check 2 - get the maximum of the confint, same as in RooStats:
    # Double_t Yat_Xmax = 0.5*ROOT::Math::chisquared_quantile(fInterval->ConfidenceLevel(),1);
    # from scipy.stats import chi2
    # df = 1 # num pars of interest, =1 if this is the only one
    # cx = np.linspace(chi2.ppf(0.85, df), chi2.ppf(0.92, df), 100)
    # cy = chi2.cdf(cx, df)
    # chi2max = cx[np.where(cy>=0.9)][0] * 0.5
    chi2max = 1.355
    print(chi2max)

    if makePlots:
        plt.close()
        plt.axhline(chi2max, c='m', lw=2, label=r"$\chi^2\mathregular{/2\ (90\%\ C.L.)}}$")

        xP0, yP0, xpb = wl.npTH1D(hP0)
        xP0, yP0 = xP0[1:] * effCorr0 - xpb/2, yP0[1:]

        # plt.plot(xP0, yP0, '-b', lw=4, label=r"w/ Measured Eff. $\mathregular{g_{ae} \leq}$ %.2e" % gae0)
        plt.plot(xP0, yP0, '-b', lw=4, label=r"Measured Eff, fitSlo95. $\mathregular{g_{ae} \leq}$ %.2e" % gae0)
        cHi0 = intHi0 * effCorr0
        cLo0 = intLo0 * effCorr0
        plt.plot([cHi0,cHi0],[0,chi2max], '-b', lw=2, alpha=0.5)
        plt.plot([cLo0,cLo0],[0,chi2max], '-b', lw=2, alpha=0.5)

        if get90:
            xP2, yP2, xpb = wl.npTH1D(hP2)
            xP2, yP2 = xP2[1:] * effCorr2 - xpb/2, yP2[1:]
            plt.plot(xP2, yP2, '-g', lw=4, label=r"Measured Eff, fitSlo90. $\mathregular{g_{ae} \leq}$ %.2e" % gae2)
            cHi2 = int(intHi2 * effCorr2)
            cLo2 = int(intLo2 * effCorr2)
            plt.plot([cHi2,cHi2],[0,chi2max], '-g', lw=2, alpha=0.5)
            plt.plot([cLo2,cLo2],[0,chi2max], '-g', lw=2, alpha=0.5)

        if getSimEff:
            xP1, yP1, xpb = wl.npTH1D(hP1)
            xP1, yP1 = xP1[1:] * effCorr1 - xpb/2, yP1[1:]
            plt.plot(xP1, yP1, '-r', lw=4, label=r"Sim-Corr. Eff, fitSlo95. $\mathregular{g_{ae} \leq}$ %.2e" % gae1)
            cHi1 = int(intHi1 * effCorr1)
            cLo1 = int(intLo1 * effCorr1)
            plt.plot([cHi1,cHi1],[0,chi2max], '-r', lw=2, alpha=0.5)
            plt.plot([cLo1,cLo1],[0,chi2max], '-r', lw=2, alpha=0.5)

        plt.xlabel(r"$\mathregular{N_{obs}}$", ha='right', x=1)
        plt.ylabel(r"-log $\mathregular{\lambda(\mu_{axion})}$", ha='right', y=1)
        plt.legend(loc=2, fontsize=14)
        plt.ylim(0, chi2max*3)
        plt.xlim(500, 1100) # DS1--5C limit w/ sim correction too
        plt.tight_layout()
        # plt.show()
        plt.savefig("%s/plots/sf-axion-profile-corr.pdf" % dsi.latSWDir)


def combineProfileLoHi():
    # plot the profiles of the lower, best fit, and upper efficiency
    # to get an idea of the systematic.  I bet it's awful

    chi2max = 1.355

    tf2 = TFile("%s/data/specPDFs.root" % dsi.latSWDir)
    hA = tf2.Get("h4")
    if res: hA = getResConv(hA, "axion")
    Nexp = hA.Integral(hA.FindBin(eLo), hA.FindBin(eHi), "width") * detExp

    # best
    tf0 = TFile("%s/data/rs-plc-0.root" % dsi.latSWDir)
    hP0 = tf0.Get("hP_0")
    hT0 = hP0.GetTitle().split()
    pars = []
    for v in hT0:
        try: pars.append(float(v))
        except: pass
    pcl0, intLo0, bestFit0, intHi0, effCorr0 = pars
    xP0, yP0, xpb0 = wl.npTH1D(hP0)
    xP0, yP0 = effCorr0 * (xP0[1:] - xpb0/2), yP0[1:]
    intHi0 *= effCorr0
    pyCts0 = xP0[np.where(yP0>=chi2max)][0]-xpb0/2
    Nobs = intHi0
    gae0 = np.power(Nobs/Nexp, 1/4)
    print("Best: Expo (kg-d) %.2f  Nobs: %.2f  Nexp %.2e  eLo %.1f  eHi %.1f  g_ae U.L. %.4e" % (detExp, Nobs, Nexp, eLo, eHi, gae0))

    # hi
    tf1 = TFile("%s/data/rs-plc-0-effHi.root" % dsi.latSWDir)
    hP1 = tf1.Get("hP_0")
    hT1 = hP1.GetTitle().split()
    pars = []
    for v in hT1:
        try: pars.append(float(v))
        except: pass
    pCL1, intLo1, bestFit1, intHi1, effCorr1 = pars
    xP1, yP1, xpb1 = wl.npTH1D(hP1)
    xP1, yP1 = effCorr1 * (xP1[1:] - xpb1/2), yP1[1:]
    intHi1 *= effCorr1
    pyCts1 = xP1[np.where(yP1>=chi2max)][0]-xpb1/2
    Nobs = intHi1
    gae1 = np.power(Nobs/Nexp, 1/4)
    print("Hi-lim: Expo (kg-d) %.2f  Nobs: %.2f  Nexp %.2e  eLo %.1f  eHi %.1f  g_ae U.L. %.4e" % (detExp, Nobs, Nexp, eLo, eHi, gae1))

    # lo
    tf2 = TFile("%s/data/rs-plc-0-effLo.root" % dsi.latSWDir)
    hP2 = tf2.Get("hP_0")
    hT2 = hP2.GetTitle().split()
    pars = []
    for v in hT2:
        try: pars.append(float(v))
        except: pass
    pCL2, intLo2, bestFit2, intHi2, effCorr2 = pars
    xP2, yP2, xpb2 = wl.npTH1D(hP2)
    xP2, yP2 = effCorr2 * (xP2[1:] - xpb2/2), yP2[1:]
    intHi2 *= effCorr2
    pyCts2 = xP2[np.where(yP2>=chi2max)][0]-xpb2/2
    Nobs = intHi2
    gae2 = np.power(Nobs/Nexp, 1/4)
    print("Lo-lim: Expo (kg-d) %.2f  Nobs: %.2f  Nexp %.2e  eLo %.1f  eHi %.1f  g_ae U.L. %.4e" % (detExp, Nobs, Nexp, eLo, eHi, gae2))

    print("Best gae: %.4e (-%.4e) (+%.4e)" % (gae0, gae0-gae1, gae2-gae0))

    plt.close()

    plt.axhline(chi2max, c='m', lw=2, label=r"$\chi^2\mathregular{/2\ (90\%\ C.L.)}}$")

    plt.plot(xP0, yP0, "-b", lw=2, label=r"Best-Fit Eff, $\mathregular{g_{ae} \leq}$ %.2e" % gae0)
    plt.plot([intHi0,intHi0],[0,chi2max], '-b', lw=2, alpha=0.5)

    plt.plot(xP1, yP1, "-g", lw=2, label=r"Upper Eff, $\mathregular{g_{ae} \leq}$ %.2e" % gae1)
    plt.plot([intHi1,intHi1],[0,chi2max], '-g', lw=2, alpha=0.5)

    plt.plot(xP2, yP2, "-r", lw=2, label=r"Lower Eff, $\mathregular{g_{ae} \leq}$ %.2e" % gae2)
    plt.plot([intHi2,intHi2],[0,chi2max], '-r', lw=2, alpha=0.5)

    plt.xlabel(r"$\mathregular{N_{obs}}$", ha='right', x=1)
    plt.ylabel(r"-log $\mathregular{\lambda(\mu_{axion})}$", ha='right', y=1)
    plt.legend(loc=1, fontsize=13)
    plt.ylim(0, chi2max*3)
    plt.xlim(700, 1500)
    plt.tight_layout()
    # plt.show()
    plt.savefig("%s/plots/sf-spec-effSystematic.pdf" % dsi.latSWDir)


def getProfileM1():
    # I do some things specific to the axion spectrum in getProfile
    # so make a separate function for the M1 peaks here.

    from ROOT import TCanvas
    from ROOT import RooStats as RS

    f = TFile("%s/data/fitWS.root" % dsi.latSWDir)
    fitWS = f.Get("fitWS")
    fData = fitWS.allData().front()
    hitE = fitWS.var("trapENFCal")
    model = fitWS.pdf("model")
    fitRes = fitWS.allGenericObjects().front()
    fPars = fitRes.floatParsFinal()
    nPars = fPars.getSize()

    # === get fit results: {name : [nCts, err]} ===
    fitVals = {}
    for i in range(nPars):
        fp = fitRes.floatParsFinal()
        name = fp.at(i).GetName()
        fitVal, fitErr = fp.at(i).getValV(), fp.at(i).getError()
        fitVals[name] = [fitVal, fitErr]

    # for f in fitVals:
        # print(f, wl.niceList(fitVals[f]))


    # === get the efficiency correction factor

    name = "axFe_M1"

    # this is if you let the peak be a RooGaussian
    mu = bkgVals[name][0]
    nCts = fitVals["amp-"+name][0]
    opt = "enr" if enr else "nat"

    sig = fitVals["sig-"+name][0] # floating width
    # sig = getSigma(mu, opt) # fixed width

    x = np.arange(eLo, eHi, ppb)
    xpb = ppb
    y = wl.gauss_function(x, 20, mu, sig)
    y /= np.sum(y)
    yc = nCts * getEffCorr(x, y, inv=True)
    nCorr = np.sum(yc)
    effCorr = nCorr/nCts

    # # this is if you force it to be a histo w/ set energy and resolution
    # pkE = bkgVals[name][0]
    # opt = "enr" if enr else "nat"
    # hAx = peakPDF(pkE, getSigma(pkE, opt), "axFe_M1", eff)
    # x, y, xpb = wl.npTH1D(hAx)
    # x, y = normPDF(x, y, eLo, eHi)
    # nCts, nErr = fitVals["amp-"+name]
    # yc = nCts * getEffCorr(x, y, inv=True)
    # nCorr = np.sum(yc)
    # effCorr = nCorr/nCts
    print("nCts %d  nCorr %.2f  effCorr %.2f" % (nCts, nCorr, effCorr))

    # plt.close()
    # plt.plot(x, y, c='b')
    # plt.plot(x, yc, c='r')
    # plt.show()

    idx = 1 if simEffCorr else 0

    # tOut = TFile("%s/data/rs-plc-%d-%s.root" % (dsi.latSWDir, idx, name), "RECREATE")
    # print("Saving fixed width file")

    # tOut = TFile("%s/data/rs-plc-%d-%s-float.root" % (dsi.latSWDir, idx, name), "RECREATE")
    # print("Saving floating width file")

    if useLoEff:
        tOut = TFile("%s/data/rs-plc-%d-effLoM1.root" % (dsi.latSWDir, idx), "RECREATE")
        print("Warning, saving the profile curve for eff=lower limit...")
    elif useHiEff:
        tOut = TFile("%s/data/rs-plc-%d-effHiM1.root" % (dsi.latSWDir, idx), "RECREATE")
        print("Warning, saving the profile curve for eff=upper limit...")
    else:
        tOut = TFile("%s/data/rs-plc-%d-floatM1.root" % (dsi.latSWDir, idx), "RECREATE")
        print("Warning, saving to the floaty file instead of default")


    start = time.clock()

    fitVal = fitVals["amp-"+name][0]
    thisVar = fitWS.var("amp-"+name)
    pCL = 0.9
    plc = RS.ProfileLikelihoodCalculator(fData, model, ROOT.RooArgSet(thisVar))
    plc.SetConfidenceLevel(pCL)
    interval = plc.GetInterval()
    lower = interval.LowerLimit(thisVar)
    upper = interval.UpperLimit(thisVar)
    plot = RS.LikelihoodIntervalPlot(interval)
    plot.SetNPoints(50)
    plot.Draw("tf1")

    c = TCanvas("c","c",800,600)
    plot.Draw("tf1")
    c.Print("%s/plots/profile-axion-test.pdf" % dsi.latSWDir)

    pName = "hP_%d" % idx if idx is not None else "hP"
    hProfile = plot.GetPlottedObject()
    hProfile.SetName(pName)
    hProfile.SetTitle("PL %.2f  %s  lo %.3f  mid %.3f  hi %.3f  eC %.3f" % (pCL, name, lower, fitVal, upper, effCorr))
    hProfile.Write()
    print(hProfile.GetTitle())

    tOut.Close()

    # print("elapsed:",time.clock()-start)


def plotProfileM1():

    name = "axFe_M1"
    idx = 1 if simEffCorr else 0

    # tf = TFile("%s/data/rs-plc-%d-%s.root" % (dsi.latSWDir, idx, name)) # fixed width file
    # print("loading fixed width file")

    tf = TFile("%s/data/rs-plc-%d-%s-float.root" % (dsi.latSWDir, idx, name)) # floating width file
    print("loading floating width file")

    hP = tf.Get("hP_%d" % idx)
    hT = hP.GetTitle().split()
    pars = []
    for v in hT:
        try: pars.append(float(v))
        except: pass
    pCL, intLo, bestFit, intHi, effCorr = pars
    # print(pars) # reminder: pars are not efficiency corrected

    xP, yP, xpb = wl.npTH1D(hP)
    xP, yP = effCorr * (xP[1:] - xpb/2), yP[1:]
    intHi *= effCorr
    chi2max = 1.355
    pyCts = xP[np.where(yP>=chi2max)][0]-xpb/2
    print("True 90CL (effcorr), %.3f" % pyCts)

    plt.close()

    plt.axhline(chi2max, c='m', lw=2, label=r"$\chi^2\mathregular{/2\ (90\%\ C.L.)}}$")

    plt.plot(xP, yP, "-b", lw=2, label="Fe-57 (M1): %.2f cts (90%%CL, eff.corr.)" % (pyCts))
    plt.plot([pyCts,pyCts],[0,chi2max], '-b', lw=2, alpha=0.5)
    plt.plot([intHi,intHi],[0,chi2max], '-r', lw=2, alpha=0.5)

    plt.xlabel(r"$\mathregular{N_{obs}}$", ha='right', x=1)
    plt.ylabel(r"-log $\mathregular{\lambda(\mu_{axion})}$", ha='right', y=1)
    plt.legend(loc=1, fontsize=14)
    plt.ylim(0, chi2max*3)
    plt.tight_layout()
    # plt.show()
    plt.savefig("%s/plots/sf-fe57-profile.pdf" % dsi.latSWDir)


def combineProfileM1():

    chi2max = 1.355

    # === 1. show effect of fixed vs floating width peak ===

    tf1 = TFile("%s/data/rs-plc-0-axFe_M1.root" % dsi.latSWDir)
    hP1 = tf1.Get("hP_0")
    hT1 = hP1.GetTitle().split()
    pars = []
    for v in hT1:
        try: pars.append(float(v))
        except: pass
    pCL1, intLo1, bestFit1, intHi1, effCorr1 = pars
    xP1, yP1, xpb1 = wl.npTH1D(hP1)
    xP1, yP1 = effCorr1 * (xP1[1:] - xpb1/2), yP1[1:]
    intHi1 *= effCorr1
    pyCts1 = xP1[np.where(yP1>=chi2max)][0]-xpb1/2
    print("fixed: True 90CL (effcorr), %.3f" % pyCts1)

    tf2 = TFile("%s/data/rs-plc-0-axFe_M1-float.root" % dsi.latSWDir)
    hP2 = tf2.Get("hP_0")
    hT2 = hP2.GetTitle().split()
    pars = []
    for v in hT2:
        try: pars.append(float(v))
        except: pass
    pCL2, intLo2, bestFit2, intHi2, effCorr2 = pars
    xP2, yP2, xpb2 = wl.npTH1D(hP2)
    xP2, yP2 = effCorr2 * (xP2[1:] - xpb2/2), yP2[1:]
    intHi2 *= effCorr2
    pyCts2 = xP2[np.where(yP2>=chi2max)][0]-xpb2/2
    print("float: True 90CL (effcorr), %.3f" % pyCts2)

    plt.close()

    plt.axhline(chi2max, c='m', lw=2, label=r"$\chi^2\mathregular{/2\ (90\%\ C.L.)}}$")

    plt.plot(xP1, yP1, "-b", lw=2, label="Fixed-width, Fe-57 (M1): %.2f cts (90%% CL, eff.corr.)" % (pyCts1))
    plt.plot([pyCts1,pyCts1],[0,chi2max], '-b', lw=2, alpha=0.5)

    plt.plot(xP2, yP2, "-r", lw=2, label="Float-width, Fe-57 (M1): %.2f cts (90%% CL, eff.corr.)" % (pyCts2))
    plt.plot([pyCts2,pyCts2],[0,chi2max], '-r', lw=2, alpha=0.5)

    plt.xlabel(r"$\mathregular{N_{obs}}$", ha='right', x=1)
    plt.ylabel(r"-log $\mathregular{\lambda(\mu_{axion})}$", ha='right', y=1)
    plt.legend(loc=1, fontsize=13)
    plt.ylim(0, chi2max*3)
    plt.tight_layout()
    # plt.show()
    plt.savefig("%s/plots/sf-fe57-systematic.pdf" % dsi.latSWDir)


    # === 2. show effect of upper and lower efficiency correction.
    tf3 = TFile("%s/data/rs-plc-0-floatM1.root" % dsi.latSWDir)
    hP3 = tf3.Get("hP_0")
    hT3 = hP3.GetTitle().split()
    pars = []
    for v in hT3:
        try: pars.append(float(v))
        except: pass
    pCL3, intLo3, bestFit3, intHi3, effCorr3 = pars
    xP3, yP3, xpb3 = wl.npTH1D(hP3)
    xP3, yP3 = effCorr3 * (xP3[1:] - xpb3/2), yP3[1:]
    intHi3 *= effCorr3
    pyCts3 = xP3[np.where(yP3>=chi2max)][0]-xpb3/2
    print("regular float: True 90CL (effcorr), %.3f" % pyCts3)

    tf4 = TFile("%s/data/rs-plc-0-effLoM1.root" % dsi.latSWDir)
    hP4 = tf4.Get("hP_0")
    hT4 = hP4.GetTitle().split()
    pars = []
    for v in hT4:
        try: pars.append(float(v))
        except: pass
    pCL4, intLo4, bestFit4, intHi4, effCorr4 = pars
    xP4, yP4, xpb4 = wl.npTH1D(hP4)
    xP4, yP4 = effCorr4 * (xP4[1:] - xpb4/2), yP4[1:]
    intHi4 *= effCorr4
    pyCts4 = xP4[np.where(yP4>=chi2max)][0]-xpb4/2
    print("lower float: True 90CL (effcorr), %.3f" % pyCts4)

    tf5 = TFile("%s/data/rs-plc-0-effHiM1.root" % dsi.latSWDir)
    hP5 = tf5.Get("hP_0")
    hT5 = hP5.GetTitle().split()
    pars = []
    for v in hT5:
        try: pars.append(float(v))
        except: pass
    pCL5, intLo5, bestFit5, intHi5, effCorr5 = pars
    xP5, yP5, xpb5 = wl.npTH1D(hP5)
    xP5, yP5 = effCorr5 * (xP5[1:] - xpb5/2), yP5[1:]
    intHi5 *= effCorr5
    pyCts5 = xP5[np.where(yP5>=chi2max)][0]-xpb5/2
    print("upper float: True 90CL (effcorr), %.3f" % pyCts5)

    plt.close()

    plt.axhline(chi2max, c='m', lw=2, label=r"$\chi^2\mathregular{/2\ (90\%\ C.L.)}}$")

    plt.plot(xP3, yP3, "-b", lw=2, label="Best-fit Eff, Fe-57 (M1): %.2f cts (90%% CL, eff.corr.)" % (pyCts3))
    plt.plot([pyCts3,pyCts3],[0,chi2max], '-b', lw=2, alpha=0.5)

    plt.plot(xP4, yP4, "-r", lw=2, label="Lower Eff, Fe-57 (M1): %.2f cts (90%% CL, eff.corr.)" % (pyCts4))
    plt.plot([pyCts4,pyCts4],[0,chi2max], '-r', lw=2, alpha=0.5)

    plt.plot(xP5, yP5, "-g", lw=2, label="Upper Eff, Fe-57 (M1): %.2f cts (90%% CL, eff.corr.)" % (pyCts5))
    plt.plot([pyCts5,pyCts5],[0,chi2max], '-g', lw=2, alpha=0.5)

    plt.xlabel(r"$\mathregular{N_{obs}}$", ha='right', x=1)
    plt.ylabel(r"-log $\mathregular{\lambda(\mu_{axion})}$", ha='right', y=1)
    plt.legend(loc=1, fontsize=13)
    plt.ylim(0, chi2max*3)
    plt.tight_layout()
    # plt.show()
    plt.savefig("%s/plots/sf-fe57-systematicEff.pdf" % dsi.latSWDir)


def hadronicCurve():
    # reproduce kris's 14.4 kev hadronic axion number and put a new one

    # edelweiss values
    massEW = np.asarray([
    3.98e-03, 8.78e-01, 1.20e+00, 1.52e+00, 1.81e+00, 2.66e+00, 3.30e+00,
    3.90e+00, 5.07e+00, 5.90e+00, 6.67e+00, 7.30e+00, 7.86e+00, 8.34e+00,
    8.79e+00, 9.18e+00, 9.52e+00, 9.85e+00, 1.02e+01, 1.08e+01, 1.12e+01,
    1.16e+01, 1.19e+01, 1.23e+01, 1.27e+01, 1.30e+01, 1.34e+01, 1.38e+01,
    1.40e+01, 1.41e+01, 1.41812e+01, 1.42393e+01, 1.42829e+01, 1.4312e+01,
    1.44e+01, 1.44e+01, 1.44e+01, 1.44e+01, 1.44e+01, 1.44e+01, 1.44e+01])
    coupEW = np.asarray([
    4.79e-17, 4.90e-17, 4.90e-17, 4.90e-17, 4.90e-17, 4.87e-17, 4.92e-17,
    5.01e-17, 5.13e-17, 5.25e-17, 5.38e-17, 5.50e-17, 5.64e-17, 5.77e-17,
    5.91e-17, 6.05e-17, 6.20e-17, 6.35e-17, 6.56e-17, 6.91e-17, 7.28e-17,
    7.60e-17, 8.04e-17, 8.63e-17, 9.45e-17, 1.03e-16, 1.20e-16, 1.46e-16,
    1.84e-16, 2.12e-16, 2.44e-16, 2.87e-16, 3.36e-16, 3.94e-16, 6.62e-16,
    9.07e-16, 1.67e-15, 3.35e-15, 6.21e-15, 1.27e-14, 2.26e-14])

    # vorren values
    massV = np.asarray([
    0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8, 8.5,
    9, 9.5, 10, 10.5, 11, 11.5, 12, 12.5, 13, 13.5, 14, 14.2, 14.3, 14.35,
    14.39, 14.399, 14.3999])

    tf = TFile("%s/data/specPDFs.root" % dsi.latSWDir)
    h1 = tf.Get("h1") # ge photoelectric XS,       [cm^2 / kg]
    h2 = tf.Get("h2") # axioelectric XS, m_A = 0   [cm^2 / kg]

    # check 1 - verify the m_a = 0 matches the histogram.  OK cgw
    E, m = 14.4, 0
    sigAe_hist = h2.GetBinContent(h2.FindBin(E))
    sigAe = wl.sig_ae(E, m) * h1.GetBinContent(h1.FindBin(E)) # allow m_a to float. [cm^2 / kg]
    # print("%.3e  %.3e" % (sigAe, sigAe_hist))

    # check 2 - verify this reproduces vorren's curve in GAT/MDStat/KVBosonic/hadronicLikelihood2.py  OK cgw
    photo144 = h1.GetBinContent(h1.FindBin(E))
    coupV = np.zeros(len(massV))
    for i in range(len(massV)):
        b = np.sqrt( (E**2 - massV[i]**2)/(E**2) )
        sigAe = wl.sig_ae(E, massV[i]) * photo144
        coupV[i] = np.sqrt(11.2 / 478 / (b**3) / 4.56e23 / 86400 / sigAe)
        # print(coupV[i])

    nObs = 17.165 # floating width, best fit << use this one
    # nObs = 17.949 # floating width, lower efficiency bound (bad)
    # nObs = 16.566 # floating width, upper efficiency bound (good)

    # check 3 - make a new curve
    photo144 = h1.GetBinContent(h1.FindBin(E))
    coupV2 = np.zeros(len(massV))
    for i in range(len(massV)):
        b = np.sqrt( (E**2 - massV[i]**2)/(E**2) )
        sigAe = wl.sig_ae(E, massV[i]) * photo144
        coupV2[i] = np.sqrt(nObs / (detExp * (b**3) * 4.56e23 * 86400 * sigAe))
        if massV[i]==0:
            print("lowest mass:",coupV2[i])

    plt.semilogy(massEW, coupEW, c='orange', label="EDELWEISS (2013)")
    plt.semilogy(massV, coupV, c='r', lw=3, label="MJD (2017): 11.2 cts (90% CL), 478 kg-d")
    plt.semilogy(massV, coupV2, c='b', lw=3, label="MJD (2018): %.2f cts (90%% CL), %.0f kg-d" % (nObs,detExp))

    plt.xlim(0, 15)
    plt.ylim(1e-17, 1e-15)
    plt.legend(loc=2, fontsize=14)
    plt.xlabel("Axion Mass (keV)", ha='right', x=1)
    plt.ylabel(r"$\mathregular{ g_{ae} \times g^{eff}_{an} }$", ha='right', y=1)
    plt.tight_layout()
    # plt.show()
    plt.savefig("%s/plots/sf-hadronic144axion-limit.pdf" % dsi.latSWDir)


def gaeProj():
    """
    Project a limit on g_ae based on the number of counts in the data (assuming c/kev/kg-d).
    This assumes that (on average) in the energy range, the axion amplitude is no higher
    than the statistical error bar.
    # g_ae = (N_observed / N_expected)**(1./4.)
    # g_ae = (B dE / Mt)**(1./8.) * (axFlux)**(-1./4.)
    """
    eLo, eHi = 1.5, 8

    f = TFile("%s/data/specPDFs.root" % dsi.latSWDir)
    hAx = f.Get("h4")
    hTr = f.Get("h5")

    def gaeBkg(expo, B, eLo, eHi):
        axFlux = hAx.Integral(hAx.FindBin(eLo), hAx.FindBin(eHi), "width") * 365.25 # [cts / kg-y]
        return ((B * (eHi-eLo) * 365.25)/(expo))**0.125 * axFlux**-0.25

    # === 1. projection based on tritium rate, to compare w/ mike marino's thesis.
    # he figured we would be tritium dominated, at a rate of 0.03 c/kkd (from 0-20, the whole trit spec)
    # so calculate what rate 0.03 corresponds to in our energy region [eLo, eHi] (it will be a little higher)

    B = 0.03 # what an optimist
    x, y, xpb = wl.npTH1D(hTr)
    x, y = normPDF(x, y, 0, 20)
    avgFac = B / np.average(y)
    y *= avgFac
    idx = np.where((x>=eLo) & (x<=eHi))
    tritB = np.average(y[idx])
    plt.plot(x, y, lw=2, c='b', label="Tritium Spectrum")
    plt.axhline(B, c='r', lw=2, label='avg, %.1f-%.1f keV: %.2f' % (0, 20, B))
    plt.axhline(tritB, c='g', lw=2, label='avg, %.1f-%.1f keV: %.3f' % (eLo, eHi, tritB))
    plt.legend(loc=3)
    plt.xlabel("Energy (keV)", ha='right', x=1)
    plt.ylabel("Counts / (keV kg-d)", ha='right', y=1)
    plt.tight_layout()
    # plt.show()
    plt.savefig("%s/plots/sf-tritBkgIdx.pdf" % (dsi.latSWDir))

    plt.close()

    # === 2. projection based on DS1--5C analysis
    # plt.figure(figsize=(6,5))

    # # Efficiency corrected rate, 1.5-20.0 keV:  0.062 c/(kev-kg-d)
    # # Efficiency corrected rate, 1.5-8.0 keV:  0.119 c/(kev-kg-d)
    B = 0.119

    xE = np.arange(1, 1000, 0.1) # 1 ton-y LEGEND exposure if it was magically a 1 kev detector

    # A key assumption here is that we are not seeing a significant number of counts in the axion pdf.
    # (we *are* in the current analysis, which is why the observed value is such an apparent outlier.)

    plt.axhline(2.6e-11, c='g', label="MALBEK (2015)")

    plt.axhline(1.92e-11, c='y', label="MJD DS1--5C (Nonzero Axion Amp.), %.2f kg-y" % (detExp/365.25))

    plt.plot(xE, gaeBkg(xE, tritB, eLo, eHi), c='b', label=r'Tritium-dominated (proj.), B=%.2f cts/(kev kg-d)' % tritB)

    plt.plot(xE, gaeBkg(xE, B, eLo, eHi), c='r', label=r'Excess-dominated (proj.), B=%.2f cts/(kev kg-d)' % B)

    plt.axhline(4.35e-12, c='m', label="PandaX (2017)")
    plt.axhline(3.5e-12, c='k', label="LUX (2017)")

    plt.xlabel("Exposure (kg-y)", ha='right', x=1)

    plt.ylabel(r"$\mathregular{g_{ae}}$ Projected UL", ha='right', y=1)

    plt.ylim(1e-12, 3e-11)

    plt.legend(loc=1, bbox_to_anchor=(0., 0.73, 1, 0.2), fontsize=13)
    plt.tight_layout()
    # plt.show()
    plt.savefig("%s/plots/gae-proj.pdf" % dsi.latSWDir)


    # === 3. insane LEGEND 10 ton-yr projection
    plt.close()
    plt.figure(figsize=(6,5))

    B = 0.119

    xE = np.arange(1, 10000, 1) # 10 ton-y LEGEND exposure if it was magically a 1 kev detector

    # A key assumption here is that we are not seeing a significant number of counts in the axion pdf.
    # (we *are* in the current analysis, which is why the observed value is such an apparent outlier.)

    # plt.axhline(2.6e-11, c='g', label="MALBEK (2015)")

    plt.axhline(1.92e-11, c='g', label="MJD DS1--5C (Nonzero Axion Amp.), %.2f kg-y" % (detExp/365.25))

    plt.loglog(xE, gaeBkg(xE, tritB, eLo, eHi), c='b', label=r'LEGEND, Tritium-dominated (proj.), B=%.2f cts/(kev kg-d)' % tritB)

    plt.loglog(xE, gaeBkg(xE, B, eLo, eHi), c='r', label=r'LEGEND, Excess-dominated (proj.), B=%.2f cts/(kev kg-d)' % B)

    plt.axhline(4.35e-12, c='m', label="PandaX (2017)")
    plt.axhline(3.5e-12, c='k', label="LUX (2017)")

    plt.xlabel("Exposure (kg-y)", ha='right', x=1)

    plt.ylabel(r"$\mathregular{g_{ae}}$ Projected UL", ha='right', y=1)

    # plt.ylim(3e-12, 1e-10)

    plt.legend(loc=1, fontsize=10, bbox_to_anchor=(0., 0.73, 1, 0.2))
    plt.tight_layout()
    # plt.show()
    plt.savefig("%s/plots/gae-proj-legend.pdf" % dsi.latSWDir)


if __name__=="__main__":
    main(sys.argv[1:])