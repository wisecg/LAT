#!/usr/bin/env python
""" Let's reproduce ALL the input PDF's for the specFit routine,
    get ALL the correct scaling factors, and make sure there are NO surprises.
    C. Wiseman 23 Oct. 2017
"""
import ROOT, warnings
from ROOT import TFile, TTree, TCanvas, TH1D, TLegend, TPad, TLine, TGraph
from ROOT import gStyle, gPad, gROOT, std
import numpy as np
from array import array
gStyle.SetOptStat(0)
gROOT.ProcessLine("gErrorIgnoreLevel = 3001;") # suppress ROOT messages

def main():
    # runAllPlots()
    generatePDFs()


def runAllPlots():
    c = TCanvas("c","c",800,600)
    c.SetLogy()
    c.SetGrid(1,1)

    # CONSTANTS
    kpb = 0.02
    eLo, eHi = 0.1, 12.
    binRange = 2.
    pks = [1.739, 1.836, 2.307, 2.464, 6.404]
    frankFlux = [4.95e+38, 4.95e+38, 3.94e+38, 2.27e+38, 4.06e+38]

    redondoScale = 1e19 * 0.511e-10**-2 # convert table to [cts / (keV cm^2 d)]
    jcapScale = 1e-13**2. * 365 * 1e4  * 1e-20 # (gae in paper * per yr * per m^2 * 10^-20 scaling)
    barnsPerAtom = 120.5
    nAvo = 6.0221409e+23
    phoScale = 72.64 / nAvo # (ge mol mass / nAvo)
    axConvScale = phoScale / 1000 # (ge mol mass / nAvo / 1000 g/kg)
    beta = 1.
    sigAeScale = beta**-1 * 3./(16. * np.pi * (1./137.) * 511.**2.) * (1 - beta**(2./3.)/3)

    # 1. Reproduce Redondo's data and Frank's axion flux points.

    axData, axEne, axFlux = [], [], []
    with open("../data/redondoFlux.txt") as f1:
        lines = f1.readlines()[11:]
        for line in lines:
            data = line.split()
            axData.append([float(data[0]),float(data[1])])
            axEne.append(float(data[0]))
            axFlux.append(float(data[1]) * redondoScale)
    axData = np.array(axData)

    g1 = TGraph(len(axEne), array('d', axEne), array('d', axFlux))
    g1.SetTitle(" ")
    g1.GetXaxis().SetTitle("Energy (kev)")
    g1.GetYaxis().SetTitle("flux (kev^{-1} cm^{-2} d^{-1})")
    g1.SetMarkerColor(ROOT.kRed)
    g1.SetMaximum(1.5e40)
    g1.SetMinimum(8e37)
    g1.Draw()
    l1 = TLine(pks[4], min(axFlux), pks[4], max(axFlux))
    l1.SetLineColor(ROOT.kRed)
    l1.Draw("same")

    # verify the algorithm to translate into a hist (used for RooFit) works correctly
    nBins = int((eHi-eLo)/kpb)
    h1 = TH1D("h1","h1",nBins, eLo, eHi)
    for i in range(nBins):
        ene = i * kpb + eLo
        eneLo, eneHi = ene - kpb/2., ene + kpb/2.
        idx = np.where((axData[:,0] >= eneLo) & (axData[:,0] <= eneHi))
        flux = 0.
        with warnings.catch_warnings(): # ignore "mean of empty slice" errors (they're harmless)
            warnings.simplefilter("ignore",category=RuntimeWarning)
            flux = np.mean(axData[idx][:,1]) * redondoScale
        if np.isnan(flux): flux = 0.
        h1.SetBinContent(i, flux)
    h1.SetLineColor(ROOT.kBlue)
    h1.Draw("same")

    g2 = TGraph(len(pks),array('d',pks),array('d',frankFlux))
    g2.SetMarkerStyle(20)
    g2.SetMarkerColor(ROOT.kRed)
    g2.Draw("P same")

    clintFlux = []
    ax = h1.GetXaxis()
    for pk in pks:
         # "* kpb" is same as using "width" option in integral
        intFlux = h1.Integral(ax.FindBin(pk-kpb*binRange), ax.FindBin(pk+kpb*binRange)) * kpb
        clintFlux.append(intFlux)
    g3 = TGraph(len(pks), array('d',pks), array('d', clintFlux))
    g3.SetMarkerStyle(20)
    g3.SetMarkerColor(ROOT.kGreen)
    g3.Draw("P same")

    leg = TLegend(0.6,0.6,0.85,0.85)
    leg.AddEntry(g1,"raw data","p")
    leg.AddEntry(h1,"histo 0.02 kev/bin", "l")
    leg.AddEntry(l1,"6.404 kev","l")
    leg.AddEntry(g2,"frank flux","p")
    leg.AddEntry(g3,"clint flux","p")
    leg.Draw("same")

    c.Print("../plots/axFluxCompare.pdf")


    # 2. Reproduce Graham's axion flux plot

    c.Clear()
    axFluxCm2Sec = [f * 86400**-1. for f in axFlux]
    g4 = TGraph(len(axEne), array('d', axEne), array('d', axFluxCm2Sec))
    g4.SetTitle(" ")
    g4.GetXaxis().SetTitle("Energy (kev)")
    g4.GetYaxis().SetTitle("flux (kev^{-1} cm^{-2} s^{-1})")
    g4.SetMarkerColor(ROOT.kBlue)
    g4.Draw()
    c.Print("../plots/axFluxGraham.pdf")


    # 3. Reproduce Redondo's JCAP figure 2 in : https://arxiv.org/pdf/1310.0823v1.pdf

    c.Clear()
    c.SetGrid(0,0)
    c.SetCanvasSize(800,800)
    c.SetLogy(0)
    axFluxJCAP = [f * jcapScale for f in axFlux]
    g5 = TGraph(len(axEne), array('d', axEne), array('d', axFluxJCAP))
    g5.SetTitle(" ")
    g5.GetXaxis().SetTitle("Energy (keV)")
    g5.GetYaxis().SetTitle("flux 10^{-20} kev^{-1} y^{-1}  m^{-2}")
    g5.SetMarkerColor(ROOT.kBlack)
    g5.SetMaximum(3.1)
    g5.SetMinimum(-0.1)
    g5.Draw()
    c.Print("../plots/axFluxJCAP.pdf")


    # 4. Plot the photoelectric cross section from mucal.

    phoData, phoEne, phoVal = [], [], []
    with open("../data/ge76peXS.txt") as f2: # 2499 entries, 0.01 kev intervals
        lines = f2.readlines()
        for line in lines:
            data = line.split()
            phoData.append([float(data[0]),float(data[1])])
            phoEne.append(float(data[0]))
            phoVal.append(float(data[1]))
    phoData, phoEne, phoVal = np.array(phoData), np.array(phoEne), np.array(phoVal)

    c.Clear()
    c.SetCanvasSize(800,600)
    c.SetLogy(1)
    phoValScaled = [v * barnsPerAtom for v in phoVal] # barns/atom
    g6 = TGraph(len(phoEne), array('d',phoEne), array('d',phoVal))
    g6.SetTitle(" ")
    g6.GetXaxis().SetTitle("Energy (keV)")
    g6.GetYaxis().SetTitle("#sigma_{pe} (cm^{2}/g)")
    g6.GetXaxis().SetRangeUser(0.1, 12.)
    g6.SetMinimum(10)
    g6.SetMaximum(3e5)
    g6.SetMarkerColor(ROOT.kBlue)
    g6.Draw()
    c.Print("../plots/mucalPho.pdf")


    # 5. Plot the axioelectric cross section (2 vals of beta, 1 and 5)

    c.Clear()
    c.SetGrid(1,1)
    axioVal = [phoVal[idx] * sigAeScale * phoEne[idx]**2. * barnsPerAtom for idx in range(len(phoVal))]

    g7 = TGraph(len(phoEne), array('d',phoEne), array('d',axioVal))
    g7.Draw("AL")

    g7.GetXaxis().SetTitle("Energy (keV)")
    g7.GetYaxis().SetTitle("#sigma_{ae} (barns/atom)")
    g7.GetXaxis().SetRangeUser(0.1, 12.)
    g7.SetTitle(" ")
    g7.SetLineColor(ROOT.kBlue)
    g7.SetLineWidth(2)
    g7.SetMinimum(1)
    g7.SetMaximum(6e2)
    c.Update()


    # 5 kev axion (new formula for sigma_ae)
    def sig_ae(E,m):
        beta = (1 - m**2./E**2.)**(1./2)
        return (1 - (1./3.)*beta**(2./3.)) * (3. * E**2.) / (16. * np.pi * (1./137.) * 511.**2. * beta)
    idx = np.where(phoEne > 5)
    axioVal2 = [phoVal[idx][i] * sig_ae(phoEne[idx][i],5.) * barnsPerAtom for i in range(len(phoVal[idx]))]
    axioVal2.insert(0,0.)
    phoEne2 = phoEne[idx].tolist()
    phoEne2.insert(0,5.)

    g8 = TGraph(len(phoEne2), array('d',phoEne2), array('d',axioVal2))
    g8.SetLineColor(ROOT.kRed)
    g8.Draw("L same")

    leg2 = TLegend(0.65,0.75,0.85,0.85)
    leg2.AddEntry(g7,"m_{a}= 0 keV","l")
    leg2.AddEntry(g8,"m_{a}= 5 keV","l")
    leg2.Draw("same")

    # convert to histogram (for RooDataHist)
    nBins = int((eHi - eLo)/kpb)
    h2 = TH1D("h2","h2",nBins, eLo, eHi)
    for i in range(nBins):
        ene = i * kpb + eLo
        eneLo, eneHi = ene - kpb/2., ene + kpb/2.
        # ignore "mean of empty slice" errors (they're harmless)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore",category=RuntimeWarning)
            idx = np.where((phoData[:,0] >= eneLo) & (phoData[:,0] <= eneHi))
            pho = np.mean(phoData[idx][:,1])
            if np.isnan(pho) or len(phoData[idx][:,1]) == 0: pho = 0.
            axio = pho * ene**2. * sigAeScale * barnsPerAtom
            h2.SetBinContent(i, axio)
    h2.SetLineColor(ROOT.kBlack)
    h2.Draw("hist same")

    c.Print("../plots/axioElectric.pdf")

    # 6. convolve the axion flux w/ the axioelectric effect.  (bins are different so use a histogram)

    nBins = int((eHi - eLo)/kpb)
    h3 = TH1D("h3","",nBins,eLo,eHi)
    for i in range(nBins):
        ene = i * kpb + eLo
        eneLo, eneHi = ene - kpb/2., ene + kpb/2.
        # ignore "mean of empty slice" errors (they're harmless)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore",category=RuntimeWarning)

            # axioelectric x-section [cm^2 / kg]
            idx = np.where((phoData[:,0] >= eneLo) & (phoData[:,0] <= eneHi))
            pho = np.mean(phoData[idx][:,1])
            if np.isnan(pho) or len(phoData[idx][:,1]) == 0: pho = 0.
            axio = pho * ene**2. * sigAeScale * 1000.

            # axion flux [cts / (cm^2 d keV)]
            idx = np.where((axData[:,0] >= eneLo) & (axData[:,0] <= eneHi))
            flux = np.mean(axData[idx][:,1]) * redondoScale
            if np.isnan(flux): flux = 0.

            # convolved [cts / (keV d kg)]
            axConv = axio * flux
            h3.SetBinContent(i, axConv)

    c.Clear()
    h3.SetLineColor(ROOT.kBlue)
    h3.SetMinimum(2e39)
    h3.SetMaximum(4e42)
    h3.GetXaxis().SetTitle("Energy (keV)")
    h3.GetYaxis().SetTitle("flux cts kev^{-1} d^{-1} kg^{-1}")
    h3.Draw("hist")
    c.Print("../plots/axionConv.pdf")


    # 7. reproduce the rest of frank's numbers in his tables from MY HISTOGRAMS
    # and calculate g_ae

    print "sig_ae factor: ",sigAeScale

    binRange = 2.
    phos, axos, cFlux, cRateHist = [], [], [], []
    ax1, ax3 = h1.GetXaxis(), h3.GetXaxis()
    for pk in pks:
        eneLo, eneHi = pk - kpb, pk + kpb

        # photoelectric
        idx = np.where((phoData[:,0] >= eneLo) & (phoData[:,0] <= eneHi))
        pho = np.mean(phoData[idx][:,1]) * phoScale # [cm^2 / atom]
        phos.append( pho )

        # axioelectric
        axos.append( pho * pk**2. * sigAeScale)

        # axion flux
        flux = h1.Integral(ax1.FindBin(pk-kpb*binRange), ax1.FindBin(pk+kpb*binRange)) * kpb
        cFlux.append(flux)

        # axion rate - integrate histogram directly (proof it's correct if it matches frank's table)
        conv = h3.Integral(ax3.FindBin(pk-kpb*binRange), ax3.FindBin(pk+kpb*binRange)) * kpb * axConvScale
        cRateHist.append(conv)

    frankRates = [axos[i] * frankFlux[i] for i in range(len(axos))]
    cRateTable = [axos[i] * cFlux[i] for i in range(len(axos))]

    # stackoverflow trick to make list printing prettier
    class prettyfloat(float):
        def __repr__(self): return "%.2e" % self

    print "T3,C2: E^2 * (sigAeScale)             ", map(prettyfloat, [np.power(ax, 2.) * sigAeScale for ax in pks])
    print "T3,C3: sig_pe (cm^2/atom)             ", map(prettyfloat, phos)
    print "T3,C4: sig_ae                         ", map(prettyfloat, axos)
    print "T4,C2 (old): Phi_a (cm^2/d)           ", map(prettyfloat, frankFlux)
    print "T4,C2 (new):                          ", map(prettyfloat, clintFlux)
    print "T4,C4 (table): Phi_a * sig_ae (cts/d) ", map(prettyfloat, cRateTable)
    print "T4,C4 (histo):                        ", map(prettyfloat, cRateHist)

    # exposure,  expected counts, upper bound on g_ae
    malbekExpo = 89.5 # kg-d
    exposure = malbekExpo * 1000 * (1./72.64) * nAvo # [atom-d] = [kg d][1000 g/kg][1/72.64 mol/g][nAvo atom/mol]

    N_obs = 10. # just a guess

    # Compare methods
    N_exp = exposure * sum(frankFlux[i] * axos[i] for i in range(4))
    g_ae = np.power(N_obs / N_exp, 1./4.) # upper bound
    print "Frank's Table g_ae: %.2e" % g_ae

    N_exp = exposure * sum(cRateTable[i] for i in range(4))
    g_ae = np.power(N_obs / N_exp, 1./4.) # upper bound
    print "Clint's Table g_ae: %.2e" % g_ae

    N_exp = exposure * sum(cRateHist[i] for i in range(4))
    g_ae = np.power(N_obs / N_exp, 1./4.) # upper bound
    print "Clint's Histo g_ae: %.2e" % g_ae

    print "Histo peaks expected counts: ",N_exp
    N_exp = h3.Integral(ax3.FindBin(1.5), ax3.FindBin(8.)) * kpb * axConvScale * exposure
    print "Continuum expected counts: ",N_exp


    # 8. do the tritium plot

    tritData, tritEne, tritVal = [], [], []
    with open("../data/TritiumSpectrum.txt") as f3: # 20000 entries
        lines = f3.readlines()[1:]
        for line in lines:
            data = line.split()
            conv = float(data[2]) # raw spectrum convolved w/ ge cross section
            if conv < 0: conv = 0.
            tritData.append([float(data[1]),conv])
            tritEne.append(float(data[1]))
            tritVal.append(conv)
    tritData = np.array(tritData)

    c.Clear()
    c.SetLogy(0)
    g8 = TGraph(len(tritEne),array('d',tritEne),array('d',tritVal))
    g8.SetTitle(" ")
    g8.GetXaxis().SetTitle("Energy (keV)")
    g8.GetYaxis().SetTitle("Counts (norm.)")
    g8.SetMarkerColor(ROOT.kRed)
    g8.Draw()

    kpb = 0.02
    eLo, eHi = 0.1, 20.
    nBins = int((eHi-eLo)/kpb)
    h4 = TH1D("h4","h4",nBins,eLo,eHi)
    for i in range(nBins):
        ene = i * kpb + eLo
        eneLo, eneHi = ene - kpb/2., ene + kpb/2.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore",category=RuntimeWarning)
            idx = np.where((tritData[:,0] >= eneLo) & (tritData[:,0] <= eneHi))
            trit = np.mean(tritData[idx][:,1])
            if np.isnan(trit): trit = 0.
            h4.SetBinContent(i, trit)
            if trit==0.: print "Trit bin %d is zero" % i

    h4.SetLineColor(ROOT.kBlue)
    h4.Draw("hist same")
    c.Print("../plots/tritSpec.pdf")


def generatePDFs(ma=0):
    """ Generate pdfs for fitting. Takes axion mass (in keV) as a parameter. """

    f = TFile("../data/inputHists.root","RECREATE")
    redondoScale = 1e19 * 0.511e-10**-2 # convert table to [cts / (keV cm^2 d)]

    axData, phoData, tritData = [], [], []
    with open("../data/redondoFlux.txt") as f1: # 23577 entries
        lines = f1.readlines()[11:]
        for line in lines:
            data = line.split()
            axData.append([float(data[0]),float(data[1])])
    with open("../data/ge76peXS.txt") as f2: # 2499 entries, 0.01 kev intervals
        lines = f2.readlines()
        for line in lines:
            data = line.split()
            phoData.append([float(data[0]),float(data[1])])
    with open("../data/TritiumSpectrum.txt") as f3: # 20000 entries
        lines = f3.readlines()[1:]
        for line in lines:
            data = line.split()
            conv = float(data[2]) # raw spectrum convolved w/ ge cross section
            if conv < 0: conv = 0.
            tritData.append([float(data[1]),conv])
    axData, phoData, tritData = np.array(axData), np.array(phoData), np.array(tritData)

    def sig_ae(E,m):
        beta = (1 - m**2./E**2.)**(1./2)
        return (1 - (1./3.)*beta**(2./3.)) * (3. * E**2.) / (16. * np.pi * (1./137.) * 511.**2. * beta)

    kpb = 0.02
    eLo, eHi = 0.1, 30.
    nBins = int((eHi-eLo)/kpb)
    h1 = TH1D("h1","photoelectric",nBins,eLo,eHi)  # [cm^2 / kg]
    h2 = TH1D("h2","axioelectric",nBins,eLo,eHi)   # [cm^2 / kg]
    h3 = TH1D("h3","axion flux",nBins,eLo,eHi)     # [cts / (keV cm^2 d)]
    h4 = TH1D("h4","convolved flux",nBins,eLo,eHi) # [cts / (keV d kg)]
    h5 = TH1D("h5","tritium",nBins,eLo,eHi)        # [cts] (normalized to 1)

    for i in range(nBins):
        ene = i * kpb + eLo
        eneLo, eneHi = ene - kpb/2., ene + kpb/2.
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

    # normalize tritium to 1 (see sandbox/gae2.py::tritBkgIndex for how to normalize to a 3H activation rate)
    normTrit = h5.Integral("width")
    h5.Scale(1./normTrit)

    h1.Write()
    h2.Write()
    h3.Write()
    h4.Write()
    h5.Write()
    f.Close()


if __name__ == "__main__":
    main()

