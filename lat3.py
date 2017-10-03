#!/usr/local/bin/python
"""
===================== LAT3.py =====================

PLACEHOLDER CODE.  THIS SHOULD BE A LAT SKIMMER,
WHICH READS DATASETINFO.PY AND PRODUCES CHANNEL-SPECIFIC
FILES WITH RUN CUTS APPLIED.

or should this be where "tuneCuts" is done?

v1: 07 Aug 2017

========= C. Wiseman (USC), B. Zhu (LANL) =========
"""
import sys, time, ROOT, glob
sys.argv += [ '-b' ] # force ROOT to be loaded in batch mode.
from ROOT import gROOT, gStyle, gPad
from ROOT import TFile, TTree, TChain, TCanvas, TH1D, TH2D, TF1, TLegend, TLine, TGraph
import numpy as np
from DataSetInfo import CalInfo
import DataSetInfo as ds
import waveLibs as wl
from decimal import Decimal

def main(argv):

    print "========================================"
    print "LAT3 started:",time.strftime('%X %x %Z')
    startT = time.clock()

    cInfo = CalInfo()
    pathToInput = "."
    dsNum, subNumm, modNum, chNum = -1, -1, -1, -1
    calTree = ROOT.TChain("skimTree")
    customPar = ""
    calList, parList, parNameList, chList = [], [], [], []
    fTune, fFor, fUpd, fastMode = False, False, False, False

    if len(argv) == 0:
        return
    for i, opt in enumerate(argv):
        # -- Cut tuning options --
        if opt == "-all":
            parList.append('bcMax'), parNameList.append('bcMax')
            parList.append('(waveS4-waveS1)/bcMax/trapENFCalC'), parNameList.append('noiseWeight')
            parList.append('(bandTime-tOffset-1100)/(matchTime-tOffset)'), parNameList.append('bcTime')
            parList.append('pol2'), parNameList.append('pol2')
            parList.append('pol3'), parNameList.append('pol3')
            parList.append('fitSlo'), parNameList.append('fitSlo')
            parList.append('riseNoise'), parNameList.append('riseNoise')
            print "Tuning all cuts"
        if opt == "-bcMax":
            parList.append('bcMax'), parNameList.append('bcMax')
            print "Tuning bcMax"
        if opt == "-noiseWeight":
            parList.append('(waveS4-waveS1)/bcMax/trapENFCalC'), parNameList.append('noiseWeight')
            print "Tuning noiseWeight ((waveS4-waveS1)/bcMax/trapENFCalC)"
        if opt == "-bcTime":
            parList.append('(bandTime-tOffset-1100)/(matchTime-tOffset)'), parNameList.append('bcTime')
            print "Tuning bcTime"
        if opt == "-tailSlope":
            parList.append('pol2'), parNameList.append('pol2')
            parList.append('pol3'), parNameList.append('pol3')
            print "Tuning tailSlope"
        if opt == "-fitSlo":
            parList.append('fitSlo'), parNameList.append('fitSlo')
            print "Tuning fitSlo"
        if opt == "-riseNoise":
            parList.append('riseNoise'), parNameList.append('riseNoise')
            print "Tuning riseNoise"
        if opt == "-Custom":
            customPar = str(argv[i+1])
            parList.append(customPar), parNameList.append('customPar')
            print "Tuning custom cut parameter: ", customPar
        # -- Input/output options --
        if opt == "-s":
            dsNum, subNum, modNum = int(argv[i+1]), int(argv[i+2]), int(argv[i+3])
        if opt == "-d":
            pathToInput = argv[i+1]
            print "Custom paths: Input %s" % (pathToInput)
        if opt == "-ch":
            chNum = int(argv[i+1])
            print "Tuning specific channel %d" % (chNum)
        # -- Database options --
        if opt == "-force":
            fFor = True
            print "Force DB update mode."
        if opt == "-tune":
            fTune = True
            print "Cut tune mode."
        if opt == "-upd":
            fUpd = True
            print "File update mode."

    # -- Load calibration files --
    if dsNum == -1 or subNum == -1 or modNum == -1:
        print "DS, subDS, or module number not set properly, exiting"
        return
    else:
        calList = cInfo.GetCalList("ds%d_m%d"%(dsNum, modNum), subNum)
        for i in calList: calTree.Add("%s/latSkimDS%d_run%d_*"%(pathToInput, dsNum, i))

    # -- Load chains for this DS --
    inPath = pathToInput + "/latSkimDS%d*.root" % dsNum
    fileList = glob.glob(inPath)
    cutFile = TFile(fileList[0])
    theCut = cutFile.Get("theCut").GetTitle()

    # -- Load channel list --
    if chNum == -1:
        chList = ds.GetGoodChanList(dsNum)
        if dsNum==5: # remove 692 and 1232
            chList = [584, 592, 598, 608, 610, 614, 624, 626, 628, 632, 640, 648, 658, 660, 662, 672, 678, 680, 688, 690, 694, 1106, 1110, 1120, 1124, 1128, 1170, 1172, 1174, 1176, 1204, 1208, 1298, 1302, 1330, 1332]
    else:
        chList = [chNum]

    print "Processing channels: ", chList
    print "Processing runs: ", calList

    # -- Tune cuts --
    tunedPars = {}
    for par, parName in zip(parList, parNameList):
        cutDict = TuneCut(dsNum, subNum, calTree, chList, par, parName, theCut, fastMode)
        key = "%s_ds%d_idx%d"%(parName,dsNum,subNum)
        print key, cutDict
        if fTune:
            wl.setDBCalRecord({"key":key,"vals":cutDict})

    stopT = time.clock()
    print "Stopped:",time.strftime('%X %x %Z'),"\nProcess time (min):",(stopT - startT)/60


def TuneCut(dsNum, subNum, cal, chList, par, parName, theCut, fastMode):
    c = TCanvas("%s"%(parName),"%s"%(parName),1600,600)
    c.Divide(3,1,0.00001,0.00001)
    cutDict = {}
    for ch in chList:
        cutDict[ch] = [0,0,0,0,0]
        eb, elo, ehi = 40,0,40
        e1dCut = 5.
        d1Cut = theCut + " && trapENFCalC > %d && trapENFCalC < %d && channel==%d" % (e1dCut,ehi,ch)
        d2Cut = theCut + " && channel==%d" % ch
        nPass = cal.Draw("trapENFCalC:%s"%(par), d1Cut, "goff")
        nEnergy = cal.GetV1()
        nCut = cal.GetV2()
        nCutList = list(float(nCut[n]) for n in xrange(nPass))
        nEnergyList = list(float(nEnergy[n]) for n in xrange(nPass))
        vb, vlo, vhi = 100000, np.amin(nCutList), np.amax(nCutList)
        outPlot = "./plots/tuneCuts/%s_ds%d_idx%d_ch%d.png" % (parName,dsNum,subNum,ch)
        cut99,cut95,cut01,cut05,cut90 = MakeCutPlot(c,cal,par,eb,elo,ehi,vb,vlo,vhi,d2Cut,d1Cut,outPlot,fastMode)
        cutDict[ch] = [cut01,cut05,cut90,cut95,cut99]
    return cutDict

def MakeCutPlot(c,cal,var,eb,elo,ehi,vb,vlo,vhi,d2Cut,d1Cut,outPlot,fastMode):
    """ Repeated code is the DEVIL.  Even if you have to pass in 1,000,000 arguments. """

    # Calculate cut vals (assumes plot range is correct)
    h1 = wl.H1D(cal,vb,vlo,vhi,var,d1Cut)
    h1Sum = h1.Integral()
    if h1Sum == 0:
        return 0,0,0,0,0
    h1.Scale(1/h1Sum)
    try:
        cut99,cut95,cut01,cut05,cut90 = wl.GetIntegralPoints(h1)
    except:
        print "Error: Failed", var, "using cut", d2Cut
        return 0,0,0,0,0
    if fastMode:
        return cut99,cut95,cut01,cut05,cut90

    # Generate the plot for inspection.
    c.cd(2)
    gPad.SetLogy(0)
    h1.GetXaxis().SetRangeUser(cut01-abs(0.25*cut01), cut99 + abs(0.25*cut99) )
    h1.Draw("hist")

    c.cd(1)
    gPad.SetLogy(0)
    cal.Draw("%s:trapENFCalC>>b(%d,%d,%d,%d,%.3E,%.3E)"%(var,eb,elo,ehi,vb,cut01-abs(0.25*cut01),cut99+abs(0.25*cut99)) ,d2Cut)

    l1, l2 = TLine(), TLine()
    l1.SetLineColor(ROOT.kGreen)
    l2.SetLineColor(ROOT.kRed)

    l1.DrawLine(elo, cut99, ehi, cut99)
    l2.DrawLine(elo, cut95, ehi, cut95)
    l2.DrawLine(elo, cut05, ehi, cut05)
    l1.DrawLine(elo, cut01, ehi, cut01)

    c.cd(3)
    x_h1, y_h1 = wl.npTH1D(h1)
    int_h1 = wl.integFunc(y_h1)
    g2 = TGraph(len(x_h1), x_h1, int_h1)
    g2.Draw("ACP")
    g2.GetXaxis().SetRangeUser(cut01-abs(0.3*cut01), cut99 + abs(0.3*cut99) )
    l1.DrawLine(cut99, 0, cut99, 1)
    l2.DrawLine(cut95, 0, cut95, 1)
    l1.DrawLine(cut01, 0, cut01, 1)
    l2.DrawLine(cut05, 0, cut05, 1)

    c.Print(outPlot)
    return cut99,cut95,cut01,cut05,cut90


if __name__ == "__main__":
    main(sys.argv[1:])
