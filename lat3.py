#!/usr/bin/env python
"""
===================== LAT3.py =====================


Tunes cut parameters in LAT, calculates the 1%, 5%, 90%, 95%, and 99%

Usage Examples:
    Tuning all cuts:
        ./lat3.py -tune /path/to/calib/files -db -s DS subDS Module -all
    Specific cut and/or channel:
        ./lat3.py -tune /path/to/calib/files -db -s DS subDS Module -ch channel -bcMax
    Custom cut:
        ./lat3.py -tune /path/to/calib/files -db -s DS subDS Module -Custom "bcMax/bcMin"


    Applying cuts:
        ./lat3.py -cut /path/to/bkg/files /path/to/output/files -db -s DS subDS Module


v1: 03 Oct 2017

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
from scipy.stats import mode

def main(argv):

    print "========================================"
    print "LAT3 started:",time.strftime('%X %x %Z')
    startT = time.clock()

    cInfo = CalInfo()
    pathToInput, pathToOutput = ".", "."
    dsNum, subNumm, modNum, chNum = -1, -1, -1, -1
    skimTree = ROOT.TChain("skimTree")
    customPar = ""
    calList, parList, parNameList, chList = [], [], [], []
    fTune, fCut, fFor, fUpd, fastMode, fDB, fCSV = False, False, False, False, False, False, False

    if len(argv) == 0:
        return
    for i, opt in enumerate(argv):
        # -- Cut tuning options --
        if opt == "-all":
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
            print "Processing DS-%d subDS-%d Module-%d"%(dsNum, subNum, modNum)
        if opt == "-d":
            pathToInput, pathToOutput = argv[i+1], argv[i+2]
            print "Custom paths: Input %s" % (pathToInput)
        if opt == "-ch":
            chNum = int(argv[i+1])
            print "Tuning specific channel %d" % (chNum)
        if opt == "-fast":
            fastMode == True
            print "Tuning cuts with fastMode ON, diagnostic plots will not be generated!"
        if opt == "-pd":
            import pandas as pd
            dfList, fCSV = [], True
            print "Saving CSV file in ./output"
        if opt == "-db":
            fDB = True
            print "DB mode"
        # -- Database options --
        #TODO -- Flesh out cut database options
        if opt == "-force":
            fFor = True
            print "Force DB update mode."
        if opt == "-tune":
            fTune = True
            pathToInput = argv[i+1]
            print "Cut tune mode."
        if opt == "-cut":
            pathToInput, pathToOutput = argv[i+1], argv[i+2]
            fCut = True
            print "Cut application mode."
        if opt == "-upd":
            fUpd = True
            print "File update mode."

    # -- Load calibration files --
    if dsNum == -1 or subNum == -1 or modNum == -1:
        print "DS, subDS, or module number not set properly, exiting"
        return
    elif fTune:
        # Limit to 10 calibration runs because that's all Clint processed!
        calList = cInfo.GetCalList("ds%d_m%d"%(dsNum, modNum), subNum, runLimit=10)
        for i in calList: skimTree.Add("%s/latSkimDS%d_run%d_*"%(pathToInput, dsNum, i))
    elif fCut:
        # Add Background Data
        skimTree.Add("%s/latSkimDS%d_%d_0.root"%(pathToInput, dsNum, subNum))
    else:
        print "Tune or Cut option not set"
        return

    # -- Load chains for this DS --
    inPath = pathToInput + "/latSkimDS%d*.root" % dsNum
    fileList = glob.glob(inPath)
    cutFile = TFile(fileList[0])
    theCut = cutFile.Get("theCut").GetTitle()

    # -- Load channel list --
    if chNum == -1:
        chList = ds.GetGoodChanList(dsNum)
        if dsNum==5 and modNum == 1: # remove 692 and 1232
            chList = [584, 592, 598, 608, 610, 614, 624, 626, 628, 632, 640, 648, 658, 660, 662, 672, 678, 680, 688, 690, 694]
        if dsNum==5 and modNum == 2:
            chList = [1106, 1110, 1120, 1124, 1128, 1170, 1172, 1174, 1176, 1204, 1208, 1298, 1302, 1330, 1332]
    else:
        chList = [chNum]

    print "Processing channels: ", chList
    print "Processing runs: ", calList

    # -- Tune cuts --
    if fTune:
        tunedPars = {}
        tuneRange = [[236, 240], [5, 50]]
        tuneNames = ["Peak", "Continuum"]
        for par, parName in zip(parList, parNameList):
            for tRange, tName in zip(tuneRange, tuneNames):
                cutDict = TuneCut(dsNum, subNum, tRange[0], tRange[1], tName, skimTree, chList, par, parName, theCut, fastMode)
                key = "%s_ds%d_idx%d_m%d_%s"%(parName,dsNum,subNum,modNum,tName)
                # print key, cutDict
                if fDB:
                    wl.setDBCalRecord({"key":key,"vals":cutDict})

                if fCSV:
                    dummyDict = {"DS":[dsNum]*5, "SubDS":[subNum]*5, "Module":[modNum]*5, "Cut":[parName]*5, "Range":[tName]*5, "Percentage":[1, 5, 90, 95, 99]}
                    dummyDict2 = dict(dummyDict.items() + cutDict.items())
                dfList.append(pd.DataFrame(dummyDict2))
        if fCSV:
            dfTot = pd.concat(dfList)
            dfTot.to_csv("./output/Cuts_ds%d_idx%d_m%d.csv"%(dsNum,subNum,modNum))

    # -- Apply cuts --
    if fCut:
        megaCut = MakeCutList(cInfo, skimTree, theCut, dsNum, modNum, chList)
        print megaCut
        for idx, ch in enumerate(chList):
            print theCut+'&&'+megaCut[ch][2:]
            outFile = TFile(pathToOutput+"latCutSkimDS%d_%d_ch%d.root"%(dsNum,subNum,ch),"RECREATE")
            outTree = TTree()
            outTree = skimTree.CopyTree(theCut+'&&'+megaCut[ch][2:])
            outTree.Write()
            # cutUsed = ROOT.TNamed("theCut_ch%d_idx%d"%(ch,subNum),theCut+'&&'+megaCut[ch][2:])
            cutUsed = ROOT.TNamed("theCut",theCut+'&&'+megaCut[ch][2:])
            cutUsed.Write()
            outFile.Close()

    stopT = time.clock()
    print "Stopped:",time.strftime('%X %x %Z'),"\nProcess time (min):",(stopT - startT)/60


def TuneCut(dsNum, subNum, tMin, tMax, tName, cal, chList, par, parName, theCut, fastMode):
    c = TCanvas("%s"%(parName),"%s"%(parName),1600,600)
    c.Divide(3,1,0.00001,0.00001)
    cutDict = {}
    for ch in chList:
        cutDict[ch] = [0,0,0,0,0]
        eb, elo, ehi = (tMax-tMin),tMin,tMax
        d1Cut = theCut + " && trapENFCalC > %d && trapENFCalC < %d && channel==%d" % (elo,ehi,ch)
        d2Cut = theCut + " && channel==%d" % ch
        nPass = cal.Draw("trapENFCalC:%s"%(par), d1Cut, "goff")
        nEnergy = cal.GetV1()
        nCut = cal.GetV2()
        nCutList = list(float(nCut[n]) for n in xrange(nPass))
        nEnergyList = list(float(nEnergy[n]) for n in xrange(nPass))
        # Error and warning messages
        if len(nCutList) == 0 or len(nEnergyList) == 0:
            print "Error: Channel %d has no entries, cut cannot be set properly, setting to [0,0,0,0,0,0,0]"%(ch)
            cutDict[ch] = [0,0,0,0,0]
            continue
        if len(nCutList) <= 1000 or len(nEnergyList) <= 1000:
            print "Warning: Channel %d has less than 1000 entries, cut values may not be accurate"%(ch)

        vb, v5, v95 = 100000, np.percentile(nCutList, 5), np.percentile(nCutList,95)
        vlo, vhi = v5-5*abs(v5), v95+5*abs(v95)
        nCutListReduced = [x for x in nCutList if x > v5 and x < v95]
        outPlot = "./plots/tuneCuts/%s_ds%d_idx%d_%s_ch%d.png" % (parName,dsNum,subNum,tName,ch)
        # cutMode, cutMedian = mode(np.round(nCutListReduced))[0][0], np.median(nCutListReduced)
        cut99,cut95,cut01,cut05,cut90 = MakeCutPlot(c,cal,par,eb,elo,ehi,vb,vlo,vhi,d2Cut,d1Cut,outPlot,fastMode)
        cutDict[ch] = [cut01,cut05,cut90,cut95,cut99]
    return cutDict

def MakeCutPlot(c,cal,var,eb,elo,ehi,vb,vlo,vhi,d2Cut,d1Cut,outPlot,fastMode):
    """ Repeated code is the DEVIL.  Even if you have to pass in 1,000,000 arguments. """

    # Calculate cut vals (assumes plot range is correct)
    h1 = wl.H1D(cal,vb,vlo,vhi,var,d1Cut)
    h1Sum = h1.Integral()
    if h1Sum == 0:
        print "Error: Failed %s, histogram sum is 0 so cannot normalize, setting to [0,0,0,0,0]"%(var)
        return 0,0,0,0,0
    h1.Scale(1/h1Sum)
    try:
        cut99,cut95,cut01,cut05,cut90 = wl.GetIntegralPoints(h1)
    except:
        print "Error: Failed %s using cut %s, setting to [0,0,0,0,0]"%(var,d1Cut)
        return 0,0,0,0,0
    if fastMode:
        return cut99,cut95,cut01,cut05,cut90

    # Generate the plot for inspection.
    c.cd(2)
    gPad.SetLogy(0)
    h1.GetXaxis().SetRangeUser(cut01-abs(0.25*cut01), cut99 + abs(0.25*cut99) )
    h1.SetTitle("")
    h1.GetXaxis().SetTitle(var)
    h1.Draw("hist")

    c.cd(1)
    gPad.SetLogy(0)
    cal.Draw("%s:trapENFCalC>>b(%d,%d,%d,%d,%.3E,%.3E)"%(var,eb+10,elo-5,ehi+5,vb,cut01-abs(0.25*cut01),cut99+abs(0.25*cut99)) ,d2Cut)

    l1, l2, l3 = TLine(), TLine(), TLine()
    l1.SetLineColor(ROOT.kGreen)
    l2.SetLineColor(ROOT.kRed)
    l3.SetLineColor(ROOT.kMagenta)

    l1.DrawLine(elo-5, cut99, ehi+5, cut99)
    l2.DrawLine(elo-5, cut95, ehi+5, cut95)
    l2.DrawLine(elo-5, cut05, ehi+5, cut05)
    l1.DrawLine(elo-5, cut01, ehi+5, cut01)

    c.cd(3)
    x_h1, y_h1 = wl.npTH1D(h1)
    int_h1 = wl.integFunc(y_h1)
    g2 = TGraph(len(x_h1), x_h1, int_h1)
    g2.GetXaxis().SetRangeUser(cut01-abs(0.3*cut01), cut99 + abs(0.3*cut99) )
    g2.SetTitle("")
    g2.GetXaxis().SetTitle(var)
    g2.GetYaxis().SetTitle("Percentile")
    g2.Draw("ACP")
    l1.DrawLine(cut99, 0, cut99, 1)
    l2.DrawLine(cut95, 0, cut95, 1)
    l1.DrawLine(cut01, 0, cut01, 1)
    l2.DrawLine(cut05, 0, cut05, 1)

    c.Print(outPlot)
    return cut99,cut95,cut01,cut05,cut90

def MakeCutList(cInfo, skimTree, basicCut, dsNum, modNum, chList=[], mode='db'):
    """ Pass in background run and it generates a dictionary of cuts for all good channels"""
    nPass = skimTree.Draw("run", basicCut, "goff")
    nRun = skimTree.GetV1()
    runList = list(set(int(nRun[n]) for n in xrange(nPass)))
    print "Processing Runs", runList
    idxMin = cInfo.GetCalIdx("ds%d_m%d"%(dsNum, modNum), runList[0])
    idxMax = cInfo.GetCalIdx("ds%d_m%d"%(dsNum, modNum), runList[-1])
    megaCut = {}
    if mode == 'db':
        for subNum in range(idxMin,idxMax+1):
            fsD = wl.getDBCalRecord("fitSlo_ds%d_idx%d_m%d_Peak"%(dsNum,subNum,modNum))
            rnD = wl.getDBCalRecord("riseNoise_ds%d_idx%d_m%d_Peak"%(dsNum,subNum,modNum))
            bcD = wl.getDBCalRecord("bcMax_ds%d_idx%d_m%d_Peak"%(dsNum,subNum,modNum))
            runMin, runMax = cInfo.master['ds%d_m%d'%(dsNum,modNum)][subNum][1], cInfo.master['ds%d_m%d'%(dsNum,modNum)][subNum][2]
            for ch in chList:
                if ch in megaCut.keys():
                    megaCut[ch] += '||(run>=%d&&run<=%d&&fitSlo<%.2f&&riseNoise<%.2f)'%(runMin,runMax,fsD[ch][2], rnD[ch][2])
                else:
                    megaCut[ch] = '||(run>=%d&&run<=%d&&fitSlo<%.2f&&riseNoise<%.2f)'%(runMin,runMax,fsD[ch][2], rnD[ch][2])
    elif mode == 'csv':
        print('Not implemented!')
        return 0

    return megaCut

if __name__ == "__main__":
    main(sys.argv[1:])
