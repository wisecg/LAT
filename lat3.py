#!/usr/bin/env python
"""
===================== LAT3.py =====================
Tunes cut parameters in LAT, calculates the 1%, 5%, 90%, 95%, and 99%.
Also applies cuts to create channel-specific files.
Usage Examples:
    Tuning all cuts:
        ./lat3.py -tune /path/to/calib/files -db -s DS subDS Module -all -Range "Min_Max"
    Specific cut and/or channel:
        ./lat3.py -tune /path/to/calib/files -db -s DS subDS Module -ch channel -bcMax
    Custom cut:
        ./lat3.py -tune /path/to/calib/files -db -s DS subDS Module -Custom "bcMax/bcMin"
    Applying cuts:
        ./lat3.py -cut [dsNum] [cutType] [dataType]
v1: 03 Oct 2017
========= C. Wiseman (USC), B. Zhu (LANL) =========
"""
import sys, time, ROOT, glob
sys.argv += [ '-b' ] # force ROOT to be loaded in batch mode.
from ROOT import gROOT, gStyle, gPad
from ROOT import TFile, TTree, TChain, TCanvas, TH1D, TH2D, TF1, TLegend, TLine, TGraph, TNamed
import numpy as np
from DataSetInfo import CalInfo
import DataSetInfo as ds
import waveLibs as wl
from scipy.stats import mode
from scipy.optimize import curve_fit
from matplotlib import pyplot as plt
import seaborn as sns
sns.set(style='whitegrid', context='talk')

def main(argv):

    print "========================================"
    print "LAT3 started:",time.strftime('%X %x %Z')
    startT = time.clock()

    cInfo = CalInfo()
    pathToInput, pathToOutput = ".", "."
    dsNum, subNumm, modNum, chNum = -1, -1, -1, -1
    skimTree = ROOT.TChain("skimTree")
    customPar = ""
    tuneNames, calList, parList, parNameList, chList = [], [], [], [], []
    fTune, fFor, fastMode, fDB, fCSV = False, False, False, False, False

    if len(argv) == 0:
        return
    for i, opt in enumerate(argv):

        # -- Cut tuning options --
        if opt == "-all":
            parList.append('bcMax'), parNameList.append('bcMax')
            parList.append('fitSlo'), parNameList.append('fitSlo')
            parList.append('riseNoise'), parNameList.append('riseNoise')
            print "Tuning all cuts"
        if opt == "-bcMax":
            parList.append('bcMax'), parNameList.append('bcMax')
            print "Tuning bcMax"
        if opt == "-noiseWeight":
            parList.append('(waveS4-waveS1)/bcMax/trapENFCal'), parNameList.append('noiseWeight')
            print "Tuning noiseWeight ((waveS4-waveS1)/bcMax/trapENFCal)"
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

        # Tune specific range -- input as string with comma separation:
        if opt == "-Range":
            rangeName = str(argv[i+1])
            tuneNames = rangeName.split(',')
            print "Tuning ranges: ", tuneNames

        # -- Input/output options --
        if opt == "-s":
            dsNum, subNum, modNum = int(argv[i+1]), int(argv[i+2]), int(argv[i+3])
            print "Processing DS-%d subDS-%d Module-%d" % (dsNum, subNum, modNum)
        if opt == "-d":
            pathToInput, pathToOutput = argv[i+1], argv[i+2]
            print "Custom paths: Input %s" % pathToInput
        if opt == "-ch":
            chNum = int(argv[i+1])
            print "Tuning specific channel %d" % (chNum)
        if opt == "-fast":
            fastMode = True
            print "Tuning cuts with fastMode ON, diagnostic plots will not be generated!"
        if opt == "-pd":
            import pandas as pd
            dfList, fCSV = [], True
            print "Saving CSV file in ./output"
        if opt == "-db":
            fDB = True
            print "DB mode"

        # -- Database options --
        if opt == "-force":
            fFor = True
            print "Force DB update mode."
        if opt == "-tune":
            fTune = True
            pathToInput = argv[i+1]
            print ("Cut tune mode. Input path for cal files: ", pathToInput)

        # -- Apply channel cuts (exits after running this) --
        if opt == "-cut":
            dsNum = int(argv[i+1])
            cutType = argv[i+2]
            dType = argv[i+3]
            print "Applying %s cut to DS-%d (%s)..." % (cutType, dsNum, dType)
            ApplyChannelCuts(dsNum, cutType, dType)
            stopT = time.clock()
            print "Stopped:",time.strftime('%X %x %Z'),"\nProcess time (min):",(stopT - startT)/60
            return

    # -- Load calibration files --
    if dsNum == -1 or subNum == -1 or modNum == -1:
        print "DS, subDS, or module number not set properly, exiting"
        return
    elif fTune:
        # Limit to 10 calibration runs because that's all Clint processed!  What a jerk.
        calList = cInfo.GetCalList("ds%d_m%d" % (dsNum, modNum), subNum, runLimit=10)
        for i in calList: skimTree.Add("%s/latSkimDS%d_run%d_*" % (pathToInput, dsNum, i))
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
        if dsNum==5 and modNum == 1: # remove 692 and 1232 (both beges, so who cares)
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
        # Default is peak only
        if not tuneNames: tuneNames.append("Peak")
        for par, parName in zip(parList, parNameList):
            for idx, tName in enumerate(tuneNames):
                tRange = []
                if tName == "Continuum": tRange = [5, 50]
                elif tName == "Peak": tRange = [236, 240]
                elif tName == "SoftPlus": tRange = [] # Maybe should be another boolean
                else: tRange = [int(tName.split("_")[0]), int(tName.split("_")[1])]

                key = "%s_ds%d_idx%d_m%d_%s"%(parName,dsNum,subNum,modNum,tName)
                if not tRange:
                    cutDict = TuneSoftPlus(dsNum, subNum, tName, skimTree, chList, par, parName, theCut, fastMode)
                else:
                    cutDict = TuneCut(dsNum, subNum, tRange[0], tRange[1], tName, skimTree, chList, par, parName, theCut, fastMode)

                if fDB: wl.setDBCalRecord({"key":key,"vals":cutDict}, forceUpdate=fFor)
                if fCSV:
                    dummyDict = {"DS":[dsNum]*5, "SubDS":[subNum]*5, "Module":[modNum]*5, "Cut":[parName]*5, "Range":[tName]*5, "Percentage":[1, 5, 90, 95, 99]}
                    dummyDict2 = dict(dummyDict.items() + cutDict.items())
                    dfList.append(pd.DataFrame(dummyDict2))
        if fCSV:
            dfTot = pd.concat(dfList)
            dfTot.to_csv("./output/Cuts_ds%d_idx%d_m%d.csv"%(dsNum,subNum,modNum))

    stopT = time.clock()
    print "Stopped:",time.strftime('%X %x %Z'),"\nProcess time (min):",(stopT - startT)/60


def TuneCut(dsNum, subNum, tMin, tMax, tName, cal, chList, par, parName, theCut, fastMode):
    c = TCanvas("%s"%(parName),"%s"%(parName),1600,600)
    c.Divide(3,1,0.00001,0.00001)
    cutDict = {}
    for ch in chList:
        cutDict[ch] = [0,0,0,0,0]
        eb, elo, ehi = (tMax-tMin),tMin,tMax
        d1Cut = theCut + " && trapENFCal > %d && trapENFCal < %d && channel==%d" % (elo,ehi,ch)
        d2Cut = theCut + " && channel==%d" % ch
        nPass = cal.Draw("trapENFCal:%s"%(par), d1Cut, "goff")
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
        cut99,cut95,cut01,cut05,cut90 = MakeCutPlot(c,cal,par,eb,elo,ehi,vb,vlo,vhi,d2Cut,d1Cut,outPlot,fastMode)
        cutDict[ch] = [cut01,cut05,cut90,cut95,cut99]
    return cutDict


def TuneSoftPlus(dsNum, subNum, tName, cal, chList, par, parName, theCut, fastMode):
    nTest = np.linspace(0, 250, 1000)
    cutDict = {}
    for ch in chList:
        cutDict[ch] = [0,0,0,0]
        nPass1 = cal.Draw("trapENFCal:%s"%(par),  theCut + "&& trapENFCal>5 && trapENFCal<50 && channel==%d"%(ch), "goff")
        nCutArray1, nCutArray2, nCutArray3 = [], [], []
        if nPass1 != 0:
            nEnergy1 = cal.GetV1()
            nCut1 = cal.GetV2()
            nCutList1 = list(float(nCut1[n]) for n in xrange(nPass1))
            nEnergyList1 = list(float(nEnergy1[n]) for n in xrange(nPass1))
            nCutArray1 = [[x,y] for x,y in zip(nCutList1, nEnergyList1) if x > np.percentile(nCutList1, 5) and x < np.percentile(nCutList1, 85)]
        nPass2 = cal.Draw("trapENFCal:%s"%(par),  theCut+"&& trapENFCal>50 && trapENFCal<150 && channel==%d"%(ch), "goff")
        if nPass2 != 0:
            nEnergy2 = cal.GetV1()
            nCut2 = cal.GetV2()
            nCutList2 = list(float(nCut2[n]) for n in xrange(nPass2))
            nEnergyList2 = list(float(nEnergy2[n]) for n in xrange(nPass2))
            nCutArray2 = [[x,y] for x,y in zip(nCutList2, nEnergyList2) if x > np.percentile(nCutList2, 5) and x < np.percentile(nCutList2, 90)]
        nPass3 = cal.Draw("trapENFCal:%s"%(par),  theCut+"&& trapENFCal>150 && trapENFCal<240 && channel==%d"%(ch), "goff")
        if nPass3 != 0:
            nEnergy3 = cal.GetV1()
            nCut3 = cal.GetV2()
            nCutList3 = list(float(nCut3[n]) for n in xrange(nPass3))
            nEnergyList3 = list(float(nEnergy3[n]) for n in xrange(nPass3))
            nCutArray3 = [[x,y] for x,y in zip(nCutList3, nEnergyList3) if x > np.percentile(nCutList3, 5) and x < np.percentile(nCutList3, 90)]
        nCutArray = np.asarray(nCutArray1 + nCutArray2 + nCutArray3)

        if len(nCutArray) == 0:
            print "No events (setting to [0,0,0,0]), skipping channel ", ch
            continue
        popt,_ = curve_fit(softplus, nCutArray[:,1], nCutArray[:,0], p0=[1., 0.005, 10, 0.5], bounds = ((-10,0,-10,0),(10,10,150,100)))

        cutDict[ch] = [popt[0], popt[1], popt[2], popt[3]]
        if fastMode: continue

        # Draw plots
        yFit = softplus(nTest, *popt)
        g2 = sns.JointGrid(x=nCutArray[:,1], y=nCutArray[:,0], size=10, space=0.2)
        g2.plot_joint(sns.kdeplot)
        plt.plot(nTest, yFit, "-", color='red')
        plt.ylabel('%s'%(parName))
        plt.xlabel('Energy')
        plt.title('Y-Offset: %.2f  Slope: %.3f  X-Shift: %.2f  Curvature: %.2f'%(popt[0], popt[1], popt[2], popt[3]))
        g2.ax_marg_y.set_axis_off()
        _ = g2.ax_marg_x.hist(np.array(nEnergyList1 + nEnergyList2 + nEnergyList3), alpha=0.8, bins=np.linspace(0, 250, 500))
        g2.savefig("./plots/tuneCuts/%s_ds%d_idx%d_%s_ch%d.png" % (parName,dsNum,subNum,tName,ch))
    return cutDict


def softplus(x, a, b, c, d):
    """ A = Y-offset, B = Slope , C = X shift for flatness, D = Curvature """
    return a + b*np.log(1+np.exp((x - c)/d) )


def MakeCutPlot(c,cal,var,eb,elo,ehi,vb,vlo,vhi,d2Cut,d1Cut,outPlot,fastMode):
    """ Creates a channel-specific energy calibration plot. """

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
        print "Returning fastMode output: ", cut99,cut95,cut01,cut05,cut90
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
    cal.Draw("%s:trapENFCal>>b(%d,%d,%d,%d,%.3E,%.3E)"%(var,eb+10,elo-5,ehi+5,vb,cut01-abs(0.25*cut01),cut99+abs(0.25*cut99)) ,d2Cut)

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


def ApplyChannelCuts(dsNum, cutType, dType):
    """ ./lat3.py -cut [dsNum] [cutType]
    This runs over whole datasets.
    cutTypes:
        fs, rn, wf, fs+rn, fs+wf, rn+wf, fs+rn+wf

    dTypes:
        bkg, cal
    """

    # setup a loop over modules and dataset ranges
    gROOT.ProcessLine("gErrorIgnoreLevel = 3001;")
    cInfo = ds.CalInfo()
    nMods = [1]
    if dsNum == 4: nMods = [2]
    if dsNum == 5: nMods = [1,2]
    # Dummy string for file writing -- adds nothing to the directories if background
    dString = ""
    if dType == "cal": dString = "cal"

    for modNum in nMods:
        # Changed so range of idx are set here to take advantage of module number
        if dType == "bkg":
            nRanges = [0, ds.dsMap[dsNum]]
            # if dsNum==5: nRanges[0] = 80 # exclude DS-5A
        elif dType == "cal":
            nRanges = [0, len(cInfo.master['ds%d_m%d'%(dsNum, modNum)])]

        # Loop over bkgIdx, even though for calibration runs this will represent calIdx
        for bkgIdx in range(nRanges[0], nRanges[1]+1):

            # load the chains and find the right calIdx's.
            skimTree = TChain("skimTree")

            # build the file list
            fRegex = ""
            if dType == "bkg":
                fRegex = "/global/homes/w/wisecg/project/bg-lat/latSkimDS%d_%d_*.root" % (dsNum, bkgIdx)
                fList = glob.glob(fRegex)
                skimTree.Add(fRegex)
            elif dType == "cal":
                calList = cInfo.GetCalList("ds%d_m%d" % (dsNum, modNum), bkgIdx, runLimit=10)
                fList = []
                for i in calList:
                    fList += glob.glob("/global/homes/w/wisecg/project/cal-lat/latSkimDS%d_run%d_*.root"%(dsNum,i))
                    skimTree.Add("/global/homes/w/wisecg/project/cal-lat/latSkimDS%d_run%d_*.root" % (dsNum, i))
            file0 = fList[0]
            print "DS-%d subset %d, Mod-%d.  N_files: %d" % (dsNum, bkgIdx, modNum, len(fList))

            # Print some basic info about files
            f = TFile(file0)
            theCut = f.Get("theCut").GetTitle()
            if dType == "bkg":
                skimTree.GetEntry(0)
                firstRun = skimTree.run
                skimTree.GetEntry(skimTree.GetEntries()-1)
                lastRun = skimTree.run
                calIdxLo = cInfo.GetCalIdx("ds%d_m%d" % (dsNum, modNum), firstRun)
                calIdxHi = cInfo.GetCalIdx("ds%d_m%d" % (dsNum, modNum), lastRun)
            elif dType == "cal":
                # All the idx are the same for calibration!
                calIdxLo = calIdxHi = bkgIdx
                firstRun, lastRun = calList[0], calList[-1]

            print "    Entries %d  firstRun %d  lastRun %d  calIdxLo %d  calIdxHi %d" % (skimTree.GetEntries(),firstRun,lastRun,calIdxLo,calIdxHi)

            # build the channel list  (remove 692 and 1232 from DS5 for now.)
            chList = ds.GetGoodChanList(dsNum)
            if dsNum==5 and modNum==1:
                chList = [ch for ch in chList if ch < 1000 and ch!=692]
            if dsNum==5 and modNum==2:
                chList = [ch for ch in chList if ch > 1000 and ch!=1232]

            # -- create a dict of cuts for each channel, covering each calIdx. --
            cutDict = {}
            for calIdx in range(calIdxLo, calIdxHi+1):

                runCovMin = cInfo.master["ds%d_m%d" % (dsNum, modNum)][calIdx][1]
                runCovMax = cInfo.master["ds%d_m%d" % (dsNum, modNum)][calIdx][2]
                runCut = "run>=%d && run<=%d" % (runCovMin, runCovMax)

                fsD = wl.getDBCalRecord("fitSlo_ds%d_idx%d_m%d_Peak" % (dsNum,calIdx,modNum))
                rnCD = wl.getDBCalRecord("riseNoise_ds%d_idx%d_m%d_Continuum" % (dsNum,calIdx,modNum))
                rnSD = wl.getDBCalRecord("riseNoise_ds%d_idx%d_m%d_SoftPlus" % (dsNum,calIdx,modNum))
                wfD = wl.getDBCalRecord("wfstd_ds%d_idx%d_mod%d" % (dsNum, calIdx, modNum)) # returns 0 for ranges where we have no data

                for ch in chList:

                    fsCut, rnCut, wfCut, chanCut = None, None, None, None

                    # print ch,":",fsD[ch][2]

                    # fitSlo: check the 90% value is positive
                    if fsD[ch][2] > 0:
                        fsCut = "fitSlo<%.2f" % fsD[ch][2]

                    # riseNoise: check the softplus curvature is positive
                    if rnSD[ch][3] > 0:
                        rnCut = "riseNoise<(%.3f+%.5f*TMath::Log(1+TMath::Exp((trapENFCal-(%.3f))/%.3f)))" % (max(rnSD[ch][0],rnCD[ch][4]), rnSD[ch][1], rnSD[ch][2], rnSD[ch][3])

                    # wfStd: check if ralph says this is ok to use
                    if wfD!=0 and ch in wfD.keys() and wfD[ch][0]==u'y':
                        wfCut = "abs(wfstd - sqrt((%.4e + %.4e*trapENFCal + %.4e*trapENFCal**2 + %.2e*pow(trapENFCal,3) + %.2e*pow(trapENFCal,4))**2 + %.4f)) < (%.2f+%.3f*trapENFCal)" % (wfD[ch][3], wfD[ch][4], wfD[ch][5], wfD[ch][6], wfD[ch][7], wfD[ch][8], wfD[ch][9], wfD[ch][10])

                    # set the combination channel cut
                    if cutType == "fs" and fsCut!=None:
                        chanCut = "(%s && %s)" % (runCut, fsCut)

                    if cutType == "rn" and rnCut!=None:
                        chanCut = "(%s && %s)" % (runCut, rnCut)

                    if cutType == "wf" and wfCut!=None:
                        chanCut = "(%s && %s)" % (runCut, wfCut)

                    if cutType == "fs+rn" and fsCut!=None and rnCut!=None:
                        chanCut = "(%s && %s && %s)" % (runCut, fsCut, rnCut)

                    if cutType == "fs+wf" and fsCut!=None and wfCut!=None:
                        chanCut = "(%s && %s && %s)" % (runCut, fsCut, wfCut)

                    if cutType == "rn+wf" and rnCut!=None and wfCut!=None:
                        chanCut = "(%s && %s && %s)" % (runCut, rnCut, wfCut)

                    if cutType == "fs+rn+wf" and fsCut!=None and rnCut!=None and wfCut!=None:
                        chanCut = "(%s && %s && %s && %s)" % (runCut, fsCut, rnCut, wfCut)

                    # create dict entry for this channel or append to existing, taking care of parentheses and OR's.
                    if ch in cutDict.keys() and chanCut!=None:
                        cutDict[ch] += " || %s" % chanCut
                    elif ch not in cutDict.keys() and chanCut!=None:
                        cutDict[ch] = "(%s" % chanCut

            # close the parens for each channel entry
            for key in cutDict:
                cutDict[key] += ")"

            # -- finally, loop over each channel we have an entry for, get its cut, and create an output file. --
            for ch in cutDict:
                # TODO: threshold cut (or at least save the value for each bkgIdx)

                chanCut = theCut + "&& gain==0 && channel==%d" % ch

                if cutType == "fs":
                    outFile = "~/project/cuts/%sfs/%sfitSlo-DS%d-%d-ch%d.root" % (dString, dString, dsNum, bkgIdx, ch)
                    chanCut += "&& fitSlo>0 && %s" % cutDict[ch]

                if cutType == "rn":
                    outFile = "~/project/cuts/%srn/%sriseNoise-DS%d-%d-ch%d.root" % (dString, dString, dsNum, bkgIdx, ch)
                    chanCut += "&& %s" % cutDict[ch]

                if cutType == "wf":
                    outFile = "~/project/cuts/%swf/%swfstd-DS%d-%d-ch%d.root" % (dString, dString, dsNum, bkgIdx, ch)
                    chanCut += "&& %s" % cutDict[ch]

                if cutType == "fs+rn":
                    outFile = "~/project/cuts/%sfs_rn/%sfs_rn-DS%d-%d-ch%d.root" % (dString, dString, dsNum, bkgIdx, ch)
                    chanCut += "&& fitSlo>0 && %s" % cutDict[ch]

                if cutType == "fs+wf":
                    outFile = "~/project/cuts/%sfs_wf/%sfs_wf-DS%d-%d-ch%d.root" % (dString, dString, dsNum, bkgIdx, ch)
                    chanCut += "&& fitSlo>0 && %s" % cutDict[ch]

                if cutType == "rn+wf":
                    outFile = "~/project/cuts/%srn_wf/%srn_wf-DS%d-%d-ch%d.root" % (dString, dString, dsNum, bkgIdx, ch)
                    chanCut += "&& %s" % cutDict[ch]

                if cutType == "fs+rn+wf":
                    outFile = "~/project/cuts/%sfs_rn_wf/%sfs_rn_wf-DS%d-%d-ch%d.root" % (dString, dString, dsNum, bkgIdx, ch)
                    chanCut += "&& fitSlo>0 && %s" % cutDict[ch]

                print "    Writing to:",outFile
                print "    Cut used:",chanCut,"\n"
                outFile = TFile(outFile,"RECREATE")
                outTree = TTree()
                outTree = skimTree.CopyTree(chanCut)
                outTree.Write()
                cutUsed = TNamed("chanCut",chanCut)
                cutUsed.Write()
                print "Wrote",outTree.GetEntries(),"entries."
                outFile.Close()
                # return


if __name__ == "__main__":
    main(sys.argv[1:])
