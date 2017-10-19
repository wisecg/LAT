#!/usr/bin/env python
import sys, time, os, glob
import ROOT
import numpy as np
import pandas as pd
import DataSetInfo as ds
import waveLibs as wl

"""
    This script is for simple Poisson counting of regions after successive cuts
    Saves the data to a CSV file
    Can also draw spectra with cuts -- but it's suuuuper slow

"""

def main(argv):

    inDir, cutDir, outDir = ".", ".", "./plots"
    specMode = False
    dsNum, modNum, chNum = -1, -1, -1
    skimTree = ROOT.TChain("skimTree")
    bins,lower,upper = 1250,0,250

    for i,opt in enumerate(argv):
        if opt == "-spec":
            specMode = True
            print "Spectrum Mode"
        if opt == "-d":
            inDir, outDir = argv[i+1], argv[i+2]
        if opt == "-ch":
            chNum = int(argv[i+1])
            print ("Drawing specific channel %d" % (chNum))
        if opt == "-s":
            dsNum, modNum = int(argv[i+1]), int(argv[i+2])
            print ("Drawing DS-%d Module-%d"%(dsNum, modNum))
        if opt == "-cal":
            calMode = True
            print ("Drawing Calibration runs")

        if opt == "-db":
            print "Loading from DB"
        if opt == "-csv":
            cutDir = argv[i+1]
            print "Loading CSV"

    cInfo = ds.CalInfo()
    EnergyList = [[1.,5.], [2., 4.], [4., 9.], [9., 12.], [12., 40.], [40., 50.], [50., 100.]]

    # -- Load channel list --
    if chNum == -1:
        chList = ds.GetGoodChanList(dsNum)
        if dsNum==5 and modNum == 1: # remove 692 and 1232
            chList = [584, 592, 598, 608, 610, 614, 624, 626, 628, 632, 640, 648, 658, 660, 662, 672, 678, 680, 688, 690, 694]
        if dsNum==5 and modNum == 2:
            chList = [1106, 1110, 1120, 1124, 1128, 1170, 1172, 1174, 1176, 1204, 1208, 1298, 1302, 1330, 1332]
    else:
        chList = [chNum]

    # -- Load calibration files --
    if dsNum == -1 or modNum == -1:
        print "DS, subDS, or module number not set properly, exiting"
        return
    else:
        # Limit to 10 calibration runs because that's all Clint processed!
        for subNum in cInfo.master["ds%d_m%d"%(dsNum,modNum)].keys():
            calList = cInfo.GetCalList("ds%d_m%d"%(dsNum, modNum), subNum, runLimit=10)
            for i in calList: skimTree.Add("%s/latSkimDS%d_run%d_*"%(inDir, dsNum, i))

    # -- Load cut files --
    cutFiles = glob.glob('%s/Cuts_ds%d_*_m%d.csv'%(cutDir, dsNum, modNum))
    cutdfList = []
    for f in cutFiles:
        try:
            cutdfList.append(pd.read_csv(f))
        except:
            print (f, "is broken, SAD")
    dfTot = pd.concat(cutdfList)

    # List of histograms (for summing together) and Dictionary of histograms (for different cuts)
    hList, hDict = [], {}

    # Create a list of DataFrames for each channel to concatenate at the end
    # cutNames = ["BasicCut", "+tailSlope", "+riseNoise", "+fitSlo"]
    cutNames = ["BasicCut", "+riseNoise", "+fitSlo"]
    # cutNames = ["BasicCut", "+tailSlope", "+bcMax", "+fitSlo"]
    dfList = []

    for idx,ch in enumerate(chList):
        # Create new key for dictionaries according to channel
        hDict[ch] = []
        # Create dummy list and series to store into DataFrame later
        cutMatrix = {}

        channelCut = "channel==%d" % (ch)

        # Get threshold info
        # goodRuns,badRuns,goodRunSigmas = ds.GetThreshDicts(dsNum)
        # threshrunCut = "&&("
        # for idx2,runRange in enumerate(goodRuns[ch]):
        #     threshrunCut += "(run>=%d&&run<=%d)||" % (runRange[0],runRange[1])
        #     # Strip all spaces to save space
        #     totalCut = megaCut.replace(" ", "") + threshrunCut[:-2] + ")"

        for subNum in cInfo.master["ds%d_m%d"%(dsNum,modNum)].keys():
            # Create new array for each subDS
            hDict[ch].append([])

            # Create separate dataframe for each subNum -- excess?
            dfSub = dfTot[dfTot.SubDS == subNum]
            # dfSub = dfTot.query("SubDS == %d")
            dfSub.set_index('Cut', inplace=True)

            # Set cuts here
            runCut = "&&run>=%d&&run<=%d" % (cInfo.master["ds%d_m%d"%(dsNum,modNum)][subNum][1], cInfo.master["ds%d_m%d"%(dsNum,modNum)][subNum][2])
            # pol2Cut = "&&pol2>%.2e&&pol2<%.2e" % (dfSub[dfSub.Range=='Peak'].loc['pol2','%d'%(ch)][0], dfSub[dfSub.Range=='Peak'].loc['pol2','%d'%(ch)][4])
            # pol3Cut = "&&pol3>%.2e&&pol3<%.2e" % (dfSub[dfSub.Range=='Peak'].loc['pol3','%d'%(ch)][0], dfSub[dfSub.Range=='Peak'].loc['pol3','%d'%(ch)][4])
            fitSloCut = "&&fitSlo<%.2f" % (0.8*dfSub[dfSub.Range=='Peak'].loc['fitSlo','%d'%(ch)][2])
            riseNoiseCut = "&&riseNoise<%.2f" % (dfSub[dfSub.Range=='Peak'].loc['riseNoise','%d'%(ch)][4])

            PSA1 = channelCut + runCut + riseNoiseCut
            PSA2 = channelCut + runCut + riseNoiseCut + fitSloCut

            # Save all cuts into a list to iterate through
            cutList = [channelCut+runCut, PSA1, PSA2]

            for idx2,cuts in enumerate(cutList):
                if specMode:
                    hDict[ch][subNum].append(ROOT.TH1D())
                    hDict[ch][subNum][idx2] = wl.H1D(skimTree,bins, lower, upper, "trapENFCal", cuts, Title="h0-Ch%d-%d-%d"%(ch,subNum,idx2))
                    print ("Drawn: h0-Ch%d-%d-%d"%(ch,subNum,idx2))
                # For each energy range and each cut, get number of events and fill to list
                # for idx3, eRange in enumerate(EnergyList):
                # cutMatrix.append([])
                #     cutMatrix[idx2].append( float(skimTree.GetEntries(cuts+"&&trapENFCal>%.1f&&trapENFCal<%.1f"%(eRange[0], eRange[1])) ) )

        # Create dataframe, add two additional columns with channel # and cut range
        # dfList.append( pd.DataFrame(np.transpose(np.array(cutMatrix)), columns = cutNames) )
        # dfList[idx].loc[:,'Channel'] = ch
        # dfList[idx].loc[:, 'Energy Range'] = pd.Series(EnergyList, index=dfList[idx].index)

    # dfTot = pd.concat(dfList)
    # dfTot.to_csv("%s/SpecList.csv"%(outDir), sep='\t')

    # Merge histograms into a list of histograms per cut
    if specMode:
        ROOT.gStyle.SetOptStat(0)
        c1 = ROOT.TCanvas("c1", "c1", 1100, 800)
        c1.SetLogy()
        leg1 = ROOT.TLegend(0.35, 0.8, 0.65, 0.89)
        leg1.SetBorderSize(0)
        for idx2,cuts in enumerate(cutList):
            hList.append(ROOT.TH1D())
            hList[idx2] = hDict[chList[0]][0][idx2]
            for subNum in cInfo.master["ds%d_m%d"%(dsNum,modNum)].keys():
                for idx, ch in enumerate(chList[1:]):
                    hList[idx2].Add(hDict[ch][subNum][idx2])

            hList[idx2].SetTitle("")
            hList[idx2].GetXaxis().SetTitle("Energy (keV)")
            hList[idx2].GetYaxis().SetTitle("Counts/ %.1f keV"%(float((upper-lower)/bins)))
            hList[idx2].SetMinimum(0.1) # Arbitrary unit right now...
            hList[idx2].SetLineColorAlpha(idx2+1, 0.75)
            hList[idx2].Draw("SAME")
            leg1.AddEntry(hList[idx2], "%s"%cutNames[idx2] , "l")

        leg1.Draw()
        c1.SaveAs("%s/Spec_ds%d_m%d.pdf"%(outDir,dsNum,modNum))
        c1.SaveAs("%s/Spec_ds%d_m%d.C"%(outDir,dsNum,modNum))

if __name__ == "__main__":
    main(sys.argv[1:])
