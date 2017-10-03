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

    inDir, outDir = ".", "./plots"
    specMode = False
    dsNum = 1
    dType = "Nat"
    bins,lower,upper = 250,0,50

    for i,opt in enumerate(argv):
		if opt == "-s":
			specMode = True
			print "Spectrum Mode"
		if opt == "-p":
			inDir, outDir = argv[i+1], argv[i+2]
			print "Custom paths: Input %s,  Output %s" % (inDir,outDir)
		if opt == "-ds":
			dsNum, dType = int(argv[i+1]), argv[i+2]
			print "Using Dataset: %d -- %s" % (dsNum,dType)


    EnergyList = [[1.,5.], [2., 4.], [4., 9.], [9., 12.], [12., 40.], [40., 50.], [50., 100.]]

    # Channel list
    # NatList0 = [646,664,608,598,600,594,592]
    chList = [600,692] # DS1
    # NatList3 = [600,692,624,694,614] #692 is bad
    # NatList4 = [1106,1144,1170,1174,1136,1232,1330]
    # NatList4 = [1144]
    # NatList5 = [614, 626, 628, 680, 688, 694, 1104, 1120, 1124, 1128, 1170, 1174, 1208, 1330]
    # NatList5 = [680]

    latPath = "/projecta/projectdirs/majorana/users/bxyzhu"
    inPath = latPath + "/lat/latSkimDS%d_*.root" % dsNum

    bkgTree = ROOT.TChain("skimTree")
    bkgTree.Add(inPath)

    # List of histograms (for summing together) and Dictionary of histograms (for different cuts)
    hList, hDict = [], {}

    cutNames = ["BasicCut", "+tailSlope", "+bandTime", "+bcTime", "+bcMax", "+noiseWeight", "+fitSlo", "+threshold"]
    # Create a list of DataFrames for each channel to concatenate at the end
    dfList = []

    for idx,ch in enumerate(chList):
        # Create new key for dictionaries according to channel
        hDict[ch] = []
        # Get threshold info
        goodRuns,badRuns,goodRunSigmas = ds.GetThreshDicts(dsNum)

        # Set cuts here
        channelCut = "channel==%d" % ch
        bandTimeCut = "&&bandTime-tOffset-1100<11000"
        bcMaxCut = "&&bcMax<%.2e" % ds.bcMax[dsNum][ch][0] # 99% value
        noiseWt = "(waveS4-waveS1)/bcMax/trapENFCal"
        noiseWtCut = " && %s < %.2e && %s > %.2e" % (noiseWt,ds.noiseWt[dsNum][ch][0],noiseWt,ds.noiseWt[dsNum][ch][1])
        tailSlopeCut = "&&pol2>%.2e&&pol2<%.2e&&pol3>%.2e&&pol3<%.2e" % (ds.pol2[dsNum][ch][1],ds.pol2[dsNum][ch][0],ds.pol3[dsNum][ch][1],ds.pol3[dsNum][ch][0])
        bcTimeCut = "&&(bandTime-tOffset-1100)/(matchTime-tOffset)>%.2e" % ds.bcTime[dsNum][ch]
        fitSloCut = "&&fitSlo<=(%.2f)" % (ds.fitSloCont[dsNum][ch][0])

        PSA1 = channelCut + tailSlopeCut
        PSA2 = channelCut + tailSlopeCut + bandTimeCut
        PSA3 = channelCut + tailSlopeCut + bandTimeCut + bcTimeCut
        PSA4 = channelCut + tailSlopeCut + bandTimeCut + bcTimeCut + bcMaxCut
        PSA5 = channelCut + tailSlopeCut + bandTimeCut + bcTimeCut + bcMaxCut + noiseWtCut
        megaCut = channelCut + bandTimeCut + noiseWtCut + tailSlopeCut + bcTimeCut + fitSloCut

        runCut = "&&("
        for idx2,runRange in enumerate(goodRuns[ch]):
            runCut += "(run>=%d&&run<=%d)||" % (runRange[0],runRange[1])
            # Strip all spaces to save space
            totalCut = megaCut.replace(" ", "") + runCut[:-2] + ")"

        # Save all cuts into a list to iterate through
        cutList = [channelCut, PSA1, PSA2, PSA3, PSA4, PSA5, megaCut, totalCut]
        # Create dummy list and series to store into DataFrame later
        cutMatrix= []
        for idx2,cuts in enumerate(cutList):
            cutMatrix.append([])
            if specMode:
                hDict[ch].append(ROOT.TH1D())
                hDict[ch][idx2] = wl.H1D(bkgTree,"h0-Ch%d-%d"%(ch,idx2),bins,lower,upper,"trapENFCal",cuts)

            # For each energy range and each cut, get number of events and fill to list
            for idx3, eRange in enumerate(EnergyList):
                cutMatrix[idx2].append( float(bkgTree.GetEntries(cuts+"&&trapENFCal>%.1f&&trapENFCal<%.1f"%(eRange[0], eRange[1])) ) )

        # Create dataframe, add two additional columns with channel # and cut range
        dfList.append( pd.DataFrame(np.transpose(np.array(cutMatrix)), columns = cutNames) )
        dfList[idx].loc[:,'Channel'] = ch
        dfList[idx].loc[:, 'Energy Range'] = pd.Series(EnergyList, index=dfList[idx].index)

    dfTot = pd.concat(dfList)
    dfTot.to_csv("%s/SpecList.csv"%(outDir), sep='\t')

    # Merge histograms into a list of histograms per cut
    if specMode:
        ROOT.gStyle.SetOptStat(0)
        c1 = ROOT.TCanvas("c1", "c1", 1100, 800)
        c1.SetLogy()
        leg1 = ROOT.TLegend(0.6, 0.6, 0.89, 0.89)
        leg1.SetBorderSize(0)
        for idx2,cuts in enumerate(cutList):
            hList.append(ROOT.TH1D())
            hList[idx2] = hDict[chList[0]][idx2]
            for idx, ch in enumerate(chList[1:]):
                hList[idx2].Add(hDict[ch][idx2])

            hList[idx2].SetTitle("")
            hList[idx2].GetXaxis().SetTitle("Energy (keV)")
            hList[idx2].GetYaxis().SetTitle("Counts/ %d keV"%((upper-lower)/bins))
            hList[idx2].SetLineColorAlpha(idx2+1, 0.5)
            hList[idx2].Draw("SAME")
            leg1.AddEntry(hList[idx2], "%s"%cutNames[idx2] , "l")

        leg1.Draw()
        c1.SaveAs("%s/SpecTest.pdf"%(outDir))


if __name__ == "__main__":
    main(sys.argv[1:])
