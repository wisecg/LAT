#!/usr/bin/env python
import sys, imp, os
sys.argv.append("-b")
ds = imp.load_source('DataSetInfo','../DataSetInfo.py')
# wl = imp.load_source('waveLibs','../waveLibs.py')
# import ROOT
from ROOT import TFile, TChain, TTree, TCanvas, TH1D, MGTWaveform
# from ROOT import gDirectory, gStyle, std
from ROOT import GATDataSet
# import numpy as np
# import matplotlib.pyplot as plt

def main():

    baselineStudy()
    # baselineDB()


def baselineStudy():
    """ 1. high-e events in skim files
        2. pulser events in gatified data
        3. calibration data

        fitBL - no better variable in skims.
        MAYBE trapENFBL in gatified data, but that might not work
        what about just loading the wf itself?
    """
    dsNum, modNum = 1, 1
    chList = ds.GetGoodChanList(dsNum)

    for bkgIdx in range(ds.dsMap[dsNum]+1):

        if bkgIdx==1: break

        runList = []
        bkgList = ds.bkgRunsDS[dsNum][bkgIdx]
        for idx in range(0,len(bkgList),2):
            [runList.append(run) for run in range(bkgList[idx], bkgList[idx+1]+1)]

        gds = GATDataSet()
        [gds.AddRunNumber(run) for run in runList]
        gatChain = gds.GetGatifiedChain()
        bltChain = gds.GetBuiltChain()
        gatChain.AddFriend(bltChain)

        theCut = "EventDC1Bits && Entry$<10000"
        gatTree.Draw(">>elist", theCut, "entrylist") # get a not insane number of pulsers
        elist = gDirectory.Get("elist")
        gatTree.SetEntryList(elist)

        for iList in range(elist.GetN()):

            entry = gatTree.GetEntryNumber(iList);
            gatTree.LoadTree(entry)
            gatTree.GetEntry(entry)
            nChans = gatTree.channel.size()
            numPass = gatTree.Draw("channel",theCut,"GOFF",1,iList)
            chans = gatTree.GetV1()
            chanList = list(set(int(chans[n]) for n in xrange(numPass)))
            event = bltTree.event

            hitList = (iH for iH in xrange(nChans) if gatTree.channel.at(iH) in chanList)
            for iH in hitList:

                if dsNum==2:
                    wf_downsampled = event.GetWaveform(iH)
                    wf_regular = event.GetAuxWaveform(iH)
                    wf = MJTMSWaveform(wf_downsampled,wf_regular)
                else:
                    wf = event.GetWaveform(iH)

                # maybe fill a buncha channel plots independently, instead of doing a buncha draw commands

                run = gatTree.run
                chan = gatTree.channel.at(iH)
                energy = gatTree.trapE.at(iH)

                signal = wl.processWaveform(wf)
                baseline,_ = signal.GetBaseNoise()














if __name__=="__main__":
    main()