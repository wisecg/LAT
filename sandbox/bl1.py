#!/usr/bin/env python
import sys, imp, time
sys.argv.append("-b")
ds = imp.load_source('DataSetInfo','../DataSetInfo.py')
wl = imp.load_source('waveLibs','../waveLibs.py')
from ROOT import gDirectory, TChain, MGTWaveform, MJTMSWaveform, GATDataSet
import numpy as np
import matplotlib.pyplot as plt

def main():
    """ baseline study:
    1. high-e events in skim files
    2. pulser events in gatified data
    3. calibration data

    fitBL - no better variable in skims.
    MAYBE trapENFBL in gatified data, but that might not work
    what about just loading the wf itself?

    Have to access individual TFile's because of this old stupid ROOT/MJSW bug:
        AttributeError: 'TChain' object has no attribute 'MGTWaveforms'
        AttributeError: 'TChain' object has no attribute 'event'
    """

    # pulserBL() # scans wf's from built files passing simple cuts
    # skimBL1() # uses waveforms in lat files
    # skimBL2() # uses "fitBL" from LAT files

    # baselineDB() # TODO: fill typical BL values for each bkgIdx


def pulserBL():
    """ For each bkgIdx, get the HG pulser events and all their baseline values.
    Fill a dict of channels w/ the values, then make plots'n stuff.
    """
    dsNum, modNum = 1, 1
    goodList = ds.GetGoodChanList(dsNum)

    # thing we want to make histos with
    baseDict = {ch:[] for ch in goodList}

    # for bkgIdx in range(ds.dsMap[dsNum]+1):
    for bkgIdx in range(1):

        print "DS-%d, bkgIdx %d" % (dsNum, bkgIdx)

        runList = []
        bkgList = ds.bkgRunsDS[dsNum][bkgIdx]
        for idx in range(0,len(bkgList),2):
            [runList.append(run) for run in range(bkgList[idx], bkgList[idx+1]+1)]

        gds = GATDataSet()
        [gds.AddRunNumber(run) for run in runList]
        gatTree = gds.GetGatifiedChain()
        bltTree = gds.GetBuiltChain()
        gatTree.AddFriend(bltTree)

        theCut = "EventDC1Bits && channel%2==0 && trapENFCal>50" # select pulser subset
        gatTree.Draw(">>elist", theCut, "entrylist")
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

            hitList = (iH for iH in xrange(nChans) if gatTree.channel.at(iH) in goodList)
            for iH in hitList:
                if dsNum==2:
                    wf_downsampled = bltTree.event.GetWaveform(iH)
                    wf_regular = bltTree.event.GetAuxWaveform(iH)
                    wf = MJTMSWaveform(wf_downsampled,wf_regular)
                else:
                    wf = bltTree.event.GetWaveform(iH)
                chan = int(gatTree.channel.at(iH))
                signal = wl.processWaveform(wf)
                baseline,_ = signal.GetBaseNoise()
                baseline = float("%.2f" % baseline) # kill precision to save on memory
                baseDict[chan].append(baseline)

    # could pickle the baseDict here, but idk, let's move on

    fig = plt.figure(figsize=(8,7),facecolor='w')
    p1 = plt.subplot(111)
    # for ch in baseDict:
    ch = 608

    arr = np.asarray(baseDict[chan])

    p1.hist(arr)#, bins='auto') # arguments are passed to np.histogram
    p1.set_title("Channel 608")

    plt.savefig("../plots/ch608-bl.pdf")


def skimBL1():
    """ Uses waveforms in LAT files to directly calculate baseline.

    """

    print "skim method"
    start = time.clock()

    dsNum, modNum = 1, 1
    goodList = ds.GetGoodChanList(dsNum)

    # thing we want to make histos with
    baseDict = {ch:[] for ch in goodList}

    skimChain = TChain("skimTree")
    skimChain.Add("~/project/bg-lat/latSkimDS%d_*.root" % dsNum)

    theCut = "channel%2==0 && trapENFCal>50" # select high-e subset
    skimChain.Draw(">>elist", theCut, "entrylist")
    elist = gDirectory.Get("elist")
    skimChain.SetEntryList(elist)

    for iList in range(elist.GetN()):

        entry = skimChain.GetEntryNumber(iList);
        skimChain.LoadTree(entry)
        skimChain.GetEntry(entry)
        nChans = skimChain.channel.size()
        numPass = skimChain.Draw("channel",theCut,"GOFF",1,iList)
        chans = skimChain.GetV1()
        chanList = list(set(int(chans[n]) for n in xrange(numPass)))

        hitList = (iH for iH in xrange(nChans) if skimChain.channel.at(iH) in goodList)
        for iH in hitList:
            wf = skimChain.MGTWaveforms.at(iH)
            chan = int(skimChain.channel.at(iH))
            signal = wl.processWaveform(wf)
            baseline,_ = signal.GetBaseNoise()
            baseline = float("%.2f" % baseline) # kill precision to save on memory
            baseDict[chan].append(baseline)

    for ch in sorted(baseDict):
        print "ch:",ch," - ",baseDict[ch]

    print "time elapsed:",time.clock()-start

    # fig = plt.figure(figsize=(8,7),facecolor='w')
    # p1 = plt.subplot(111)
    # # for ch in baseDict:
    # ch = 608
    #
    # arr = np.asarray(baseDict[chan])
    # print arr
    # return
    #
    # p1.hist(arr, bins=50, range=(0,100)) # arguments are passed to np.histogram
    # p1.set_title("Channel 608")
    #
    # plt.savefig("../plots/ch608-skim.pdf")


def skimBL2():
    """ Uses "fitBL" variable in LAT files to do a draw command """
    print "fitBL method"
    start = time.clock()

    dsNum, modNum = 1, 1
    goodList = ds.GetGoodChanList(dsNum)

    # thing we want to make histos with
    baseDict = {ch:[] for ch in goodList}

    skimChain = TChain("skimTree")
    skimChain.Add("~/project/bg-lat/latSkimDS%d_*.root" % dsNum)

    theCut = "channel%2==0 && trapENFCal>50" # select high-e subset

    n = skimChain.Draw("fitBL:channel",theCut,"goff")
    if n == 0: return
    arrBL = skimChain.GetV1()
    arrCH = skimChain.GetV2()
    arrBL = [float("%.2f" % arrBL[idx]) for idx in range(n)]
    arrCH = [float("%.2f" % arrCH[idx]) for idx in range(n)]

    for idx in range(n):
        bl, ch = arrBL[idx], arrCH[idx]
        baseDict[ch].append(bl)

    for ch in sorted(baseDict):
        print "ch:",ch," - ",baseDict[ch]

    print "time elapsed:",time.clock()-start




if __name__=="__main__":
    main()
    print '\a'