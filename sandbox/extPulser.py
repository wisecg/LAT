#!/usr/bin/env python3
import sys, os, imp
sys.argv.append("-b")
ds = imp.load_source('DataSetInfo',os.environ['LATDIR']+'/DataSetInfo.py')
import matplotlib.pyplot as plt
import numpy as np

def main():
    # getSettings()
    range1()


def range1():
    from ROOT import TFile, GATDataSet, MJTChannelMap, MJTChannelSettings

    calInfo = ds.CalInfo()
    runList = calInfo.GetSpecialRuns("extPulser",1)
    run = 4549

    # P6D3 HV pulser, sync on GRETINA Card 10, ch 5
    # Have to load the gat tree b/c chan 677 would have been cut out by the wave-skim DC cut
    hvpsChan = 624
    syncChan = 677

    dsNum = ds.GetDSNum(run)
    runPath = "/global/project/projectdirs/majorana/users/wisecg/special/lat/latSkimDS%d_run%d.root" % (dsNum,run)

    gds = GATDataSet()
    bltPath = gds.GetPathToRun(run,GATDataSet.kBuilt)
    bltFile = TFile(bltPath)
    chMap = bltFile.Get("ChannelMap")
    chSet = bltFile.Get("ChannelSettings")
    enabledIDs = chSet.GetEnabledIDList()
    enabledIDs = [enabledIDs[idx] for idx in range(enabledIDs.size())]
    # chMap.DumpChannelMap()
    detPos = {}
    for ch in enabledIDs:
        if ch%2==1: continue
        pos = "%sD%d" % (chMap.GetString(ch,"kStringName"), chMap.GetInt(ch,"kDetectorPosition"))
        detPos[ch] = pos
    # for k,v in sorted(detPos.items()): print(k,v)

    gatPath = gds.GetPathToRun(run,GATDataSet.kGatified)
    gatFile = TFile(gatPath)
    gatTree = gatFile.Get("mjdTree")

    latFile = TFile(runPath)
    latTree = latFile.Get("skimTree")

    # try to compare the events of the sync channel with the ext pulser channel

    nLat = latTree.Draw("iEvent","channel==624","GOFF")
    hvpsEvts = latTree.GetV1()
    hvpsEvts = [int(hvpsEvts[idx]) for idx in range(nLat)]

    nGat = gatTree.Draw("Entry$","channel==677","GOFF")
    syncEvts = gatTree.GetV1()
    syncEvts = [int(syncEvts[idx]) for idx in range(nGat)]

    nGat2 = gatTree.Draw("Entry$","channel==624","GOFF")
    hvpsEvts2 = gatTree.GetV1()
    hvpsEvts2 = [int(hvpsEvts2[idx]) for idx in range(nGat)]


    print(hvpsEvts[:10])
    print(syncEvts[:10])
    print(hvpsEvts2[:10])
    print(len(hvpsEvts),len(syncEvts),len(hvpsEvts2))


def getSettings():
    from ROOT import TFile, GATDataSet, MJTChannelMap, MJTChannelSettings

    calInfo = ds.CalInfo()
    runList = calInfo.GetSpecialRuns("extPulser",1)

    for run in runList:

        dsNum = ds.GetDSNum(run)
        runPath = "/global/project/projectdirs/majorana/users/wisecg/special/lat/latSkimDS%d_run%d.root" % (dsNum,run)

        gds = GATDataSet()
        bltPath = gds.GetPathToRun(run,GATDataSet.kBuilt)
        bltFile = TFile(bltPath)
        chMap = bltFile.Get("ChannelMap")
        chSet = bltFile.Get("ChannelSettings")
        enabledIDs = chSet.GetEnabledIDList()
        enabledIDs = [enabledIDs[idx] for idx in range(enabledIDs.size())]

        detPos = {}
        for ch in enabledIDs:
            if ch%2==1: continue
            pos = "%sD%d" % (chMap.GetString(ch,"kStringName"), chMap.GetInt(ch,"kDetectorPosition"))
            detPos[ch] = pos
        # for k,v in sorted(detPos.items()): print(k,v)

        latFile = TFile(runPath)
        latTree = latFile.Get("skimTree")

        nPass = latTree.Draw("channel:trapENFCal","","GOFF")
        arrCh = latTree.GetV1()
        arrEn = latTree.GetV2()
        arrCh = [int(arrCh[idx]) for idx in range(nPass)]
        arrEn = [float("%.3f" % arrEn[idx]) for idx in range(nPass)]

        minCh, maxCh = min(arrCh), max(arrCh)
        cts, chans = np.histogram(arrCh, bins=maxCh-minCh+1, range=(minCh,maxCh+1))
        maxCh = int(chans[np.argmax(cts)])
        arrEnMaxCh = [arrEn[idx] for idx in range(nPass) if arrCh[idx]==maxCh]
        loE, hiE = min(arrEnMaxCh), max(arrEnMaxCh)
        print("Run: %d  maxCh %d  %s  loE %.2f  hiE %.2f  cts %d" % (run,maxCh,detPos[maxCh],loE,hiE,np.amax(cts)))

        continue

        from matplotlib.colors import LogNorm
        fig = plt.figure(figsize=(10,5),facecolor='w')
        minCh, maxCh = min(arrCh), max(arrCh)
        minEn, maxEn = min(arrEn), max(arrEn)
        bpC, bpE = 1., 0.2
        nBinsC, nBinsE = int((maxCh-minCh)/bpC), int((maxEn-minEn)/bpE)
        plt.hist2d(arrEn,arrCh, bins=[nBinsE,nBinsC], range=[[minEn,maxEn],[minCh,maxCh]], norm=LogNorm())
        plt.colorbar()
        plt.xlabel("trapENFCal (keV)",horizontalalignment='right',x=1.0)
        plt.ylabel("channel",horizontalalignment='right',y=1.0)
        plt.savefig("../plots/extPulserTest.pdf")


if __name__=="__main__":
    main()