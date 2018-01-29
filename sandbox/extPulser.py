#!/usr/bin/env python3
import sys, os, imp, glob
sys.argv.append("-b")
ds = imp.load_source('DataSetInfo',os.environ['LATDIR']+'/DataSetInfo.py')
wl = imp.load_source('waveLibs',os.environ['LATDIR']+'/waveLibs.py')
import matplotlib.pyplot as plt
import numpy as np

def main():

    # checkSyncChannel()
    # plotHits()
    runByRun()

def runByRun():
    """ Directly confirm settings of ext pulser scripts. """
    import time
    from ROOT import TChain

    calInfo = ds.CalInfo()

    runList = calInfo.GetSpecialRuns("extPulser",7)
    for run in runList:
        if run != 5942:
            continue

        fileList = []
        subFiles = glob.glob("%s/lat/latSkimDS%d_run%d_*.root" % (ds.specialDir, ds.GetDSNum(run), run))
        for idx in range(len(subFiles)):
            thisFile = "%s/lat/latSkimDS%d_run%d_%d.root" % (ds.specialDir, ds.GetDSNum(run), run, idx)
            if not os.path.isfile(thisFile):
                print("File doesn't exist: ",thisFile)
            else:
                fileList.append(thisFile)
        latChain = TChain("skimTree")
        for f in fileList: latChain.Add(f)

        syncChan = 677
        extPChan = 674
        # theCut = "channel==%d || channel==%d" % (syncChan,extPChan)
        # theCut = "mH==2 && Entry$ < 50"
        # theCut = "mH==2 && Entry$ < 50 && !muVeto"
        # theCut = "Entry$ < 100"
        # theCut = "Entry$ < 50 && trapENFCal > 1"
        # theCut = "mH==4 && (channel==640 || channel==646)"
        theCut = "Entry$ < 100"

        start = time.time()
        tNames = ["Entry$","run","channel","mH","trapENFCal","den90","den10","fitSlo","localTime_s","tOffset"]
        tVals = wl.GetVX(latChain,tNames,theCut)
        print("took",time.time()-start)

        for idx in range(tVals["run"].size):
            ent    = tVals["Entry$"][idx]
            run    = tVals["run"][idx]
            chan   = tVals["channel"][idx]
            mH     = tVals["mH"][idx]
            enf    = tVals["trapENFCal"][idx]
            d90    = tVals["den90"][idx]
            d10    = tVals["den10"][idx]
            fitSlo = tVals["fitSlo"][idx]
            gt     = tVals["localTime_s"][idx]
            print("%d  e%d  m%d  t%.10f  %-4d  %-9.2f  %-8.2f  %.2f" % (run,ent,mH,gt,chan,enf,d90-d10,fitSlo))




def plotHits():
    from ROOT import TChain, TFile, GATDataSet, MJTChannelMap, MJTChannelSettings

    calInfo = ds.CalInfo()

    for sIdx in range(7,23):
        runList = calInfo.GetSpecialRuns("extPulser",sIdx)
        fileList = []
        for run in runList:
            subFiles = glob.glob("%s/lat/latSkimDS%d_run%d_*.root" % (ds.specialDir, ds.GetDSNum(run), run))
            for idx in range(len(subFiles)):
                thisFile = "%s/lat/latSkimDS%d_run%d_%d.root" % (ds.specialDir, ds.GetDSNum(run), run, idx)
                if not os.path.isfile(thisFile):
                    print("File doesn't exist: ",thisFile)
                else:
                    fileList.append(thisFile)

        gds = GATDataSet()
        bltPath = gds.GetPathToRun(runList[0],GATDataSet.kBuilt)
        bltFile = TFile(bltPath)
        chMap = bltFile.Get("ChannelMap")
        chSet = bltFile.Get("ChannelSettings")
        enabledIDs = chSet.GetEnabledIDList()
        enabledIDs = [enabledIDs[idx] for idx in range(enabledIDs.size())]
        detPos, detID = {}, {}
        for ch in enabledIDs:
            if ch%2==1: continue
            detPos[ch] = "%sD%d" % (chMap.GetString(ch,"kStringName"), chMap.GetInt(ch,"kDetectorPosition"))
            detID[ch] = chMap.GetString(ch,"kDetectorName")

        latChain = TChain("skimTree")
        for f in fileList: latChain.Add(f)

        nPass = latChain.Draw("channel:trapENFCal","","GOFF")
        arrCh = latChain.GetV1()
        arrEn = latChain.GetV2()
        arrCh = [int(arrCh[idx]) for idx in range(nPass)]
        arrEn = [float("%.3f" % arrEn[idx]) for idx in range(nPass)]

        minCh, maxCh = min(arrCh), max(arrCh)
        cts, chans = np.histogram(arrCh, bins=maxCh-minCh+1, range=(minCh,maxCh+1))
        maxCh = int(chans[np.argmax(cts)])
        arrEnMaxCh = [arrEn[idx] for idx in range(nPass) if arrCh[idx]==maxCh]
        loE, hiE = min(arrEnMaxCh), max(arrEnMaxCh)
        print("sIdx %d  maxCh %d  %s  %s loE %.2f  hiE %.2f  cts %d" % (sIdx, maxCh,detPos[maxCh],detID[maxCh],loE,hiE,np.amax(cts)))

        from matplotlib.colors import LogNorm
        fig = plt.figure(figsize=(10,5),facecolor='w')
        minCh, maxCh = min(arrCh), max(arrCh)
        # minEn, maxEn = min(arrEn), max(arrEn) # be careful, the max energy is probably 99999
        minEn, maxEn = -10., 300.
        bpC, bpE = 1., 1.
        nBinsC, nBinsE = int((maxCh-minCh)/bpC), int((maxEn-minEn)/bpE)
        plt.hist2d(arrEn,arrCh, bins=[nBinsE,nBinsC], range=[[minEn,maxEn],[minCh,maxCh]], norm=LogNorm())
        plt.colorbar()
        plt.xlabel("trapENFCal (keV)",horizontalalignment='right',x=1.0)
        plt.ylabel("channel",horizontalalignment='right',y=1.0)
        plt.savefig("../plots/extPulser_idx%d.pdf" % idx)


def checkSyncChannel():
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


if __name__=="__main__":
    main()