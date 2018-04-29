#!/usr/bin/env python3
""" slo-cut.py
Sandbox code to set the fitSlo cut
and produce verification plots.
"""
import sys, imp, os
import tinydb as db
import numpy as np

import matplotlib as mpl
mpl.use('Agg')
sys.argv.append("-b")
import matplotlib.pyplot as plt
# plt.style.use('../pltReports.mplstyle')
from matplotlib.colors import LogNorm, Normalize

dsi = imp.load_source('dsi', '../dsi.py')
bkg = dsi.BkgInfo()
cal = dsi.CalInfo()
det = dsi.DetInfo()
skipDS6Cal=True
import waveLibs as wl


def main():

    # testStats()
    # plotStats()
    # getCalRunTime()
    # plotEff()
    # dumpCutVals()
    # checkWidth()
    getShift()


def testStats():

    # load the last calibration run set in DS1 and figure out how many
    # counts we have in the m=2 s=238 population to work with.

    ds, calIdx = 1, 33
    calLo, calHi = 12726, 12733 # this is probably a lunchtime cal

    calDB = db.TinyDB("%s/calDB-v2.json" % dsi.latSWDir)
    pars = db.Query()

    # trap and HV thresholds for this calidx
    trapKey = "trapThr_ds1_m1_c33"
    trapVal = dsi.getDBRecord(trapKey,calDB=calDB,pars=pars)
    hvKey = "hvBias_ds1_m1_c33"
    hvVal = dsi.getDBRecord(hvKey,calDB=calDB,pars=pars)

    # pull thresh (keV) values for the bkgIdx closest to this calibration
    cLo, cHi = cal.GetCalRunCoverage("ds1_m1",calIdx)
    bkgRuns = bkg.getRunList(ds)
    bkgRanges = set()
    for run in bkgRuns:
        if cLo <= run <= cHi:
            bkgRanges.add( bkg.GetBkgIdx(ds, run) )
    bkgIdx = list(bkgRanges)[0] # it's 35

    # account for sub-ranges
    bkgRuns = bkg.getRunList(ds,bkgIdx)
    subRanges = bkg.GetSubRanges(ds,bkgIdx)
    if len(subRanges) == 0: subRanges.append((bkgRuns[0], bkgRuns[-1]))
    for subIdx, (runLo, runHi) in enumerate(subRanges):
        threshKey = "thresh_ds%d_bkg%d_sub%d" % (ds, bkgIdx, subIdx) # returns "thresh_ds1_bkg35_sub0"

    # load threshKeV values from bkg/auto-thrsh/db
    threshVal = dsi.getDBRecord(threshKey,calDB=calDB,pars=pars)
    chList = []
    print("DB results")
    for chan in sorted(threshVal):
        thrBad = threshVal[chan][2]
        if thrBad: continue
        thrMu = threshVal[chan][0]
        thrSig = threshVal[chan][1]
        thrKeV = thrMu + 3*thrSig
        print("%d  %.3f  %.3f  %d: %.3f keV" % (chan,thrMu,thrSig,thrBad,thrKeV))
        chList.append(chan)

    # ok, now let's load the cal runs themselves
    calRuns = cal.GetCalList("ds1_m1",calIdx)
    fileList = []
    for run in calRuns:
        latList = dsi.getSplitList("%s/latSkimDS%d_run%d*" % (dsi.calLatDir, ds, run), run)
        tmpList = [f for idx, f in sorted(latList.items())]
        fileList.extend(tmpList)

    # declare the output stuff
    evtIdx, evtSumET, evtHitE, evtChans = [], [], [], []
    thrCal = {ch:[] for ch in chList}

    # loop over LAT cal files
    from ROOT import TFile, TTree
    prevRun = 0
    evtCtr, totCtr, runTime = 0, 0, 0
    for iF, f in enumerate(fileList):

        print("%d/%d %s" % (iF, len(fileList), f))
        tf = TFile(f)
        tt = tf.Get("skimTree")

        # increment the run time and fill the output dict of thresholds
        tt.GetEntry(0)
        run = tt.run
        if run!=prevRun:
            start = tt.startTime_s
            stop = tt.stopTime_s
            runTime += stop-start

            # before applying thresholds (and getting sumET and mHT)
            # save them into the output dict (so we can compare w/ DB later).
            n = tt.Draw("channel:threshKeV:threshSigma","","goff")
            chan, thrM, thrS = tt.GetV1(), tt.GetV2(), tt.GetV3()
            tmpThresh = {}
            for i in range(n):
                if chan[i] not in chList:
                    continue
                if chan[i] in tmpThresh.keys():
                    continue
                thrK = thrM[i] + 3*thrS[i]
                tmpThresh[chan[i]] = [run,thrM[i],thrS[i],thrK]

            # fill the output dict
            for ch in tmpThresh:
                thrCal[ch].append(tmpThresh[ch]) # [run, thrM, thrS, thrK]

        prevRun = run

        # loop over tree
        for iE in range(tt.GetEntries()):
            tt.GetEntry(iE)
            if tt.EventDC1Bits != 0: continue
            totCtr += 1

            # calculate mHT and sumET

            n = tt.channel.size()
            chTmp = np.asarray([tt.channel.at(i) for i in range(n)])
            idxRaw = [i for i in range(tt.channel.size()) if tt.channel.at(i) in chList]
            hitERaw = np.asarray([tt.trapENFCal.at(i) for i in idxRaw])

            # get indexes of hits above threshold
            idxList = [i for i in range(tt.channel.size())
                if tt.channel.at(i) in chList
                and tt.trapENFCal.at(i) > threshVal[tt.channel.at(i)][0] + 3*threshVal[tt.channel.at(i)][1]
                and 0.7 < tt.trapENFCal.at(i) < 9999
                ]
            hitE = np.asarray([tt.trapENFCal.at(i) for i in idxList])

            mHT, sumET = len(hitE), sum(hitE)

            # for now, let's just grab m=2 s=238 evts.
            if mHT!=2: continue
            if not 237.28 < sumET < 239.46: continue

            hitChans = np.asarray([tt.channel.at(i) for i in idxList])

            # save event pairs to output
            evtIdx.append([run,iE])
            evtSumET.append(sumET)
            evtHitE.append(hitE)
            evtChans.append(hitChans)
            evtCtr += 1

    # output stats we got
    print("m2s238 evts:",evtCtr, "total evts:",totCtr, "runTime:",runTime)

    # save output
    np.savez("../plots/slo-m2s238-test.npz", evtIdx, evtSumET, evtHitE, evtChans, thrCal, evtCtr, totCtr, runTime)


def plotStats():

    # load data from testStats
    f = np.load('../plots/slo-m2s238-test.npz')
    evtIdx, evtSumET, evtHitE, evtChans = f['arr_0'], f['arr_1'], f['arr_2'], f['arr_3']
    thrCal = f['arr_4'].item()
    evtCtr, totCtr, runTime = f['arr_5'], f['arr_6'], f['arr_7']

    # load threshKeV values from bkg/auto-thrsh/db
    calDB = db.TinyDB("%s/calDB-v2.json" % dsi.latSWDir)
    pars = db.Query()
    threshDB = dsi.getDBRecord("thresh_ds1_bkg35_sub0",calDB=calDB,pars=pars)

    # throw a threshold warning if any det is above 1 keV (and by how much)
    for ch in thrCal:
        thrChan = np.asarray([val[3] for val in thrCal[ch]])
        thrMean, thrStd = np.mean(thrChan), np.std(thrChan)
        thrDB = threshDB[ch][0] + 3*threshDB[ch][1]
        errString = "Above 1" if thrMean > 1.0 else ""
        # print("ch %d  DB %.3f  CAL %.3f keV (%.3f), nRuns %d  %s" % (ch, thrDB, thrMean, thrStd, len(thrChan), errString))

    # fill hit arrays
    hitE, chan = [], []
    for iE in range(len(evtHitE)):
        hitE.extend(evtHitE[iE])
        chan.extend(evtChans[iE])

    # map channels
    chMap = list(sorted(set(chan)))
    chDict = {chMap[i]:i for i in range(len(chMap))}
    chan = [chDict[chan] for chan in chan]


    # -- plot 1 - hit E spectrum
    fig = plt.figure()

    xLo, xHi, xpb = 0, 250, 1
    x, hE = wl.GetHisto(hitE, xLo, xHi, xpb)

    plt.plot(x, hE, ls='steps', lw=1.5, c='b', label='m=2,s=238 hits')
    plt.xlabel("Energy (keV)", ha='right', x=1.)
    plt.ylabel("Counts", ha='right', y=1.)
    plt.legend(loc=1)
    plt.savefig("../plots/slo-hitE-test.png")


    # -- plot 2 - counts per channel vs E (2d), low-E region
    plt.cla()

    xLo, xHi, xpb = 0.5, 5, 0.2
    yLo, yHi = 0, len(chMap)
    nbx, nby = int((xHi-xLo)/xpb), len(chMap)

    h1,_,_ = np.histogram2d(hitE,chan,bins=[nbx,nby], range=[[xLo,xHi],[yLo,yHi]])
    h1 = h1.T
    im1 = plt.imshow(h1,cmap='jet')#,aspect='auto')#),vmin=hMin,vmax=hMax)#,norm=LogNorm())

    xticklabels = ["%.1f" % t for t in np.arange(0, 5.5, 0.5)]
    yticks = np.arange(0, len(chMap))
    plt.xlabel("Energy (keV)", ha='right', x=1.)
    plt.gca().set_xticklabels(xticklabels)

    plt.ylabel("channel", ha='right', y=1.)
    plt.yticks(yticks)
    plt.gca().set_yticklabels(chMap, fontsize=12)

    # note: can control z axis limits w/ code in LAT/sandbox/sea-plot.py
    fig.colorbar(im1, ax=plt.gca(), fraction=len(chMap)/941, pad=0.04)

    plt.tight_layout()
    plt.savefig("../plots/slo-fsVsHitE-test.png")


    # -- output: counts in each detector under 5 keV

    cLo, cHi, nbx = 0, len(chMap), len(chMap)
    x, hC = wl.GetHisto(chan, cLo, cHi, 1)

    hLow = [0]
    for idx,ch in enumerate(chMap):
        nTot = hC[idx+1] # 0-250 kev
        nLow = np.sum(h1[idx,:]) # 0-5 keV
        hLow.append(nLow)
        nCPB = nLow/(xHi-xLo)/xpb # avg counts per bin, assume flat for now.
        rTot = nTot/runTime
        rLow = nLow/runTime
        rCPB = nCPB/nbx/runTime   # counts/bin/runTime
        rt100Cts = (100/rCPB)/3600. if rCPB !=0 else -1
        print("rt %d  ch %d  rTot %.2f  rLow %.4f  rCPB %.4f / %.1f keV  need RT:%d hrs" % (runTime, ch, rTot, rLow, rCPB, xpb, rt100Cts))


    # -- plot 3 - counts per channel (1d), and a few different energy regions
    plt.cla()

    plt.bar(x-0.5, hC, 0.95, color='b', label='all hits %d-%d' % (0, 250))
    plt.bar(x-0.5, hLow, 0.95, color='r', label='hits %d-%d' % (xLo, xHi))

    plt.xlabel("channel", ha='right', x=1.)
    xticks = np.arange(0, len(chMap))
    plt.xticks(xticks)
    plt.gca().set_xticklabels(chMap, fontsize=12, rotation=90)

    plt.ylabel("Counts, mHT=2, sumET=238 hits", ha='right', x=1.)

    plt.legend(loc=1)
    plt.savefig("../plots/slo-chans-test.png")


def getStats():
    # calculate mu, sig of fitSlo for each channel in each slice

    makePlots = False

    # fname = "/global/projecta/projectdirs/majorana/users/wisecg/cal/eff/eff_ds1_m1_c1.npz"
    fname = "/global/projecta/projectdirs/majorana/users/wisecg/cal/eff/eff_ds0_m1_c16.npz"
    key = fname.split('/')[-1].split(".")[0]
    tmp = key.split("_")
    ds, mod, cIdx = tmp[1], tmp[2], tmp[3]
    print("Scanning:",key)

    f1 = np.load(fname)
    fSloSpec = f1['arr_9'].item()
    x = f1['arr_10']

    fig = plt.figure(figsize=(18,6))
    p1 = plt.subplot(131)
    p2 = plt.subplot(132)
    p3 = plt.subplot(133)

    chList = sorted(list(fSloSpec.keys()))
    for ch in chList:

        h1 = fSloSpec[ch][0] # 0-10 keV
        h2 = fSloSpec[ch][1] # 10-200 keV
        h3 = fSloSpec[ch][2] # 236-240 keV

        max1, avg1, std1, pct1, wid1 = wl.getHistInfo(x,h1)
        max2, avg2, std2, pct2, wid2 = wl.getHistInfo(x,h2)
        max3, avg3, std3, pct3, wid3 = wl.getHistInfo(x,h3)

        print("channel",ch)
        print("0-10:    %-6.2f  %-6.2f  %-6.2f  %-6.2f " % (max1, avg1, std1, wid1), wl.niceList(pct1))
        print("10-200:  %-6.2f  %-6.2f  %-6.2f  %-6.2f " % (max2, avg2, std2, wid2), wl.niceList(pct2))
        print("236-240: %-6.2f  %-6.2f  %-6.2f  %-6.2f " % (max3, avg3, std3, wid3), wl.niceList(pct3))

        if not makePlots: continue

        # save a diagnostic plot
        p1.cla()
        p1.plot(x, h1,'b',lw=2.,ls='steps',label='ch %d' % ch)
        p1.plot(np.nan,np.nan,'.w',label='0-10 keV')
        p1.plot(np.nan,np.nan,'.w',label='max %.2f' % max1)
        p1.plot(np.nan,np.nan,'.w',label='avg %.2f' % avg1)
        p1.plot(np.nan,np.nan,'.w',label='wid %.2f' % wid1)
        p1.axvline(pct1[0], c='g', label='pct5 %.2f' % pct1[0])
        p1.axvline(pct1[2], c='r', label='pct90 %.2f' % pct1[2])
        p1.set_xlabel("fitSlo",ha='right',x=1)
        p1.legend(loc=2, fontsize=12)

        p2.cla()
        p2.plot(x, h2,'b',lw=2.,ls='steps',label='ch %d' % ch)
        p2.plot(np.nan,np.nan,'.w',label='10-200 keV')
        p2.plot(np.nan,np.nan,'.w',label='max %.2f' % max2)
        p2.plot(np.nan,np.nan,'.w',label='avg %.2f' % avg2)
        p2.plot(np.nan,np.nan,'.w',label='wid %.2f' % wid2)
        p2.axvline(pct2[0], c='g', label='pct5 %.2f' % pct2[0])
        p2.axvline(pct2[2], c='r', label='pct90 %.2f' % pct2[2])
        p2.set_xlabel("fitSlo",ha='right',x=1)
        p2.legend(loc=2, fontsize=12)

        p3.cla()
        p3.plot(x, h3,'b',lw=2.,ls='steps',label='ch %d' % ch)
        p3.plot(np.nan,np.nan,'.w',label='236-240 keV')
        p3.plot(np.nan,np.nan,'.w',label='max %.2f' % max3)
        p3.plot(np.nan,np.nan,'.w',label='avg %.2f' % avg3)
        p3.plot(np.nan,np.nan,'.w',label='wid %.2f' % wid3)
        p3.axvline(pct3[0], c='g', label='pct5 %.2f' % pct3[0])
        p3.axvline(pct3[2], c='r', label='pct90 %.2f' % pct3[2])
        p3.set_xlabel("fitSlo",ha='right',x=1)
        p3.legend(loc=2, fontsize=12)

        plt.tight_layout()
        plt.savefig('./plots/lat2-diag-ch%d.png' % ch)

        # return


def getCalRunTime():
    """
    Need to know the total run time of all cal runs in each DS.
    that's how we can predict the statistics before going through
    the trouble of a full scan over calibration data

    Rough prediction from plotStats:
    Need ~200 hours to get to 100 cts in every 0.2 keV bin.
    """
    from ROOT import GATDataSet, TFile, TTree, MJTRun

    for ds in [0,1,2,3,4,5]:

        runList = []

        # load standard cals
        for key in cal.GetKeys(ds):
            for sub in range(cal.GetIdxs(key)):
                runList.extend(cal.GetCalList(key,sub))
        print("DS",ds,"num standard cals:",len(runList))

        # load long cals
        lIdx = {0:[0], 1:[1], 2:[], 3:[2], 4:[3], 5:[5,6]}
        for l in lIdx[ds]:
            runList.extend(cal.GetSpecialRuns("longCal",l))
        runList = sorted(list(set(runList)))
        print("DS",ds,"num adding longcals:",len(runList))

        # use GDS once just to pull out the path.
        gds = GATDataSet()
        runPath = gds.GetPathToRun(runList[0],GATDataSet.kBuilt)
        filePath = '/'.join(runPath.split('/')[:-1])

        totCalRunTime = 0

        # get run time from built files (no tree loading)
        for iR, run in enumerate(runList):

            # print progress
            # if np.fabs(100*iR/len(runList) % 10) < 0.1:
                # print("%d/%d  run %d  RT %.2f hrs" % (iR, len(runList), run, totCalRunTime/3600))

            f = filePath+"/OR_run%d.root" % run
            tf = TFile(f)
            rInfo = tf.Get("run")
            start = rInfo.GetStartTime()
            stop = rInfo.GetStopTime()
            runTime = stop-start
            if runTime < 0 or runTime > 9999:
                print("error, run",run,"start",start,"stop")
                continue
            totCalRunTime += stop-start
            tf.Close()

        print("Total cal run time, DS%d: %.2f hrs." % (ds, totCalRunTime/3600))


def plotEff():

    # arrays to plot m2s238 data
    effHitE = []  # [hitE1, hitE2 , ...] (remove sub-list of input format)
    effChan = []  # [chan1, chan2 , ...]
    effSlo = []   # [fSlo1, fSlo2, ...]
    effRise = []  # [rise1, rise2, ...]
    effRun = []   # [run1, run1, ...]

    sloSpec = [] # array of fitSlo histo dicts (i should have used pandas probably)

    # load efficiency files
    fList = []
    for ds in [4]:
        print("Loading DS-%d" % ds)
        for key in cal.GetKeys(ds):
            mod = -1
            if "m1" in key: mod = 1
            if "m2" in key: mod = 2
            for cIdx in range(cal.GetIdxs(key)):
                eFile = "%s/eff_%s_c%d.npz" % (dsi.effDir, key, cIdx)
                if os.path.isfile(eFile):
                    fList.append([ds,cIdx,mod,eFile])
                else:
                    print("File not found:",eFile)
                    continue
    for ds,ci,mod,ef in fList:
        # print(ds,ci,mod,ef)
        f = np.load(ef)
        evtIdx = f['arr_0']          # m2s238 event [[run,iE] , ...]
        evtSumET = f['arr_1']        # m2s238 event [sumET , ...]
        evtHitE = f['arr_2']         # m2s238 event [[hitE1, hitE2] , ...]
        evtChans = f['arr_3']        # m2s238 event [[chan1, chan2] , ...]
        thrCal = f['arr_4'].item()   # {ch : [run,thrM,thrS,thrK] for ch in goodList(ds)}
        thrFinal = f['arr_5'].item() # {ch : [thrAvg, thrDev] for ch in goodList(ds)}
        evtCtr = f['arr_6']          # num m2s238 evts
        totCtr = f['arr_7']          # num total evts
        runTime = f['arr_8']         # cal run time
        fSloSpec = f['arr_9'].item() # fitSlo histos (all hits) {ch:[h10, h200, h238] for ch in chList}
        fSloX = f['arr_10']          # xVals for fitSlo histos
        evtSlo = f['arr_11']         # m2s238 event [[fSlo1, fSlo2], ...]
        evtRise = f['arr_12']        # m2s238 event [[rise1, rise2], ...]

        sloSpec.append(fSloSpec)

        # remove the hit pair
        for i in range(len(evtHitE)):
            effHitE.extend(evtHitE[i])
            effChan.extend(evtChans[i])
            effSlo.extend(evtSlo[i])
            effRise.extend(evtRise[i])
            effRun.extend([evtIdx[i][0],evtIdx[i][0]])

    effHitE = np.asarray(effHitE)
    effChan = np.asarray(effChan)
    effSlo = np.asarray(effSlo)
    effRise = np.asarray(effRise)
    effRun = np.asarray(effRun)

    chList = det.getGoodChanList(ds)

    # -- MAKE PLOTS --
    fig = plt.figure(figsize=(9,7))

    # # -- 1. hit spectrum, all channels, 0-250
    # plt.cla()
    # xLo, xHi, xpb = 0, 250, 0.2
    # x1, h1 = wl.GetHisto(effHitE,xLo,xHi,xpb)
    # plt.plot(x1, h1, ls='steps', label='ds%d' % ds)
    # plt.xlabel("Energy (keV)",ha='right',x=1.)
    # plt.ylabel("Counts",ha='right',y=1.)
    # plt.legend()
    # plt.savefig("../plots/slo-specTest.png")
    #
    # # -- 2. bar plot, hits in all channels
    # plt.cla()
    # x, yAll = wl.GetHisto(effChan, chList[0], chList[-1], 1)
    # idx = np.where(yAll!=0)
    # x, yAll = x[idx]-0.5, yAll[idx]
    # x = [int(ch) for ch in x]
    # xb = np.arange(0,len(x),1)
    # plt.bar(xb, yAll, 0.95, color='b', label='ds%d all hits' % ds)
    #
    # # hits under 20 keV
    # idx = np.where(effHitE < 20)
    # x, yLow = wl.GetHisto(effChan[idx], chList[0], chList[-1]+1, 1)
    # idx = np.where(yLow!=0)
    # x, yLow = x[idx]-0.5, yLow[idx]
    # x = [int(t) for t in x]
    # print("x:",x)
    # plt.bar(xb, yLow, 0.95, color='r', label='ds%d < 20 keV' % ds)
    #
    # # plt.gca().set_ylim(1)
    # plt.gca().set_yscale('log')
    #
    # plt.xticks(xb)
    # plt.gca().set_xticklabels(x, fontsize=12, rotation='vertical')
    # plt.xlabel("channel", ha='right', x=1.)
    # leg = plt.legend(fontsize=14, ncol=2)
    # leg.get_frame().set_alpha(0.6)
    # plt.tight_layout()
    # plt.savefig("../plots/slo-chans.png")
    #
    # # -- 3. 2d hits vs channels, 0-20 keV
    # plt.cla()
    #
    # chDict = {chList[i]:i for i in range(len(chList))}
    # chan = [chDict[chan] for chan in effChan]
    # xLo, xHi, xpb = 0, 20., 0.2
    # yLo, yHi = 0, len(chList)
    # nbx, nby = int((xHi-xLo)/xpb), len(chList)
    #
    # plt.hist2d(effHitE, chan, bins=[nbx, nby], range=[[xLo,xHi],[yLo,yHi]], cmap='jet')
    # plt.xlabel("Energy (keV)", ha='right', x=1.)
    # plt.xticks(np.arange(xLo, xHi+1, 1.0))
    # plt.ylabel("channel", ha='right', y=1.)
    # yticks = np.arange(0, len(chList))
    # plt.yticks(yticks+0.5)
    # plt.gca().set_yticklabels(chList, fontsize=10)
    # plt.tight_layout()
    # plt.savefig("../plots/slo-hist2d.png")

    # -- 4. typical fitSlo values
    # for ch in chList[:]:
    #     plt.cla()
    #     plt.semilogy(fSloX, fSloSpec[ch][0], 'r', ls='steps', label='ds%d ch%d 0-10' % (ds,ch))
    #     plt.semilogy(fSloX, fSloSpec[ch][1], 'g', ls='steps', label='ds%dch%d 10-200' % (ds,ch))
    #     plt.semilogy(fSloX, fSloSpec[ch][2], 'b', ls='steps', label='ds%dch%d 236-240' % (ds,ch))
    #
    #     maxLo = fSloX[np.argmax(fSloSpec[ch][0])]
    #     maxHi = fSloX[np.argmax(fSloSpec[ch][1])]
    #     plt.axvline(maxLo, c='k', label='max 10-200: %.1f' % maxLo)
    #     plt.axvline(maxHi, c='m', label='max 0-10: %.1f' % maxHi)
    #
    #     plt.xlabel("fitSlo",ha='right',x=1)
    #     plt.legend()
    #     plt.tight_layout()
    #     plt.savefig("../plots/slo-spec-ds%d-%d.png" % (ds,ch))

    # -- 5. m2s238 slowness
    # for ch in chList[:]:
    #
    #     idx = np.where(effChan==ch)
    #     tmpE = effHitE[idx]
    #     tmpS = effSlo[idx]
    #
    #     idx2 = np.where(tmpE < 10)
    #     nCtsLo = len(idx2[0])
    #
    #     plt.cla()

        # 1d energy
        # xLo, xHi, xpb = 0, 250, 1
        # plt.plot(*(wl.GetHisto(tmpE, xLo, xHi, xpb)),c='r',ls='steps')

        # 2d energy vs slowness
        # xLo, xHi, xpb = 0, 250, 1
        # nbx = int((xHi-xLo)/xpb)
        # yLo, yHi, ypb = -50, 400, 1
        # nby = int((yHi-yLo)/ypb)
        # _,_,_,im = plt.hist2d(tmpE, tmpS, bins=[nbx, nby], range=[[xLo,xHi],[yLo,yHi]], norm=LogNorm(), cmap='jet')
        # plt.plot(np.nan, np.nan, c='w', label='m2s238 ch%d  nCts %d  nCts0-10 %d' % (ch, len(tmpE), nCtsLo))
        # # cb = plt.colorbar()
        # plt.xlabel("Energy (keV)", ha='right', x=1.)
        # plt.ylabel("fitSlo", ha='right', y=1.)
        # plt.legend()
        # plt.tight_layout()
        # plt.savefig("../plots/slo-ch%d-tmp.png" % ch)

        # # 1d slowness (with 90% cut value)
        # yLo, yHi, ypb = -50, 400, 1
        # x, hSlo = wl.GetHisto(tmpS, yLo, yHi, ypb)
        #
        # max1, avg1, std1, pct1, wid1 = wl.getHistInfo(x,hSlo)
        #
        # plt.plot(x, hSlo, c='b', ls='steps', label="m2s238 ch %d" % ch)
        # plt.axvline(pct1[2], c='r', label='90%% value: %.1f' % pct1[2])
        # plt.legend()
        #
        # plt.savefig("../plots/slo-ds%d-ch%d-m2s238.png" % (ds, ch))


    # -- 5. difference between m2s238 and pk238 90pct cut values
    # for ch in chList[:]:
    #
    #     # compare the 10-200 max w. the m2s238 max from 10-200
    #
    #     maxLo = fSloX[np.argmax(fSloSpec[ch][0])] # 0-10
    #     maxHi = fSloX[np.argmax(fSloSpec[ch][1])] # 10-200
    #     maxPk = fSloX[np.argmax(fSloSpec[ch][2])] # 236-240
    #
    #     idx = np.where(effChan==ch)
    #     tmpS = effSlo[idx]
    #     yLo, yHi, ypb = -50, 400, 1
    #     x, hSlo = wl.GetHisto(tmpS, yLo, yHi, ypb)
    #     maxEff = x[np.argmax(hSlo)]
    #
    #     max1, avg1, std1, pct1, wid1 = wl.getHistInfo(x,hSlo)
    #     pct90 = pct1[2]
    #
    #     max2, avg2, std2, pct2, wid2 = wl.getHistInfo(fSloX,fSloSpec[ch][2])
    #     pk90 = pct2[2]
    #
    #     diff = pk90-pct90
    #
    #     # print(maxLo, maxHi, maxPk, maxEff, "m2s238 90",pct90, "pk90",pk90,"diff",diff)
    #
    #     print(ch,diff)


    # -- 6. stability of fitSlo vs run number (calIdx)
    # and plot as a function of run number
    # gonna also need to pull in HV changes from the DB

    # sweep over values
    # for ch in chList:
    #     fs10, fs200 = [], []
    #     for ci in range(len(sloSpec)):
    #         fs10.append(fSloX[np.argmax(sloSpec[ci][ch][0])])
    #         fs200.append(fSloX[np.argmax(sloSpec[ci][ch][1])])
    #     fs10, fs200 = np.asarray(fs10), np.asarray(fs200)
    #     print("%d  fs10 %.2f pm %.2f  fs200 %.2f pm %.2f" % (ch, np.mean(fs10), np.std(fs10), np.mean(fs200), np.std(fs200)))

    # plot vals by calIdx
    # nCal = np.arange(len(sloSpec))
    cmap = plt.cm.get_cmap('hsv', len(chList)+1)
    plt.cla()

    fig2 = plt.figure(figsize=(10,8))
    p1 = plt.subplot(211)
    p2 = plt.subplot(212)

    for i, ch in enumerate(chList[:]):

        # TODO: compute avg num counts,
        # then throw a warning if the avg counts are low

        fs200, x200, fsm2s238 = [], [], []
        for ci in range(len(sloSpec)):

            # only save the value if we have a nonzero number of counts
            spec = sloSpec[ci][ch][1]
            nCts = np.sum(spec)
            if nCts < 2: continue
            # print(ds,ch,ci,nCts)

            # get the width
            max, avg, std, pct, wid = wl.getHistInfo(fSloX, sloSpec[ci][ch][1])

            # TODO: smarter way to get the width
            # like a FWHM.  find the max, then find the point of 50% reduction on either side

            fs200.append(fSloX[np.argmax(sloSpec[ci][ch][1])])
            x200.append(ci)

            # get m2s238 events from this calIdx and find the typical value
            idx = np.where(effChan==ch)
            tmpS = effSlo[idx]
            tmpC = effChan[idx]
            tmpR = effRun[idx]
            thisFS = []
            for j in range(len(tmpR)):
                key = "ds%d_m1" % ds if ch < 1000 else "ds%d_m2" % ds
                if ci == cal.GetCalIdx(key,tmpR[j]):
                    thisFS.append(tmpS[j])
            nEff = len(thisFS)

            yLo, yHi, ypb = -50, 400, 1
            x, hSlo = wl.GetHisto(thisFS, yLo, yHi, ypb)
            maxEff = np.nan if len(thisFS)==0 else x[np.argmax(hSlo)]

            # NOTE: the diff is NEVER more than 1.
            print("%d  %-3d  nTot %-8d  nEff %-5d  wid %-4.0f  fs200 %-4.0f  fsEff %-4.0f  diff %.0f" % (ch, ci, nCts, nEff, wid, fs200[-1], maxEff, fs200[-1]-maxEff))

            fsm2s238.append(maxEff)


        # plot the raw value (stability)
        p1.plot(x200, fs200, ".", c=cmap(i))
        p1.axhline(np.mean(fs200), c=cmap(i), linewidth=0.5, label="ch%d: %.2f" % (ch, np.mean(fs200)))
        p1.set_ylim(-50,400)

        # plot the difference from the average (deviation)
        fAvg = np.mean(fs200)
        fDev = [(f-fAvg) for f in fs200]
        p2.plot(x200, fDev, ".", c=cmap(i), label="ch%d  fAvg %.0f" % (ch, fAvg))

        # plot the difference between the raw value and the m2s238 value
        # man, i shoulda just added the calIdx of the m2s238 hits



    p1.set_xlabel("calIdx", ha='right', x=1)
    p1.set_ylabel("fitSlo", ha='right', y=1)
    if ds!=5: p1.legend(ncol=3)
    else: p1.legend(ncol=6, fontsize=8)
    p2.set_ylabel("fitSlo Deviation from avg", ha='right', y=1)
    plt.tight_layout()
    plt.savefig("../plots/slo-stability-ds%d.png" % (ds))


def dumpCutVals():

    # compare to this cut value
    # slo-ds1-ch608-m2s238: 88.5 NICE, THIS IS LOWER

    ds, ch, mod = 1, 608, 1

    calDB = db.TinyDB('../calDB.json')
    pars = db.Query()

    for cIdx in range(cal.GetNCalIdxs(ds,mod)):
        fsD = dsi.getDBRecord("fitSlo_ds%d_idx%d_m%d_Peak" % (ds, cIdx, mod), False, calDB, pars)

        tmpCut = fsD[ch] # [1%,5%,90%,95%,99%] v1 used 90%
        db90 = tmpCut[2]
        print(ch, cIdx, db90)


def checkWidth():

    # arrays to plot m2s238 data
    effHitE = []  # [hitE1, hitE2 , ...] (remove sub-list of input format)
    effChan = []  # [chan1, chan2 , ...]
    effSlo = []   # [fSlo1, fSlo2, ...]
    effRise = []  # [rise1, rise2, ...]
    effRun = []   # [run1, run1, ...]
    sloSpec = [] # array of fitSlo histo dicts (i should have used pandas probably)

    # load efficiency files
    fList = []
    for ds in [4]:
        print("Loading DS-%d" % ds)
        for key in cal.GetKeys(ds):
            mod = -1
            if "m1" in key: mod = 1
            if "m2" in key: mod = 2
            for cIdx in range(cal.GetIdxs(key)):
                eFile = "%s/eff_%s_c%d.npz" % (dsi.effDir, key, cIdx)
                if os.path.isfile(eFile):
                    fList.append([ds,cIdx,mod,eFile])
                else:
                    print("File not found:",eFile)
                    continue
    for ds,ci,mod,ef in fList:
        # print(ds,ci,mod,ef)
        f = np.load(ef)
        evtIdx = f['arr_0']          # m2s238 event [[run,iE] , ...]
        evtSumET = f['arr_1']        # m2s238 event [sumET , ...]
        evtHitE = f['arr_2']         # m2s238 event [[hitE1, hitE2] , ...]
        evtChans = f['arr_3']        # m2s238 event [[chan1, chan2] , ...]
        thrCal = f['arr_4'].item()   # {ch : [run,thrM,thrS,thrK] for ch in goodList(ds)}
        thrFinal = f['arr_5'].item() # {ch : [thrAvg, thrDev] for ch in goodList(ds)}
        evtCtr = f['arr_6']          # num m2s238 evts
        totCtr = f['arr_7']          # num total evts
        runTime = f['arr_8']         # cal run time
        fSloSpec = f['arr_9'].item() # fitSlo histos (all hits) {ch:[h10, h200, h238] for ch in chList}
        fSloX = f['arr_10']          # xVals for fitSlo histos
        evtSlo = f['arr_11']         # m2s238 event [[fSlo1, fSlo2], ...]
        evtRise = f['arr_12']        # m2s238 event [[rise1, rise2], ...]

        sloSpec.append(fSloSpec)

        # remove the hit pair
        for i in range(len(evtHitE)):
            effHitE.extend(evtHitE[i])
            effChan.extend(evtChans[i])
            effSlo.extend(evtSlo[i])
            effRise.extend(evtRise[i])
            effRun.extend([evtIdx[i][0],evtIdx[i][0]])

    effHitE = np.asarray(effHitE)
    effChan = np.asarray(effChan)
    effSlo = np.asarray(effSlo)
    effRise = np.asarray(effRise)
    effRun = np.asarray(effRun)
    chList = det.getGoodChanList(ds)

    for ci in range(len(sloSpec)):
        for ch in chList:

            # Get mode (maximum of hist) and 50% width of the 10-200 hits.
            # This is what we use to shift m2s238.

            h200 = sloSpec[ci][ch][1]
            if np.sum(h200)==0:
                print("ci %d  ch %d  no counts" % (ci, ch))
            h200Bin = np.argmax(h200)
            h200Max = fSloX[h200Bin]
            h200BinLo, h200BinHi = -1, -1
            for j in range(len(h200)):
                if h200BinLo==-1 and h200[j] >= h200[h200Bin]/2.:
                    h200BinLo = j
                if j > h200Bin and h200[j] <= h200[h200Bin]/2.:
                    h200BinHi = j
                    break
            h200Lo, h200Hi = fSloX[h200BinLo], fSloX[h200BinHi]
            h200Wid = h200Hi - h200Lo

            # get the maximum of the m2s238 hits in this range (limited stats)
            idx = np.where(effChan==ch)
            tmpS = effSlo[idx]
            tmpC = effChan[idx]
            tmpR = effRun[idx]
            thisFS = []
            for j in range(len(tmpR)):
                key = "ds%d_m1" % ds if ch < 1000 else "ds%d_m2" % ds
                if ci == cal.GetCalIdx(key,tmpR[j]):
                    thisFS.append(tmpS[j])
            nEff = len(thisFS)
            yLo, yHi, ypb = -50, 400, 1
            x, hSlo = wl.GetHisto(thisFS, yLo, yHi, ypb)
            effMax = np.nan if len(thisFS)==0 else x[np.argmax(hSlo)]

            modeDiff = h200Max - effMax

            print("%d  %-4d  %-4d  %-4d  wid %d  h200-eff %.1f" % (ch, h200Lo, h200Max, h200Hi, h200Wid, modeDiff))

            # plot the fitSlo 10-200 distibution, mode, and the width
            # fig = plt.figure()
            # plt.plot(fSloX, h200, ls='steps', c='b')
            # plt.xlim(50,120)
            # plt.axvline(h200Max-0.5,c='g')
            # plt.axvline(h200Lo-0.5, c='r')
            # plt.axvline(h200Hi-0.5, c='r')
            # plt.axhline(h200[h200Bin]/2.)
            # plt.savefig("../plots/slo-width-ch%d.png" % ch)
            # return


def getShift():
    # Brian says shifting might introduce systematic error
    # and I should just throw away the calIdx's that deviate from the mean for the DS.
    # if i do that, then we lose huge chunks of data.

    # arrays to plot m2s238 data
    effHitE = []  # [hitE1, hitE2 , ...] (remove sub-list of input format)
    effChan = []  # [chan1, chan2 , ...]
    effSlo = []   # [fSlo1, fSlo2, ...]
    effRise = []  # [rise1, rise2, ...]
    effRun = []   # [run1, run1, ...]
    sloSpec = [] # array of fitSlo histo dicts (i should have used pandas probably)

    # load efficiency files
    fList = []
    for ds in [1]:
        print("Loading DS-%d" % ds)
        for key in cal.GetKeys(ds):
            mod = -1
            if "m1" in key: mod = 1
            if "m2" in key: mod = 2
            for cIdx in range(cal.GetIdxs(key)):
                eFile = "%s/eff_%s_c%d.npz" % (dsi.effDir, key, cIdx)
                if os.path.isfile(eFile):
                    fList.append([ds,cIdx,mod,eFile])
                else:
                    print("File not found:",eFile)
                    continue
    for ds,ci,mod,ef in fList:
        # print(ds,ci,mod,ef)
        f = np.load(ef)
        evtIdx = f['arr_0']          # m2s238 event [[run,iE] , ...]
        evtSumET = f['arr_1']        # m2s238 event [sumET , ...]
        evtHitE = f['arr_2']         # m2s238 event [[hitE1, hitE2] , ...]
        evtChans = f['arr_3']        # m2s238 event [[chan1, chan2] , ...]
        thrCal = f['arr_4'].item()   # {ch : [run,thrM,thrS,thrK] for ch in goodList(ds)}
        thrFinal = f['arr_5'].item() # {ch : [thrAvg, thrDev] for ch in goodList(ds)}
        evtCtr = f['arr_6']          # num m2s238 evts
        totCtr = f['arr_7']          # num total evts
        runTime = f['arr_8']         # cal run time
        fSloSpec = f['arr_9'].item() # fitSlo histos (all hits) {ch:[h10, h200, h238] for ch in chList}
        fSloX = f['arr_10']          # xVals for fitSlo histos
        evtSlo = f['arr_11']         # m2s238 event [[fSlo1, fSlo2], ...]
        evtRise = f['arr_12']        # m2s238 event [[rise1, rise2], ...]

        sloSpec.append(fSloSpec)

        # remove the hit pair
        for i in range(len(evtHitE)):
            effHitE.extend(evtHitE[i])
            effChan.extend(evtChans[i])
            effSlo.extend(evtSlo[i])
            effRise.extend(evtRise[i])
            effRun.extend([evtIdx[i][0],evtIdx[i][0]])

    effHitE = np.asarray(effHitE)
    effChan = np.asarray(effChan)
    effSlo = np.asarray(effSlo)
    effRise = np.asarray(effRise)
    effRun = np.asarray(effRun)

    chList = det.getGoodChanList(ds)

    # for every calIdx:
    # get the avg value of fs200 and the width for every channel
    # fsAvg = {ch : [avg, width]}
    # then shift the m2s238 hits for that channel by the avg value
    # fsShift = []
    # then at the end of the DS, compute the 90% value for the shifted fitSlo of every channel
    # then we save into the DB the 90% value and the mean, for every calIdx

    # this stores the shifted m2s238 spectra for each channel in the DS
    yLo, yHi, ypb = -50, 50, 1
    nby = int((yHi-yLo)/ypb)
    shiftSpec = {ch:np.zeros(nby+1) for ch in chList}

    shiftDict = {ci:None for ci in range(len(sloSpec))}

    for ci in range(len(sloSpec)):

        shiftDict[ci] = {ch:[] for ch in chList}

        for ch in chList:

            # Get mode (maximum) and 50% width of the 10-200 hits.
            # This is what we use to shift m2s238.
            h200 = sloSpec[ci][ch][1]
            if np.sum(h200)==0:
                print("ci %d  ch %d  no counts" % (ci, ch))
            fsBin = np.argmax(h200)
            fsMax = fSloX[fsBin]
            fsBinLo, fsBinHi = -1, -1
            for j in range(len(h200)):
                if fsBinLo==-1 and h200[j] >= h200[fsBin]/2.:
                    fsBinLo = j
                if j > fsBin and h200[j] <= h200[fsBin]/2.:
                    fsBinHi = j
                    break
            fsLo, fsHi = fSloX[fsBinLo], fSloX[fsBinHi]
            fsWid = fsHi - fsLo

            # save the max and width of h200
            shiftDict[ci][ch].extend([fsMax, fsWid])

            # now histogram the shifted m2s238 vals
            idx = np.where(effChan==ch)
            tmpS = effSlo[idx]
            tmpC = effChan[idx]
            tmpR = effRun[idx]
            thisFS = []
            for j in range(len(tmpR)):
                key = "ds%d_m1" % ds if ch < 1000 else "ds%d_m2" % ds
                if ci == cal.GetCalIdx(key,tmpR[j]):
                    thisFS.append(tmpS[j] - fsMax) # ** apply the shift **
            if len(thisFS)==0: continue
            x, hSlo = wl.GetHisto(thisFS, yLo, yHi, ypb)

            # add to the histogram for this ch in this DS
            shiftSpec[ch] = np.add(shiftSpec[ch], hSlo)

            # print("%d  %d  %-4d  %-4d  %-4d  wid %d  n238 %d" % (ci, ch, fsLo, fsMax, fsHi, fsWid, np.sum(shiftSpec[ch])))

    # now plot the shifted m2s238 spectra (could also plot against the unshifted)
    for ch in shiftSpec:

        # find the 90% cut value (shifted)
        max1, avg1, std1, pct1, wid1 = wl.getHistInfo(x,shiftSpec[ch])
        pct90 = pct1[2]

        # add it to shiftDict
        for ci in shiftDict:
            shiftDict[ci][ch].extend([pct90])

        plt.cla()
        plt.plot(x, shiftSpec[ch], c='b', ls='steps', label='ds %d ch %d' % (ds,ch))
        plt.axvline(pct90, c='r', label='90pct cut: %d' % (pct90))
        plt.xlabel("fitSlo", ha='right', x=1)
        plt.legend(loc=1)
        plt.savefig("../plots/slo-m2s238shift-ds%d-ch%d.png" % (ds,ch))
        # return

    # now print all the shifted values and compare to the previous DB value
    calDB = db.TinyDB('../calDB.json')
    pars = db.Query()

    for ci in shiftDict:
        print("cIdx",ci)

        for ch in shiftDict[ci]:

            v2Cut90 = shiftDict[ci][ch][0] + shiftDict[ci][ch][2] # fs max + m2s238 90% val

            mod = 1 if ch < 1000 else 2
            fsD = dsi.getDBRecord("fitSlo_ds%d_idx%d_m%d_Peak" % (ds, ci, mod), False, calDB, pars)
            if ch in fsD.keys():
                v1Cut90 = fsD[ch][2] # [1%,5%,90%,95%,99%] v1 used 90%
                print("ds %d  cIdx %d  ch%d  v1Cut90 %-6.1f  v2Cut90 %-6.1f  diff %.1f" % (ds, ci, ch, v1Cut90, v2Cut90, v2Cut90-v1Cut90))
            else:
                print("ds %d  cIdx %d  ch%d  v1Cut90 %-6.1f  v2Cut90 %-6.1f  diff %.1f" % (ds, ci, ch, np.nan, v2Cut90, np.nan))



if __name__=="__main__":
    main()