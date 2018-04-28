#!/usr/bin/env python3
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
    plotEff()
    # dumpCutVals()


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
    for ds in [5]:
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

    for i, ch in enumerate(chList[:]):

        fs200, x200 = [], []
        for ci in range(len(sloSpec)):

            # diagnostic spectrum
            # if ch==594:
            #     print(ci, fs200[-1])
            #     if fs200[-1] < 0:
            #         plt.cla()
            #         plt.plot(fSloX, sloSpec[ci][ch][1], ls='steps')
            #         plt.xlabel("fitSlo", ha='right', x=1)
            #         plt.savefig("../plots/slo-ds%d-ch%d.png" % (ds, ch))

            # only save the value if we have a nonzero number of counts
            spec = sloSpec[ci][ch][1]
            nCts = np.sum(spec)
            if nCts < 2: continue
            fs200.append(fSloX[np.argmax(sloSpec[ci][ch][1])])
            x200.append(ci)

        # plot the raw value
        plt.plot(x200, fs200, ".", c=cmap(i))
        plt.axhline(np.mean(fs200), c=cmap(i), linewidth=0.5, label="ch%d: %.2f" % (ch, np.mean(fs200)))
        plt.ylim(-50,400)

    plt.xlabel("calIdx", ha='right', x=1)
    plt.ylabel("fitSlo", ha='right', y=1)

    plt.legend(ncol=3)
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



if __name__=="__main__":
    main()