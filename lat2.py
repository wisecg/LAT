#!/usr/bin/env python3
"""
====== LAT2.py =======
Tunes cut parameters,
calculates cut efficiencies.
C. Wiseman, USC
v1. 17 Apr 2018
======================
"""
import sys, time
import numpy as np
import tinydb as db

import matplotlib as mpl
mpl.use('Agg')
sys.argv.append("-b")
import matplotlib.pyplot as plt
plt.style.use('pltReports.mplstyle')
from matplotlib.colors import LogNorm, Normalize

import waveLibs as wl
import dsi
bkg = dsi.BkgInfo()
cal = dsi.CalInfo()
det = dsi.DetInfo()
skipDS6Cal = True # ignore DS6 cal runs until they're processed


def main(argv):

    ds, cIdx, mod = None, None, None

    for i, opt in enumerate(argv):

        # set dataset and options
        if opt=="-ds":
            ds = int(argv[i+1])
        if opt=="-ci":
            ds = int(argv[i+1])
            cIdx = int(argv[i+2])
        if opt=="-m":
            mod = int(argv[i+1])

        # wrapper function for scanRuns
        if opt=="-load":
            loadRuns(ds,cIdx,mod)

        # call scanRuns directly (used by job-panda)
        if opt=="-scan":
            ds, key, mod, cIdx = int(argv[i+1]), argv[i+2], int(argv[i+3]), int(argv[i+4])
            scanRuns(ds,key,mod,cIdx)

        if opt=="-g":
            getStats()


def loadRuns(dsIn=None,subIn=None,modIn=None):

    # loop over datasets, skipping DS6 cal runs till they're processed
    for ds in [0,1,2,3,4,5,6]:
        if skipDS6Cal is True and ds==6:
            continue

        if dsIn is not None and ds!=dsIn:
            continue

        # loop over keys in this DS
        for key in cal.GetKeys(ds):

            mod = -1
            if "m1" in key: mod = 1
            if "m2" in key: mod = 2

            # loop over cIdx's for this key
            for cIdx in range(cal.GetIdxs(key)):
                if subIn is not None and cIdx!=subIn:
                    continue
                if modIn is not None and mod!=modIn:
                    continue

                # now that ds, cIdx, and module are determined, we can call scanRuns
                scanRuns(ds, key, mod, cIdx)


def scanRuns(ds, key, mod, cIdx):
    from ROOT import TFile, TTree

    # load file and channel list
    fileList = []
    calRuns = cal.GetCalList(key,cIdx)
    for run in calRuns:
        latList = dsi.getSplitList("%s/latSkimDS%d_run%d*" % (dsi.calLatDir, ds, run), run)
        tmpList = [f for idx, f in sorted(latList.items())]
        fileList.extend(tmpList)
    chList = det.getGoodChanList(ds)

    print("Scanning DS:%d  calIdx %d  mod %d  key %s  nFiles:%d" % (ds, cIdx, mod, key, len(fileList)), time.strftime('%X %x %Z'))
    outFile = "%s/eff_%s_c%d.npz" % (dsi.effDir, key, cIdx)
    print("Saving output in:",outFile)

    # declare the output stuff
    evtIdx, evtSumET, evtHitE, evtChans = [], [], [], []
    thrCal = {ch:[] for ch in chList}
    fLo, fHi, fpb = -200, 400, 1
    nbf = int((fHi-fLo)/fpb)+1
    fSloSpec = {ch:[np.zeros(nbf) for i in range(3)] for ch in chList} # 0-10, 10-200, 236-240

    # loop over LAT cal files
    scanStart = time.time()
    prevRun = 0
    evtCtr, totCtr, runTime = 0, 0, 0
    for iF, f in enumerate(fileList):

        print("%d/%d %s" % (iF, len(fileList), f))
        tf = TFile(f)
        tt = tf.Get("skimTree")

        # histogram these for each file (and then add to total)
        fs10 = {ch:[] for ch in chList}
        fs200 = {ch:[] for ch in chList}
        fs238 = {ch:[] for ch in chList}

        # increment the run time and fill the output dict of thresholds
        tt.GetEntry(0)
        run = tt.run
        if run!=prevRun:
            start = tt.startTime_s
            stop = tt.stopTime_s
            runTime += stop-start
            if runTime < 0 or runTime > 9999:
                print("run time error, run",run,"start",start,"stop")
            else:
                totCalRunTime += runTime

            # find thresholds for this run s/t we can apply them
            # and calculate sumET and mHT in the loop.
            # save them into the output dict (so we can compare w/ DB later).

            n = tt.Draw("channel:threshKeV:threshSigma","","goff")
            chan, thrM, thrS = tt.GetV1(), tt.GetV2(), tt.GetV3()
            tmpThresh = {}
            for i in range(n):
                if chan[i] not in chList:
                    continue
                if chan[i] in tmpThresh.keys():
                    continue
                if thrM[i] < 9999:
                    thrK = thrM[i] + 3*thrS[i]
                    tmpThresh[chan[i]] = [run,thrM[i],thrS[i],thrK]
            for ch in chList:
                if ch not in tmpThresh.keys():
                    tmpThresh[ch] = [-1,-1,-1,-1]

            # fill the output dict
            for ch in tmpThresh:
                thrCal[ch].append(tmpThresh[ch]) # [run, thrM, thrS, thrK]

        prevRun = run
        # continue

        # loop over tree
        for iE in range(tt.GetEntries()):
            tt.GetEntry(iE)
            if tt.EventDC1Bits != 0: continue
            totCtr += 1

            n = tt.channel.size()
            chTmp = np.asarray([tt.channel.at(i) for i in range(n)])
            idxRaw = [i for i in range(tt.channel.size()) if tt.channel.at(i) in chList]
            hitERaw = np.asarray([tt.trapENFCal.at(i) for i in idxRaw])

            # get indexes of hits above threshold (use thresholds from THIS CAL RUN)
            idxList = [i for i in range(tt.channel.size())
                if tt.channel.at(i) in chList
                and tt.trapENFCal.at(i) > tmpThresh[tt.channel.at(i)][3]
                and 0.7 < tt.trapENFCal.at(i) < 9999
                ]
            hitE = np.asarray([tt.trapENFCal.at(i) for i in idxList])

            # calculate mHT and sumET
            mHT, sumET = len(hitE), sum(hitE)

            # save fitSlo data for 0-10 and 10-200 kev ranges for each channel
            for i in idxList:
                en = tt.trapENFCal.at(i)
                ch = tt.channel.at(i)
                fs = tt.fitSlo.at(i)
                if fLo < fs < fHi:
                    if 0 < en < 10:
                        fs10[ch].append(fs)
                    if 10 < en < 200:
                        fs200[ch].append(fs)
                    if 236 < en < 240:
                        fs238[ch].append(fs)

            # Save m2s238 events to output, skip everything else
            if mHT!=2: continue
            if not 237.28 < sumET < 239.46: continue
            hitChans = np.asarray([tt.channel.at(i) for i in idxList])
            evtIdx.append([run,iE])
            evtSumET.append(sumET)
            evtHitE.append(hitE)
            evtChans.append(hitChans)
            evtCtr += 1

        # fill the fitSlo histograms w/ the events from this file
        for ch in chList:
            x, h1 = wl.GetHisto(fs10[ch],fLo,fHi,fpb)
            x, h2 = wl.GetHisto(fs200[ch],fLo,fHi,fpb)
            x, h3 = wl.GetHisto(fs238[ch],fLo,fHi,fpb)
            fSloSpec[ch][0] = np.sum([fSloSpec[ch][0], h1], axis=0)
            fSloSpec[ch][1] = np.sum([fSloSpec[ch][1], h2], axis=0)
            fSloSpec[ch][2] = np.sum([fSloSpec[ch][2], h3], axis=0)
            # n1 = np.sum(fSloSpec[ch][0])
            # n2 = np.sum(fSloSpec[ch][1])
            # n3 = np.sum(fSloSpec[ch][2])
            # print("ch:%d  n10 %d  n200 %d  n238 %d" % (ch, n1, n2, n3))

    # get average threshold for each channel in this file list
    thrFinal = {chan:[] for chan in thrCal}
    for chan in thrCal:
        thrVals = []
        for iT in range(len(thrCal[chan])):
            run, thrM, thrS, thrK = thrCal[chan][iT]
            # print("%d  %d  %.3f  %.3f  %.3f" % (chan,run,thrM,thrS,thrK))
            if thrK > -1:
                thrVals.append(thrK)
        thrVals = np.asarray(thrVals)
        thrAvg = np.mean(thrVals)
        thrDev = np.std(thrVals)
        # print("%d  %.3f  %.3f" % (chan, thrAvg, thrDev))
        thrFinal[chan] = [thrAvg,thrDev]

    # print to screen the final thresholds, stdev, and an error message if necessary
    print("Detector Thresholds:")
    for chan in sorted(thrFinal):
        thKeV = thrFinal[chan][0]
        thE = thrFinal[chan][1]
        errString = ""
        if thE/thKeV > 0.5:
            errString = ">50pct error:  thE/thKeV=%.2f" % (thE/thKeV)
        if thKeV > 2:
            errString = ">2kev"
        print("%d  %.3f  %.3f  %s" % (chan,thKeV,thE,errString))

    # save output
    np.savez(outFile, evtIdx, evtSumET, evtHitE, evtChans, thrCal, thrFinal, evtCtr, totCtr, totCalRunTime, fSloSpec, x)

    # output stats
    print("Done:",time.strftime('%X %x %Z'),", %.2f sec/file." % ((time.time()-scanStart)/len(fileList)))
    print("  m2s238 evts:",evtCtr, "total evts:",totCtr, "runTime:",totCalRunTime)


def getHistInfo(x,h):
    """ Computes max, mean, width, percentiles of a numpy
    array based histogram , w/ x values 'x' and counts 'h'. """
    if np.sum(h)==0:
        return 0, 0, 0, [0,0,0,0], 0

    max = x[np.argmax(h)]
    avg = np.average(x, weights=h/np.sum(h))
    std = np.sqrt(np.average((h-max)**2, weights=h)/np.sum(h))
    pct = []
    for p in [5, 10, 90, 95]:
        tmp = np.cumsum(h)/np.sum(h)*100
        idx = np.where(tmp > p)
        pct.append(x[idx][0])
    wid = pct[2]-pct[0]
    return max, avg, std, pct, wid


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

        max1, avg1, std1, pct1, wid1 = getHistInfo(x,h1)
        max2, avg2, std2, pct2, wid2 = getHistInfo(x,h2)
        max3, avg3, std3, pct3, wid3 = getHistInfo(x,h3)

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


if __name__=="__main__":
    main(sys.argv[1:])