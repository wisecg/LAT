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

        # scan cal runs and make some outputs
        if opt=="-scan":
            loadRuns(ds,cIdx,mod)

        if opt=="-p":
            plotHist()


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

    fileList = []
    calRuns = cal.GetCalList(key,cIdx)
    for run in calRuns:
        latList = dsi.getSplitList("%s/latSkimDS%d_run%d*" % (dsi.calLatDir, ds, run), run)
        tmpList = [f for idx, f in sorted(latList.items())]
        fileList.extend(tmpList)

    print("Scanning DS:%d  calIdx %d  mod %d  key %s  nFiles:%d" % (ds, cIdx, mod, key, len(fileList)), time.strftime('%X %x %Z'))

    # load good channel list for this DS
    chList = det.getGoodChanList(ds)

    # declare the output stuff
    evtIdx, evtSumET, evtHitE, evtChans = [], [], [], []
    thrCal = {ch:[] for ch in chList}

    # store histograms for fitSlo in every channel,
    # slicing on energy ranges 0-10 and 10-200 keV
    # {channel : [hist 0-10, hist 10-200] }
    fLo, fHi, fpb = -200, 400, 10
    nbf = int((fHi-fLo)/fpb)+1
    fSloSpec = {ch:[np.zeros(nbf) for i in range(2)] for ch in chList}
    # fSloSpec = {ch:[] for ch in chList}
    # fSloX = np.arange(fLo-fpb/2.,fHi,fpb)

    print(nbf, len(fSloX))
    # return

    # debug: fs and ene arrays, to verify the histo slice is working
    fsTot, enTot = [], []


    # loop over LAT cal files
    scanStart = time.time()
    prevRun = 0
    evtCtr, totCtr, runTime = 0, 0, 0
    for iF, f in enumerate(fileList[:3]):

        print("%d/%d %s" % (iF, len(fileList), f))
        tf = TFile(f)
        tt = tf.Get("skimTree")

        # histogram these for each file
        fs10 = {ch:[] for ch in chList}
        # en10 = {ch:[] for ch in chList}
        fs200 = {ch:[] for ch in chList}
        # en200 = {ch:[] for ch in chList}

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

            # get indexes of hits above threshold (use thresholds from THIS CAL RUN)
            idxList = [i for i in range(tt.channel.size())
                if tt.channel.at(i) in chList
                and tt.trapENFCal.at(i) > tmpThresh[tt.channel.at(i)][3]
                and 0.7 < tt.trapENFCal.at(i) < 9999
                ]
            hitE = np.asarray([tt.trapENFCal.at(i) for i in idxList])
            mHT, sumET = len(hitE), sum(hitE)

            # increment histogram of fitSlo for this channel
            for i in idxList:
                en = tt.trapENFCal.at(i)
                ch = tt.channel.at(i)
                fs = tt.fitSlo.at(i)
                if fLo < fs < fHi:
                    if ch==610:
                        fsTot.append(fs) # TODO: remove this
                        enTot.append(en)
                    if 0 < en < 10:
                        fs10[ch].append(fs)
                        en10[ch].append(en)
                    elif 10 < en < 200:
                        fs200[ch].append(fs)
                        en200[ch].append(en)

            # Save m2s238 events to output, skip everything else
            if mHT!=2: continue
            if not 237.28 < sumET < 239.46: continue
            hitChans = np.asarray([tt.channel.at(i) for i in idxList])
            evtIdx.append([run,iE])
            evtSumET.append(sumET)
            evtHitE.append(hitE)
            evtChans.append(hitChans)
            evtCtr += 1

        # fill the fitSlo histograms w/ the events from this run
        for ch in chList:
            x1, h10 = wl.GetHisto(fs10[ch],fLo,fHi,fpb)
            x2, h200 = wl.GetHisto(fs200[ch],fLo,fHi,fpb)
            fSloSpec[ch][0] = np.sum([fSloSpec[ch][0], h10], axis=0)
            fSloSpec[ch][1] = np.sum([fSloSpec[ch][1], h200], axis=0)

    # get average threshold for each channel
    thrFinal = {cpd:[] for cpd in thrCal}
    for cpd in thrCal:
        thrVals = []
        for iT in range(len(thrCal[cpd])):
            run, thrM, thrS, thrK = thrCal[cpd][iT]
            # print("%d  %d  %.3f  %.3f  %.3f" % (cpd,run,thrM,thrS,thrK))
            thrVals.append(thrK)
        thrVals = np.asarray(thrVals)
        thrAvg = np.mean(thrVals)
        thrDev = np.std(thrVals)
        # print("%d  %.3f  %.3f" % (cpd, thrAvg, thrDev))
        thrFinal[cpd] = [thrAvg,thrDev]

    print("Detector Thresholds:")
    for cpd in thrFinal:
        thKeV = thrFinal[cpd][0]
        thE = thrFinal[cpd][1]

        errString = ""
        if thE/thKeV > 0.5:
            errString = ">50pct error:  thE/thKeV=%.2f" % (thE/thKeV)
        if thKeV > 2:
            errString = ">2kev threshold"

        print("%d  %.3f  %.3f  %s" % (cpd,thKeV,thE,errString))

    # TODO: calculate mu, sig of fitSlo for each channel in each slice

    # output stats
    print("Done:",time.strftime('%X %x %Z'),", %.2f sec/file." % ((time.time()-scanStart)/len(fileList)))
    print("  m2s238 evts:",evtCtr, "total evts:",totCtr, "runTime:",runTime)

    # save output
    outFile = "%s/eff_%s_c%d.npz" % (dsi.effDir, key, cIdx)
    print("Saving output in:",outFile)
    np.savez(outFile, evtIdx, evtSumET, evtHitE, evtChans, thrCal, thrFinal, evtCtr, totCtr, runTime, fSloSpec, fSloX)

    # debug output
    np.savez("./plots/histo-tmp.npz",fsTot,enTot)


def plotHist():

    f1 = np.load("/global/projecta/projectdirs/majorana/users/wisecg/cal/eff/eff_ds1_m1_c1.npz")
    fSloSpec = f1['arr_9'].item()
    fSloX = f1['arr_10']

    print(list(fSloSpec.keys()))

    myHist1 = fSloSpec[610][0] # 0-10 keV
    myHist2 = fSloSpec[610][1] # 10-200 keV

    # myHist1 = np.insert(myHist1, 0, 0, axis=0)
    # myHist2 = np.insert(myHist2, 0, 0, axis=0)

    # if shift: x = x-xpb/2.

    f2 = np.load("./plots/histo-tmp.npz")
    fsTot, enTot = f2['arr_0'], f2['arr_1']

    fig = plt.figure()
    plt.plot(enTot, fsTot, ".k", ms=1)
    plt.xlim(0,250)
    plt.ylim(-50, 400)
    plt.xlabel("Energy",ha='right',x=1)
    plt.ylabel("fitSlo",ha='right',y=1)
    plt.savefig("./plots/lat2-histoTest.png")

    plt.cla()

    fLo, fHi, fpb = -200, 400, 10
    idx = np.where((enTot>0) & (enTot<10))

    print(len(fsTot[idx]), np.sum(myHist1))

    x, histLo = wl.GetHisto(fsTot[idx],fLo,fHi,fpb)
    plt.plot(x, histLo,'b',lw=2.,ls='steps',label='fsTot 0-10')
    plt.plot(fSloX, myHist1,'r',lw=2.,ls='steps',label='myHist 0-10')
    plt.xlabel("fitSlo",ha='right',x=1)
    plt.xlim(50,150)
    plt.legend(loc=1)
    plt.savefig('./plots/lat2-histoTest2.png')

    plt.cla()
    idx = np.where((enTot>10) & (enTot<200))
    x, histHi = wl.GetHisto(fsTot[idx],fLo,fHi,fpb)
    plt.plot(x, histHi,'b',lw=2.,ls='steps',label='fsTot 10-200')
    plt.plot(fSloX, myHist2,'r',lw=2.,ls='steps',label='myHist 10-200')
    plt.xlim(50,150)
    plt.xlabel("fitSlo",ha='right',x=1)
    plt.legend(loc=1)
    plt.savefig('./plots/lat2-histoTest3.png')

if __name__=="__main__":
    main(sys.argv[1:])