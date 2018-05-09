#!/usr/bin/env python
import sys, os, imp, glob
import numpy as np
from scipy.optimize import curve_fit

import matplotlib as mpl
# mpl.use('Agg')
# sys.argv.append("-b")
import matplotlib.pyplot as plt
plt.style.use('../pltReports.mplstyle')
from matplotlib.colors import LogNorm, Normalize

# load LAT libraries
ds = imp.load_source('DataSetInfo',os.environ['LATDIR']+'/DataSetInfo.py')
wl = imp.load_source('waveLibs',os.environ['LATDIR']+'/waveLibs.py')
calInfo = ds.CalInfo()

# load threshold data
import tinydb as db
dsNum, bkgIdx = 5, 83
calDB = db.TinyDB('../calDB.json')
pars = db.Query()
thD = ds.getDBRecord("thresh_ds%d_bkgidx%d" % (dsNum, bkgIdx), False, calDB, pars)

def main():

    # getSpec()
    plotSumHitCut()
    plotSpec()

    # plotHitSpec() < move to mult4.  everything should use mHT and sumET!
    # findPeaks,
    # plotPeaks,
    # plotBkgPeaks
    # retuneFitSlo()


def getSpec():
    """ Get sum and hit spectra w/ threshold cut, without ch. 598 (it's noisy.)
    Also w/ !EventDC1Bits and trapENFCal > 0.7, all hits.  (Skim file was generated w/ 'dontSkipAnything')
    Save detailed info for events with hits < 100 keV.
    """
    from ROOT import TFile, TTree

    fileList = []
    calInfo = ds.CalInfo()
    runList = calInfo.GetSpecialRuns("longCal",5) # 54 total
    # runList = runList[:3]
    runList = runList[:15] # this gets ~10k mH=4 events
    # runList = runList[:] # histats file, see filename below
    for run in runList: fileList.extend(ds.getLATRunList([run],"%s/lat" % (ds.specialDir)))
    nFiles = len(fileList)
    runTime = getRunTime(runList)

    chList = ds.GetGoodChanListNew(dsNum=5)

    eCut = 250
    hitList, hitData = [], []
    xLo, xHi, xpb = 0, 4000, 1
    nb = int((xHi-xLo)/xpb)+1
    sumTSpec = {mHT:np.zeros(nb) for mHT in range(7)}
    sumSpec = {mH:np.zeros(nb) for mH in range(7)}
    dtVals = {mHT:[] for mHT in range(7)}

    for iFile, f in enumerate(fileList):
        print("%d/%d %s" % (iFile,nFiles,f))
        # tf = TFile("%s/lat/%s" % (ds.specialDir,f))
        tf = TFile("%s/%s" % (os.environ['SLURM_TMP'], f))
        lTree = tf.Get("skimTree")

        bnames = ["startClockTime_s","clockTime_s","EventDC1Bits","sumEH",
        "mH","channel","trapENFCal","dtPulserCard","fitSlo","riseNoise"]
        lTree.SetBranchStatus('*',0)
        for name in bnames: lTree.SetBranchStatus(name,1)

        lTree.GetEntry(0)
        startTime = lTree.startClockTime_s
        dt, prev = {}, {}
        sumTVals = {mHT:[] for mHT in range(7)}
        sumVals = {mH:[] for mH in range(7)}

        for iEnt in range(lTree.GetEntries()):
            lTree.GetEntry(iEnt)
            if lTree.EventDC1Bits != 0: continue

            # get iteration numbers for hits passing cuts
            idxList = [i for i in range(lTree.channel.size())
                if lTree.channel.at(i) in chList
                and lTree.channel.at(i) != 598
                and lTree.channel.at(i) in thD.keys()
                and lTree.trapENFCal.at(i) > thD[lTree.channel.at(i)][0] + 3*thD[lTree.channel.at(i)][1]
                and lTree.trapENFCal.at(i) < 9999
                and lTree.trapENFCal.at(i) > 0.7
                ]
            if len(idxList) == 0: continue
            hitE = np.asarray([lTree.trapENFCal.at(i) for i in idxList])

            # save multiplicity and sum energy
            mH = lTree.mH
            mHT = len(hitE)
            sumE = lTree.sumEH
            sumET = sum(hitE)
            if mH < 7: sumVals[mH].append(sumE)
            if mHT < 7: sumTVals[mHT].append(sumET)

            # calculate dt since last m=N event
            # Tried w/ mH - fitted rates were lower than raw rates (could be from removing 592)
            # Now trying w/ mHT ...
            clockTime = lTree.clockTime_s
            evtTime = clockTime-startTime
            if mHT not in prev:
                prev[mHT], dt[mHT] = 0, -1
            if prev[mHT] != 0:
                dt[mHT] = evtTime - prev[mHT]
            prev[mHT] = evtTime
            if mHT < 7: dtVals[mHT].append(dt[mHT])

            # save all information for events that have at least one hit under 'eCut'.
            if mHT > 1 and len(hitE[np.where(hitE < eCut)]) > 0:
                chan = np.asarray([int(lTree.channel.at(i)) for i in idxList])
                fSlo = np.asarray([lTree.fitSlo.at(i) for i in idxList])
                rise = np.asarray([lTree.riseNoise.at(i) for i in idxList])
                dtpc = np.asarray([lTree.dtPulserCard.at(i) for i in idxList])

                hitData.append([mH, mHT, sumE, sumET, dt[mHT]])#, iFile, iEnt])
                hitList.append([hitE, chan, fSlo, rise, dtpc])

        # save sum spectra
        for i in range(7):
            x, y = wl.GetHisto(sumTVals[i], xLo, xHi, xpb)
            sumTSpec[i] = np.add(sumTSpec[i], y)
            x, y = wl.GetHisto(sumVals[i], xLo, xHi, xpb)
            sumSpec[i] = np.add(sumSpec[i], y)

    np.savez("../data/mult3-sumE.npz", runTime, x, sumSpec, sumTSpec)
    np.savez("../data/mult3-dtmHT.npz", runTime, dtVals)
    np.savez("../data/mult3-hitE.npz", runTime, hitList, hitData, eCut)


def getRunTime(runList):
    from ROOT import TFile, TTree, GATDataSet
    runTime = 0
    for run in runList:
        gds = GATDataSet()
        gatPath = gds.GetPathToRun(run, GATDataSet.kGatified)
        tf = TFile(gatPath)
        gatTree = tf.Get("mjdTree")
        gatTree.GetEntry(0)
        runTime += gatTree.timeinfo.stopTime - gatTree.timeinfo.startTime
        tf.Close()
        print(run, runTime)
    return runTime


def expoDist(t, rate, amp):
    # from scipy.stats import poisson
    # return amp * poisson.pmf(0, rate * t)
    return amp * np.exp(-1*rate*t)


def plotSumHitCut():
    """ Keeping only events w/ hits > some threshold has some strange effects on the higher sum peaks.
    Just making sure I didn't mess up any code.

    To illustrate: plot the hitE spec below with smaller values of eCut.

    sumSpec : {mH: histos}
    hitData : [mH, mHT, sumE, sumET, dt[mH], iFile, iEnt]
    hitList : [hitE, chan, fSlo, rise, dtpc]  (same length as hitData)
    """
    f1 = np.load("../data/mult3-sumE.npz")
    f2 = np.load("../data/mult3-hitE.npz")
    runTime, x, sumSpec, sumSpecT = f1['arr_0'], f1['arr_1'], f1['arr_2'].item(), f1['arr_3'].item()
    hitList, hitData, eCut = f2['arr_1'], f2['arr_2'], f2['arr_3']

    mH = 2
    # cut = eCut
    cut = 50

    xLo, xHi, xpb = 0, 4000, 1
    plt.semilogy(x, sumSpec[mH],'r',lw=1.,ls='steps',label='sumEH')
    plt.semilogy(x, sumSpecT[mH],'b',lw=1.,ls='steps',label='sumET')

    n = len(hitList)
    hitE = [sum(hitList[i][0]) for i in range(n) if hitData[i][1]==mH and (any(i < cut for i in hitList[i][0]))]
    plt.semilogy(*wl.GetHisto(hitE,xLo,xHi,xpb),c='m',lw=1.,ls='steps',label='hitList %d' % cut)

    plt.xlabel("sumE (keV)", ha='right', x=1.)
    plt.legend(loc=1, fontsize=14)
    plt.savefig("../plots/mult3-sumECut.png")


def plotSpec():
    """ Plot sum spectrum, multiplicity counts, and dt[mHT] plots.
    dtVals  : {mH : [dt vals]}
    hitData : [mH, mHT, sumE, sumET, dt[mH], iFile, iEnt]
    hitList : [hitE, chan, fSlo, rise, dtpc]  (same length as hitData)
    """
    f1 = np.load("../data/mult3-hitE.npz")
    f2 = np.load("../data/mult3-dtmHT.npz")
    f3 = np.load("../data/mult3-sumE.npz")
    runTime, hitList, hitData, eCut = f1['arr_0'], f1['arr_1'], f1['arr_2'], f1['arr_3']
    runTime, dtVals = f2['arr_0'], f2['arr_1'].item()
    runTime, x, sumSpec, sumSpecT = f3['arr_0'], f3['arr_1'], f3['arr_2'].item(), f3['arr_3'].item()

    # counts & rates for each multiplicity (use mHT)
    nHits = [sum(sumSpecT[i]) for i in range(7)]
    nErr = [np.sqrt(nHits[i]) for i in range(7)]
    nPct = [100 / nErr[i] if nHits[i] > 0 else 0 for i in range(7)]
    rate = [nHits[i] / runTime for i in range(7)]
    rErr = [nErr[i] / runTime for i in range(7)]

    # sum spectrum
    fig = plt.figure()
    xLo, xHi, xpb = 0, 4000, 2
    cols = [0,'r','b','m','c','g','k']
    for mH in range(1,7):
        pLabel = r'mHT=%d %.2f $\pm$ %.3f Hz' % (mH,rate[mH],rErr[mH])
        plt.semilogy(x, sumSpecT[mH], ls='steps', lw=1.5, c=cols[mH],label=pLabel)
    plt.xlabel("sumE (keV)", ha='right', x=1.)
    plt.ylabel("Counts / %.1f keV" % (xpb), ha='right', y=1.)
    plt.legend(fontsize=12)
    plt.savefig("../plots/mult3-sumSpec.png")

    # delta-t plots.  some question as to use mHT or mH ...

    fig2 = plt.figure(figsize=(18,12))
    p = [0,0,0,0,0]
    p[1] = plt.subplot(221) # two columns (slides)
    p[2] = plt.subplot(222)
    p[3] = plt.subplot(223)
    p[4] = plt.subplot(224)
    ranges = [0, (0,0.005,0.000005), (0,0.05,0.00005), (0,0.5,0.0005), (0,5.,0.005)]

    for i in range(1,5):
        xLo, xHi, xpb = ranges[i][0], ranges[i][1], ranges[i][2]
        x, y = wl.GetHisto(dtVals[i], xLo, xHi, xpb)
        # y = np.divide(y, np.sum(y))

        p[i].plot(x, y, c=cols[i], lw=1., ls='steps', label='mH=%d, Raw %.2f Hz' % (i, rate[i]))

        g = (rate[i], np.mean(y[2:5])*1.3)
        popt, pcov = curve_fit(expoDist, x, y, p0=g)
        perr = np.sqrt(np.diag(pcov))
        rateFit, rateErr = popt[0], perr[0]
        p[i].plot(x, expoDist(x, popt[0], popt[1]), '-b', lw=1.5, label=r'Fit: %.2f$\pm$%.3f Hz' % (rateFit, rateErr))

        p[i].set_ylabel("Cts / %.0e sec" % xpb, ha='right', y=1.)
        p[i].legend()
        if i==1:
            p[i].ticklabel_format(axis='x',style='sci',scilimits=(1,2))
            p[i].set_xlabel("sec         ", ha='right', x=1.)
        else:
            p[i].set_xlabel("sec", ha='right', x=1.)

    plt.savefig("../plots/mult3-expoFit-mH.png")




if __name__=="__main__":
    main()