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
    # plotSpec()
    # plotHitSpec()
    plotTest()

    # findPeaks, < use sumList
    # plotPeaks,
    # plotBkgPeaks


def getSpec():
    """ Get sum and hit spectra w/ threshold cut, minus ch. 598 (it's noisy.)
    Also w/ !EventDC1Bits and trapENFCal > 0.7, all hits.  (Skim file was generated w/ 'dontSkipAnything')
    """
    from ROOT import TFile, TTree

    fileList = []
    calInfo = ds.CalInfo()
    runList = calInfo.GetSpecialRuns("longCal",5) # 54 total
    runList = runList[:3]
    # runList = runList[:15] # this gets ~10k mH=4 events
    # runList = runList[:] # histats file, see filename below
    for run in runList: fileList.extend(ds.getLATRunList([run],"%s/lat" % (ds.specialDir)))
    nFiles = len(fileList)
    runTime = getRunTime(runList)

    chList = ds.GetGoodChanListNew(dsNum=5)

    sumList, sumList50, hitList50 = [], [], []

    for iFile, f in enumerate(fileList):
        print("%d/%d %s" % (iFile,nFiles,f))
        # tf = TFile("%s/lat/%s" % (ds.specialDir,f))
        tf = TFile("%s/%s" % (os.environ['SLURM_TMP'], f))
        lTree = tf.Get("skimTree")

        bnames = ["startClockTime_s","clockTime_s","EventDC1Bits","isGood",
        "gain","mH","channel","trapENFCal","dtPulserCard","fitSlo","riseNoise"]
        lTree.SetBranchStatus('*',0)
        for name in bnames: lTree.SetBranchStatus(name,1)

        lTree.GetEntry(0)
        startTime = lTree.startClockTime_s
        dt, prev = {}, {}

        for iEnt in range(lTree.GetEntries()):
            lTree.GetEntry(iEnt)
            if lTree.EventDC1Bits != 0: continue

            # get iteration numbers for hits passing cuts
            hitList = [i for i in range(lTree.channel.size())
                if lTree.channel.at(i) in chList
                and lTree.channel.at(i) != 598
                and lTree.channel.at(i) in thD.keys()
                and lTree.trapENFCal.at(i) > thD[lTree.channel.at(i)][0] + 3*thD[lTree.channel.at(i)][1]
                and lTree.trapENFCal.at(i) < 9999
                and lTree.trapENFCal.at(i) > 0.7
                ]
            if len(hitList) == 0: continue

            hitE = np.asarray([lTree.trapENFCal.at(i) for i in hitList])
            chan = np.asarray([lTree.channel.at(i) for i in hitList])
            fSlo = np.asarray([lTree.fitSlo.at(i) for i in hitList])
            rise = np.asarray([lTree.riseNoise.at(i) for i in hitList])
            dtpc = np.asarray([lTree.dtPulserCard.at(i) for i in hitList])

            sumEThr = sum(hitE)
            mHT = len(hitE)
            mH = lTree.mH
            sumE = lTree.sumEH

            # calculate dtmH (based on mH. could change to mHT)
            clockTime = lTree.clockTime_s
            evtTime = clockTime-startTime
            if mH not in prev:
                prev[mH], dt[mH] = 0, -1
            if prev[mH] != 0:
                dt[mH] = evtTime - prev[mH]
            prev[mH] = evtTime

            sumList.append([mH, mHT, sumE, sumEThr, dt[mH], iFile, iEnt])

            # idx = np.where(hitE < 50)
            # if len(hitE[idx]) > 0:
            sumList50.append([mH, mHT, sumE, sumEThr, dt[mH], iFile, iEnt])
            hitList50.append([hitE, chan, fSlo, rise, dtpc])

    np.savez("../plots/mult3-sumE-test.npz", runTime, sumList)
    np.savez("../plots/mult3-hit50List-test.npz", runTime, sumList50, hitList50)

    # np.savez("../plots/mult3-sumE.npz", runTime, sumList)
    # np.savez("../plots/mult3-sumE-histats.npz", runTime, sumList)


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


def plotSpec():
    """ sumList: (mHT, sumEThr, dt[mHT])
    NEW FORMAT: [mH, mHT, sumE, sumEThr, dt[mH], iFile, iEnt]
    Plot sum spectrum, multiplicity counts, and dt[mHT] plots.
    """
    f = np.load("../plots/mult3-sumE.npz")
    # f = np.load("../plots/mult3-sumE-histats.npz")
    runTime, sumList = f['arr_0'], f['arr_1']

    # # counts & rates for each multiplicity
    # multHits = [0]
    # for mH in range(1,7):
    #     nM = sum([evt[0] for evt in sumList if evt[0]==mH])
    #     multHits.append(nM)
    # multErr = [np.sqrt(multHits[mH]) for mH in range(7)]
    # multPct = [100/multErr[mH] if multHits[mH] > 0 else 0 for mH in range(7)]
    # rMult = [multHits[i]/runTime for i in range(7)]
    # rErr = [multErr[i]/runTime for i in range(7)]
    # # for i in range(1,7):
    #     # print("mH %d  cts %d  err %.2f  pct %.2f  rate %.3f  err %.3f" % (i, multHits[i], multErr[i], multPct[i], multHits[i]/runTime, multErr[i]/runTime))
    #
    # # sum spectrum
    # fig = plt.figure()
    # xLo, xHi, xpb = 0, 4000, 2
    # cols = [0,'r','b','m','c','g','k']
    # for mH in range(1,7):
    #     sumE = [evt[1] for evt in sumList if evt[0]==mH]
    #     plt.semilogy(*wl.GetHisto(sumE, xLo, xHi, xpb), ls='steps', lw=1.5, c=cols[mH],
    #         label=r'mH=%d %.2f $\pm$ %.3f Hz' % (mH, rMult[mH], rErr[mH]))
    # plt.xlabel("sumE (keV)", ha='right', x=1.)
    # plt.ylabel("Counts / %.1f keV" % (xpb), ha='right', y=1.)
    # plt.legend(fontsize=12)
    # plt.savefig("../plots/mult3-sumSpec.png")
    #
    # # delta-t plots
    # plt.cla()
    #
    # fig2 = plt.figure(figsize=(18,12))
    # p = [0,0,0,0,0]
    # p[1] = plt.subplot(221) # two columns (slides)
    # p[2] = plt.subplot(222)
    # p[3] = plt.subplot(223)
    # p[4] = plt.subplot(224)
    # ranges = [0, (0.00001,0.004,0.00001), (0.0001,0.03,0.0001), (0.001,0.3,0.001), (0.01,3.,0.01)]
    #
    # for i in range(1,5):
    #     xLo, xHi, xpb = ranges[i][0], ranges[i][1], ranges[i][2]
    #     dtv = [evt[2] for evt in sumList if evt[0]==i]
    #
    #     x, y = wl.GetHisto(dtv, xLo, xHi, xpb)
    #     # norm = np.sum(y)
    #     # y = np.divide(y, norm)
    #     p[i].plot(x, y, c=cols[i], lw=1., ls='steps', label='mH=%d, Raw: %.2f Hz' % (i, rMult[i]))
    #     popt, pcov = curve_fit(expoDist, x, y, p0=(rMult[i],1))
    #     perr = np.sqrt(np.diag(pcov))
    #     rateFit, rateErr = popt[0], perr[0]
    #     p[i].plot(x, expoDist(x, *popt), '-b', lw=1.5, label=r'Fit: %.2f$\pm$%.2f Hz ' % (rateFit, rateErr))
    #     p[i].set_ylabel("Cts / %.0e sec" % xpb, ha='right', y=1.)
    #     p[i].legend()
    #     if i==1:
    #         p[i].ticklabel_format(axis='x',style='sci',scilimits=(1,2))
    #         p[i].set_xlabel("sec       ", ha='right', x=1.)
    #     else:
    #         p[i].set_xlabel("sec", ha='right', x=1.)
    #
    # plt.savefig("../plots/mult3-expoFit.png")


def plotTest():
    """
    sumList   : [mH, mHT, sumE, sumEThr, dt[mH], iFile, iEnt]
    sumList50 : [mH, mHT, sumE, sumEThr, dt[mH], iFile, iEnt]
    hitList : [hitE, chan, fSlo, rise, dtpc]  (same length as sumList50)
    """
    f1 = np.load("../plots/mult3-sumE-test.npz")
    f2 = np.load("../plots/mult3-hit50List-test.npz")
    runTime, sumList = f1['arr_0'], f1['arr_1']
    sumList50, hitList = f2['arr_1'], f2['arr_2']

    mH = 2

    sumE = [evt[3] for evt in sumList if evt[0]==mH]
    sumE50 = [evt[3] for evt in sumList50 if evt[0]==mH]

    xLo, xHi, xpb = 0, 3000, 20
    plt.semilogy(*wl.GetHisto(sumE,xLo,xHi,xpb),'r',lw=2.,ls='steps',label='sumList')

    for cut, col in [(500,'b'), (250,'g'), (100,'m'), (50,'k')]:
        hitECut = [sum(hitList[i][0]) for i in range(len(hitList)) if sumList50[i][0]==mH and (any(i < cut for i in hitList[i][0]))]
        plt.semilogy(*wl.GetHisto(hitECut,xLo,xHi,xpb),c=col,lw=1.,ls='steps',label='hitList %d' % cut)

    # for thresh in [500, 200, 100, 50]:
    #     print(thresh)
    #     for i in range(len(sumList50)):
    #         hitE50 = []
    #         if sumList50[i][0] == mH:
    #             iFile, iEnt = sumList50[i][5], sumList50[i][6]
    #             sumE1 = sumList50[i][3]
    #             hitE = hitList[i][0]
    #             idx = np.where(hitE < thresh)
    #             if len(hitE[idx]) > 0:
    #                 hitE50.append(sum(hitE))
    #             # print("iE %d  iF %d  mH %d  sumE %-8.2f  hSum %-8.2f  u50 %d  hits" % (iEnt, iFile, mH, sumE1, hitSum, under50), hitE)
    #             # if i > 100: break
    #             plt.semilogy(*wl.GetHisto(hitE50,xLo,xHi,xpb),lw=1.,ls='steps',label='hitList')

    plt.legend()
    plt.savefig("../plots/mult3-sumSpec-test.png")



# def retuneFitSlo():



if __name__=="__main__":
    main()