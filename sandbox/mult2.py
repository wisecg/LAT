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
ds = imp.load_source('DataSetInfo',os.environ['LATDIR']+'/sandbox/DataSetInfo.py')
wl = imp.load_source('waveLibs',os.environ['LATDIR']+'/waveLibs.py')
calInfo = ds.CalInfo()

# load threshold data
import tinydb as db
dsNum, bkgIdx = 5, 83
calDB = db.TinyDB('../calDB.json')
pars = db.Query()
thD = ds.getDBRecord("thresh_ds%d_bkgidx%d" % (dsNum, bkgIdx), False, calDB, pars)

def main():

    # getSumSpec()
    # plotSumSpec()
    # plotHitSpec()
    # plotMultipRates()

    # getLoHitSpec()
    # plotLoHitSpec()

    # getHiDT() # <- deprecated
    # plotHiDT()

    # getDT() # <- use this one
    # plotFitRates()
    # plotDT_mH1Ene()
    # plotDT_mH1Ene_wThr()
    # plotDT_mH1Low()
    # plotThreshFunc()
    # plotDT_mH1Low_wThr()
    # plotDT_wThr()
    # plotChannelRates_mH1()
    # plotDTCut()
    # plotPCut()

    # noiseProbability()
    # noiseProbability2D()

    # getPeaks()
    # plotPeaks()

    # getSumEvents()
    # plotComptonEdge()
    # plotLoHits()
    plotSloHits()
    # plotChannelEff()

    # getExtPulser()
    # plotExtPulser()


# def unpackData():
#     cd $SLURM_TMP
#     cp /global/project/projectdirs/majorana/users/wisecg/special/lat/latLongCal.tar.gz .
#     tar -zxvf latLongCal.tar.gz


def getSumE(tree, tNames, theCut):
    """ Get sum energy for all events passing cuts """
    sumArr = []
    tvals = wl.GetVX(tree, tNames, theCut, False)
    nPass = len(tvals["Entry$"])
    prevEnt, sumE = -1, 0.

    for idx in range(nPass):
    # for idx in range(20):
        ent = tvals["Entry$"][idx]
        mH = tvals["mH"][idx]
        chan = tvals["channel"][idx]
        enf = tvals["trapENFCal"][idx]
        gatSumE = tvals["sumEH"][idx] # equals sumEHClean when isGood is applied in cut
        if enf > 99999: enf = 0
        if ent!=prevEnt and idx!=0:
            # print("app %.2f" % sumE)
            sumArr.append(sumE)
            sumE = 0
        sumE += enf
        prevEnt = ent
        # print("%d  %d  %d  %.2f  lat %2f  gat %.2f" % (ent,chan,mH,enf,sumE,gatSumE))

    sumArr.append(sumE) # append last one
    # print("app %.2f" % sumE)

    sumArr = np.asarray(sumArr)
    return sumArr


def getSumSpec():
    from ROOT import TFile, TTree, GATDataSet

    fileList = []
    calInfo = ds.CalInfo()
    runList = calInfo.GetSpecialRuns("longCal",5) # 54 total
    # runList = runList[:1]
    # runList = runList[:15] # this gets ~10k mH=4 events
    for run in runList: fileList.extend(ds.getLATRunList([run],"%s/lat" % (ds.specialDir)))
    nFiles = len(fileList)

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

    eLo, eHi, epb = 0., 4000., 2.
    nbe = int((eHi-eLo)/epb) + 1
    sumSpec = [np.zeros(nbe) for i in range(0,6)]
    hitSpec = [np.zeros(nbe) for i in range(0,6)]

    mLo, mHi = 1, 39
    mult = np.zeros(mHi)

    for iFile, f in enumerate(fileList):
        print("%d/%d %s" % (iFile,nFiles,f))
        # tf = TFile("%s/lat/%s" % (ds.specialDir,f))
        tf = TFile("%s/%s" % (os.environ['SLURM_TMP'], f))
        lTree = tf.Get("skimTree")

        # activate only branches we want
        bNames = ["Entry$", "mH", "gain", "isGood", "EventDC1Bits", "trapENFCal","sumEH","channel"]
        lTree.SetBranchStatus('*',0)
        for name in bNames[1:]:
            lTree.SetBranchStatus(name,1)

        # hit and sum spectra, multiplicity
        for i in range(1,6):
            theCut = "mH==%d && gain==0 && isGood && !EventDC1Bits" % i

            n = lTree.Draw("trapENFCal", theCut, "GOFF")
            hitE = lTree.GetV1()
            hitE = [hitE[i] for i in range(n)]
            x, y = wl.GetHisto(hitE, eLo, eHi, epb)
            hitSpec[i] = np.add(hitSpec[i], y)

            sumArr = getSumE(lTree, bNames, theCut)
            x, y = wl.GetHisto(sumArr, eLo, eHi, epb)
            sumSpec[i] = np.add(sumSpec[i], y)

            mult[i] += len(sumArr) # bug: used an = instead of += when i generated the npz

        tf.Close()

    # note: "mult" was wrong, do sum(sumSpec[i]) to get the num. cts
    np.savez("../data/mult2-sumSpec.npz", x, sumSpec, hitSpec, mult, runTime)


def plotSumSpec():

    f = np.load("../data/mult2-sumSpec.npz")
    x, sumSpec, hitSpec, runTime = f['arr_0'], f['arr_1'], f['arr_2'], f['arr_4']
    fig = plt.figure()

    # calculate rates
    rm = [(0,0)]
    for i in range(1,6):
        sumCts = sum(sumSpec[i])
        hitCts = sum(hitSpec[i])
        rate = sumCts/runTime
        error = np.sqrt(sumCts)/runTime
        rm.append((rate,error))
        # print("rt %d sec  mH %d  sumCts %d  hitCts %d  rate %.2f pm %.3f" % (runTime, i, sumCts, hitCts, rate, error))

    plt.semilogy(x, sumSpec[1], ls='steps', c='r', lw=1.0, alpha=0.6, label=r"mH=1 %.2f$\pm$%.3f $s^{-1}$ " % (rm[1][0],rm[1][1]))
    plt.semilogy(x, sumSpec[2], ls='steps', c='b', lw=0.9, alpha=0.7, label=r"mH=2 %.2f$\pm$%.3f $s^{-1}$ " % (rm[2][0],rm[2][1]))
    plt.semilogy(x, sumSpec[3], ls='steps', c='m', lw=0.8, alpha=0.8, label=r"mH=3 %.2f$\pm$%.3f $s^{-1}$ " % (rm[3][0],rm[3][1]))
    plt.semilogy(x, sumSpec[4], ls='steps', c='c', lw=0.7, alpha=0.9, label=r"mH=4 %.2f$\pm$%.3f $s^{-1}$ " % (rm[4][0],rm[4][1]))
    plt.semilogy(x, sumSpec[5], ls='steps', c='g', lw=0.6, alpha=1.0, label=r"mH=5 %.2f$\pm$%.3f $s^{-1}$ " % (rm[5][0],rm[5][1]))

    plt.xlabel("sumE (keV)", ha='right', x=1.)
    plt.ylabel("Counts / 2 keV", ha='right', y=1.)
    plt.legend(fontsize=10)
    plt.savefig("../plots/mult2-sumSpec.png")


def plotHitSpec():
    f = np.load("../data/mult2-sumSpec.npz")
    x, hitSpec = f['arr_0'], f['arr_2']
    fig = plt.figure()

    plt.semilogy(x, hitSpec[1], ls='steps', c='r', lw=1.0, alpha=0.6, label="mH=1")
    plt.semilogy(x, hitSpec[2], ls='steps', c='b', lw=0.9, alpha=0.7, label="mH=2")
    plt.semilogy(x, hitSpec[3], ls='steps', c='m', lw=0.8, alpha=0.8, label="mH=3")
    plt.semilogy(x, hitSpec[4], ls='steps', c='c', lw=0.7, alpha=0.9, label="mH=4")
    plt.semilogy(x, hitSpec[5], ls='steps', c='g', lw=0.6, alpha=1.0, label="mH=5")

    plt.xlabel("hitE (keV)", ha='right', x=1.)
    plt.ylabel("Counts", ha='right', y=1.)
    plt.legend()
    plt.savefig("../plots/mult2-hitSpec.png")


def plotMultipRates():
    # 'mult' in the orig. npz file was busted.  use sum(sumSpec[i]) to get the num counts.

    f = np.load("../data/mult2-sumSpec.npz")
    sumSpec, runTime = f['arr_1'], f['arr_4']
    fig = plt.figure()

    # calculate rates
    mH = np.arange(1,6,1)
    rates, errors = [], []
    for i in range(1,6):
        sumCts = sum(sumSpec[i])
        rates.append(sumCts/runTime)
        errors.append(np.sqrt(sumCts)/runTime)

    plt.bar(mH, rates, 0.95, color='b', log=True, yerr=errors, label='runTime: %.0f sec' % (runTime))
    plt.xlim(0, 7)
    plt.xlabel("mH", ha='right', x=1.)
    plt.ylabel("Rate (Hz)", ha='right', y=1.)
    plt.legend()
    plt.savefig("../plots/mult2-multip.png")


def getLoHitSpec():
    from ROOT import TFile, TTree, GATDataSet
    fileList = []
    calInfo = ds.CalInfo()
    runList = calInfo.GetSpecialRuns("longCal",5) # 54 total
    # runList = runList[:5]
    runList = runList[:] # hi stats - different filename, see below
    for run in runList: fileList.extend(ds.getLATRunList([run],"%s/lat" % (ds.specialDir)))
    nFiles = len(fileList)

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

    eLo, eHi, epb = 0., 50., 0.1
    nbe = int((eHi-eLo)/epb) + 1
    hitSpec = [np.zeros(nbe) for i in range(0,6)]

    for iFile, f in enumerate(fileList):
        print("%d/%d %s" % (iFile,nFiles,f))
        # tf = TFile("%s/lat/%s" % (ds.specialDir,f))
        tf = TFile("%s/%s" % (os.environ['SLURM_TMP'], f))
        lTree = tf.Get("skimTree")

        # activate only branches we want
        bNames = ["Entry$", "mH", "gain", "isGood", "EventDC1Bits", "trapENFCal"]
        lTree.SetBranchStatus('*',0)
        for name in bNames[1:]:
            lTree.SetBranchStatus(name,1)

        # low-e hit spectra
        for i in range(1,6):
            theCut = "mH==%d && gain==0 && isGood && !EventDC1Bits && trapENFCal < %.1f" % (i, eHi)
            n = lTree.Draw("trapENFCal", theCut, "GOFF")
            hitE = lTree.GetV1()
            hitE = [hitE[i] for i in range(n)]
            x, y = wl.GetHisto(hitE, eLo, eHi, epb)
            hitSpec[i] = np.add(hitSpec[i], y)

    # np.savez("../data/mult2-lowHitSpec.npz", x, hitSpec, runTime)
    np.savez("../data/mult2-lowHitSpec-histats.npz", x, hitSpec, runTime)


def plotLoHitSpec():
    f = np.load("../data/mult2-lowHitSpec-histats.npz")
    x, hitSpec, runTime = f['arr_0'], f['arr_1'], f['arr_2']
    fig = plt.figure()

    plt.semilogy(x, hitSpec[1], ls='steps', c='r', lw=1.0, alpha=0.6, label="mH=1")
    plt.semilogy(x, hitSpec[2], ls='steps', c='b', lw=0.9, alpha=0.7, label="mH=2")
    plt.semilogy(x, hitSpec[3], ls='steps', c='m', lw=0.8, alpha=0.8, label="mH=3")
    plt.semilogy(x, hitSpec[4], ls='steps', c='c', lw=0.7, alpha=0.9, label="mH=4")
    plt.semilogy(x, hitSpec[5], ls='steps', c='g', lw=0.6, alpha=1.0, label="mH=5")

    plt.xlabel("hitE (keV)", ha='right', x=1.)
    plt.ylabel("Counts / 0.1 keV", ha='right', y=1.)
    plt.legend()
    plt.savefig("../plots/mult2-hitSpecLow.png")

    plt.cla()
    idx = np.where(x < 11)
    plt.semilogy(x[idx], hitSpec[1][idx], ls='steps', c='r', lw=1.0, alpha=0.6, label="mH=1")
    plt.semilogy(x[idx], hitSpec[2][idx], ls='steps', c='b', lw=0.9, alpha=0.7, label="mH=2")
    plt.semilogy(x[idx], hitSpec[3][idx], ls='steps', c='m', lw=0.8, alpha=0.8, label="mH=3")
    plt.semilogy(x[idx], hitSpec[4][idx], ls='steps', c='c', lw=0.7, alpha=0.9, label="mH=4")
    plt.semilogy(x[idx], hitSpec[5][idx], ls='steps', c='g', lw=0.6, alpha=1.0, label="mH=5")

    plt.xlabel("hitE (keV)", ha='right', x=1.)
    plt.ylabel("Counts / 0.1 keV", ha='right', y=1.)
    plt.legend()
    plt.savefig("../plots/mult2-hitSpecLow10.png")


def getHiDT():
    """ Cut used in the hit spectrum generation.  Reproduced here.
    theCut = "mH==%d && gain==0 && isGood && !EventDC1Bits && trapENFCal < %.1f" % (i, eHi)
    There are untagged pulsers w/ high dt.  Save data for the energy & dt of these events.

    This one is less important b/c I realized the longCal file was generated w/ dontSkipAnything
    """
    from ROOT import TFile, TTree, GATDataSet
    fileList = []
    calInfo = ds.CalInfo()
    runList = calInfo.GetSpecialRuns("longCal",5) # 54 total
    # runList = runList[:6]
    # runList = runList[:15] # this gets ~10k mH=4 events
    runList = runList[:30] # get big stats (saved as different filename below)
    for run in runList: fileList.extend(ds.getLATRunList([run],"%s/lat" % (ds.specialDir)))
    nFiles = len(fileList)

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

    dtVals = {mH:[] for mH in range(30)}
    highDtVals = {mH:[] for mH in range(30)}

    for iFile, f in enumerate(fileList):
        print("%d/%d %s" % (iFile,nFiles,f))
        # tf = TFile("%s/lat/%s" % (ds.specialDir,f))
        tf = TFile("%s/%s" % (os.environ['SLURM_TMP'], f)) # moved my shit to here
        lTree = tf.Get("skimTree")

        bNames = ["startClockTime_s","clockTime_s","mH","channel","gain","isGood","EventDC1Bits","trapENFCal"]
        lTree.SetBranchStatus('*',0)
        for name in bNames: lTree.SetBranchStatus(name,1)

        lTree.GetEntry(0)
        startTime = lTree.startClockTime_s

        dt, prev = {}, {}
        for iEnt in range(lTree.GetEntries()):
            lTree.GetEntry(iEnt)
            clockTime = lTree.clockTime_s
            nHit = lTree.channel.size()
            mH = lTree.mH
            evtTime = clockTime-startTime
            if lTree.EventDC1Bits != 0: continue

            # calculate dt
            if mH not in prev:
                prev[mH] = 0.
            if prev[mH] != 0:
                dt[mH] = evtTime - prev[mH]
                # print("%-4d  %-3d  evt %-10.9f  dt[%d] %-10.9f  pr[%d] %-10.9f  DCB %d" % (iEnt, mH, evtTime, mH, dt[mH], mH, prev[mH], lTree.EventDC1Bits))
                dtVals[mH].append(dt[mH])

                if dt[mH] > 2:
                    chans = [lTree.channel.at(i) for i in range(nHit) if lTree.gain.at(i)%2==0]
                    hitE = [lTree.trapENFCal.at(i) for i in range(nHit) if lTree.gain.at(i)%2==0]
                    # lTree.isGood.at(0): <ROOT._Bit_reference object at 0x6eee3c0> ... hard to access.
                    highDtVals[mH].append((dt[mH],chans,hitE))

            # else:
                # print("%-4d  %-3d  evt %-10.9f  pr[%d] %-10.9f" % (iEnt, mH, evtTime, mH, prev[mH]))
            prev[mH] = evtTime

    np.savez("../data/mult2-dtVals.npz", dtVals, highDtVals, runTime)
    np.savez("../data/mult2-dtVals-histats.npz", dtVals, highDtVals, runTime)


def expoDist(t, rate, amp):
    # from scipy.stats import poisson
    # return amp * poisson.pmf(0, rate * t)
    return amp * rate * np.exp(-1*rate*t)


def plotHiDT():
    """ This one is less important b/c I realized the longCal file was generated w/ dontSkipAnything """

    # f = np.load("../data/mult2-dtVals.npz")
    f = np.load("../data/mult2-dtVals-histats.npz")
    dtVals = f['arr_0'].item()
    highDtVals = f['arr_1'].item()

    fig = plt.figure()

    xLo, xHi, xpb = 0., 10., 0.02
    cols = [0,'r','b','m','c','g']
    for i in range(1,6):
        plt.semilogy(*wl.GetHisto(dtVals[i], xLo, xHi, xpb), c=cols[i], ls='steps', label='mH=%d'%i)

    plt.legend()
    plt.xlabel("delta-t since last m=N event (sec)", horizontalalignment='right', x=1.0)
    plt.ylabel("Counts (arb)", horizontalalignment='right', y=1.0)
    plt.savefig("../plots/mult2-dt10sec.png")


    plt.cla()
    for mH in range(1,6):
        sumE = [sum(highDtVals[mH][i][2]) for i in range(len(highDtVals[mH]))]
        hiDt = [highDtVals[mH][i][0] for i in range(len(highDtVals[mH]))]

        # get percentage of high-dt events
        x, y = wl.GetHisto(dtVals[mH], xLo, xHi, xpb)
        numTot = np.sum(y)
        hiDtPct = 100*len(hiDt)/numTot
        # print("mH %d  total %d  hiDt %d  %.2f%%" % (mH, numTot, len(hiDt), hiDtPct))
        plt.loglog(sumE, hiDt, '.', c=cols[mH],label='mH=%d  hi-dt %.2f%%' % (mH, hiDtPct))

    plt.xlabel("sumE", ha='right' ,x=1.)
    plt.ylabel("delta-t since last m=N event (sec)", ha='right',y=1.)
    plt.legend()
    plt.savefig("../plots/mult2-highDtVals.png")

    plt.cla()


def getDT():
    from ROOT import TFile, TTree, GATDataSet
    fileList = []
    calInfo = ds.CalInfo()
    runList = calInfo.GetSpecialRuns("longCal",5) # 54 total

    # runList = runList[:5] # make a lo-stats file (check filename below)
    # runList = runList[:15] # this gets ~10k mH=4 events
    runList = runList[:30] # get big stats (saved as different filename below)

    for run in runList: fileList.extend(ds.getLATRunList([run],"%s/lat" % (ds.specialDir)))
    nFiles = len(fileList)

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

    dtVals = {mH:[] for mH in range(30)}

    for iFile, f in enumerate(fileList):
        print("%d/%d %s" % (iFile,nFiles,f))
        # tf = TFile("%s/lat/%s" % (ds.specialDir,f))
        tf = TFile("%s/%s" % (os.environ['SLURM_TMP'], f)) # moved my shit to here
        lTree = tf.Get("skimTree")

        bNames = ["startClockTime_s","clockTime_s","mH","channel","gain","EventDC1Bits","trapENFCal","dtPulserCard"]
        lTree.SetBranchStatus('*',0)
        for name in bNames: lTree.SetBranchStatus(name,1)

        lTree.GetEntry(0)
        startTime = lTree.startClockTime_s

        dt, prev = {}, {}
        for iEnt in range(lTree.GetEntries()):
            lTree.GetEntry(iEnt)
            clockTime = lTree.clockTime_s
            nHit = lTree.channel.size()
            mH = lTree.mH
            evtTime = clockTime-startTime
            if lTree.EventDC1Bits != 0: continue

            # calculate dt
            if mH not in prev:
                prev[mH] = 0.
            if prev[mH] != 0:
                dt[mH] = evtTime - prev[mH]
                # print("%-4d  %-3d  evt %-10.9f  dt[%d] %-10.9f  pr[%d] %-10.9f  DCB %d" % (iEnt, mH, evtTime, mH, dt[mH], mH, prev[mH], lTree.EventDC1Bits))

                chans = [lTree.channel.at(i) for i in range(nHit) if lTree.gain.at(i)%2==0]
                hitE = [lTree.trapENFCal.at(i) for i in range(nHit) if lTree.gain.at(i)%2==0]
                dtpc = [lTree.dtPulserCard.at(i) for i in range(nHit) if lTree.gain.at(i)%2==0]

                dtVals[mH].append((dt[mH],chans,hitE,dtpc))

            # else:
                # print("%-4d  %-3d  evt %-10.9f  pr[%d] %-10.9f" % (iEnt, mH, evtTime, mH, prev[mH]))
            prev[mH] = evtTime

    # np.savez("../data/mult2-dtVals-ene.npz", dtVals, runTime)
    np.savez("../data/mult2-dtVals-ene-histats.npz", dtVals, runTime)


def plotFitRates():

    # f = np.load("../data/mult2-dtVals.npz")
    f = np.load("../data/mult2-dtVals-ene-histats.npz")
    dtVals, runTime = f['arr_0'].item(), f['arr_1']

    # get raw rate
    f = np.load("../data/mult2-sumSpec.npz") # generated from full long cal
    sumSpec, runTime = f['arr_1'], f['arr_4']
    rates = [sum(sumSpec[i])/runTime for i in range(0,6)]

    # one column (report)
    # fig2 = plt.figure(figsize=(9,24))
    # p = [0,0,0,0,0]
    # p[1] = plt.subplot(411)
    # p[2] = plt.subplot(412)
    # p[3] = plt.subplot(413)
    # p[4] = plt.subplot(414)

    # two columns (slides)
    fig2 = plt.figure(figsize=(18,12))
    p = [0,0,0,0,0]
    p[1] = plt.subplot(221)
    p[2] = plt.subplot(222)
    p[3] = plt.subplot(223)
    p[4] = plt.subplot(224)

    ranges = [0, (0,0.002,0.000005), (0,0.02,0.00005), (0,0.2,0.0005), (0,2.,0.005)]
    cols = [0,'r','b','m','c','g']

    for i in range(1,5):
        xLo, xHi, xpb = ranges[i][0], ranges[i][1], ranges[i][2]
        dtv = [dtVals[i][j][0] for j in range(len(dtVals[i]))]
        x, y = wl.GetHisto(dtv, xLo, xHi, xpb)
        p[i].plot(x, y, c=cols[i], lw=1., ls='steps', label='mH=%d, Raw: %.2f Hz' % (i, rates[i]))
        popt, pcov = curve_fit(expoDist, x, y, p0=(rates[i], 1))
        perr = np.sqrt(np.diag(pcov))
        rateFit, rateErr = popt[0], perr[0]
        p[i].plot(x, expoDist(x, *popt), '-b', lw=1.5, label=r'Fit: %.2f$\pm$%.2f Hz ' % (rateFit, rateErr))
        p[i].set_ylabel("Counts / %.0e sec" % xpb, ha='right', y=1.)
        p[i].legend()
        if i==1:
            p[i].ticklabel_format(axis='x',style='sci',scilimits=(1,2))
            p[i].set_xlabel("sec       ", ha='right', x=1.)
        else:
            p[i].set_xlabel("sec", ha='right', x=1.)


    plt.savefig("../plots/mult2-expoFit.png")


def plotDT_mH1Ene():

    f = np.load("../data/mult2-dtVals-ene.npz")
    dtVals = f['arr_0'].item()

    mH = 1
    n = len(dtVals[mH])
    delt = [dtVals[mH][i][0] for i in range(n)]
    chan = [dtVals[mH][i][1][0] for i in range(n)]
    hitE = [dtVals[mH][i][2][0] for i in range(n)]

    delt10 = [dtVals[mH][i][0] for i in range(n) if hitE[i] < 10]
    chan10 = [dtVals[mH][i][1][0] for i in range(n) if hitE[i] < 10]
    hitE10 = [dtVals[mH][i][2][0] for i in range(n) if hitE[i] < 10]

    # see if we have any spiky behavior
    # plot by energy, and by channel

    f = plt.figure(figsize=(9,12))
    p1 = plt.subplot(211)
    p2 = plt.subplot(212)

    chMap = list(sorted(set(chan10)))
    chDict = {chMap[i]:i for i in range(len(chMap))}
    chan10 = [chDict[chan] for chan in chan10]

    xLo, xHi, xpb = 0, 0.001, 0.000005
    yLo, yHi = 0, len(chMap)
    nbx, nby = int((xHi-xLo)/xpb), len(chMap)
    _,_,_,im = p1.hist2d(delt10, chan10, bins=[nbx, nby], range=[[xLo,xHi],[yLo,yHi]], norm=LogNorm(), cmap='jet')
    f.colorbar(im, ax=p1)

    p1.set_xlabel(r'$\Delta t$ (sec)', ha='right', x=1.)
    xLabels = ["%.1e" % x for x in np.arange(xLo, xHi, xpb*10)]
    p1.set_xticklabels(xLabels, fontsize=12)
    p1.set_ylabel("channel", ha='right', y=1.)
    yticks = np.arange(0, len(chMap))+0.5
    p1.set_yticks(yticks)
    p1.set_yticklabels(chMap, fontsize=12)

    xLo, xHi, xpb = 0.5, 5., 0.1
    nbx = int((xHi-xLo)/xpb)
    _,_,_,im2 = p2.hist2d(hitE10, chan10, bins=[nbx,nby], range=[[xLo,xHi],[yLo,yHi]], norm=LogNorm(), cmap='jet')
    f.colorbar(im2, ax=p2)

    p2.set_xlabel("hitE (keV)", ha='right', x=1.)
    p2.set_ylabel("channel", ha='right', y=1.)
    p2.set_yticks(yticks)
    p2.set_yticklabels(chMap, fontsize=12)

    plt.savefig("../plots/mult2-dtm1-ene.png")


def plotDT_mH1Ene_wThr():

    # same as above, but only plot hits if they are above the detector threshold.

    f = np.load("../data/mult2-dtVals-ene.npz")
    dtVals = f['arr_0'].item()

    mH = 1
    n = len(dtVals[mH])
    delt = [dtVals[mH][i][0] for i in range(n)]
    chan = [dtVals[mH][i][1][0] for i in range(n)]
    hitE = [dtVals[mH][i][2][0] for i in range(n)]

    deltTh, chanTh, hitETh = [], [], []
    for i in range(n):
        if chan[i] in thD.keys() and hitE[i] > thD[chan[i]][0] + 3 * thD[chan[i]][1]:
            deltTh.append(delt[i])
            chanTh.append(chan[i])
            hitETh.append(hitE[i])

    deltTh10 = [deltTh[i] for i in range(len(chanTh)) if hitETh[i] < 10]
    chanTh10 = [chanTh[i] for i in range(len(chanTh)) if hitETh[i] < 10]
    hitETh10 = [hitETh[i] for i in range(len(chanTh)) if hitETh[i] < 10]

    # rename back to the names used in the previous function
    chan, hitE, chan10, hitE10, delt10 = chanTh, hitETh, chanTh10, hitETh10, deltTh10

    # see if we have any spiky behavior
    # plot by energy, and by channel

    f = plt.figure(figsize=(9,12))
    p1 = plt.subplot(211)
    p2 = plt.subplot(212)

    chMap = list(sorted(set(chan10)))
    chDict = {chMap[i]:i for i in range(len(chMap))}
    chan10 = [chDict[chan] for chan in chan10]

    xLo, xHi, xpb = 0, 0.001, 0.000005
    yLo, yHi = 0, len(chMap)
    nbx, nby = int((xHi-xLo)/xpb), len(chMap)
    _,_,_,im = p1.hist2d(delt10, chan10, bins=[nbx, nby], range=[[xLo,xHi],[yLo,yHi]], norm=LogNorm(), cmap='jet')
    f.colorbar(im, ax=p1)

    p1.set_xlabel(r'$\Delta t$ (sec)', ha='right', x=1.)
    xLabels = ["%.1e" % x for x in np.arange(xLo, xHi, xpb*10)]
    p1.set_xticklabels(xLabels, fontsize=12)
    p1.set_ylabel("channel", ha='right', y=1.)
    yticks = np.arange(0, len(chMap))+0.5
    p1.set_yticks(yticks)
    p1.set_yticklabels(chMap, fontsize=12)

    xLo, xHi, xpb = 0.5, 5., 0.1
    nbx = int((xHi-xLo)/xpb)
    _,_,_,im2 = p2.hist2d(hitE10, chan10, bins=[nbx,nby], range=[[xLo,xHi],[yLo,yHi]], norm=LogNorm(), cmap='jet')
    f.colorbar(im2, ax=p2)

    p2.set_xlabel("hitE (keV)", ha='right', x=1.)
    p2.set_ylabel("channel", ha='right', y=1.)
    p2.set_yticks(yticks)
    p2.set_yticklabels(chMap, fontsize=12)

    plt.savefig("../plots/mult2-dtm1-thresh.png")


def plotDT_mH1Low():

    f = np.load("../data/mult2-dtVals-ene.npz")
    dtVals = f['arr_0'].item()

    mH = 1
    n = len(dtVals[mH])
    delt = [dtVals[mH][i][0] for i in range(n)]
    chan = [dtVals[mH][i][1][0] for i in range(n)]
    hitE = [dtVals[mH][i][2][0] for i in range(n)]
    delt10 = [dtVals[mH][i][0] for i in range(n) if hitE[i] < 10]
    delt_no598 = [dtVals[mH][i][0] for i in range(n) if chan[i]!=598]
    delt10_no598 = [dtVals[mH][i][0] for i in range(n) if hitE[i] < 10 and chan[i]!=598]

    xLo, xHi, xpb = 0, 0.005, 0.000005

    plt.semilogy(*wl.GetHisto(delt, xLo, xHi, xpb), c='r', lw=1., ls='steps', label='mH=1')

    x, y = wl.GetHisto(delt_no598, xLo, xHi, xpb)
    plt.semilogy(x, y, c='g', lw=1., ls='steps', label='mH=1, no598')

    popt, pcov = curve_fit(expoDist, x, y, p0=(300, 1))
    perr = np.sqrt(np.diag(pcov))
    rateFit, rateErr = popt[0], perr[0]
    plt.semilogy(x, expoDist(x, *popt), '-b', lw=1.5, label=r"Fit, no598 %.2f$\pm$%.2f Hz " % (rateFit, rateErr))

    x, y = wl.GetHisto(delt10_no598, xLo, xHi, xpb)
    plt.semilogy(x, y, c='m', lw=1., ls='steps', label='mH=1, E<10, no598')

    idx = np.where(x > 0.0002)
    popt, pcov = curve_fit(expoDist, x[idx], y[idx], p0=(30., 1))
    perr = np.sqrt(np.diag(pcov))
    rateFit, rateErr = popt[0], perr[0]
    plt.semilogy(x, expoDist(x, *popt), '-b', lw=1.5, label=r"Fit, no598, E<10 %.2f$\pm$%.2f Hz " % (rateFit, rateErr))

    plt.xlabel(r'$\Delta t$ (sec)', ha='right', x=1.)
    plt.ylabel("Counts / %.0e sec" % xpb, ha='right', y=1.)
    plt.legend(fontsize=12)
    plt.savefig("../plots/mult2-dtm1.png")


def threshFunc(x,mu,sig):
    from scipy.special import erf
    # return erf((x-mu)/sig)
    return 0.5*(1 + erf( (x - mu)/(sig*np.sqrt(2) ) ))
    # from scipy.special import expit
    # return expit((x-mu)/sig)


def plotThreshFunc():

    f = plt.figure(figsize=(4,2))
    mu = 0.5073093022
    sig = 0.1092046244
    x = np.arange(0, 1, 0.001)
    y = threshFunc(x, mu, sig)

    s1 = plt.axvline(mu+sig, c='b')
    s2 = plt.axvline(mu+2*sig, c='g')
    s3 = plt.axvline(mu+3*sig, c='m')

    plt.plot(x, y, '-r')
    plt.show()


def plotDT_mH1Low_wThr():

    # same as above, but only plot hits if they are above the detector threshold.

    f = np.load("../data/mult2-dtVals-ene.npz")
    dtVals = f['arr_0'].item()

    mH = 1
    n = len(dtVals[mH])
    delt = [dtVals[mH][i][0] for i in range(n)]
    chan = [dtVals[mH][i][1][0] for i in range(n)]
    hitE = [dtVals[mH][i][2][0] for i in range(n)]

    hitE10 = [hitE[i] for i in range(n) if hitE[i] < 10.]
    delt10 = [delt[i] for i in range(n) if hitE[i] < 10.]

    deltTh, chanTh, hitETh = [], [], []
    for i in range(n):
        if chan[i] in thD.keys() and hitE[i] > thD[chan[i]][0] + 3 * thD[chan[i]][1]:
            deltTh.append(delt[i])
            chanTh.append(chan[i])
            hitETh.append(hitE[i])

    deltTh_no598 = [deltTh[i] for i in range(len(chanTh)) if chanTh[i] != 598]

    deltTh10 = [deltTh[i] for i in range(len(chanTh)) if hitETh[i] < 10]
    chanTh10 = [chanTh[i] for i in range(len(chanTh)) if hitETh[i] < 10]
    hitETh10 = [hitETh[i] for i in range(len(chanTh)) if hitETh[i] < 10]

    deltTh10_no598 = [deltTh[i] for i in range(len(chanTh)) if hitETh[i] < 10 and chanTh[i] != 598]
    chanTh10_no598 = [chanTh[i] for i in range(len(chanTh)) if hitETh[i] < 10 and chanTh[i] != 598]
    hitETh10_no598 = [hitETh[i] for i in range(len(chanTh)) if hitETh[i] < 10 and chanTh[i] != 598]

    # try with a floor
    floor = 1.
    deltTh10_no598_f = [deltTh10_no598[i] for i in range(len(deltTh10_no598)) if hitETh10_no598[i] > floor]
    hitETh10_no598_f = [hitETh10_no598[i] for i in range(len(deltTh10_no598)) if hitETh10_no598[i] > floor]

    # now do plots

    f = plt.figure(figsize=(9,12))
    p0 = plt.subplot(211)
    p1 = plt.subplot(212)

    # xLo, xHi, xpb = 0, 0.005, 0.000005
    # p0.semilogy(*wl.GetHisto(delt, xLo, xHi, xpb), c='r', lw=1., ls='steps', label='mH=1, no thresh')
    # p0.semilogy(*wl.GetHisto(deltTh, xLo, xHi, xpb), c='g', lw=1., ls='steps', label='mH=1, w/ thresh')
    # p0.semilogy(*wl.GetHisto(deltTh_no598, xLo, xHi, xpb), c='b', lw=1., ls='steps', label='mH=1, w/ thresh, no 598')

    xLo, xHi, xpb = 0, 0.001, 0.000005
    p0.plot(*wl.GetHisto(delt10, xLo, xHi, xpb), c='r', lw=1., ls='steps', label='mH=1, no thresh')
    p0.plot(*wl.GetHisto(deltTh10, xLo, xHi, xpb), c='g', lw=1., ls='steps', label='w/ thresh')
    p0.plot(*wl.GetHisto(deltTh10_no598, xLo, xHi, xpb), c='b', lw=1., ls='steps', label='no 598')
    p0.plot(*wl.GetHisto(deltTh10_no598_f, xLo, xHi, xpb), c='m', lw=1., ls='steps', label='hits > %.1f kev' % floor)

    p0.set_xlabel(r'$\Delta t$ (sec)', ha='right', x=1.)
    p0.set_ylabel("Counts / %.0e sec" % xpb, ha='right', y=1.)
    p0.legend(fontsize=12)

    xLo, xHi, xpb = 0.5, 5, 0.05
    p1.semilogy(*wl.GetHisto(hitE10, xLo, xHi, xpb), c='r', lw=1., ls='steps', label='mH=1, no thresh')
    p1.semilogy(*wl.GetHisto(hitETh10, xLo, xHi, xpb), c='g', lw=1., ls='steps', label='w/thresh')
    p1.semilogy(*wl.GetHisto(hitETh10_no598, xLo, xHi, xpb), c='b', lw=1., ls='steps', label='no 598')
    p1.semilogy(*wl.GetHisto(hitETh10_no598_f, xLo, xHi, xpb), c='m', lw=1., ls='steps', label='hits > %.1f kev' % floor)

    p1.set_xlabel("hitE", ha='right', x=1.)
    p1.set_ylabel("Counts", ha='right', y=1.)
    p1.legend(fontsize=12)

    plt.savefig("../plots/mult2-lowE-thresh.png")


def plotDT_wThr():

    # same as above, but show fewer steps and do it for all multiplicities

    # f = np.load("../data/mult2-dtVals-ene.npz")
    f = np.load("../data/mult2-dtVals-ene-histats.npz")
    dtVals = f['arr_0'].item()

    for mH in [1,2,3,4]:
    # for mH in [2]:

        n = len(dtVals[mH])

        deltAll = [dtVals[mH][i][0] for i in range(n)]

        delt, chan, hitE, dtpc = [], [], [], []

        # no e cut
        # for i in range(n):
        #     chan.extend( dtVals[mH][i][1] )
        #     hitE.extend( dtVals[mH][i][2] )
        #     dtpc.extend( dtVals[mH][i][3] )
        #     delt.extend( [deltAll[i] for n in range(len(dtVals[mH][i][2]))] ) # creates copies

        # with e cut
        for i in range(n):
            hTmp = dtVals[mH][i][2] # mH values
            dTmp = dtVals[mH][i][0] # single value
            chan.extend( [dtVals[mH][i][1][j] for j in range(len(hTmp)) if 0.7 < hTmp[j] < 10.] )
            hitE.extend( [dtVals[mH][i][2][j] for j in range(len(hTmp)) if 0.7 < hTmp[j] < 10.] )
            dtpc.extend( [dtVals[mH][i][3][j] for j in range(len(hTmp)) if 0.7 < hTmp[j] < 10.] )
            delt.extend( [dTmp for j in range(len(hTmp)) if 0.7 < hTmp[j] < 10.] )

        n = len(hitE)
        delt_no598 = [delt[i] for i in range(n) if chan[i]!=598]
        hitE_no598 = [hitE[i] for i in range(n) if chan[i]!=598]

        deltTh, hitETh, dtpcTh = [], [], []
        for i in range(len(hitE)):
            if chan[i] not in thD.keys():
                continue
            mu, sig = thD[chan[i]][0], thD[chan[i]][1]
            if hitE[i] > (mu + 3*sig) and chan[i]!=598:
                deltTh.append(delt[i])
                hitETh.append(hitE[i])
                dtpcTh.append(dtpc[i])

        pT =  8.388608
        eps = 0.01
        n = len(hitETh)
        delt_dtp = [deltTh[i] for i in range(n) if abs(dtpcTh[i]-pT/2.) > eps] # could add 8.388 condition also
        hitE_dtp = [hitETh[i] for i in range(n) if abs(dtpcTh[i]-pT/2.) > eps]

        plt.cla()
        f = plt.figure(figsize=(9,12))
        p0 = plt.subplot(211)
        p1 = plt.subplot(212)

        ranges = [0, (0,0.002,0.000005), (0,0.02,0.0001), (0,0.2,0.001), (0,2.,0.01)]

        xLo, xHi, xpb = ranges[mH][0], ranges[mH][1], ranges[mH][2]

        # 1 for each event
        p0.semilogy(*wl.GetHisto(deltAll, xLo, xHi, xpb), c='k', lw=1., ls='steps', label='mH=%d, all' % mH)
        p0.semilogy(*wl.GetHisto(delt, xLo, xHi, xpb), c='r', lw=1., ls='steps', label='E < 10')
        p0.semilogy(*wl.GetHisto(delt_no598, xLo, xHi, xpb), c='b', lw=1., ls='steps', label='no 598')
        p0.semilogy(*wl.GetHisto(deltTh, xLo, xHi, xpb), c='g', lw=1., ls='steps', label='w/ thresh')
        p0.semilogy(*wl.GetHisto(delt_dtp, xLo, xHi, xpb), c='m', lw=1., ls='steps', label='w/ dtpc')
        # delt_dtp

        p0.set_xlabel(r'$\Delta t$ (sec)', ha='right', x=1.)
        p0.set_ylabel("Counts / %.0e sec" % xpb, ha='right', y=1.)
        p0.legend(fontsize=12)

        # mH=N for each event
        xLo, xHi, xpb = 0.5, 5, 0.05
        p1.semilogy(*wl.GetHisto(hitE, xLo, xHi, xpb), c='r', lw=1., ls='steps', label='mH=%d, no thresh' % mH)
        p1.semilogy(*wl.GetHisto(hitE_no598, xLo, xHi, xpb), c='b', lw=1., ls='steps', label='no 598')
        p1.semilogy(*wl.GetHisto(hitETh, xLo, xHi, xpb), c='g', lw=1., ls='steps', label='w/ thresh')
        p1.semilogy(*wl.GetHisto(hitE_dtp, xLo, xHi, xpb), c='m', lw=1., ls='steps', label='w/ dtpc')

        p1.set_xlabel("hitE", ha='right', x=1.)
        p1.set_ylabel("Counts", ha='right', y=1.)
        p1.legend(fontsize=12)

        plt.savefig("../plots/mult2-lowE-thresh-mh%d.png" % mH)


def plotChannelRates_mH1():

    f = np.load("../data/mult2-dtVals-ene.npz")
    dtVals = f['arr_0'].item()
    runTime = 2473.0 # kludge - this is for runList[7:14] in getDT

    f = plt.figure()

    mH = 1
    n = len(dtVals[mH])

    plt.cla()
    chan = [dtVals[mH][i][1][0] for i in range(n)]
    hitE = [dtVals[mH][i][2][0] for i in range(n)]
    chan10 = [chan[i] for i in range(n) if chan[i]!=598 and hitE[i]<10]
    chan = [chan[i] for i in range(n) if chan[i]!=598]
    chMap = list(sorted(set(chan)))

    x, y = wl.GetHisto(chan, chMap[0], chMap[-1], 1)
    idx = np.where(y!=0)
    x, y = x[idx]-0.5, y[idx]
    x = [int(ch) for ch in x]
    xb = np.arange(0,len(x),1)
    rateY = y/runTime
    rateFull = sum(rateY)
    plt.bar(xb, rateY, 0.95, color='blue', log=True, label='All mH=1 (no ch598), %.2f Hz' % (rateFull))

    x10, y10 = wl.GetHisto(chan10, chMap[0], chMap[-1], 1)
    idx = np.where(y10!=0)
    rate10Y = y10[idx]/runTime
    rate10 = sum(rate10Y)
    plt.bar(xb, rate10Y, 0.95, color='red', log=True, label='E < 10, %.2f Hz' % (rate10))

    plt.xticks(xb)
    plt.gca().set_xticklabels(x, fontsize=12, rotation='vertical')
    plt.xlabel("channel", ha='right', x=1.)
    plt.ylabel("Rate (Hz)", ha='right', y=1.)
    plt.legend(fontsize=14, ncol=2)
    plt.savefig("../plots/mult2-chanCts.png")


def noiseProbability():
    from math import factorial

    from scipy.stats import poisson
    # poisson.pmf(n, r * dt) # this is prob. of getting n noise events w/ mean rate mu = r*dt
    # def pois(n, mu): return (mu)**n * np.exp(-1*mu) / factorial(n) # equivalent.

    rNoise = 40
    rmH = [0, 344.73, 69.60, 12.14, 1.96, 0.58]  # observed rates

    # observed physics probabilities
    pPhys = [ poisson.pmf(1, rmH[i] * i * 4e-6) for i in range(0,6) ]

    # multi-hit noise probabilities
    pNoise = [ poisson.pmf(i, rNoise * i * 4e-6) for i in range(0,6) ]

    # mixed physics/noise probabilieis
    pMix = [0, 0]
    for N in range(2, 6):
        pN = 0
        pairs = [ (i, N-i) for i in range(1, N) ]
        for i, j in pairs:
            pN += pPhys[i] * pNoise[j]
        pMix.append(pN)

    for i in range(1, 6):
        pctMix = 100 * pMix[i] / pPhys[i]
        pctNoise = 100 * pNoise[i] / pPhys[i]
        print("N %d  %-3dus  Pp %.2e  Pn %.2e  Pm %.2e || %%Mix %.2e  %%Noise %.2e" % (i,i*4,pPhys[i], pNoise[i], pMix[i], pctMix, pctNoise))


def noiseProbability2D():

    f = np.load("../data/mult2-dtVals-ene.npz")
    dtVals = f['arr_0'].item()

    # recalculated the run time, i forgot to save it
    runTime = 2473.0

    mH = 1
    n = len(dtVals[mH])
    chan = [dtVals[mH][i][1][0] for i in range(n)]
    hitE = [dtVals[mH][i][2][0] for i in range(n)]

    deltTh, chanTh, hitETh = [], [], []
    for i in range(n):
        if chan[i] in thD.keys() and hitE[i] > thD[chan[i]][0] + 3 * thD[chan[i]][1]:
            chanTh.append(chan[i])
            hitETh.append(hitE[i])

    chanTh10 = [chanTh[i] for i in range(len(chanTh)) if hitETh[i] < 10]
    hitETh10 = [hitETh[i] for i in range(len(chanTh)) if hitETh[i] < 10]

    # rename back to the names used in the previous function
    chan, hitE, chan10, hitE10 = chanTh, hitETh, chanTh10, hitETh10

    fig = plt.figure()

    chMap = list(sorted(set(chan10)))
    chDict = {chMap[i]:i for i in range(len(chMap))}
    chan10 = [chDict[chan] for chan in chan10]

    xLo, xHi, xpb = 0.5, 5., 0.1
    yLo, yHi = 0, len(chMap)
    nbx, nby = int((xHi-xLo)/xpb), len(chMap)
    cts,_,_,_ = plt.hist2d(hitE10, chan10, bins=[nbx,nby], range=[[xLo,xHi],[yLo,yHi]], norm=LogNorm(), cmap='jet')
    plt.cla()

    def pois(n, mu):
        from math import factorial
        return (mu)**n * np.exp(-1*mu) / factorial(n) # equivalent to scipy.stats.poisson.pmf

    rates = [0,385.,69.2,11.95,1.82] # fitted mH rates
    mH = 2
    totRate = 0

    nx, ny = cts.shape[0], cts.shape[1]
    pMat = np.zeros((nx,ny))
    for iX in range(nx):
        for iY in range(ny):

            binCts = cts[iX,iY]
            if binCts == 0:
                pMat[iX,iY] = 0
                continue

            binRate = binCts/runTime
            dt = mH * 4e-6
            pCoin = pois(1, binRate*dt) * pois(1, rates[mH]*dt)

            pMat[iX,iY] = pCoin
            totRate += binRate

    idx = np.where(pMat>0)
    pMax, pMin = pMat.max(), pMat[idx].min()
    print("mH = %d  pMax %.2e  pMin %.2e" % (mH, pMat.max(), pMat[idx].min()))

    img2 = plt.imshow(pMat.T, extent=[xLo,xHi,yLo,yHi], origin='lower', cmap='jet', norm=LogNorm(), aspect='auto')

    # this fking wins the title of world's most complicated colorbar
    _, e1 = '{:.2e}'.format(pMin).split('e')
    _, e2 = '{:.2e}'.format(pMax).split('e')
    e1, e2 = int(e1), int(e2)
    pTicks = np.asarray([10**b for b in range(e1,e2+2)])
    cbar = fig.colorbar(img2, ax=plt.gca(), fraction=0.046, pad=0.04, norm=Normalize(), ticks=pTicks, format="%.0e")
    cbar.ax.set_title("Prob, mH=%d" % mH, fontsize=14)

    plt.xlabel("hitE (keV)", ha='right', x=1.)
    plt.ylabel("channel", ha='right', y=1.)

    yticks = np.arange(0, len(chMap))+0.5
    plt.yticks(yticks)
    plt.gca().set_yticklabels(chMap, fontsize=12)

    plt.savefig("../plots/mult2-dtm1-prob.png")


def plotDTCut():

    hitESpec = {mH:[] for mH in [1,2,3,4]}      # E cut, 598 cut, thresh cut
    hitESpec_dtp = {mH:[] for mH in [1,2,3,4]}  # dtPulser cut

    # f = np.load("../data/mult2-dtVals-ene.npz")
    f = np.load("../data/mult2-dtVals-ene-histats.npz")
    dtVals, runTime = f['arr_0'].item(), f['arr_1']
    for mH in [1,2,3,4]:
        n = len(dtVals[mH])

        # ecut, 598 cut
        delt, chan, hitE, dtpc = [], [], [], []
        for i in range(n):
            hTmp = dtVals[mH][i][2] # mH values
            dTmp = dtVals[mH][i][0] # single value
            chan.extend( [dtVals[mH][i][1][j] for j in range(len(hTmp)) if 0.7 < hTmp[j] < 10. and dtVals[mH][i][1][j]!=598] )
            hitE.extend( [dtVals[mH][i][2][j] for j in range(len(hTmp)) if 0.7 < hTmp[j] < 10. and dtVals[mH][i][1][j]!=598] )
            dtpc.extend( [dtVals[mH][i][3][j] for j in range(len(hTmp)) if 0.7 < hTmp[j] < 10. and dtVals[mH][i][1][j]!=598] )
            delt.extend( [dTmp for j in range(len(hTmp)) if 0.7 < hTmp[j] < 10.] )

        # thresh cut
        deltTh, hitETh, dtpcTh = [], [], []
        for i in range(len(hitE)):
            if chan[i] not in thD.keys(): continue
            mu, sig = thD[chan[i]][0], thD[chan[i]][1]
            if hitE[i] > (mu + 3*sig):
                # deltTh.append(delt[i])
                hitETh.append(hitE[i])
                dtpcTh.append(dtpc[i])

        # pulser retrigger cut
        pT =  8.388608
        eps = 0.01
        n = len(hitETh)
        # delt_dtp = [deltTh[i] for i in range(n) if abs(dtpcTh[i]-pT/2.) > eps] # could add 8.388 condition also
        hitE_dtp = [hitETh[i] for i in range(n) if abs(dtpcTh[i]-pT/2.) > eps]

        # save spectra
        hitESpec[mH].extend(hitETh)
        hitESpec_dtp[mH].extend(hitE_dtp)

    # save to file for comparison in plotPCut()
    np.savez("../data/mult2-dtp-cutSpec.npz", hitESpec_dtp)

    f, ax = plt.subplots()

    cols = [0, 'r','b','m','c']

    xLo, xHi, xpb = 0.7, 10, 0.1
    for mH in [1,2,3,4]:

        ax.semilogy(*wl.GetHisto(hitESpec[mH], xLo, xHi, xpb), c=cols[mH], lw=1., ls='steps', label='mH=%d' % mH)
        # line1, = ax.semilogy(*wl.GetHisto(hitESpec_dtp[mH], xLo, xHi, xpb), '--', c=cols[mH], lw=1., ls='steps', label='mH=%d w/dtp' % mH)
        # line1.set_dashes([1, 1])

    plt.xlabel("hitE (keV)", ha='right', x=1.)
    plt.legend(loc=1, fontsize=12)
    plt.savefig("../plots/mult2-dtcut.png")


def plotPCut():

    # use the 2d histo to cut on probability

    f = np.load("../data/mult2-dtVals-ene.npz")
    # f = np.load("../data/mult2-dtVals-ene-histats.npz")
    dtVals, runTime = f['arr_0'].item(), f['arr_1']

    mH = 1
    n = len(dtVals[mH])

    # ecut only
    chan, hitE = [], []
    for i in range(n):
        hTmp = dtVals[mH][i][2] # mH values
        chan.extend( [dtVals[mH][i][1][j] for j in range(len(hTmp)) if 0.7 < hTmp[j] < 10.] )
        hitE.extend( [dtVals[mH][i][2][j] for j in range(len(hTmp)) if 0.7 < hTmp[j] < 10.] )

    fig = plt.figure()

    chMap = list(sorted(set(chan)))
    chDict = {chMap[i]:i for i in range(len(chMap))}
    chan = [chDict[chan] for chan in chan]

    xLo, xHi, xpb = 0.5, 5., 0.1
    yLo, yHi = 0, len(chMap)
    nbx, nby = int((xHi-xLo)/xpb), len(chMap)
    cts,_,_,_ = plt.hist2d(hitE, chan, bins=[nbx,nby], range=[[xLo,xHi],[yLo,yHi]], norm=LogNorm(), cmap='jet')
    plt.cla()

    def pois(n, mu):
        from math import factorial
        return (mu)**n * np.exp(-1*mu) / factorial(n) # equivalent to scipy.stats.poisson.pmf

    rates = [0,385.,69.2,11.95,1.82] # fitted mH rates
    mH = 2
    totRate = 0

    nx, ny = cts.shape[0], cts.shape[1] # nx: num e bins,  ny: num detectors
    pMat = np.zeros((nx,ny))
    for iX in range(nx):
        for iY in range(ny):
            binCts = cts[iX,iY]
            if binCts == 0:
                pMat[iX,iY] = 0
                continue
            binRate = binCts/runTime
            dt = mH * 4e-6
            pCoin = pois(1, binRate*dt) * pois(1, rates[mH]*dt) # raw probability
            pMat[iX,iY] = pCoin
            totRate += binRate

    idx = np.where(pMat>0)
    pMax, pMin = pMat.max(), pMat[idx].min()
    print("mH = %d  pMax %.2e  pMin %.2e" % (mH, pMat.max(), pMat[idx].min()))

    img2 = plt.imshow(pMat.T, extent=[xLo,xHi,yLo,yHi], origin='lower', cmap='jet', norm=LogNorm(), aspect='auto')

    # this fking wins the title of world's most complicated colorbar
    _, e1 = '{:.2e}'.format(pMin).split('e')
    _, e2 = '{:.2e}'.format(pMax).split('e')
    e1, e2 = int(e1), int(e2)
    pTicks = np.asarray([10**b for b in range(e1,e2+2)])
    cbar = fig.colorbar(img2, ax=plt.gca(), fraction=0.046, pad=0.04, norm=Normalize(), ticks=pTicks, format="%.0e")
    cbar.ax.set_title("Prob, mH=%d" % mH, fontsize=14)

    plt.xlabel("hitE (keV)", ha='right', x=1.)
    plt.ylabel("channel", ha='right', y=1.)

    yticks = np.arange(0, len(chMap))+0.5
    plt.yticks(yticks)
    plt.gca().set_yticklabels(chMap, fontsize=12)

    plt.savefig("../plots/mult2-dtm1-prob.png")


    plt.cla()

    pCutVal = 1e-10
    pCut = np.zeros((nx,ny)) # fill w/ unscaled counts
    for iX in range(nx):
        for iY in range(ny):
            binCts = cts[iX,iY]
            if binCts == 0:
                pCut[iX,iY] = 0
                continue
            binRate = binCts/runTime
            dt = mH * 4e-6
            pCoin = pois(1, binRate*dt) * pois(1, rates[mH]*dt)
            if pCoin < pCutVal:
                pCut[iX, iY] = binCts
    # img2 = plt.imshow(pCut.T, extent=[xLo,xHi,yLo,yHi], origin='lower', cmap='jet', norm=LogNorm(), aspect='auto')

    hitESpec = pCut.sum(axis=1) # sum over rows (channels)
    hitESpec = np.insert(hitESpec, 0, 0, axis=0) # prepend a 0

    x, y = wl.GetHisto(hitE, xLo, xHi, xpb)

    plt.semilogy(x, y, ls='steps', c='b', lw=1., label='all mH=%d cts' % mH)
    plt.semilogy(x, hitESpec, ls='steps', c='r', lw=1., label='prob < %.1e' % pCutVal)

    plt.xlabel("hitE (keV)", ha='right', x=1.)
    plt.ylabel("Counts", ha='right', y=1.)
    plt.legend()
    plt.savefig("../plots/mult2-dtm1-pcut.png")


    plt.cla()
    fig2 = plt.figure()

    f2 = np.load("../data/mult2-dtp-cutSpec.npz")
    hitESpec_dtp = f2['arr_0'].item()

    plt.plot(*wl.GetHisto(hitESpec_dtp[mH], xLo, xHi, xpb), ls='steps', lw=1., c='r', label='mH=%d dtpcut' % mH)
    plt.plot(x, hitESpec, ls='steps', c='g', lw=1., label='prob < %.1e' % pCutVal)

    plt.xlabel("hitE (keV)", ha='right', x=1)
    plt.legend()
    plt.savefig("../plots/mult2-pcut-vs-dtpcut.png")


def getPeaks():

    from ROOT import TFile, TTree, GATDataSet
    fileList = []
    calInfo = ds.CalInfo()
    runList = calInfo.GetSpecialRuns("longCal",5) # 54 total
    # runList = runList[:5]
    # runList = runList[:15] # this gets ~10k mH=4 events
    runList = runList[:] # hi stats, using different file name below
    for run in runList: fileList.extend(ds.getLATRunList([run],"%s/lat" % (ds.specialDir)))
    nFiles = len(fileList)

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

    # 6 peaks total, for m=2 and m=3:  238, 583, 2615
    pks = [(235,242,0.1), (580,586,0.1), (2605,2625,0.2)]

    # note: to access x and y histogram values:
    # x, y = pkHist[mH][iPk][0], pkHist[mH][iPk][1]
    pkHist = {2:[], 3:[]}
    for pk in pks:
        x, y = wl.GetHisto([0], *pk)
        pkHist[2].append([x, y])
        pkHist[3].append([x, y])

    bNames = ["Entry$", "mH", "gain", "isGood", "EventDC1Bits", "trapENFCal","sumEH","channel"]

    for iFile, f in enumerate(fileList):
        print("%d/%d %s" % (iFile,nFiles,f))
        # tf = TFile("%s/lat/%s" % (ds.specialDir,f))
        tf = TFile("%s/%s" % (os.environ['SLURM_TMP'],f)) # I moved my shit here
        lTree = tf.Get("skimTree")

        for mH in range(2,4):
            theCut = "mH==%d && gain==0 && isGood && !EventDC1Bits" % mH
            sumE = getSumE(lTree, bNames, theCut)
            for iPk, pk in enumerate(pks):
                idx = np.where((sumE > pk[0]) & (sumE < pk[1]))
                x, y = wl.GetHisto(sumE[idx], *pk)
                pkHist[mH][iPk][0] = x
                pkHist[mH][iPk][1] = np.add(pkHist[mH][iPk][1], y)

        tf.Close()

    # np.savez("../data/mult2-peaks.npz", runTime, pks, pkHist)
    np.savez("../data/mult2-peaks-histats.npz", runTime, pks, pkHist)


def roughSigma(ene):
    """ from gpxFitter, just used to set initial guesses """
    p0, p1, p2 = 0.2, 0.02, 0.0003
    return np.sqrt(p0**2. + p1**2. * ene + p2**2. * ene**2.)


def gaus(x, b, a, mu, sig):
    """ gaussian + flat bg """
    return b + a * np.exp(-(x-mu)**2. / (2. * sig**2.))


def plotPeaks():

    # f = np.load("../data/mult2-peaks.npz")
    f = np.load("../data/mult2-peaks-histats.npz") # full cal run
    runTime, pks, pkHist = f['arr_0'], f['arr_1'], f['arr_2'].item()

    # note: to access x and y histogram values:
    # x, y = pkHist[mH][iPk][0], pkHist[mH][iPk][1]

    pkE = [238, 583, 2615]
    for mH in [2,3]:
        for iPk in range(len(pks)):
            x, y = pkHist[mH][iPk][0], pkHist[mH][iPk][1]

            plt.cla()
            plt.plot(x, y, c='b', ls='steps', label="m=%d" % (mH))

            p0 = (np.mean(y[:5]), max(y), pkE[iPk], roughSigma(pkE[iPk]))
            popt,_ = curve_fit(gaus, x, y, p0=p0)
            bpx = x[1] - x[0]
            bgRate = popt[0] * bpx
            mu, sig = popt[2], popt[3]
            idx = np.where((x > mu-3*sig) & (x < mu+3*sig))
            totCts = np.sum(y[idx])
            bgCts = bgRate * len(y[idx])
            pkCts = totCts - bgCts
            pbr = pkCts / bgCts
            print("%d  %.2f  %.2f  %.2f  %.2f  %d  %d  %d  P/B: %.2f" % (pkE[iPk], bpx, bgRate, mu, sig, totCts, bgCts, pkCts, pbr))
            print("sumLo, sumHi = %.2f, %.2f" % (mu-3*sig, mu+3*sig))

            plt.plot(x, gaus(x, *popt), 'r-', label='fit.  mu %.2f  sig %.2f\nP/B %.2f' % (popt[2], popt[3], pbr))

            plt.title("Peak-to-bkg ratio: %.3f" % pbr)
            plt.xlabel("sumE (keV)", ha='right', x=1.)
            plt.ylabel("Counts", ha='right', y=1.)
            plt.legend(loc=1, fontsize=14)
            plt.savefig("../plots/mult2-m%d-pk%d.png" % (mH,iPk))

            # return


def getSumEvents():

    from ROOT import gROOT, TFile, TTree, GATDataSet
    gROOT.ProcessLine("gErrorIgnoreLevel = 3001;") # suppress ROOT error messages

    fileList = []
    calInfo = ds.CalInfo()
    runList = calInfo.GetSpecialRuns("longCal",5) # 54 total

    # runList = runList[:2] # make a lo-stats file (check filename below)
    # runList = runList[:15] # this gets ~10k mH=4 events
    runList = runList[:] # get big stats (saved as different filename below)

    for run in runList: fileList.extend(ds.getLATRunList([run],"%s/lat" % (ds.specialDir)))
    nFiles = len(fileList)

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

    # each entry: (chan, sumE, hitE, dtpc, fitSlo) - only keep the lo hits.
    sumEvts = {
        2:{238:[], 583:[], 2615:[]},
        3:{238:[], 583:[], 2615:[]}
        }
    sumPks = [
        (2,238,237.28,239.46), (2,583,581.26,584.46), (2,2615,2610.57,2618.01),
        (3,238,237.13,239.43), (3,583,581.04,584.36), (3,2615,2610.10,2617.92)
        ]

    # save a dict of just hit groups
    comptonEvts2 = { 238:[], 583:[], 2615:[] }
    comptonEvts3 = { 238:[], 583:[], 2615:[] }

    bNames = ["Entry$","mH","gain","isGood","EventDC1Bits","sumEH","channel","trapENFCal","dtPulserCard","fitSlo"]

    for iFile, f in enumerate(fileList):
        print("%d/%d %s" % (iFile,nFiles,f))
        # tf = TFile("%s/lat/%s" % (ds.specialDir,f))
        tf = TFile("%s/%s" % (os.environ['SLURM_TMP'], f)) # moved my shit to here
        lTree = tf.Get("skimTree")

        for sp in sumPks:
            mH, pk, sLo, sHi = sp

            # this draw trusts sumEH to be true.  it's ok - i've already confirmed it matches the GetSumE function.
            theCut = "mH==%d && gain==0 && isGood && !EventDC1Bits && sumEH > %.2f && sumEH < %.2f" % (mH, sLo, sHi)
            tvals = wl.GetVX(lTree, bNames, theCut)

            # create a set of entry numbers - not so bad for single files
            eList = sorted(set(tvals["Entry$"]))
            for iEnt in eList:
                idx = np.where(tvals["Entry$"]==iEnt)
                if len(idx[0]) != mH: continue

                hitE = tvals["trapENFCal"][idx]
                chan = tvals["channel"][idx]
                fSlo = tvals["fitSlo"][idx]
                dtpc = tvals["dtPulserCard"][idx]

                sumE = tvals["sumEH"][idx][0]
                if not sLo < sumE < sHi:
                    print("ERROR: some stuff went funny. sumEH: %.2f  sLo, sHi: %.2f %.2f" % (sumE, sLo, sHi))
                    return

                # save the hit array
                if mH==2: comptonEvts2[pk].append(hitE)
                if mH==3: comptonEvts3[pk].append(hitE)

                # save the lo hit only (floats)
                iLo = np.argmin(hitE)
                sumEvts[mH][pk].append((chan[iLo],sumE,hitE[iLo],dtpc[iLo],fSlo[iLo]))

        tf.Close()

    # np.savez("../data/mult2-comptonEdgeEvts.npz", runTime, comptonEvts2, comptonEvts3)
    # np.savez("../data/mult2-sumLoEvts.npz", runTime, sumEvts)

    np.savez("../data/mult2-comptonEdgeEvts-histats.npz", runTime, comptonEvts2, comptonEvts3)
    np.savez("../data/mult2-sumLoEvts-histats.npz", runTime, sumEvts)


def plotComptonEdge():

    # f = np.load("../data/mult2-comptonEdgeEvts.npz")
    f = np.load("../data/mult2-comptonEdgeEvts-histats.npz")
    runTime, comptonEvts2, comptonEvts3 = f['arr_0'], f['arr_1'].item(), f['arr_2'].item()

    f = plt.figure()

    for i, pk in enumerate(comptonEvts2):

        evts = comptonEvts2[pk]
        eLo, eHi = [], []
        for evt in evts:
            iLo, iHi = np.argmin(evt), np.argmax(evt)
            eLo.append(evt[iLo])
            eHi.append(evt[iHi])
        eLo, eHi = np.asarray(eLo), np.asarray(eHi)

        if i==0: xLo, xHi, xpb, pkE = 0., 240., 1., 238.
        if i==1: xLo, xHi, xpb, pkE = 0., 590., 1., 583.
        if i==2: xLo, xHi, xpb, pkE = 0., 2620., 1., 2615.

        plt.cla()
        plt.plot(*wl.GetHisto(eLo, xLo, xHi, xpb), c='r', lw=2., ls='steps',label='lo hits')
        plt.plot(*wl.GetHisto(eHi, xLo, xHi, xpb), c='b', lw=2., ls='steps',label='hi hits')

        # eCrit = pkE * (1 - 1 / (1 + 2*pkE/511))
        eCrit = pkE / (1 + 2*pkE/511)
        plt.axvline(eCrit, color='g', alpha=0.3, label=r'$E_C=$%.1f keV' % eCrit)
        plt.axvline(pkE, color='g', alpha=0.3, label=r'$E=$%.0f keV' % pkE)

        plt.xlabel("hitE (keV)", ha='right', x=1.)
        plt.ylabel("Counts / %.1f keV" % xpb, ha='right', y=1.)
        plt.legend(loc=2, fontsize=14)
        plt.savefig("../plots/mult2-comptonEdge-%.0f.png" % (pkE))

    for i, pk in enumerate(comptonEvts3):

        evts = comptonEvts3[pk]
        eLo, eMid, eHi = [], [], []
        for evt in evts:
            iLo, iHi = np.argmin(evt), np.argmax(evt)
            iMid = list(set([0,1,2]) - set([iLo,iHi]))[0]
            eLo.append(evt[iLo])
            eMid.append(evt[iMid])
            eHi.append(evt[iHi])
        eLo, eMid, eHi = np.asarray(eLo), np.asarray(eMid), np.asarray(eHi)

        if i==0: xLo, xHi, xpb, pkE = 0., 240., 1., 238.
        if i==1: xLo, xHi, xpb, pkE = 0., 590., 1., 583.
        if i==2: xLo, xHi, xpb, pkE = 0., 2620., 1., 2615.

        plt.cla()
        plt.plot(*wl.GetHisto(eLo, xLo, xHi, xpb), c='r', lw=2., ls='steps',label='lo hits')
        plt.plot(*wl.GetHisto(eMid, xLo, xHi, xpb), c='g', lw=2., ls='steps',label='mid hits')
        plt.plot(*wl.GetHisto(eHi, xLo, xHi, xpb), c='b', lw=2., ls='steps',label='hi hits')

        # eCrit = pkE * (1 - 1 / (1 + 2*pkE/511))
        eCrit = pkE / (1 + 2*pkE/511)
        plt.axvline(eCrit, color='g', alpha=0.3, label=r'$E_C=$%.1f keV' % eCrit)
        plt.axvline(pkE, color='g', alpha=0.3, label=r'$E=$%.0f keV' % pkE)

        plt.xlabel("hitE (keV)", ha='right', x=1.)
        plt.ylabel("Counts / %.1f keV" % xpb, ha='right', y=1.)
        plt.legend(loc=2, fontsize=14)
        plt.savefig("../plots/mult2-comptonEdge-m3-%.0f.png" % (pkE))


def plotLoHits():

    f = np.load("../data/mult2-sumLoEvts.npz")
    # f = np.load("../data/mult2-sumLoEvts-histats.npz")
    runTime, sumEvts = f['arr_0'], f['arr_1'].item()

    f = plt.figure()

    cols = ['r','b','m','c','k','g','o','y']
    for iP, mH, pk in [(0,2,238),(1,2,583),(2,2,2615),(3,3,238),(4,3,583),(5,3,2615)]:

        evts = sumEvts[mH][pk]

        # chan = [evt[0] for evt in evts]
        # sumE = [evt[1] for evt in evts]
        hitE = [evt[2] for evt in evts]
        # dtpc = [evt[3] for evt in evts]
        # fSlo = [evt[4] for evt in evts]

        # plt.cla()
        xLo, xHi, xpb = 0., 50., 0.2

        plt.plot(*wl.GetHisto(hitE,xLo,xHi,xpb), ls='steps', c=cols[iP], lw=2., label='mH=%d E=%d' % (mH,pk))

    plt.xlabel("hitE (keV)", ha='right', x=1.)
    plt.ylabel("Counts / %.1f keV" % xpb, ha='right', y=1.)
    plt.legend()
    plt.savefig("../plots/mult2-loHitSpec.png")


def plotSloHits():
    """ This is the original plot that got us excited about 85% efficiency of slowness.
        **** BUG, DO NOT USE ***
    """

    # f = np.load("../data/mult2-sumLoEvts.npz")
    f = np.load("../data/mult2-sumLoEvts-histats.npz")
    runTime, sumEvts = f['arr_0'], f['arr_1'].item()

    mH, pk = 2, 238

    # load fitSlo constants for cal run range closest to the run range [22513, 22566]
    dsNum, modNum, calIdx = 5, 1, 11  # calIdx 11: [[22568,22635],22568,22841],
    calDB = db.TinyDB('../calDB.json')
    pars = db.Query()
    fsDict = ds.getDBRecord("fitSlo_ds%d_idx%d_m%d_Peak" % (dsNum, calIdx, modNum), False, calDB, pars)

    evts = sumEvts[mH][pk]
    chan = [evt[0] for evt in evts if evt[0] in fsDict.keys() and evt[2] > 0.7] # only plot detectors we have fs values for
    hitE = [evt[2] for evt in evts if evt[0] in fsDict.keys() and evt[2] > 0.7]
    fSlo = [evt[4] for evt in evts if evt[0] in fsDict.keys() and evt[2] > 0.7]
    n = len(hitE)

    # # apply threshold cut
    # **** BUG, DO NOT USE ***
    # I forgot to copy over the slowness values here, so the stuff below is looking at the wrong index.
    chanTh, hitETh = [], []
    for i in range(n):
        if chan[i] in thD.keys() and hitE[i] > thD[chan[i]][0] + 3 * thD[chan[i]][1]:
            chanTh.append(chan[i])
            hitETh.append(hitE[i])
    chan, hitE = chanTh, hitETh
    nTot = len(hitE)

    n = len(hitE)
    hitEPass = [hitE[i] for i in range(nTot) if fSlo[i] < fsDict[chan[i]][2]]
    fSloPass = [fSlo[i] for i in range(nTot) if fSlo[i] < fsDict[chan[i]][2]]
    nPass = len(hitEPass)

    hitECut = [hitE[i] for i in range(nTot) if fSlo[i] > fsDict[chan[i]][2]]
    fSloCut = [fSlo[i] for i in range(nTot) if fSlo[i] > fsDict[chan[i]][2]]
    nCut = len(hitECut)

    # -- plot 1: 2d scatter plot, fSlo vs hitE, hits passing/failing v1 of fitSlo cut.
    fig = plt.figure()
    plt.cla()

    plt.plot(hitEPass, fSloPass, ".", ms=0.5, c='k', label='pass')
    plt.plot(hitECut, fSloCut, '.', ms=1.0, c='r', label='cut')

    plt.ylim(-10, 600)
    plt.xlim(0, 30)

    plt.xlabel("hitE (keV)", ha='right', x=1.)
    plt.ylabel("fitSlo", ha='right', y=1.)
    plt.legend(loc=1)
    plt.savefig("../plots/mult2-fitSloScatter.png")

    # -- plot 2: 1-d energy histo, hits passing/failing
    plt.cla()

    xLo, xHi, xpb = 0, 15, 0.2
    x, hTot = wl.GetHisto(hitE, xLo, xHi, xpb)
    x, hPass = wl.GetHisto(hitEPass, xLo, xHi, xpb)
    x, hCut = wl.GetHisto(hitECut, xLo, xHi, xpb)

    plt.plot(x, hTot, c='b', lw=1.5, ls='steps', label='m=%d, pk %d' % (mH, pk))
    plt.plot(x, hPass, c='g', lw=1.5, ls='steps', label='fitSlo pass')
    plt.plot(x, hCut, c='r', lw=1.5, ls='steps', label='fitSlo cut')

    plt.xlabel("hitE (keV)", ha='right', x=1.)
    plt.ylabel("Counts / %.1f kev" % xpb, ha='right', y=1.)
    plt.legend(loc=1)
    plt.savefig("../plots/mult2-m2s238-hitE.png")

    # -- plot 3: efficiency of v1 slowness cut, assuming all m2s238 hits are fast.
    plt.cla()

    from statsmodels.stats import proportion
    idx = np.where((hTot > 0) & (hPass > 0))
    ci_low, ci_upp = proportion.proportion_confint(hPass[idx], hTot[idx], alpha=0.1, method='beta')
    sloEff = hPass[idx] / hTot[idx]

    nPad = len(hPass)-len(hPass[idx])
    sloEff = np.pad(sloEff, (nPad,0), 'constant', constant_values=0)
    ci_low = np.pad(ci_low, (nPad,0), 'constant', constant_values=0)
    ci_upp = np.pad(ci_upp, (nPad,0), 'constant', constant_values=0)

    plt.plot(x, sloEff, '.b', ms=10., label='efficiency')
    plt.errorbar(x, sloEff, yerr=[sloEff - ci_low, ci_upp - sloEff], color='k', linewidth=0.8, fmt='none')

    idx = np.where(np.abs(ci_upp) > 0)
    xCut = x[6]
    plt.axvline(xCut,color='g',label='cutoff: %.2f keV' % xCut)

    plt.xlabel("hitE (keV)", ha='right', x=1.)
    plt.ylabel("Efficiency", ha='right', y=1.)
    plt.legend(loc=1)
    plt.savefig("../plots/mult2-m2s238-eff.png")


def plotChannelEff():

    # Idea: show a bar plot for efficiency of each detector (integrated over energy)

    # f = np.load("../data/mult2-sumLoEvts.npz")
    f = np.load("../data/mult2-sumLoEvts-histats.npz")
    runTime, sumEvts = f['arr_0'], f['arr_1'].item()

    # load fitSlo constants for cal run range closest to the run range [22513, 22566]
    dsNum, modNum, calIdx = 5, 1, 11  # calIdx 11: [[22568,22635],22568,22841],
    calDB = db.TinyDB('../calDB.json')
    pars = db.Query()
    fsDict = ds.getDBRecord("fitSlo_ds%d_idx%d_m%d_Peak" % (dsNum, calIdx, modNum), False, calDB, pars)

    mH, pk = 2, 238
    evts = sumEvts[mH][pk]
    chan = [evt[0] for evt in evts if evt[0] in fsDict.keys() and evt[2] > 0.7] # only plot detectors we have fs values for
    hitE = [evt[2] for evt in evts if evt[0] in fsDict.keys() and evt[2] > 0.7]
    fSlo = [evt[4] for evt in evts if evt[0] in fsDict.keys() and evt[2] > 0.7]
    n = len(hitE)

    chanTh = [chan[i] for i in range(n) if chan[i] in thD.keys() and hitE[i] > thD[chan[i]][0] + 3*thD[chan[i]][1]]
    hitETh = [hitE[i] for i in range(n) if chan[i] in thD.keys() and hitE[i] > thD[chan[i]][0] + 3*thD[chan[i]][1]]
    fSloTh = [fSlo[i] for i in range(n) if chan[i] in thD.keys() and hitE[i] > thD[chan[i]][0] + 3*thD[chan[i]][1]]
    n = len(hitETh)

    chanPass = [chanTh[i] for i in range(n) if fSloTh[i] < fsDict[chan[i]][2]]
    hitEPass = [hitETh[i] for i in range(n) if fSloTh[i] < fsDict[chan[i]][2]]
    fSloPass = [fSloTh[i] for i in range(n) if fSloTh[i] < fsDict[chan[i]][2]]

    chanMap = sorted(list(set(chanTh)))
    chanIdx = np.arange(0, len(chanMap), 1)
    chanDict = {chanMap[i]:i for i in range(len(chanIdx))}

    print(len(chanTh), len(chanPass))

    chanTh = [chanDict[ch] for ch in chanTh]
    chanPass = [chanDict[ch] for ch in chanPass]

    f = plt.figure()
    p1 = plt.subplot(111)

    x, hThr = wl.GetHisto(chanTh, 0, len(chanMap), 1.)
    p1.barh(x, hThr, 0.95, color='r', label='w/thresh')

    x, hPass = wl.GetHisto(chanPass, 0, len(chanMap), 1.)
    p1.barh(x, hPass, 0.95, color='b', label='+fitSlo')

    p1.set_xlabel("Counts", ha='right', x=1.)
    yticks = np.arange(0, len(chanMap)) + 0.5
    p1.set_yticks(yticks)
    p1.set_yticklabels(chanMap, fontsize=12)
    p1.set_ylabel("Channel", ha='right', y=1.)
    p1.set_ylim(-0.5)
    p1.legend()
    plt.savefig("../plots/mult2-channelEff.png")


def getExtPulser():
    """ Get ext pulser data and save it for comparison w/ the Compton events. """
    from ROOT import TFile, TChain, GATDataSet

    effs = [] # run, runTime, nExp, trigEff, cutEff
    evts = [] # chan, hitE, run, fSlo

    extPulserInfo = calInfo.GetSpecialList()["extPulserInfo"]
    syncChan = wl.getChan(0,10,0) # 672
    dsNum, modNum, calIdx = 0, 1, 33
    fsDXP = ds.getDBRecord("fitSlo_ds%d_idx%d_m%d_Peak" % (dsNum, calIdx, modNum), False, calDB, pars)

    for pIdx in [19, 20, 21]:

        runList = calInfo.GetSpecialRuns("extPulser", pIdx)
        extChan = extPulserInfo[pIdx][-1]
        fsCut = fsDXP[extChan][2]

        for i, run in enumerate(runList):
            if run in [7225, 7233]: continue

            gds = GATDataSet()
            gatPath = gds.GetPathToRun(run, GATDataSet.kGatified)
            tf = TFile(gatPath)
            gatTree = tf.Get("mjdTree")
            gatTree.GetEntry(0)
            runTime = gatTree.timeinfo.stopTime - gatTree.timeinfo.startTime
            tf.Close()

            pulseRate = 20 # Hz
            nExp = runTime * pulseRate

            fileList = ds.getLATRunList([run],"%s/lat" % (ds.specialDir))
            lTree = TChain("skimTree")
            # for f in fileList: lTree.Add("%s/%s" % (os.environ['SLURM_TMP'], f))
            for f in fileList: lTree.Add("%s/lat/%s" % (ds.specialDir,f))

            bNames = ["Entry$","mH","channel","trapENFCal","fitSlo"]
            theCut = "(channel==%d || channel==%d) && mH==2" % (syncChan, extChan) # enforce correct sync
            tvals = wl.GetVX(lTree, bNames, theCut)
            nPass = len(tvals["Entry$"])

            hitE = [tvals["trapENFCal"][i] for i in range(nPass) if tvals["channel"][i] == extChan]
            fSlo = [tvals["fitSlo"][i] for i in range(nPass) if tvals["channel"][i] == extChan]
            chan = [tvals["channel"][i] for i in range(nPass) if tvals["channel"][i] == extChan] # redundant, but who cares

            if len(hitE)==0:
                print("Run %d, no hits in channel %d found.  Continuing ..." % (run, extChan))
                continue

            for iEvt in range(len(hitE)):
                hit, fit, ch = hitE[iEvt], fSlo[iEvt], chan[iEvt]
                evts.append((hit,fit,run,ch))

            fPass = [fSlo[i] for i in range(len(fSlo)) if fSlo[i] < fsDXP[chan[i]][2]]
            cutEff = len(fPass)/len(fSlo)
            effs.append((run,runTime,nExp,cutEff,extChan))

            print(run,runTime,np.mean(hitE),nExp,cutEff,extChan)

    np.savez("../data/mult2-extPulser.npz", effs, evts)


def plotExtPulser():

    f = np.load("../data/mult2-extPulser.npz")
    effs, evts = f['arr_0'], f['arr_1']

    # effs: (run,runTime,nExp,cutEff,extChan)
    # evts: (hit,fit,run,ch)


if __name__=="__main__":
    main()
