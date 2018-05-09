#!/usr/bin/env python
import sys, os, imp, glob
import numpy as np
from scipy.optimize import curve_fit

import matplotlib as mpl
mpl.use('Agg')
sys.argv.append("-b")
import matplotlib.pyplot as plt
plt.style.use('../pltReports.mplstyle')
from matplotlib.colors import LogNorm, Normalize

# load LAT libraries
ds = imp.load_source('dsi',os.environ['LATDIR']+'/dsi.py')
wl = imp.load_source('waveLibs',os.environ['LATDIR']+'/waveLibs.py')
calInfo = ds.CalInfo()

# load threshold data
import tinydb as db
dsNum, bkgIdx = 5, 83
calDB = db.TinyDB('../calDB.json')
pars = db.Query()
thD = ds.getDBRecord("thresh_ds%d_bkgidx%d" % (dsNum, bkgIdx), False, calDB, pars)

# load fitSlo vals for cal run range closest to the run range [22513, 22566]
dsNum, modNum, calIdx = 5, 1, 11  # calIdx 11: [[22568,22635],22568,22841],
fsD = ds.getDBRecord("fitSlo_ds%d_idx%d_m%d_Peak" % (dsNum, calIdx, modNum), False, calDB, pars)

# load riseNoise vals
rnSD = ds.getDBRecord("riseNoise_ds%d_idx%d_m%d_SoftPlus" % (dsNum, calIdx, modNum), False, calDB, pars)
rnCD = ds.getDBRecord("riseNoise_ds%d_idx%d_m%d_Continuum" % (dsNum, calIdx, modNum), False, calDB, pars)


def main():

    # getSpec()
    # plotSpec()

    # plotRiseNoise()
    # tuneFitSlo()

    # testSimData()
    # plotSimTest()

    # getSimData()
    # plotSimData()

    # getSpecTest() # a test version of getSpec
    # plotSpecTest()

    # plotFitSlo()
    plotFitSloHist()

    # getSimDataEvtLoop()
    # plotSimDataLoop()

    # compareDataSim()


def getSpecTest():
    from ROOT import TFile, TTree
    fileList = []
    calInfo = ds.CalInfo()
    runList = calInfo.GetSpecialRuns("longCal",5) # 54 total
    runList = runList[:5]
    # runList = runList[:] # histats file, see filename below
    for run in runList: fileList.extend(ds.getLATRunList([run],"%s/lat" % (ds.specialDir)))
    nFiles = len(fileList)
    runTime = getRunTime(runList)
    chList = ds.GetGoodChanListNew(dsNum=5)

    eCut = 120
    hitList, hitData = [], []
    xLo, xHi, xpb = 0, 4000, 1
    nb = int((xHi-xLo)/xpb)+1
    sumSpec = {mHT:np.zeros(nb) for mHT in range(7)}
    hitSpec = {mHT:np.zeros(nb) for mHT in range(7)}

    hLo, hHi, hpb = 0, 20, 0.1
    nbh = int((hHi-hLo)/hpb)+1
    hitSpecLo = {mHT:np.zeros(nbh) for mHT in range(7)}

    for iFile, f in enumerate(fileList):
        print("%d/%d %s" % (iFile,nFiles,f))
        tf = TFile("%s/%s" % (os.environ['SLURM_TMP'], f))

        lTree = tf.Get("skimTree")
        bnames = ["startClockTime_s","clockTime_s","EventDC1Bits","sumEH",
        "mH","channel","trapENFCal","dtPulserCard","fitSlo","riseNoise"]
        lTree.SetBranchStatus('*',0)
        for name in bnames: lTree.SetBranchStatus(name,1)

        lTree.GetEntry(0)
        sumVals = {mHT:[] for mHT in range(7)}
        hitVals = {mHT:[] for mHT in range(7)}

        for iEnt in range(lTree.GetEntries()):
            lTree.GetEntry(iEnt)
            if lTree.EventDC1Bits != 0: continue

            n = lTree.channel.size()
            chTmp = np.asarray([lTree.channel.at(i) for i in range(n)])

            idxRaw = [i for i in range(lTree.channel.size())
                if lTree.channel.at(i) in chList
            ]
            hitERaw = np.asarray([lTree.trapENFCal.at(i) for i in idxRaw])

            idxList = [i for i in range(lTree.channel.size())
                if lTree.channel.at(i) in chList
                and lTree.channel.at(i) != 598
                and lTree.channel.at(i) in thD.keys()
                and lTree.trapENFCal.at(i) > thD[lTree.channel.at(i)][0] + 3*thD[lTree.channel.at(i)][1]
                and 0.7 < lTree.trapENFCal.at(i) < 9999
                ]
            hitE = np.asarray([lTree.trapENFCal.at(i) for i in idxList])

            mH, sumE = lTree.mH, lTree.sumEH

            mHT, sumET = len(hitE), sum(hitE)

            if 237.28 < sumET < 239.46 and mHT > 1:

                print("%-7d m %d %d  s %.2f %.2f  diff %.2f" % (iEnt, mH, mHT, sumE, sumET, sumET-sumE))
                print("  chan:",chTmp)
                print("   raw:",idxRaw," ".join("%.2f" % f for f in hitERaw))
                print("  pass:",idxList," ".join("%.2f" % f for f in hitE))

            if iEnt > 10000: exit(1)
            continue


            if len(idxList) == 0: continue

            if sumET > 250 or mHT < 2: continue

            if mHT < 7:
                sumVals[mHT].append(sumET)
                hitVals[mHT].extend(hitE)

            # if mHT > 1 and len(hitE[np.where(hitE < eCut)]) > 0:
            chan = np.asarray([int(lTree.channel.at(i)) for i in idxList])
            fSlo = np.asarray([lTree.fitSlo.at(i) for i in idxList])

            hitData.append([mHT, sumET])
            hitList.append([hitE, chan, fSlo])


        # save sum and hit spectra
        for i in range(7):
            x, y = wl.GetHisto(sumVals[i], xLo, xHi, xpb)
            sumSpec[i] = np.add(sumSpec[i], y)

            x, y = wl.GetHisto(hitVals[i], xLo, xHi, xpb)
            hitSpec[i] = np.add(hitSpec[i], y)

            xHit, y = wl.GetHisto(hitVals[i], hLo, hHi, hpb)
            hitSpecLo[i] = np.add(hitSpecLo[i], y)

    np.savez("../data/mult4-sumE-test.npz", runTime, x, sumSpec, hitSpec, xHit, hitSpecLo)
    np.savez("../data/mult4-hitE-test.npz", runTime, hitList, hitData, eCut)


def plotSpecTest():
    """
    sumSpec : {mHT: histos}
    hitData : [mHT, sumET, dt[mHT]]
    hitList : [hitE, chan, fSlo]  (same length as hitData)
    """
    f1 = np.load("../data/mult4-sumE-histats.npz")
    f2 = np.load("../data/mult4-hitE-histats.npz")
    runTime, x, sumSpec = f1['arr_0'], f1['arr_1'], f1['arr_2'].item()
    hitSpec, xHit, hitSpecLo = f1['arr_3'].item(), f1['arr_4'], f1['arr_5'].item()
    hitList, hitData, eCut = f2['arr_1'], f2['arr_2'], f2['arr_3']

    # fitSlo results from tuneFitSlo.
    fsVals = {
        584: 102.5, 592: 75.5, 608: 73.5, 610: 76.5, 614: 94.5, 624: 69.5,
        626: 81.5, 628: 102.5, 632: 81.5, 640: 73.5, 648: 74.5, 658: 75.5,
        660: 127.5, 662: 84.5, 672: 80.5, 678: 82.5, 680: 86.5, 688: 77.5,
        690: 80.5, 694: 80.5
        }
    chList = fsVals.keys()

    mHT = 2

    hitE, chan, fSlo = [], [], []
    for i in range(len(hitData)):
        if hitData[i][0]==mHT and 237.28 < hitData[i][1] < 239.46:
        # if hitData[i][0]==mHT and 235 < hitData[i][1] < 240:
            hitE.extend(hitList[i][0])
            chan.extend(hitList[i][1])
            fSlo.extend(hitList[i][2])
    n = len(hitE)
    hitE = [hitE[i] for i in range(n) if chan[i] in fsVals.keys()]
    fSloShift = [fSlo[i]-fsVals[chan[i]] for i in range(n) if chan[i] in chList]

    hitESlow = [hitE[i] for i in range(len(hitE)) if fSloShift[i] > 30]
    hitEFast = [hitE[i] for i in range(len(hitE)) if fSloShift[i] < 30]

    xLo, xHi, xpb = 0, 250, 1
    x, hSlo = wl.GetHisto(hitESlow,xLo,xHi,xpb)
    x, hFast = wl.GetHisto(hitEFast,xLo,xHi,xpb)

    # load sim data
    fs = np.load("../data/mult4-evtTrans.npz")
    hits = fs['arr_0'].item()
    evtTotal, evtBulk, evtTrans = hits[0], hits[1], hits[2]
    eneTotal, eneBulk, eneTrans = hits[3], hits[4], hits[5]
    x, hTotal = wl.GetHisto(evtTotal, xLo, xHi, xpb)
    x, hBulk = wl.GetHisto(evtBulk, xLo, xHi, xpb)
    x, hTrans = wl.GetHisto(evtTrans, xLo, xHi, xpb)
    x, hETotal = wl.GetHisto(eneTotal, xLo, xHi, xpb)
    x, hEBulk = wl.GetHisto(eneBulk, xLo, xHi, xpb)
    x, hETrans = wl.GetHisto(eneTrans, xLo, xHi, xpb)

    f = plt.figure()

    plt.plot(x, hEBulk/np.sum(hEBulk), ls='steps', c='g', lw=2., label='sim bulk')
    plt.plot(x, hFast/np.sum(hFast), ls='steps', c='b', lw=2., label='data fast')

    plt.legend(loc=3,fontsize=12)
    plt.savefig("../plots/mult4-238test.png")


def getSpec():
    """ Get sum and hit spectra w/ threshold cut, without ch. 598 (it's noisy.)
        Here we use "mHT" and "sumET" exclusively.
        Use !EventDC1Bits and trapENFCal > 0.7, all hits.  (Skim file was generated w/ 'dontSkipAnything')
        Save detailed info for events with hits < 100 keV.
    """
    from ROOT import TFile, TTree
    fileList = []
    calInfo = ds.CalInfo()
    runList = calInfo.GetSpecialRuns("longCal",5) # 54 total
    # runList = runList[:10]
    runList = runList[:] # histats file, see filename below
    for run in runList: fileList.extend(ds.getLATRunList([run],"%s/lat" % (ds.specialDir)))
    nFiles = len(fileList)
    runTime = getRunTime(runList)
    chList = ds.GetGoodChanListNew(dsNum=5)

    eCut = 120
    hitList, hitData = [], []
    xLo, xHi, xpb = 0, 4000, 1
    nb = int((xHi-xLo)/xpb)+1
    sumSpec = {mHT:np.zeros(nb) for mHT in range(7)}
    hitSpec = {mHT:np.zeros(nb) for mHT in range(7)}

    hLo, hHi, hpb = 0, 20, 0.1
    nbh = int((hHi-hLo)/hpb)+1
    hitSpecLo = {mHT:np.zeros(nbh) for mHT in range(7)}

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
        sumVals = {mHT:[] for mHT in range(7)}
        hitVals = {mHT:[] for mHT in range(7)}

        for iEnt in range(lTree.GetEntries()):
            lTree.GetEntry(iEnt)
            if lTree.EventDC1Bits != 0: continue

            # get iteration numbers for hits passing cuts
            idxList = [i for i in range(lTree.channel.size())
                if lTree.channel.at(i) in chList
                and lTree.channel.at(i) != 598
                and lTree.channel.at(i) in thD.keys()
                and lTree.trapENFCal.at(i) > thD[lTree.channel.at(i)][0] + 3*thD[lTree.channel.at(i)][1]
                and 0.7 < lTree.trapENFCal.at(i) < 9999
                ]
            if len(idxList) == 0: continue
            hitE = np.asarray([lTree.trapENFCal.at(i) for i in idxList])

            # save multiplicity, sum energy (using mHT and sumET only!)
            # mH, sumE = lTree.mH, lTree.sumEH <-- do not use
            mHT, sumET = len(hitE), sum(hitE)
            if mHT < 7:
                sumVals[mHT].append(sumET)
                hitVals[mHT].extend(hitE)

            # calculate dt since last mHT=N event
            clockTime = lTree.clockTime_s
            evtTime = clockTime - startTime
            if mHT not in prev:
                prev[mHT], dt[mHT] = 0, -1
            if prev[mHT] != 0:
                dt[mHT] = evtTime - prev[mHT]
            prev[mHT] = evtTime

            # save all information for events that have at least one hit under 'eCut'.
            if mHT > 1 and len(hitE[np.where(hitE < eCut)]) > 0:
                chan = np.asarray([int(lTree.channel.at(i)) for i in idxList])
                fSlo = np.asarray([lTree.fitSlo.at(i) for i in idxList])
                rise = np.asarray([lTree.riseNoise.at(i) for i in idxList])
                dtpc = np.asarray([lTree.dtPulserCard.at(i) for i in idxList])

                hitData.append([mHT, sumET, dt[mHT]])#, iFile, iEnt])
                hitList.append([hitE, chan, fSlo, rise, dtpc])

        # save sum and hit spectra
        for i in range(7):
            x, y = wl.GetHisto(sumVals[i], xLo, xHi, xpb)
            sumSpec[i] = np.add(sumSpec[i], y)

            x, y = wl.GetHisto(hitVals[i], xLo, xHi, xpb)
            hitSpec[i] = np.add(hitSpec[i], y)

            xHit, y = wl.GetHisto(hitVals[i], hLo, hHi, hpb)
            hitSpecLo[i] = np.add(hitSpecLo[i], y)

    # np.savez("../data/mult4-sumE.npz", runTime, x, sumSpec, hitSpec, xHit, hitSpecLo)
    # np.savez("../data/mult4-hitE.npz", runTime, hitList, hitData, eCut)

    np.savez("../data/mult4-sumE-histats.npz", runTime, x, sumSpec, hitSpec, xHit, hitSpecLo)
    np.savez("../data/mult4-hitE-histats.npz", runTime, hitList, hitData, eCut)


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


def plotSpec():
    """
    sumSpec : {mHT: histos}
    hitData : [mHT, sumET, dt[mHT]]
    hitList : [hitE, chan, fSlo, rise, dtpc]  (same length as hitData)
    """
    # f1 = np.load("../data/mult4-sumE.npz")
    # f2 = np.load("../data/mult4-hitE.npz")
    f1 = np.load("../data/mult4-sumE-histats.npz")
    f2 = np.load("../data/mult4-hitE-histats.npz")
    runTime, x, sumSpec = f1['arr_0'], f1['arr_1'], f1['arr_2'].item()
    hitSpec, xHit, hitSpecLo = f1['arr_3'].item(), f1['arr_4'], f1['arr_5'].item()
    hitList, hitData, eCut = f2['arr_1'], f2['arr_2'], f2['arr_3']

    xLo, xHi, xpb = 0, 4000, 1

    # counts & rates for each multiplicity (use mHT)
    nHits = [sum(sumSpec[i]) for i in range(7)]
    nErr = [np.sqrt(nHits[i]) for i in range(7)]
    nPct = [100 / nErr[i] if nHits[i] > 0 else 0 for i in range(7)]
    rate = [nHits[i] / runTime for i in range(7)]
    rErr = [nErr[i] / runTime for i in range(7)]

    # sum spectrum
    fig = plt.figure()
    xLo, xHi, xpb = 0, 4000, 2
    cols = [0,'r','b','m','c','g','k']
    for mHT in range(1,7):
        pLabel = r'mHT=%d %.2f $\pm$ %.3f Hz' % (mHT,rate[mHT],rErr[mHT])
        plt.semilogy(x, sumSpec[mHT], ls='steps', lw=1.5, c=cols[mHT],label=pLabel)
    plt.xlabel("sumET (keV)", ha='right', x=1.)
    plt.ylabel("Counts / %.1f keV" % (xpb), ha='right', y=1.)
    plt.legend(fontsize=12)
    plt.savefig("../plots/mult4-sumSpec.png")

    # sum spectrum, events w/ 1 or more hit under 'eCut' keV
    plt.cla()
    n = len(hitList)
    for mHT in range(2,7):
        evts = [hitData[i][1] for i in range(n) if hitData[i][0] == mHT]
        plt.semilogy(*wl.GetHisto(evts,xLo,xHi,xpb), ls='steps', lw=1.5, c=cols[mHT], label='mHT=%d, eCut:%.0f keV' % (mHT, eCut))
    plt.xlabel("sumE (keV)", ha='right', x=1.)
    plt.ylabel("Counts / %.1f keV" % (xpb), ha='right', y=1.)
    plt.legend(fontsize=12)
    plt.savefig("../plots/mult4-selectSpec.png")

    # hit spectrum
    plt.cla()
    for mHT in range(1,7):
        plt.semilogy(x, hitSpec[mHT], ls='steps', lw=1.5, c=cols[mHT], label='mHT=%d' % mHT)
    plt.xlabel("hitE (keV)", ha='right', x=1.)
    plt.ylabel("Counts / %.1f keV" % (xpb), ha='right', y=1.)
    plt.legend(fontsize=12)
    plt.savefig("../plots/mult4-hitSpec.png")

    # hit spectrum, under 20 kev
    plt.cla()
    for mHT in range(1,7):
        plt.semilogy(xHit, hitSpecLo[mHT], ls='steps', lw=1.5, c=cols[mHT], label='mHT=%d' % mHT)
    plt.xlabel("hitE (keV)", ha='right', x=1.)
    plt.ylabel("Counts / %.1f keV" % (xpb), ha='right', y=1.)
    plt.legend(fontsize=12)
    plt.savefig("../plots/mult4-hitSpecLow.png")


def plotRiseNoise():
    """ riseNoise removes electronics noise (pulsers, HF behavior).
    Its efficiency can be measured w/o a slow pulse cut.
    hitData : [mHT, sumET, dt[mHT]]
    hitList : [hitE, chan, fSlo, rise, dtpc]  (same length as hitData)
    """
    # f1 = np.load("../data/mult4-hitE.npz")
    f1 = np.load("../data/mult4-hitE-histats.npz")
    hitList, hitData, eCut = f1['arr_1'], f1['arr_2'], f1['arr_3']

    # this is pretty arbitrary, but it has high stats and no e-noise
    mHT = 2

    hitE, rise, chan = [], [], []
    for i in range(len(hitList)):
        if hitData[i][0] == mHT:
            idx = np.where(hitList[i][0] < 50)
            hitE.extend(hitList[i][0][idx])
            chan.extend(hitList[i][1][idx])
            rise.extend(hitList[i][3][idx])

    hitEFail, riseFail = [], []
    hitEPass, risePass = [], []
    for i in range(len(hitE)):
        if chan[i] not in rnSD.keys(): continue
        if rise[i] < 0 : continue

        a = max(rnSD[chan[i]][0], rnCD[chan[i]][4])
        b = rnSD[chan[i]][1]
        c = rnSD[chan[i]][2]
        d = rnSD[chan[i]][3]
        if d == 0: continue

        rnCut = a + b * np.log(1 + np.exp(hitE[i] - c/d))
        if rise[i] < rnCut:
            hitEPass.append(hitE[i])
            risePass.append(rise[i])
        else:
            hitEFail.append(hitE[i])
            riseFail.append(rise[i])

    rnEff = 100 * len(hitEPass) / (len(hitEPass) + len(hitEFail))

    plt.plot(hitEPass, risePass, '.', c='k', ms=1., label='pass')
    plt.plot(hitEFail, riseFail, '.', c='r', ms=5., label='fail')
    plt.plot(np.nan, np.nan, c='w', label='%.2f efficient' % rnEff )

    plt.xlabel("hitE (keV)", ha='right', x=1.)
    plt.ylabel("riseNoise", ha='right', y=1.)
    plt.legend(loc=1)
    plt.savefig("../plots/mult4-riseNoise.png")


def tuneFitSlo():
    """
    hitData : [mHT, sumET, dt[mHT]]
    hitList : [hitE, chan, fSlo, rise, dtpc]  (same length as hitData)
    """
    f1 = np.load("../data/mult4-hitE.npz")
    # f1 = np.load("../data/mult4-hitE-histats.npz")
    runTime, hitList, hitData, eCut = f1['arr_0'], f1['arr_1'], f1['arr_2'], f1['arr_3']

    # peak fit results from mult2.  could retune, but it shouldn't matter much
    mHT = 2
    sumPks = [
        (2,238,237.28,239.46), (2,583,581.26,584.46), (2,2615,2610.57,2618.01),
        (3,238,237.13,239.43), (3,583,581.04,584.36), (3,2615,2610.10,2617.92)
        ]

    hitE, chan, fSlo = [], [], []
    for i in range(len(hitData)):
        if hitData[i][0]==mHT and 237.28 < hitData[i][1] < 239.46:
            hitE.extend(hitList[i][0])
            chan.extend(hitList[i][1])
            fSlo.extend(hitList[i][2])

    plt.ylim(-5, 1000)
    plt.plot(hitE, fSlo, ".", ms=1., c='k')
    plt.savefig("../plots/mult4-fitSlo.png")

    plt.cla()
    chList = list(fsD.keys())

    xLo, xHi, xpb = -50, 150, 1
    fsTot = np.zeros(int((xHi-xLo)/xpb))
    xTot = np.arange(xLo, xHi, xpb)

    # save a dict of typical fast fitSlo values for each channel
    maxVals = {}

    n = len(hitE)
    cmap = plt.cm.get_cmap('hsv',len(chList)+1)
    for idx, ch in enumerate(chList):

        fsChan15 = [fSlo[i] for i in range(n) if chan[i]==ch and 30 < hitE[i] < 150 and 0 < fSlo[i] < 10000]
        if len(fsChan15)==0: continue

        fLo, fHi, fpb = 50, 200, 1
        x, y = wl.GetHisto(fsChan15,fLo,fHi,fpb)

        fMax = x[np.argmax(y)]
        avg = np.average(x, weights=y)
        stdv = np.sqrt(np.average((x-avg)**2, weights=y))
        maxVals[ch] = (fMax,avg,stdv)

        plt.plot(x-fMax,y,c=cmap(idx),lw=1.,ls='steps')

        xNew = x-fMax
        idxTot = np.where((xNew > xLo) & (xNew < xHi))
        yNew = y[idxTot]
        nPadLo = int(xNew[idxTot][0] - xTot[0])
        nPadHi = int(xTot[-1] - xNew[idxTot][-1])
        yPad = np.pad(yNew, (nPadLo,nPadHi), 'constant', constant_values=0)
        # print(ch, len(fsTot), len(yPad))

        fsTot = np.add(fsTot, yPad)
        maxVals[ch] = fMax

    plt.plot(xTot, fsTot, ls='steps')

    plt.xlabel("fitSlo",ha='right',x=1.)
    # plt.legend(loc=1,fontsize=14)
    plt.savefig("../plots/mult4-fitSloChan.png")

    print(maxVals)


def testSimData():
    from ROOT import TChain, TFile
    inDir = os.environ['SLURM_TMP']
    fileList = sorted(glob.glob("%s/DS3processed*.root" % inDir))
    fileList = fileList[:40]
    # fileList = fileList[:]

    hits = {i:[] for i in range(7)}
    actList = []

    totESpec = []
    sumHitESpec = []

    nPk, nCont = 0,0

    for i, f in enumerate(fileList):
        # print("%d/%d %s" % (i,len(fileList),f))
        print("%d/%d" % (i,len(fileList)))
        tf = TFile(f)
        simTree = tf.Get("simTree")

        for iEnt in range(simTree.GetEntries()):
            simTree.GetEntry(iEnt)

            ae = simTree.fAnalysisEvent
            sumE = ae.GetTotalEnergy()

            mH = ae.GetNElements()

            if 0.2 < sumE < 0.25 and mH==2:
                totESpec.append(sumE*1000)

                sumHitE = 0
                for iH in range(mH):
                    ele = ae.GetElement(iH)
                    ene = ele.GetEnergy()*1000
                    act = ele.GetActiveness()
                    if ene < 0.7 or np.isnan(ene) or np.isnan(act):
                        continue
                    sumHitE += ene
                sumHitESpec.append(sumHitE)

            if 0.23775 < sumE < 0.23924 and mH==2 : # sims window
            # if 0.23728 < sumE < 0.23946 and mH==2 : # data window

                hitE = []
                for iH in range(mH):
                    ele = ae.GetElement(iH)
                    ene = ele.GetEnergy()*1000
                    act = ele.GetActiveness()
                    # if ene < 0.7 or np.isnan(ene) or np.isnan(act):
                        # continue
                    actList.append(act)
                    hits[mH].append(ene/act)
                    hitE.append(ene)

                hitE = np.asarray(hitE)
                idx = np.where((hitE > 237.5) & (hitE < 239.2))
                if len(idx[0]) > 0:
                    nPk += 1
                    print("%-8d  pk   sumE %.2f  hits:[" % (iEnt, sumE*1000), "  ".join("%.2f" % e for e in hitE),"]")

                idx2 = np.where((hitE > 233.6) & (hitE < 235.55))
                if len(idx2[0]) > 0:
                    nCont += 1
                    print("%-8d  con  sumE %.2f  hits:[" % (iEnt, sumE*1000), "  ".join("%.2f" % e for e in hitE),"]")

    print("Peak cts: %d  Cont cts: %d" % (nPk,nCont))

    np.savez("../data/mult4-simTest.npz", hits, actList, totESpec, sumHitESpec)


def roughSigma(ene):
    """ from gpxFitter, just used to set initial guesses """
    p0, p1, p2 = 0.2, 0.02, 0.0003
    return np.sqrt(p0**2. + p1**2. * ene + p2**2. * ene**2.)


def gaus(x, b, a, mu, sig):
    """ gaussian + flat bg """
    return b + a * np.exp(-(x-mu)**2. / (2. * sig**2.))


def plotSimTest():

    f = np.load("../data/mult4-simTest.npz")
    hits = f['arr_0'].item()

    evts = hits[2]
    plt.plot(*wl.GetHisto(evts,0,250,1),ls='steps')
    plt.savefig("../plots/mult4-testPkSpec.png")

    # plt.cla()
    # act = f['arr_1']
    # idx = np.where(act < 0)
    # print(len(idx[0]))
    # plt.hist(act)
    # plt.savefig("../plots/mult4-activeness.png")

    plt.cla()
    totESpec, sumHitESpec = f['arr_2'], f['arr_3']
    plt.plot(*wl.GetHisto(totESpec,235,242,0.1),c='r',ls='steps',label="fTotalEnergy")
    plt.plot(*wl.GetHisto(sumHitESpec,235,242,0.1),c='b',ls='steps',label="fSumHitE")

    x, hSumE = wl.GetHisto(totESpec,235,242,0.1)
    popt,_ = curve_fit(gaus, x, hSumE, p0=(10,100,238,roughSigma(238.)))

    bgRate, mu, sig = popt[0] * (x[1]-x[0]), popt[2], popt[3]
    plt.plot(x, gaus(x, *popt), 'g-', label='fit mu %.1f sig %.1f' % (mu, sig))

    print(mu-2*sig, mu+2*sig)

    plt.legend()
    plt.savefig("../plots/mult4-simTotE.png")


def getSimData():
    from ROOT import TChain, TFile
    # inDir = "/global/projecta/projectdirs/majorana/users/mbuuck/sim/MJDWithGrahamTDL_byDetESmearing/5.0_TDL/0.75_transition_point/0.50_transition_level/MJDemonstrator/linesource/M1CalSource/A224_Z88"
    inDir = os.environ['SLURM_TMP']
    fileList = sorted(glob.glob("%s/DS3processed*.root" % inDir))
    # fileList = fileList[:10]
    fileList = fileList[:]

    hits = {i:[] for i in range(7)}

    for i, f in enumerate(fileList):
        print("%d/%d %s" % (i,len(fileList),f))
        tf = TFile(f)
        simTree = tf.Get("simTree")

        # use same peak parameters as above
        # theCut = "fNWaveforms==2 && fTotalEnergy/fActiveness > 0.23728 && fTotalEnergy/fActiveness < 0.23946"
        theCut = "fNWaveforms==2 && fTotalEnergy > 0.23728 && fTotalEnergy < 0.23946 && fEnergy > 0.0001"

        # these are the 'true' (simulated) energies

        n0 = simTree.Draw('fEnergy*1000/fActiveness', theCut, 'goff')
        eTot = simTree.GetV1()
        hits[0].extend([eTot[i] for i in range(n0)])

        n1 = simTree.Draw('fEnergy*1000/fActiveness', theCut + '&& fActiveness == 1', 'goff')
        eBulk = simTree.GetV1()
        hits[1].extend([eBulk[i] for i in range(n1)])

        n2 = simTree.Draw('fEnergy*1000/fActiveness', theCut + '&& fActiveness < 1', 'goff')
        eTrans = simTree.GetV1()
        hits[2].extend([eTrans[i] for i in range(n2)])

        # these are the 'degraded' (observed) energies

        n3 = simTree.Draw('fEnergy*1000', theCut, 'goff')
        enTot = simTree.GetV1()
        hits[3].extend([enTot[i] for i in range(n3)])

        n4 = simTree.Draw('fEnergy*1000', theCut + '&& fActiveness == 1', 'goff')
        enBulk = simTree.GetV1()
        hits[4].extend([enBulk[i] for i in range(n4)])

        n5 = simTree.Draw('fEnergy*1000', theCut + '&& fActiveness < 1', 'goff')
        enTrans = simTree.GetV1()
        hits[5].extend([enTrans[i] for i in range(n5)])

        tf.Close()

    for i in hits:
        hits[i] = np.asarray(hits[i])
        idx = np.where((~np.isnan(hits[i])) & (~np.isinf(hits[i])))
        if len(idx[0]) != len(hits[i]):
            print("hit set %d, found %d bad vals in %d entries." % (i, len(hits[i])-len(idx[0]),len(hits[i])))
            hits[i] = hits[i][idx]
        print("hit set %d, %d entries" % (i, len(hits[i])))

    np.savez("../data/mult4-evtTrans.npz", hits)


def plotSimData():

    f = np.load("../data/mult4-evtTrans.npz")
    hits = f['arr_0'].item()
    evtTotal, evtBulk, evtTrans = hits[0], hits[1], hits[2]
    eneTotal, eneBulk, eneTrans = hits[3], hits[4], hits[5]

    print(len(evtTotal),len(evtBulk),len(evtTrans))
    print(len(eneTotal),len(eneBulk),len(eneTrans))

    xLo, xHi, xpb = 0, 250, 1
    x, hTotal = wl.GetHisto(evtTotal, xLo, xHi, xpb)
    x, hBulk = wl.GetHisto(evtBulk, xLo, xHi, xpb)
    # x, hTrans = wl.GetHisto(evtTrans, xLo, xHi, xpb)
    x, eTotal = wl.GetHisto(eneTotal, xLo, xHi, xpb)
    x, eBulk = wl.GetHisto(eneBulk, xLo, xHi, xpb)
    # x, eTrans = wl.GetHisto(eneTrans, xLo, xHi, xpb)

    plt.plot(x, hTotal, ls='steps', c='r', lw=2., label='h total')
    plt.plot(x, hBulk, ls='steps', c='g', lw=2., label='h bulk')

    # plt.plot(x, hTrans, ls='steps', c='b', lw=2., label='transition')

    plt.plot(x, hTotal, ls='steps', c='b', lw=2., label='e total')
    plt.plot(x, hBulk, ls='steps', c='k', lw=2., label='e bulk')

    plt.xlabel("fEnergy (keV)", ha='right', x=1.)
    plt.ylabel("Simulated Counts", ha='right', y=1.)
    plt.legend(loc=4)
    plt.savefig("../plots/mult4-transSpec.png")

    # plt.cla()
    # import seaborn as sns
    # transPct = 100 * (np.divide(hTrans, hTotal, dtype=float))
    # sns.regplot(x=x, y = transPct, scatter_kws={'s':20})
    # plt.xlabel("fEnergy (keV)", ha='right', x=1.)
    # plt.ylabel("% of Transition Layer Events", ha='right', y=1.)
    # plt.savefig("../plots/mult4-transPct.png")


def plotFitSlo():
    """ hitData : [mHT, sumET, dt[mHT]]
        hitList : [hitE, chan, fSlo, rise, dtpc]  (same length as hitData)
    """

    # fitSlo results from tuneFitSlo.
    fsVals = {
        584: 102.5, 592: 75.5, 608: 73.5, 610: 76.5, 614: 94.5, 624: 69.5,
        626: 81.5, 628: 102.5, 632: 81.5, 640: 73.5, 648: 74.5, 658: 75.5,
        660: 127.5, 662: 84.5, 672: 80.5, 678: 82.5, 680: 86.5, 688: 77.5,
        690: 80.5, 694: 80.5
        }

    chList = list(fsVals.keys())

    # peak fit results from mult2.  could retune, but it shouldn't matter much
    mHT = 2
    sumPks = [
        (2,238,237.28,239.46), (2,583,581.26,584.46), (2,2615,2610.57,2618.01),
        (3,238,237.13,239.43), (3,583,581.04,584.36), (3,2615,2610.10,2617.92)
        ]

    # f1 = np.load("../data/mult4-hitE.npz")
    f1 = np.load("../data/mult4-hitE-histats.npz")
    runTime, hitList, hitData, eCut = f1['arr_0'], f1['arr_1'], f1['arr_2'], f1['arr_3']

    hitE, chan, fSlo = [], [], []
    for i in range(len(hitData)):
        if hitData[i][0]==mHT and 237.28 < hitData[i][1] < 239.46:
            hitE.extend(hitList[i][0])
            chan.extend(hitList[i][1])
            fSlo.extend(hitList[i][2])
    n = len(hitE)
    hitE = [hitE[i] for i in range(n) if chan[i] in fsVals.keys()]
    fSloShift = [fSlo[i]-fsVals[chan[i]] for i in range(n) if chan[i] in chList]
    print("peak evts:",n)

    plt.cla()
    plt.plot(hitE, fSloShift, '.', c='k', ms=1.)
    plt.ylim(-50, 400)
    plt.savefig("../plots/mult4-fitSlo-shift.png")

    cLo, cHi = 0, 230

    hitE2, chan2, fSlo2 = [], [], []
    for i in range(len(hitData)):
        if hitData[i][0]==mHT and cLo < hitData[i][1] < cHi:
            hitE2.extend(hitList[i][0])
            chan2.extend(hitList[i][1])
            fSlo2.extend(hitList[i][2])
    n = len(hitE2)
    hitE2 = [hitE2[i] for i in range(n) if chan2[i] in fsVals.keys()]
    fSloShift2 = [fSlo2[i]-fsVals[chan2[i]] for i in range(n) if chan2[i] in chList]
    print("cont evts:",n)

    plt.cla()
    plt.plot(hitE2, fSloShift2, '.', c='k', ms=1.)
    plt.ylim(-50, 400)
    plt.savefig("../plots/mult4-fitSlo-shift2.png")


    plt.cla()
    xLo, xHi, xpb = -50, 300, 1

    x, y238 = wl.GetHisto(fSloShift,xLo,xHi,xpb)
    x, yCon = wl.GetHisto(fSloShift2,xLo,xHi,xpb)

    plt.semilogy(x, y238/np.sum(y238),'b',lw=2.,ls='steps',label='m=2 238 pk')
    plt.semilogy(x, yCon/np.sum(yCon),'r',lw=2.,ls='steps',label='%d - %d kev' % (cLo, cHi))
    plt.legend()
    plt.savefig("../plots/mult4-fitSlo-hist.png")


    # plt.cla()
    # xLo, xHi, xpb = 0, 240, 1
    #
    # hitESlow = [hitE[i] for i in range(len(hitE)) if fSloShift[i] > 30]
    # x, hSlo = wl.GetHisto(hitESlow,xLo,xHi,xpb)
    #
    # hitEFast = [hitE[i] for i in range(len(hitE)) if fSloShift[i] < 30]
    # x, hFast = wl.GetHisto(hitEFast,xLo,xHi,xpb)
    #
    # fs = np.load("../data/mult4-evtTrans.npz")
    # evtTotal, evtBulk, evtTrans = fs['arr_0'], fs['arr_1'], fs['arr_2']
    # eneTotal, eneBulk, eneTrans = fs['arr_3'], fs['arr_4'], fs['arr_5']
    #
    # x, hTotal = wl.GetHisto(evtTotal, xLo, xHi, xpb)
    # x, hBulk = wl.GetHisto(evtBulk, xLo, xHi, xpb)
    # x, hTrans = wl.GetHisto(evtTrans, xLo, xHi, xpb)
    #
    # x, hETotal = wl.GetHisto(eneTotal, xLo, xHi, xpb)
    # x, hEBulk = wl.GetHisto(eneBulk, xLo, xHi, xpb)
    # x, hETrans = wl.GetHisto(eneTrans, xLo, xHi, xpb)
    #
    # plt.plot(x, hETotal/np.sum(hETotal), ls='steps', c='r', lw=2., label='sim total')
    # plt.plot(x, hEBulk/np.sum(hEBulk), ls='steps', c='g', lw=2., label='sim bulk')
    # plt.plot(x, hFast/np.sum(hFast), ls='steps', c='b', lw=2., label='data fast')
    #
    # plt.legend(loc=1,fontsize=12)
    # plt.savefig("../plots/mult4-hitSlow.png")


def plotFitSloHist():
    """ hitData : [mHT, sumET, dt[mHT]]
        hitList : [hitE, chan, fSlo, rise, dtpc]  (same length as hitData)
    """

    # fitSlo results from tuneFitSlo.
    fsVals = {
        584: 102.5, 592: 75.5, 608: 73.5, 610: 76.5, 614: 94.5, 624: 69.5,
        626: 81.5, 628: 102.5, 632: 81.5, 640: 73.5, 648: 74.5, 658: 75.5,
        660: 127.5, 662: 84.5, 672: 80.5, 678: 82.5, 680: 86.5, 688: 77.5,
        690: 80.5, 694: 80.5
        }

    chList = list(fsVals.keys())

    # peak fit results from mult2.  could retune, but it shouldn't matter much
    mHT = 2
    sumPks = [
        (2,238,237.28,239.46), (2,583,581.26,584.46), (2,2615,2610.57,2618.01),
        (3,238,237.13,239.43), (3,583,581.04,584.36), (3,2615,2610.10,2617.92)
        ]

    # f1 = np.load("../data/mult4-hitE.npz")
    f1 = np.load("../data/mult4-hitE-histats.npz")
    runTime, hitList, hitData, eCut = f1['arr_0'], f1['arr_1'], f1['arr_2'], f1['arr_3']

    hitE, chan, fSlo = [], [], []
    for i in range(len(hitData)):
        if hitData[i][0]==mHT and 237.28 < hitData[i][1] < 239.46:
            hitE.extend(hitList[i][0])
            chan.extend(hitList[i][1])
            fSlo.extend(hitList[i][2])
    n = len(hitE)
    hitE = [hitE[i] for i in range(n) if chan[i] in fsVals.keys()]
    fSloShift = [fSlo[i]-fsVals[chan[i]] for i in range(n) if chan[i] in chList]
    print("peak evts:",n)

    # plt.cla()
    # plt.plot(hitE, fSloShift, '.', c='k', ms=1.)

    fig = plt.figure()

    xLo, xHi, xpb = 0, 250, 1
    nbx = int((xHi-xLo)/xpb)
    yLo, yHi, ypb = -50, 400, 1
    nby = int((yHi-yLo)/ypb)
    _,_,_,im = plt.hist2d(hitE, fSloShift, bins=[nbx, nby], range=[[xLo,xHi],[yLo,yHi]], norm=LogNorm(), cmap='jet', label='m=2, sumE=238 hits')
    cb = plt.colorbar()
    plt.xlabel("Energy (keV)", ha='right', x=1.)
    plt.ylabel("fitSlo", ha='right', y=1.)
    plt.tight_layout()
    plt.savefig("../plots/mult4-fitSlo-shift-hist.png")

    cLo, cHi = 0, 230

    hitE2, chan2, fSlo2 = [], [], []
    for i in range(len(hitData)):
        if hitData[i][0]==mHT and cLo < hitData[i][1] < cHi:
            hitE2.extend(hitList[i][0])
            chan2.extend(hitList[i][1])
            fSlo2.extend(hitList[i][2])
    n = len(hitE2)
    hitE2 = [hitE2[i] for i in range(n) if chan2[i] in fsVals.keys()]
    fSloShift2 = [fSlo2[i]-fsVals[chan2[i]] for i in range(n) if chan2[i] in chList]
    print("cont evts:",n)

    cb.remove()
    plt.cla()

    # plt.plot(hitE2, fSloShift2, '.', c='k', ms=1.)
    _,_,_,im = plt.hist2d(hitE2, fSloShift2, bins=[nbx, nby], range=[[xLo,xHi],[yLo,yHi]], norm=LogNorm(), cmap='jet', label='m=2, hits 0-230 keV')
    cb = plt.colorbar()
    plt.xlabel("Energy (keV)", ha='right', x=1.)
    plt.ylabel("fitSlo", ha='right', y=1.)
    plt.tight_layout()
    plt.savefig("../plots/mult4-fitSlo-shift2-hist.png")


    cb.remove()
    plt.cla()
    xLo, xHi, xpb = -50, 300, 1

    x, y238 = wl.GetHisto(fSloShift,xLo,xHi,xpb)
    x, yCon = wl.GetHisto(fSloShift2,xLo,xHi,xpb)

    # integral238
    tot238 = np.sum(y238)
    int238, x238 = 0, 0
    for i in range(len(y238)):
        int238 += y238[i]
        if int238/tot238 > 0.95:
            x238 = x[i]
            break

    # integralCon
    totCon = np.sum(yCon)
    intCon, xCon = 0, 0
    for i in range(len(yCon)):
        intCon += yCon[i]
        if intCon/totCon > 0.95:
            xCon = x[i]
            break

    plt.semilogy(x, y238/np.sum(y238),'b',lw=2.,ls='steps',label='m=2 238 pk')
    plt.semilogy(x, yCon/np.sum(yCon),'r',lw=2.,ls='steps',label='%d - %d kev' % (cLo, cHi))

    plt.axvline(xCon, c='m', label='95% m=2')
    plt.axvline(x238, c='g', label='95% m=2,s=238')

    plt.xlabel("fitSlo", ha='right', x=1.)
    plt.ylabel("Counts", ha='right', y=1.)

    plt.legend(loc=1)
    plt.savefig("../plots/mult4-fitSlo-hist.png")



def getSimDataEvtLoop():
    from ROOT import TChain, TFile
    inDir = os.environ['SLURM_TMP']
    fileList = sorted(glob.glob("%s/DS3processed*.root" % inDir))
    # fileList = fileList[:3]
    fileList = fileList[:100]

    eTotHits = {mH:[] for mH in range(7)} # plot fEnergy (this is 'incident' energy)
    eActHits = {mH:[] for mH in range(7)} # plot fEnergy/fActiveness (this is 'observed' energy)

    nPk, nCont = 0,0

    for i, f in enumerate(fileList):
        # print("%d/%d %s" % (i,len(fileList),f))
        print("%d/%d" % (i,len(fileList)))
        tf = TFile(f)
        simTree = tf.Get("simTree")

        for iEnt in range(simTree.GetEntries()):
            simTree.GetEntry(iEnt)
            ae = simTree.fAnalysisEvent
            sumE = ae.GetTotalEnergy()*1000
            mH = ae.GetNElements()

            if not (mH==2 and 237.28 < sumE < 239.46):
                continue

            for iH in range(mH):
                evt = ae.GetElement(iH)
                hitE = evt.GetEnergy()*1000
                fAct = evt.GetActiveness()

                if hitE < 0.7 or np.isnan(hitE) or np.isnan(fAct): continue
                eTotHits[mH].append(hitE)
                eActHits[mH].append(fAct)

    np.savez("../data/mult4-simDataLoop.npz", eTotHits, eActHits)


def plotSimDataLoop():

    f = np.load("../data/mult4-simDataLoop.npz")
    eTotHits, eActHits = f['arr_0'].item(), f['arr_1'].item()
    mH = 2
    n = len(eTotHits[mH])

    eTotal = [eTotHits[mH][i]/eActHits[mH][i] for i in range(n)]
    eBulk = [eTotHits[mH][i]/eActHits[mH][i] for i in range(n) if abs(eActHits[mH][i]-1) < 0.001]
    eDead = [eTotHits[mH][i]/eActHits[mH][i] for i in range(n) if abs(eActHits[mH][i]) < 1.]

    xLo, xHi, xpb = 0, 250, 1
    x, hTotal = wl.GetHisto(eTotal,xLo,xHi,xpb)
    x, hBulk = wl.GetHisto(eBulk,xLo,xHi,xpb)
    x, hDead = wl.GetHisto(eDead,xLo,xHi,xpb)

    print(np.mean(hTotal), np.mean(hDead))

    fig = plt.figure()
    plt.plot(x, hTotal, ls='steps', c='b', label='all sim hits, sumE=238')
    plt.plot(x, hBulk, ls='steps', c='g', label='sim bulk hits')
    plt.plot(x, hDead, ls='steps', c='r', label='incident E, transition evts')

    plt.xlabel("hitE (keV)", ha='right', x=1.)
    plt.ylabel("Counts (A.U.)", ha='right', y=1.)
    plt.legend(loc=1)
    plt.savefig("../plots/mult4-sim238hitSpec.png")


def compareDataSim():

    xLo, xHi, xpb = 0, 240, 1

    # sim
    f = np.load("../data/mult4-simDataLoop.npz")
    eTotHits, eActHits = f['arr_0'].item(), f['arr_1'].item()
    mH = 2
    n = len(eTotHits[mH])
    # simTotal = [eTotHits[mH][i]/eActHits[mH][i] for i in range(n)]
    # simBulk = [eTotHits[mH][i]/eActHits[mH][i] for i in range(n) if abs(eActHits[mH][i]-1) < 0.001]
    # simDead = [eTotHits[mH][i]/eActHits[mH][i] for i in range(n) if abs(eActHits[mH][i]) < 1.]

    simTotal = [eTotHits[mH][i] for i in range(n)]
    simBulk = [eTotHits[mH][i] for i in range(n) if abs(eActHits[mH][i]-1) < 0.001]
    simDead = [eTotHits[mH][i] for i in range(n) if abs(eActHits[mH][i]) < 1.]

    x, sTotal = wl.GetHisto(simTotal,xLo,xHi,xpb)
    x, sBulk = wl.GetHisto(simBulk,xLo,xHi,xpb)
    x, sDead = wl.GetHisto(simDead,xLo,xHi,xpb)

    # data
    mHT = 2
    fsVals = {
        584: 102.5, 592: 75.5, 608: 73.5, 610: 76.5, 614: 94.5, 624: 69.5,
        626: 81.5, 628: 102.5, 632: 81.5, 640: 73.5, 648: 74.5, 658: 75.5,
        660: 127.5, 662: 84.5, 672: 80.5, 678: 82.5, 680: 86.5, 688: 77.5,
        690: 80.5, 694: 80.5
        }
    chList = list(fsVals.keys())
    # f1 = np.load("../data/mult4-hitE.npz")
    f1 = np.load("../data/mult4-hitE-histats.npz")
    runTime, hitList, hitData, eCut = f1['arr_0'], f1['arr_1'], f1['arr_2'], f1['arr_3']
    hitE, chan, fSlo = [], [], []
    for i in range(len(hitData)):
        if hitData[i][0]==mHT and 237.28 < hitData[i][1] < 239.46:
            hitE.extend(hitList[i][0])
            chan.extend(hitList[i][1])
            fSlo.extend(hitList[i][2])
    n = len(hitE)
    hitE = [hitE[i] for i in range(n) if chan[i] in fsVals.keys()]
    fSloShift = [fSlo[i]-fsVals[chan[i]] for i in range(n) if chan[i] in chList]
    print("peak evts:",n)

    hitESlow = [hitE[i] for i in range(len(hitE)) if fSloShift[i] > 30]
    hitEFast = [hitE[i] for i in range(len(hitE)) if fSloShift[i] < 30]

    x, hSlow = wl.GetHisto(hitESlow,xLo,xHi,xpb)
    x, hFast = wl.GetHisto(hitEFast,xLo,xHi,xpb)

    # -----------------------------------------------

    # plt.plot(x, sBulk/np.sum(sBulk), ls='steps', c='r', lw=2., label='sim bulk')
    # plt.plot(x, hFast/np.sum(hFast), ls='steps', c='b', lw=2., label='data fast')
    # plt.plot(x, sDead/np.sum(sBulk), ls='steps', c='m', lw=2., label='sim dead')
    # plt.plot(x, hSlow/np.sum(hFast), ls='steps', c='g', lw=2., label='data slow')

    # plt.plot(x, sBulk/np.sum(sBulk), ls='steps', c='r', lw=2., label='sim bulk')
    # plt.plot(x, hFast/np.sum(hFast), ls='steps', c='b', lw=2., label='data fast')
    plt.plot(x, sDead/np.sum(sDead), ls='steps', c='m', lw=2., label='sim dead')
    plt.plot(x, hSlow/np.sum(hSlow), ls='steps', c='g', lw=2., label='slow data, m=2,s=238')

    plt.xlabel("hitE (keV)", ha='right', x=1.)
    plt.ylabel("Normalized Cts", ha='right', y=1.)
    plt.legend(loc=1,fontsize=12)
    plt.savefig("../plots/mult4-dataSimSpec.png")




if __name__=="__main__":
    main()