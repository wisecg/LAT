#!/usr/bin/env python3
import sys, os
import numpy as np
import matplotlib.pyplot as plt
plt.style.use('../clint.mpl')

import waveLibs as wl
import dsi
bkg = dsi.BkgInfo()
det = dsi.DetInfo()

def main(argv):

    # checkLTOutput()
    getRates()
    # printRates()
    # plotRates()
    # getOutliers()


def checkLTOutput():
    """ NOTE: Livetime for low energy uses HG channels only, so the results
    are slightly (O(1 minute)) different from the OR mode used in the 0nbb analysis.
    """
    from ROOT import TFile, TTree
    ds = 2
    nBkg = bkg.dsMap()[ds]
    chList = det.getGoodChanList(ds)
    runRanges = bkg.getRanges(ds)

    tl = TFile("../data/ds_%s_livetime.root" % str(ds))
    lt = tl.Get("dsTree")
    livetimeTot = {ch:0 for ch in chList}

    for bIdx in range(nBkg+1):

        rLo, rHi = runRanges[bIdx][0], runRanges[bIdx][-1]
        print("DS-%d  bIdx %d  runLo %d  runHi %d" % (ds, bIdx, rLo, rHi))

        # get the livetime of each channel in this bkgIdx
        livetime = {ch:0 for ch in chList}
        n = lt.Draw("run:channel:livetime","run>=%d && run<=%d" % (rLo, rHi), 'goff')
        ltRun, ltChan, ltLive = lt.GetV1(), lt.GetV2(), lt.GetV3()
        for i in range(n):
            livetime[ltChan[i]] += ltLive[i]
            livetimeTot[ltChan[i]] += ltLive[i]
            # print("run %d  chan %d  lt %.1f  tot %.1f" % (ltRun[i], ltChan[i], ltLive[i], livetime[ltChan[i]]))

    print("summary")
    for ch in chList:
        print("chan %d  lt %.4f" % (ch, livetimeTot[ch]/86400))


def getRates():

    from ROOT import TFile, TTree, gROOT
    gROOT.ProcessLine("gErrorIgnoreLevel = 3001;") # suppress ROOT error messages

    # for ds in [0,1,2,3,4,"5A","5B","5C"]:
    for ds in ["5B"]:

        dsNum = int(ds[0]) if isinstance(ds,str) else ds
        nBkg = bkg.dsMap()[dsNum]
        bLo, bHi = 0, nBkg
        if ds=="5A": bLo, bHi = 0, 79
        if ds=="5B": bLo, bHi = 80, 112
        if ds=="5C": bLo, bHi = 113, 121
        chList = det.getGoodChanList(dsNum)
        runRanges = bkg.getRanges(dsNum)
        tl = TFile("../data/ds_%s_livetime.root" % str(ds))
        lt = tl.Get("dsTree")

        dsRate = {det.getChanCPD(dsNum,ch):[] for ch in chList}

        for bIdx in range(bLo, bHi+1):

            rLo, rHi = runRanges[bIdx][0], runRanges[bIdx][-1]
            # print("DS-%s  bIdx %d  runLo %d  runHi %d" % (ds, bIdx, rLo, rHi))

            # get the livetime
            livetime = {ch:0 for ch in chList}
            n = lt.Draw("run:channel:livetime","run>=%d && run<=%d" % (rLo, rHi), 'goff')
            ltRun, ltChan, ltLive = lt.GetV1(), lt.GetV2(), lt.GetV3()
            for i in range(n):
                livetime[ltChan[i]] += ltLive[i]

            # get the exposure
            exposure = {}
            for ch in chList:
                detID = det.getDetIDChan(dsNum,ch)
                aMass = det.allActiveMasses[detID]
                exposure[ch] = livetime[ch]*aMass/86400/1000
                # print(ch, exposure[ch])

            # now access the cut files and look at the rates
            for ch in sorted(chList):
                cpd = det.getChanCPD(dsNum,ch)

                fName = "%s/bkg/cut/fr/fr_ds%d_%d_ch%d.root" % (dsi.dataDir, dsNum, bIdx, ch)

                # skip nonexistent files.  other parts of chsel should tell us why these aren't here.
                if not os.path.isfile(fName):
                    # print("nope",fName)
                    continue

                tf = TFile(fName)
                tt = tf.Get("skimTree")
                nEvt = tt.GetEntries()
                live = livetime[ch]/86400
                if exposure[ch] == 0:
                    print("'Dead' detector: DS-%s  bIdx %d  cpd %s  ch %d exposure 0, hits %d" % (ds, bIdx, cpd, ch, nEvt))
                    continue

                n10 = tt.Draw("Entry$:Iteration$","trapENFCal>=0 && trapENFCal<=10","goff")
                n250 = tt.Draw("Entry$:Iteration$","trapENFCal>=10 && trapENFCal<=250","goff")
                r10 = n10/exposure[ch]/10
                r250 = n250/exposure[ch]/240

                print("DS %s  bIdx %d  ch %d  cpd %s  exp %.1f - %d  %d  %.1f  %.1f" % (ds, bIdx, ch, cpd, exposure[ch], n10,n250,r10,r250))

                dsRate[cpd].append([r10,r250])
                # print("b %-2d  %s  live %-8.2f  expo %-8.2f  nEvt %-4d  r10 %-8.2f  rAll %-8.2f" % (bIdx, cpd, live, expo, nEvt, r10, r250))

        np.savez('../data/cs3-rates-ds%s.npz' % ds, dsRate, ds)

        printRates(ds)


def printRates(ds):

    dsNum = int(ds[0]) if isinstance(ds,str) else ds

    f = np.load('../data/cs3-rates-ds%s.npz' % ds)
    dsRate = f['arr_0'].item()

    enr10, enrAll = [], []
    nat10, natAll = [], []

    print("Detector summary")
    for cpd in sorted(dsRate):

        rate10 = [r[0] for r in dsRate[cpd]]
        mu10, std10 = 0, 0
        if len(rate10) > 0:
            mu10, std10 = np.mean(rate10), np.std(rate10)

        rate250 = [r[1] for r in dsRate[cpd]]
        muAll, stdAll = 0, 0
        if len(rate250) > 0:
            muAll, stdAll = np.mean(rate250), np.std(rate250)

        detID = det.allDetIDs[cpd]
        pct10 = std10/mu10 if mu10 !=0 else -1
        pctAll = stdAll/muAll if muAll !=0 else -1
        mu10All = mu10/muAll if muAll !=0 else -1
        print("cpd %s  ch %-4d  %-8s  r10: %.2f ± %.2f (%.2f) %-1s  rAll: %.2f ± %.2f (%.2f) %-1s  %.2f" % (cpd, det.getCPDChan(dsNum,cpd), detID, mu10, std10, pct10, " ", muAll, stdAll, pctAll, " ", mu10All))

        isEnr = True if detID > 100000 else False
        if isEnr:
            enr10.append(mu10)
            enrAll.append(muAll)
        else:
            nat10.append(mu10)
            natAll.append(muAll)

    print("Totals: ")
    enrMean10, enrStd10 = np.mean(enr10), np.std(enr10)
    enrMeanAll, enrStdAll = np.mean(enrAll), np.std(enrAll)
    natMean10, natStd10 = np.mean(nat10), np.std(nat10)
    natMeanAll, natStdAll = np.mean(natAll), np.std(natAll)
    print("Enriched:  rate10 - mu %-6.2f pm %-6.2f  rate250 - mu %-6.2f pm %-6.2f" % (enrMean10,enrStd10,enrMeanAll,enrStdAll))
    print("Natural:   rate10 - mu %-6.2f pm %-6.2f  rate250 - mu %-6.2f pm %-6.2f" % (natMean10,natStd10,natMeanAll,natStdAll))


def plotRates():
    """
    https://www.itl.nist.gov/div898/handbook/prc/section1/prc16.htm
    https://matplotlib.org/2.0.1/examples/pylab_examples/boxplot_demo.html
    tukey boxplot of rates in enr/nat vs dataset
    """
    ds = "5A"

    for i, ds in enumerate([0,1,2]):

        f = np.load('../data/cs3-rates-ds%s.npz' % ds)
        dsRate = f['arr_0'].item()

        def jitter(vals,idx,sc):
            return [idx + np.random.randn()*sc for v in vals]

        enrRatio, natRatio = [], []

        for cpd in sorted(dsRate):

            enrRatio.extend([r[0]/r[1] for r in dsRate[cpd] if det.allDetIDs[cpd] > 100000 and r[1]!=0])
            natRatio.extend([r[0]/r[1] for r in dsRate[cpd] if det.allDetIDs[cpd] < 100000 and r[1]!=0])

        plt.boxplot([enrRatio, natRatio], positions=[i,i+1], vert=False, widths=0.6)

        # plt.plot(enrRatio, jitter(enrRatio,1,0.02), '.r', ms=5, label='enr')
        # plt.plot(natRatio, jitter(natRatio,2,0.02), '.g', ms=5, label='nat')

    # ax.set_xticklabels(['A', 'B', 'C'])
    # ax.set_xticks([1.5, 4.5, 7.5])
    plt.xscale('log')
    plt.legend()
    plt.savefig('../plots/cs3-rates.pdf')


def outliers_iqr(ys):
    """ http://colingorrie.github.io/outlier-detection.html """
    quartile_1, quartile_3 = np.percentile(ys, [25, 75])
    iqr = quartile_3 - quartile_1
    lower_bound = quartile_1 - (iqr * 1.5)
    upper_bound = quartile_3 + (iqr * 1.5)
    return np.where((ys > upper_bound) | (ys < lower_bound))


def getOutliers():

    e10, e250, n10, n250 = [], [], [], []

    for i, ds in enumerate([0,1,2]):
        f = np.load('../data/cs3-rates-ds%s.npz' % ds)
        dsRate = f['arr_0'].item()
        for cpd in sorted(dsRate):
            e10.extend([r[0] for r in dsRate[cpd] if det.allDetIDs[cpd] > 100000 and r[1]!=0])
            n10.extend([r[0] for r in dsRate[cpd] if det.allDetIDs[cpd] < 100000 and r[1]!=0])
            e250.extend([r[1] for r in dsRate[cpd] if det.allDetIDs[cpd] > 100000 and r[1]!=0])
            n250.extend([r[1] for r in dsRate[cpd] if det.allDetIDs[cpd] < 100000 and r[1]!=0])

    e10, n10, e250, n250 = np.asarray(e10), np.asarray(n10), np.asarray(e250), np.asarray(n250)

    # there are a lot of ranges with rate10 == 0.  that's concerning.
    idx = np.where(e10>0)
    print("len e10 %d  num zeros %d" % (len(e10), len(e10[idx])))
    print("  outliers:", outliers_iqr(e10[idx]))
    # print("  outliers:",wl.niceList(outliers_iqr(e10[idx])))



if __name__=="__main__":
    main(sys.argv[1:])