#!/usr/bin/env python3
import sys, os, math
import numpy as np
import matplotlib.pyplot as plt
plt.style.use('../pltReports.mplstyle')

import waveLibs as wl
import dsi
bkg = dsi.BkgInfo()
det = dsi.DetInfo()

def main(argv):

    # getRates()
    printRates(0)
    getOutliers()

def getRates():

    from ROOT import TFile, TTree, gROOT
    gROOT.ProcessLine("gErrorIgnoreLevel = 3001;") # suppress ROOT error messages

    opt = ""
    # opt = "verbose"  # print results for every channel, every bIdx
    # opt = "zombie" # hits from a detector declared dead by livetime calc
    # opt = "dead"   # detectors w/ no hits at all

    for ds in [0,1,2,3,4,"5A","5B","5C"]:
    # for ds in ["5B"]:

        dsNum = int(ds[0]) if isinstance(ds,str) else ds
        nBkg = bkg.dsMap()[dsNum]
        bLo, bHi = 0, nBkg
        if ds=="5A": bLo, bHi = 0, 79
        if ds=="5B": bLo, bHi = 80, 112
        if ds=="5C": bLo, bHi = 113, 121
        runRanges = bkg.getRanges(dsNum)
        tl = TFile("../data/ds_%s_livetime.root" % str(ds))
        lt = tl.Get("dsTree")
        chList = det.getGoodChanList(dsNum)
        dsRate = {det.getChanCPD(dsNum,ch):[] for ch in chList} # output dict

        for bIdx in range(bLo, bHi+1):
            rLo, rHi = runRanges[bIdx][0], runRanges[bIdx][-1]
            # print("DS-%s  bIdx %d  runLo %d  runHi %d" % (ds, bIdx, rLo, rHi))

            # get the livetime
            livetime = {ch:0 for ch in chList}
            n = lt.Draw("run:channel:livetime","run>=%d && run<=%d" % (rLo, rHi), 'goff')
            ltRun, ltChan, ltLive = lt.GetV1(), lt.GetV2(), lt.GetV3()
            for i in range(n):
                livetime[ltChan[i]] += ltLive[i]

            # now access the cut files and look at the rates
            for ch in sorted(chList):

                cpd = det.getChanCPD(dsNum,ch)
                detID = det.getDetIDChan(dsNum,ch)
                aMass = det.allActiveMasses[detID]

                # skip nonexistent files.  other parts of chsel should tell us why these aren't here.
                fName = "%s/bkg/cut/fr/fr_ds%d_%d_ch%d.root" % (dsi.dataDir, dsNum, bIdx, ch)
                if not os.path.isfile(fName): continue

                tf = TFile(fName)
                tt = tf.Get("skimTree")
                nEvt = tt.GetEntries()
                live = livetime[ch]/86400
                expo = livetime[ch]*aMass/86400/1000
                if live == 0:
                    if opt == "zombie":
                        print("'Zombie' detector: DS-%s  bIdx %d  cpd %s  ch %d exposure 0, hits %d" % (ds, bIdx, cpd, ch, nEvt))
                    continue

                n10 = tt.Draw("Entry$:Iteration$","trapENFCal >= 0 && trapENFCal <= 10","goff")
                n250 = tt.Draw("Entry$:Iteration$","trapENFCal >= 10 && trapENFCal <= 250","goff")
                r10 = n10/expo/10
                r250 = n250/expo/240

                if opt == "verbose":
                    print("DS %s  bIdx %d  ch %d  cpd %s  exp %.1f - %d  %d  %.1f  %.1f" % (ds, bIdx, ch, cpd, exposure[ch], n10,n250,r10,r250))

                dsRate[cpd].append([r10,r250,expo,bIdx])

        np.savez('../data/cs4-rates-ds%s.npz' % ds, dsRate, ds)
        printRates(ds)


def printRates(ds):

    f = np.load('../data/cs4-rates-ds%s.npz' % ds)
    dsRate = f['arr_0'].item()
    dsNum = int(ds[0]) if isinstance(ds,str) else ds

    # average rate in each detector
    avgRateDS = {}
    for cpd in sorted(dsRate):

        r10  = [r[0] for r in dsRate[cpd]]
        r250 = [r[1] for r in dsRate[cpd]]
        expo = [r[2] for r in dsRate[cpd]]
        bIdx = [r[3] for r in dsRate[cpd]]

        dType = "enr" if det.allDetIDs[cpd] > 100000 else "nat"

        if sum(expo)!=0:
            ar10 = np.average(r10, weights=expo)
            sr10 = math.sqrt( np.average((r10 - ar10)**2, weights=expo) )

            ar250 = np.average(r250, weights=expo)
            sr250 = math.sqrt( np.average((r250 - ar250)**2, weights=expo) )

            # print("%s  %s  r10 %.3f ± %.3f  r250 %.3f ± %.3f" % (cpd, dType, ar10, sr10, ar250, sr250))
        else:
            ar10, sr10, ar250, sr250 = 0, 0, 0, 0
            # print("%s  %s  exp = 0" % (cpd, dType))

        avgRateDS[cpd] = [ar10,sr10,ar250,sr250]

    # check output
    # for cpd in avgRateDS:
        # print(cpd, avgRateDS[cpd])

    # average enriched rate, average natural rate
    enr10, enr250, enrExp, nat10, nat250, natExp = [], [], [], [], [], []

    for cpd in sorted(dsRate):
        enr10.extend( [r[0] for r in dsRate[cpd] if det.allDetIDs[cpd] > 100000] )
        enr250.extend( [r[1] for r in dsRate[cpd] if det.allDetIDs[cpd] > 100000] )
        enrExp.extend( [r[2] for r in dsRate[cpd] if det.allDetIDs[cpd] > 100000] )

        nat10.extend( [r[0] for r in dsRate[cpd] if det.allDetIDs[cpd] < 100000] )
        nat250.extend( [r[1] for r in dsRate[cpd] if det.allDetIDs[cpd] < 100000] )
        natExp.extend( [r[2] for r in dsRate[cpd] if det.allDetIDs[cpd] < 100000] )

    ear10 = np.average(enr10, weights=enrExp)
    esr10 = math.sqrt( np.average((enr10 - ear10)**2, weights=enrExp) )
    ear250 = np.average(enr250, weights=enrExp)
    esr250 = math.sqrt( np.average((enr250 - ear250)**2, weights=enrExp) )

    nar10 = np.average(nat10, weights=natExp)
    nsr10 = math.sqrt( np.average((nat10 - nar10)**2, weights=natExp) )
    nar250 = np.average(nat250, weights=natExp)
    nsr250 = math.sqrt( np.average((nat250 - nar250)**2, weights=natExp) )

    print("DS-%s Enriched: r10 %.3f ± %.3f  r250 %.3f ± %.3f   Natural: r10 %.3f ± %.3f  r250 %.3f ± %.3f" % (ds,ear10,esr10,ear250,esr250,nar10,nsr10,nar250,nsr250))


def outliers_iqr(ys,k=1.5):
    """ http://colingorrie.github.io/outlier-detection.html """
    quartile_1, quartile_3 = np.percentile(ys, [25, 75])
    iqr = quartile_3 - quartile_1
    lower_bound = quartile_1 - (iqr * k)
    upper_bound = quartile_3 + (iqr * k)
    return np.where((ys > upper_bound) | (ys < lower_bound))


def getOutliers():
    """ https://www.itl.nist.gov/div898/handbook/prc/section1/prc16.htm
    1. Find outliers by cpd
    2. Find outliers within a given cpd
    3. Find runs causing the outliers
    """
    ds, dType = 0, "enr"

    f = np.load('../data/cs4-rates-ds%s.npz' % ds)
    dsRate = f['arr_0'].item()

    # exposure-weighted rate for each detector
    avgRate = {}
    for cpd in sorted(dsRate):

        isEnr = True if det.allDetIDs[cpd] > 100000 else False
        if dType == "enr" and not isEnr: continue
        elif dType == "nat" and isEnr: continue

        r10 = np.asarray([r[0] for r in dsRate[cpd]])
        r250 = np.asarray([r[1] for r in dsRate[cpd]])
        expo = np.asarray([r[2] for r in dsRate[cpd]])
        bIdx = np.asarray([r[3] for r in dsRate[cpd]])

        a10 = np.average(r10, weights=expo)
        s10 = math.sqrt( np.average((r10 - a10)**2, weights=expo) )

        idxOut10 = outliers_iqr(r10)
        nOut = len(idxOut10[0])

        print("DS-%s  CPD %s  nTot %d  nOut %d  r10: %.2f ± %.2f" % (ds, cpd, len(r10), nOut, a10, s10))
        for i in idxOut10[0]:
            print("    bIdx %-3d  r10 %.2f" % (bIdx[i], r10[i]))

        avgRate[cpd] = [a10,s10]



if __name__=="__main__":
    main(sys.argv[1:])