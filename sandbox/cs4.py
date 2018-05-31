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
        rateData = {det.getChanCPD(dsNum,ch):[] for ch in chList} # output dict

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

                rateData[cpd].append([r10,r250,expo,bIdx,int(cpd)])

        np.savez('../data/cs4-rates-ds%s.npz' % ds, rateData, ds)


def getMuStd(vals, wts):
    mu = np.average(vals, weights=wts)
    std = math.sqrt( np.average((vals-mu)**2, weights=wts) )
    return mu, std


def outliers_iqr(ys, k=1.5):
    """ k is the Tukey fence.  1.5:"far out", 3:"way out"
    https://www.itl.nist.gov/div898/handbook/prc/section1/prc16.htm
    Returns indexes of outliers.
    Thx Colin http://colingorrie.github.io/outlier-detection.html
    """
    quartile_1, quartile_3 = np.percentile(ys, [25, 75])
    iqr = quartile_3 - quartile_1
    lower_bound = quartile_1 - (iqr * k)
    upper_bound = quartile_3 + (iqr * k)
    return np.where((ys > upper_bound) | (ys < lower_bound))


def closeFence(ds, rates, kList, verbose=False):
    """ Applies the IQR / Tukey fence method to reject outliers and update
    the 'rates' object, repeating for the values in kList.
    Returns a list [[cpd1,bIdx1],[cpd2,bIdx2],...] and the cleaned 'rates' object
    """

    excList = []

    for k in kList:

        # tag outliers from both 0-10 and 10-250
        iE1 = outliers_iqr(rates[:,0], k)[0]
        iE2 = outliers_iqr(rates[:,1], k)[0]
        iE = np.append(iE1, iE2)

        if verbose:
            ea10, es10 = getMuStd(rates[:,0], rates[:,2])
            print("DS-%s  Enr  r10 %.3f ± %.3f  %d/%d outliers, k=%.1f" % (ds, ea10, es10, len(iE), len(rates[:,0]), k))

        for i in iE:

            # don't add zeros to the exclude list
            if rates[:,0][i] == 0: continue

            excList.append( [int(rates[:,4][i]), int(rates[:,3][i]) ]) # cpd, bkgIdx

            if verbose:
                print("iE %-4d  %-5d  bIdx %-3d  r10 %-8.2f" % (i, rates[:,4][i], rates[:,3][i], rates[:,0][i]))

        # exclude iE from rates, like the opposite of np.where
        rates = rates[~np.in1d(range(len(rates)),iE)]

    return excList, rates


def getOutliers():
    """
    1. Calculate typical rates in a DS
    2. Find outliers by cpd
    3. Find outliers within a given cpd
    4. Find runs causing the outliers
    5. Recalculate typical rates after rejecting outliers

    https://math.stackexchange.com/questions/966331/why-john-tukey-set-1-5-iqr-to-detect-outliers-instead-of-1-or-2
    """

    for ds in [0,1,2,3,4,"5A","5B","5C"]:
    # for ds in [1]:

        f = np.load('../data/cs4-rates-ds%s.npz' % ds)
        rateData = f['arr_0'].item()
        dsNum = int(ds[0]) if isinstance(ds,str) else ds

        # split rate data into enr/nat groups
        enr, nat = [], []
        for i, cpd in enumerate(sorted(rateData)):
            if len(rateData[cpd])==0:
                continue
            isEnr = True if det.allDetIDs[cpd] > 100000 else False
            for v in np.asarray(rateData[cpd]):
                if isEnr: enr.append(v)
                else: nat.append(v)
        enr, nat = np.asarray(enr), np.asarray(nat) # [:,0]=rate10, [:,1]=rate250, [:,2]=expo, [:,3}=bkgIdx, [:,4]=cpd

        excList, exc = closeFence(ds, enr, [5,3], False)
        # for cpd,bIdx in excList:
            # print(cpd, bIdx)

        ea10, es10 = getMuStd(exc[:,0], exc[:,2])
        print("DS-%s  Enr  r10 %.3f ± %.3f" % (ds, ea10, es10))



if __name__=="__main__":
    main(sys.argv[1:])