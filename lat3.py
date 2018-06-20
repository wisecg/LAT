#!/usr/bin/env python3
"""
===================== lat3.py ======================
Apply a final run/channel selection "burst cut" to
data passing the PSA cuts from LAT2.
Use the "Tukey fence" method to detect anomalous rates
in cpd/bIdx's, and flag them for removal.
Generates "cut files" from LAT2 data.
=================== C. Wiseman (USC) ===================
"""
import sys, os, math, glob
import numpy as np
import matplotlib.pyplot as plt
plt.style.use('./pltReports.mplstyle')

import waveLibs as wl
import dsi
bkg = dsi.BkgInfo()
det = dsi.DetInfo()

rateWin1 = [0, 5] # < -- use this one
rateWin2 = [5, 20]
# kList = [5, 1.5]
kList = [5, 2] # < -- use this one

# a very important parameter, use 90 or 95
global pctTot
# pctTot = 90
pctTot = 95 # < -- use this one
print("Using pctTot ==",pctTot)

def main(argv):
    """
    1. Calculate typical enr/nat rates in a DS
    2. Find outliers in specific cpd/bIdx combinations
    3. Recalculate typical rates after rejecting outliers
    4? Find runs causing the outliers
    """

    # these can all be run sequentially
    # getRates()
    # getOutliers(True,usePass2=False)
    # plotRates()
    makeCutFiles()
    # plotSpecBeforeAfter()
    # plotSpectraAfter()
    # combineSpectra()


def getRates():

    from ROOT import TFile, TTree, gROOT
    gROOT.ProcessLine("gErrorIgnoreLevel = 3001;") # suppress ROOT error messages

    cutType = "fr"

    tOffCut = "tOffset < 100"

    opt = ""
    # opt = "verbose"  # print results for every channel, every bIdx
    # opt = "zombie" # hits from a detector declared dead by livetime calc
    # opt = "dead"   # detectors w/ no hits at all

    for ds in [0,1,2,3,4,"5A","5B","5C"]:
    # for ds in [0]:
        print("Getting DS-%s rates ..." % str(ds))

        dsNum = int(ds[0]) if isinstance(ds,str) else ds
        nBkg = bkg.dsMap()[dsNum]
        bLo, bHi = 0, nBkg
        if ds=="5A": bLo, bHi = 0, 79
        if ds=="5B": bLo, bHi = 80, 112
        if ds=="5C": bLo, bHi = 113, 121
        runRanges = bkg.getRanges(dsNum)
        tl = TFile("./data/ds_%s_livetime.root" % str(ds))
        lt = tl.Get("dsTree")
        chList = det.getGoodChanList(dsNum)

        rateData = {det.getChanCPD(dsNum,ch):[] for ch in chList} # output dict

        dsTmp = ds # this is a hack to keep the dataset in the numpy array
        if ds=="5A": dsTmp=50
        if ds=="5B": dsTmp=51
        if ds=="5C": dsTmp=52

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
                fName = "%s/bkg/cut/%s%d/%s_ds%d_%d_ch%d.root" % (dsi.dataDir, cutType, pctTot, cutType, dsNum, bIdx, ch)
                if not os.path.isfile(fName): continue

                tf = TFile(fName)
                tt = tf.Get("skimTree")
                nEvt = tt.GetEntries()
                live = livetime[ch]/86400
                expo = livetime[ch]*aMass/86400/1000
                if live == 0:
                    # if opt == "zombie":
                    # OK, we remove these hits/files when we run makeCutFiles.
                    # print("'Zombie' detector: DS-%s  bIdx %d  cpd %s  ch %d exposure 0, hits %d" % (ds, bIdx, cpd, ch, nEvt))
                    continue

                n1 = tt.Draw("run:channel","trapENFCal >= %.1f && trapENFCal <= %.1f && %s" % (rateWin1[0], rateWin1[1], tOffCut), "goff")
                hitC = tt.GetV2()
                hitC = list(set([hitC[i] for i in range(n1)]))
                if len(hitC) > 1:
                    print("ERROR: Looking for channel %d, found channels" % ch, hitC)
                    exit(1)

                n2 = tt.Draw("Entry$:Iteration$","trapENFCal >= %.1f && trapENFCal <= %.1f && %s" % (rateWin2[0], rateWin2[1], tOffCut), "goff")
                r1 = n1/expo/(rateWin1[1]-rateWin1[0]) # cts/kg-d-kev
                r2 = n2/expo/(rateWin2[1]-rateWin2[0])

                if opt == "verbose":
                    if ch==690:
                        # print("DS %s  bIdx %d  ch %d  cpd %s  exp %-4.1f - %d  %d  %.1f  %.1f" % (ds, bIdx, ch, cpd, expo, n1, n2, r1, r2))
                        print("bIdx %d  %d/%s  exp %-4.1f - %-4d  %.1f" % (bIdx, ch, cpd, expo, n1, r1))

                rateData[cpd].append([r1,r2,expo,bIdx,int(cpd),dsTmp,n1])

        np.savez('./data/lat3-rates-ds%s-e%d.npz' % (ds, pctTot), rateData, ds)


def getOutliers(verbose=False, usePass2=False, noSkip=False):
    """ Apply the "closeFence" method to enriched and natural detectors in each data set separately.
    Return a list of excluded [ds,cpd,bIdx]'s.
    https://math.stackexchange.com/questions/966331/why-john-tukey-set-1-5-iqr-to-detect-outliers-instead-of-1-or-2
    """
    if verbose:
        print("Rate window1: %.1f - %.1f keV" % (rateWin1[0], rateWin1[1]))
        # print("rateWin2: %.1f - %.1f keV" % (rateWin2[0], rateWin2[1])) # not using this, only window 1
        print("Fence values:", kList)
        print("Using Pass 2:", usePass2)

    enrExc, natExc = [], []  # return [ds,cpd,bIdx] to be excluded

    # 'rates' objects: [:,0]=rate1, [:,1]=rate2, [:,2]=expo, [:,3]=bkgIdx, [:,4]=cpd, [:,5]=ds, [:,6]=n1
    enrRates, natRates = [], []

    for ds in [0,1,2,3,4,"5A","5B","5C"]:
    # for ds in ["5A"]:

        dsNum = int(ds[0]) if isinstance(ds,str) else ds # this is used for getGoodChanList

        f = np.load('./data/lat3-rates-ds%s-e%d.npz' % (ds, pctTot))
        rateData = f['arr_0'].item()

        dsTmp = ds # this is a hack to keep the dataset in the numpy array
        if ds=="5A": dsTmp=50
        if ds=="5B": dsTmp=51
        if ds=="5C": dsTmp=52

        # split rate data into enr/nat groups
        enr, nat = [], []
        for i, cpd in enumerate(sorted(rateData)):
            if len(rateData[cpd])==0:
                continue
            isEnr = True if det.allDetIDs[cpd] > 100000 else False
            for v in np.asarray(rateData[cpd]):
                if isEnr: enr.append(v)
                else: nat.append(v)
        enr, nat = np.asarray(enr), np.asarray(nat) # rates objects

        # this is where the magic happens

        # 1. find cpd/bIdx outliers for all enr and nat detectors in the DS
        tmpEnrExc, tmpEnr = closeFence(ds, enr, kList, "Enr", noSkip=noSkip)
        tmpNatExc, tmpNat = closeFence(ds, nat, kList, "Nat", noSkip=noSkip)

        # 2. find bIdx outliers for each detector separately (try to remove thresholds noise)
        #    ugg, this didn't really work for DS0 (still see peaking)
        #    let's not use it, it just unnecessarily removes exposure.
        # if usePass2:
        #     k = 1.5
        #
        #     chList = det.getGoodChanList(dsNum)
        #     for ch in chList:
        #         cpd = int(det.getChanCPD(dsNum,ch))
        #         isEnr = True if det.getDetIDChan(dsNum,ch) > 100000 else False
        #         if isEnr:
        #             idx = np.where(tmpEnr[:,4]==cpd)
        #             if len(idx[0])==0 or len(tmpEnr[idx]) < 2: continue
        #             detEnrExc,_ = closeFence(ds, tmpEnr[idx], [k], cpd, 0, iZ=True)
        #
        #             # add to the exclude list and delete from the rate object
        #             if len(detEnrExc)>0:
        #
        #                 if len(tmpEnrExc)==0: tmpEnrExc = detEnrExc
        #                 tmpEnrExc = np.append(tmpEnrExc,detEnrExc,axis=0)
        #
        #                 for d in detEnrExc:
        #                     iR = np.where((tmpEnr[:,3]==d[2]) & (tmpEnr[:,5]==dsTmp) & (tmpEnr[:,4]==cpd))
        #                     if len(iR[0])!=1:
        #                         print("ERROR, found more than one index")
        #                         exit(1)
        #                     tmpEnr = np.delete(tmpEnr, iR[0], 0)
        #         else:
        #             idx = np.where(tmpNat[:,4]==cpd)
        #             if len(idx[0])==0 or len(tmpNat[idx]) < 2: continue
        #             detNatExc,_ = closeFence(ds, tmpNat[idx], [k], cpd, 0, iZ=True)
        #             if len(detNatExc)>0:
        #
        #                 if len(tmpNatExc)==0: tmpNatExc = detNatExc
        #                 tmpNatExc = np.append(tmpNatExc,detNatExc,axis=0)
        #
        #                 for d in detNatExc:
        #                     iR = np.where((tmpNat[:,3]==d[2]) & (tmpNat[:,5]==dsTmp) & (tmpNat[:,4]==cpd))
        #                     if len(iR[0])!=1:
        #                         print("ERROR, found more than one index")
        #                         exit(1)
        #                     tmpNat = np.delete(tmpNat, iR[0], 0)

        ae1, se1 = getMuStd(tmpEnr[:,0])
        ae2, se2 = getMuStd(tmpEnr[:,1])
        an1, sn1 = getMuStd(tmpNat[:,0])
        an2, sn2 = getMuStd(tmpNat[:,1])

        # print("DS-%-3s  Enr, kept %d/%d  rLo %.3f ± %.3f  rHi %.3f ± %.3f" % (ds,len(tmpEnr),len(enr),ae1,se1,ae2,se2))
        # print("DS-%-3s  Nat, kept %d/%d  rLo %.3f ± %.3f  rHi %.3f ± %.3f" % (ds,len(tmpNat),len(nat),an1,sn1,an2,sn2))
        if verbose:
            print("DS-%-2s  Enr %-4d/%-4d  rLo %.3f ± %.3f  rHi %.3f ± %.3f    Nat %-4d/%-4d  rLo %.3f ± %.3f  rHi %.3f ± %.3f" % (ds,len(tmpEnr),len(enr),ae1,se1,ae2,se2,  len(tmpNat),len(nat),an1,sn1,an2,sn2))

        enrExc.extend(tmpEnrExc)
        natExc.extend(tmpNatExc)

        enrRates.append([v for v in tmpEnr])
        natRates.append([v for v in tmpNat])

    enrRates = np.vstack(enrRates)
    natRates = np.vstack(natRates)

    enrExc = np.asarray(enrExc)
    natExc = np.asarray(natExc)

    return enrExc, natExc, enrRates, natRates


def getMuStd(vals, wts=None):
    """ Do a weighted average, and a weighted standard deviation. """
    mu = np.average(vals, weights=wts)
    std = math.sqrt( np.average((vals-mu)**2, weights=wts) )
    return mu, std


def closeFence(ds, rates, kList, name="", verbose=0, iZ=True, noSkip=False):
    """ Recursively applies the IQR / Tukey fence method to reject upper outliers, for the values in kList.
    Returns a list [[ds,cpd1,bIdx1],[ds,cpd2,bIdx2],...] and the updated averages.
    """
    excList = []

    dsNum = ds # this is a hack to keep the dataset in the numpy array
    if ds=="5A": dsNum=50
    if ds=="5B": dsNum=51
    if ds=="5C": dsNum=52

    for k in kList:

        # tag outliers from just the lower rate window
        iE1 = outliersIQR(rates[:,0], k, "hi", ignoreZeros=True)[0]
        # iE2 = outliersIQR(rates[:,1], k, "hi", ignoreZeros=True)[0]
        # iE = sorted(list(set(np.append(iE1, iE2))))
        iE = iE1

        if noSkip:
            print("Warning, no-skip mode active.")
            iE = np.asarray([])

        # get the average before excluding
        avg1, std1 = getMuStd(rates[:,0], None)
        avg2, std2 = getMuStd(rates[:,1], None)
        if verbose != 0:
            print("DS-%s %s rates:" % (ds, name))
            print("    rLo %.3f ± %-10.3f  %d/%d outliers" % (avg1, std1, len(iE1), len(rates)))
            # print("    rHi %.3f ± %-10.3f  %d/%d outliers" % (avg2, std2, len(iE2), len(rates)))
            print("    Applying outlier cut, k=%.1f, %d/%d total outliers" % (k, len(iE), len(rates)))

        for i in iE:
            excList.append( np.asarray([dsNum, int(rates[:,4][i]), int(rates[:,3][i]) ])) # dsNum, cpd, bkgIdx

            if verbose == 1:
                msg = ""
                if i in iE1: msg += "lo "
                # if i in iE2: msg += "hi "
                print("    iE %-4d  %-5d  bIdx %-3d  rLo %-8.2f  rHi %-8.2f  nLo %-3d  %s" % (i, rates[:,4][i], rates[:,3][i], rates[:,0][i], rates[:,1][i], rates[:,6][i], msg))

        # super verbose output
        if verbose == 2:
            print("k value, ",k)
            for i in range(len(rates[:,0])):
                msg = ""
                if i in iE1: msg += "lo "
                # if i in iE2: msg += "hi "
                print("i %-4d  %-5d  bIdx %-3d  rLo %-8.2f  rHi %-8.2f  nLo %-3d  %s" % (i, rates[:,4][i], rates[:,3][i], rates[:,0][i], rates[:,1][i], rates[:,6][i], msg))

        # exclude iE from rates, like the opposite of np.where
        rates = rates[~np.in1d(range(len(rates)),iE)]

    excList = np.asarray(excList)

    return excList, rates


def outliersIQR(vals, k=1.5, opt="", ignoreZeros=False):
    """ k is the Tukey fence.  1.5:"far out", 3:"way out", 5:"our detector sucks"
    https://www.itl.nist.gov/div898/handbook/prc/section1/prc16.htm
    Returns numpy indexes of outlying points in array 'vals'.
    Thx Colin http://colingorrie.github.io/outlier-detection.html
    """
    if ignoreZeros:
        qt1, qt3 = np.percentile(vals[np.where(vals>0)], [25,75])
    else:
        qt1, qt3 = np.percentile(vals,[25,75])

    iqr = qt3 - qt1
    lower_fence = qt1 - (iqr * k)
    upper_fence = qt3 + (iqr * k)

    if opt == "hi":
        return np.where(vals > upper_fence)  # hi outliers only
    else:
        return np.where((vals > upper_fence) | (vals < lower_fence))


def plotRates():
    """ Make box-and-whisker plots of the rates after the chan sel cut.
    rate objects: [:,0]=rate1, [:,1]=rate2, [:,2]=expo, [:,3]=bkgIdx, [:,4]=cpd, [:,5]=ds
    dumb DS5 trick: 5A==50, 5B==51, 5C==52
    """

    # plotName = "./plots/lat3-rates-before-burst-withzeros-e%d.pdf" % (pctTot)
    # enrExc, natExc, enrRates, natRates = getOutliers(False, noSkip=True)

    plotName = "./plots/lat3-rates-after-burst-withzeros-e%d.pdf" % (pctTot)
    enrExc, natExc, enrRates, natRates = getOutliers(False)

    def jitter(vals,idx,sc):
        return [idx + np.random.randn()*sc for v in vals]

    def set_box_color(bp, color):
        plt.setp(bp['boxes'], color=color)
        plt.setp(bp['whiskers'], color=color)
        plt.setp(bp['caps'], color=color)
        plt.setp(bp['medians'], color=color)

    fig = plt.figure()

    dsList = [0,1,2,3,4,"5A","5B","5C"]

    dataBox, xVals = [], []
    for i, ds in enumerate(dsList):
        dsNum = ds
        if ds=="5A": dsNum=50
        if ds=="5B": dsNum=51
        if ds=="5C": dsNum=52
        # idx = np.where((enrRates[:,5]==dsNum) & (enrRates[:,0] > 0))
        idx = np.where((enrRates[:,5]==dsNum))
        yvals = enrRates[:,0][idx]
        xvals = jitter(yvals,i,0.05)
        dataBox.append(yvals)
        xVals.append(xvals)

    plt.boxplot(dataBox, positions=np.array(range(len(dataBox))), sym='', widths=0.6)
    # plt.boxplot(natBox, positions=np.array(range(len(natBox)))*2.0-0.4, sym='', widths=0.6)

    for i in range(len(dataBox)):
        # plt.semilogy(xVals[i], dataBox[i], ".b", ms=5) # before
        plt.plot(xVals[i], dataBox[i], ".b", ms=3) # after

    plt.gca().set_xticklabels(dsList)

    plt.xlabel("Data Set", ha='right', x=1)
    plt.ylabel("Counts/kg-d-keV", ha='right', y=1)

    plt.tight_layout()
    # plt.show()
    plt.savefig(plotName)

    # Brian says stop messing w/ the box plot and fit a histogram of rates in each DS to a Poisson distribution
    # he also says try a violin plot https://seaborn.pydata.org/generated/seaborn.violinplot.html


def getBadRuns():
    """ For each bad [ds,cpd,bIdx], figure out if it's localized to just a few runs.
        Maybe apply the analysis threshold ??
        Ehhh, if we can live with the exposure reduction from the full cpd/bIdx cut, it's not worth it.
    """
    enrExc, natExc, enrRates, natRates = getOutliers(False)

    print("enr")
    for ds,cpd,bIdx in enrExc:
        print(ds,cpd,bIdx)

    # print("nat")
    # for ds,cpd,bIdx in natExc:
    #     print(ds,cpd,bIdx)


def makeCutFiles():

    from ROOT import TFile, TTree, MGTWaveform

    enrExc, natExc, enrRates, natRates = getOutliers(True)
    # enrExc, natExc: [:,0]=dsNum, [:,1]=cpd, [:,2]=bkgIdx
    # enrRates, natRates: [:,0]=rate1, [:,1]=rate2, [:,2]=expo, [:,3]=bkgIdx, [:,4]=cpd, [:,5]=ds

    cutType = "fr"

    # additional DC cuts can go here
    tOffCut = "tOffset < 100"

    # which burst cut do we want?
    pass2 = False
    outType = "frb2%d" % pctTot if pass2 else "frb%d" % pctTot

    for ds in [0,1,2,3,4,"5A","5B","5C"]:

        dsNum = int(ds[0]) if isinstance(ds,str) else ds
        nBkg = bkg.dsMap()[dsNum]
        bLo, bHi = 0, nBkg
        if ds=="5A": bLo, bHi = 0, 79
        if ds=="5B": bLo, bHi = 80, 112
        if ds=="5C": bLo, bHi = 113, 121
        runRanges = bkg.getRanges(dsNum)
        chList = det.getGoodChanList(dsNum)

        # clear out any files from a previous attempt
        fList = ["%s/bkg/cut/%s/%s_ds%d_%d_*.root" % (dsi.dataDir, outType, outType, dsNum, bIdx) for bIdx in range(bLo, bHi+1)]
        for f in fList:
            fTmp = glob.glob(f)
            for f in fTmp:
                os.remove(f)

        # build skip list
        dsTmp = ds
        if ds=="5A": dsTmp=50
        if ds=="5B": dsTmp=51
        if ds=="5C": dsTmp=52
        iE = np.where(enrExc[:,0]==dsTmp)
        iN = np.where(natExc[:,0]==dsTmp)
        skipList = np.vstack((enrExc[iE], natExc[iN]))
        print("DS-%s, skipList:" % ds)
        print(skipList)
        continue # debug to print the skip list

        # load ds_livetime output
        tl = TFile("./data/ds_%s_livetime.root" % str(ds))
        lt = tl.Get("dsTree")

        for bIdx in range(bLo, bHi+1):
            print("DS-%s  bIdx %d" % (ds, bIdx))

            # get channel livetimes (to identify 'zombie' channels declared dead that have hits)
            live = {ch:0 for ch in chList}
            n = lt.Draw("run:channel:livetime","run>=%d && run<=%d" % (runRanges[bIdx][0], runRanges[bIdx][-1]), 'goff')
            ltRun, ltChan, ltLive = lt.GetV1(), lt.GetV2(), lt.GetV3()
            for i in range(n):
                ch = ltChan[i]
                cpd = det.getChanCPD(dsNum,ch)
                detID = det.getDetIDChan(dsNum,ch)
                aMass = det.allActiveMasses[detID]
                live[ch] += ltLive[i]

            # load the cut files
            for ch in sorted(chList):

                cpd = det.getChanCPD(dsNum,ch)

                if len(np.where((skipList == (dsTmp, int(cpd), bIdx)).all(axis=1))[0]) > 0:
                    # print("skipping det %s in bkgIdx %d" % (cpd, bIdx))
                    continue

                # skip nonexistent files.  other parts of chsel should tell us why these aren't here.
                fName = "%s/bkg/cut/%s%d/%s_ds%d_%d_ch%d.root" % (dsi.dataDir, cutType, pctTot, cutType, dsNum, bIdx, ch)
                if not os.path.isfile(fName):
                    # print("no file for det %s in bkgIdx %d" % (cpd, bIdx))
                    continue

                tf = TFile(fName)
                tt = tf.Get("skimTree")
                nEvt = tt.GetEntries()

                # skip zombie detectors
                if live[ch]==0:
                    print("Zombie detector, cpd %s  ch %d, lt=0, nHits %d.  Excluding..." % (cpd, ch, nEvt))
                    continue

                outName = "%s/bkg/cut/%s/%s_ds%d_%d_ch%d.root" % (dsi.dataDir, outType, outType, dsNum, bIdx, ch)
                outFile = TFile(outName, "RECREATE")
                outTree = TTree()
                outTree = tt.CopyTree(tOffCut)
                # print("Wrote %d entries." % outTree.GetEntries())

                # if nEvt != outTree.GetEntries():
                    # print("ERROR, number of entries don't match: input %d  output %d" % (nEvt, outTree.GetEntries()))
                    # return

                outTree.Write()
                outFile.Close()
                tf.Close()


def plotSpecBeforeAfter():

    from ROOT import TChain

    dsList = [0,1,2,3,4,"5A","5B","5C"]
    # dsList = ["5C"]

    plotName = "./plots/lat3-comp-allDS-e%d.pdf" % pctTot

    tB = TChain("skimTree") # before
    tA = TChain("skimTree") # after

    for ds in dsList:

        # plotName = "./plots/lat3-comp-DS%s.pdf" % ds
        # tB = TChain("skimTree") # before
        # tA = TChain("skimTree") # after

        # for ds in dsList:
        dsNum = int(ds[0]) if isinstance(ds,str) else ds
        nBkg = bkg.dsMap()[dsNum]
        bLo, bHi = 0, nBkg
        if ds=="5A": bLo, bHi = 0, 79
        if ds=="5B": bLo, bHi = 80, 112
        if ds=="5C": bLo, bHi = 113, 121
        fB = ["%s/bkg/cut/fr%d/fr_ds%d_%d_*.root" % (dsi.dataDir, pctTot, dsNum, bIdx) for bIdx in range(bLo, bHi+1)]
        for f in fB: tB.Add(f)
        fA = ["%s/bkg/cut/frb%d/frb%d_ds%d_%d_*.root" % (dsi.dataDir, pctTot, pctTot, dsNum, bIdx) for bIdx in range(bLo, bHi+1)]
        for f in fA: tA.Add(f)

    print("before",tB.GetEntries(),"after",tA.GetEntries())

    xLo, xHi, xpb = 0, 20, 0.1

    fig = plt.figure(figsize=(8,7))
    p0 = plt.subplot(211)
    p1 = plt.subplot(212)

    tCut = "!isEnr"

    n = tB.Draw("trapENFCal",tCut,"goff")
    hitE = tB.GetV1()
    hitE = [hitE[i] for i in range(n)]
    x, hB = wl.GetHisto(hitE, xLo, xHi, xpb, shift=False)

    n = tA.Draw("trapENFCal",tCut,"goff")
    hitE = tA.GetV1()
    hitE = [hitE[i] for i in range(n)]
    x, hA = wl.GetHisto(hitE, xLo, xHi, xpb, shift=False)

    p0.semilogy(x, hB, lw=2, ls='steps', c='r', label='Natural, Before')
    p0.semilogy(x, hA, lw=2, ls='steps', c='b', label='Natural, After')
    p0.axvline(1., c='g', lw=2, alpha=0.5, label="1.0 keV")
    p0.set_ylabel("Counts", ha='right', y=1)
    p0.legend()

    tCut = "isEnr"

    n = tB.Draw("trapENFCal",tCut,"goff")
    hitE = tB.GetV1()
    hitE = [hitE[i] for i in range(n)]
    x, hB = wl.GetHisto(hitE, xLo, xHi, xpb, shift=False)

    n = tA.Draw("trapENFCal",tCut,"goff")
    hitE = tA.GetV1()
    hitE = [hitE[i] for i in range(n)]
    x, hA = wl.GetHisto(hitE, xLo, xHi, xpb, shift=False)

    p1.semilogy(x, hB, lw=2, ls='steps', c='r', label='Enriched, Before')
    p1.semilogy(x, hA, lw=2, ls='steps', c='b', label='Enriched, After')
    p1.axvline(1., c='g', lw=2, alpha=0.5, label="1.0 keV")
    p1.set_xlabel("Energy (keV)", ha='right', x=1)
    p1.legend()

    plt.tight_layout()
    # plt.show()

    plt.savefig(plotName)

    tB.Reset()
    tA.Reset()


def plotSpectraAfter():

    from ROOT import TChain, TFile, TTree
    from matplotlib import colors

    dsList = [0,1,2,3,4,"5A","5B","5C"]
    # dsList = ["5A","5B","5C"]

    # cType = "fr"             # before burst cut
    cType = "frb%d" % pctTot   # after burst cut <-- use this one
    # cType = "frb2"           # after pass 2 burst cut

    for ds in dsList:

        plotName1 = "./plots/lat3-%s-DS%s-spec.pdf" % (cType,ds)
        plotName2 = "./plots/lat3-%s-DS%s-vsch.pdf" % (cType,ds)
        plotName3 = "./plots/lat3-%s-DS%s-rvsc.pdf" % (cType,ds)

        tt = TChain("skimTree")

        # for ds in dsList: # indent this block
        dsNum = int(ds[0]) if isinstance(ds,str) else ds
        nBkg = bkg.dsMap()[dsNum]
        bLo, bHi = 0, nBkg
        if ds=="5A": bLo, bHi = 0, 79
        if ds=="5B": bLo, bHi = 80, 112
        if ds=="5C": bLo, bHi = 113, 121
        fList = ["%s/bkg/cut/%s/%s_ds%d_%d_*.root" % (dsi.dataDir, cType, cType, dsNum, bIdx) for bIdx in range(bLo, bHi+1)]
        for f in fList: tt.Add(f)

        fig = plt.figure(figsize=(8,7))
        p0 = plt.subplot(211)
        p1 = plt.subplot(212)

        eSpec = True
        eVsCh = True
        rVsCh = True

        # ====================== 1. energy spectra ======================
        if eSpec:

            xLo, xHi, xpb = 0, 20, 0.1
            # xLo, xHi, xpb = 0, 50, 0.1

            tCut = "!isEnr"

            n = tt.Draw("trapENFCal",tCut,"goff")
            hitE = tt.GetV1()
            hitE = [hitE[i] for i in range(n)]
            x, hNat = wl.GetHisto(hitE, xLo, xHi, xpb, shift=False)

            p0.plot(x, hNat, lw=2, ls='steps', c='b', label='Natural')
            p0.axvline(1., c='g', lw=2, alpha=0.5, label="1.0 keV")
            p0.set_ylabel("Counts", ha='right', y=1)
            p0.legend()

            tCut = "isEnr"

            n = tt.Draw("trapENFCal",tCut,"goff")
            hitE = tt.GetV1()
            hitE = [hitE[i] for i in range(n)]
            x, hEnr = wl.GetHisto(hitE, xLo, xHi, xpb, shift=False)

            n = tt.Draw("trapENFCal",tCut,"goff")
            hitE = tt.GetV1()
            hitE = [hitE[i] for i in range(n)]
            x, hA = wl.GetHisto(hitE, xLo, xHi, xpb, shift=False)

            p1.plot(x, hEnr, lw=2, ls='steps', c='b', label='Enriched')
            p1.axvline(1., c='g', lw=2, alpha=0.5, label="1.0 keV")
            p1.set_xlabel("Energy (keV)", ha='right', x=1)
            p1.legend()

            plt.tight_layout()
            # plt.show()
            plt.savefig(plotName1)


        # ====================== 2. energy vs channel counts (2d) ======================

        if eVsCh:

            dsTmp = int(ds[0]) if isinstance(ds,str) else ds # this is used for getGoodChanList
            chList = det.getGoodChanList(dsTmp)

            p0.cla()
            tCut = "!isEnr"

            natList = [c for c in chList if det.getDetIDChan(dsTmp,c) < 100000] # natural
            chLabels = ["%s/%s" % (det.getChanCPD(dsTmp,ch),ch) for ch in natList]
            chMap = {natList[i]:i for i in range(len(natList))}

            xLo, xHi, xpb = 0, 20., 0.2
            yLo, yHi = 0, len(natList)
            nbx, nby = int((xHi-xLo)/xpb), len(natList)

            n = tt.Draw("trapENFCal:channel",tCut,"goff")
            hitE, hitC = tt.GetV1(), tt.GetV2()
            hitE = [hitE[i] for i in range(n)]
            hitC = [chMap[hitC[i]] for i in range(n)]

            p0.hist2d(hitE, hitC, bins=[nbx, nby], range=[[xLo,xHi],[yLo,yHi]], cmap='jet')
            # p0.set_xlabel("Energy (keV)", ha='right', x=1.)
            p0.set_xticks(np.arange(xLo, xHi+1, 1.0))
            p0.set_ylabel("cpd/channel", ha='right', y=1.)
            p0.set_yticks(np.arange(0, len(natList))+0.5)
            p0.set_yticklabels(chLabels, fontsize=8)


            p1.cla()

            tCut = "isEnr"

            enrList = [c for c in chList if det.getDetIDChan(dsTmp,c) > 100000] # enriched
            chLabels = ["%s/%s" % (det.getChanCPD(dsTmp,ch),ch) for ch in enrList]
            chMap = {enrList[i]:i for i in range(len(enrList))}

            xLo, xHi, xpb = 0, 20., 0.2
            yLo, yHi = 0, len(enrList)
            nbx, nby = int((xHi-xLo)/xpb), len(enrList)

            n = tt.Draw("trapENFCal:channel",tCut,"goff")
            hitE, hitC = tt.GetV1(), tt.GetV2()
            hitE = [hitE[i] for i in range(n)]
            hitC = [chMap[hitC[i]] for i in range(n)]

            p1.hist2d(hitE, hitC, bins=[nbx, nby], range=[[xLo,xHi],[yLo,yHi]], cmap='jet')
            p1.set_xlabel("Energy (keV)", ha='right', x=1.)
            p1.set_xticks(np.arange(xLo, xHi+1, 1.0))
            # p1.set_ylabel("cpd/channel", ha='right', y=1.)
            p1.set_yticks(np.arange(0, len(enrList))+0.5)
            p1.set_yticklabels(chLabels, fontsize=8)

            plt.tight_layout()
            # plt.show()
            plt.savefig(plotName2)


        # ====================== 3. rate @ 0-5 kev, channel vs. bIdx (2d) ======================

        if rVsCh:

            dsTmp = int(ds[0]) if isinstance(ds,str) else ds # this is used for getGoodChanList
            nBkg = bkg.dsMap()[dsNum]
            bLo, bHi = 0, nBkg
            if ds=="5A": bLo, bHi = 0, 79
            if ds=="5B": bLo, bHi = 80, 112
            if ds=="5C": bLo, bHi = 113, 121
            chList = det.getGoodChanList(dsTmp)
            runRanges = bkg.getRanges(dsNum)
            tl = TFile("./data/ds_%s_livetime.root" % str(ds))
            lt = tl.Get("dsTree")

            # calculate weights for each chan/bkgIdx
            expo = {}
            for bIdx in range(bLo, bHi+1):
                livetime = {ch:0 for ch in chList}
                rLo, rHi = runRanges[bIdx][0], runRanges[bIdx][-1]
                m = lt.Draw("run:channel:livetime","run>=%d && run<=%d" % (rLo, rHi), 'goff')
                ltRun, ltChan, ltLive = lt.GetV1(), lt.GetV2(), lt.GetV3()
                for i in range(m): livetime[ltChan[i]] += ltLive[i]
                for ch in chList:
                    cpd = det.getChanCPD(dsNum,ch)
                    detID = det.getDetIDChan(dsNum,ch)
                    aMass = det.allActiveMasses[detID]
                    expo["%d/%d" % (ch,bIdx)] = livetime[ch]*aMass/86400/1000
                    # expo["%d/%d" % (ch,bIdx)] = livetime[ch]*aMass/86400/1000

            p0.cla()
            tCut = "!isEnr && trapENFCal >= %.1f && trapENFCal <= %.1f" % (rateWin1[0], rateWin1[1])

            chNat = [c for c in chList if det.getDetIDChan(dsTmp,c) < 100000] # natural
            chLabels = ["%s/%s" % (det.getChanCPD(dsTmp,ch),ch) for ch in chNat]
            chMap = {chNat[i]:i for i in range(len(chNat))}

            xLo, xHi, xpb = bLo, bHi+1, 1
            yLo, yHi = 0, len(chNat)
            nbx, nby = int((xHi-xLo)/xpb), len(chNat)

            n = tt.Draw("run:channel",tCut,"goff")
            hitR, hitC = tt.GetV1(), tt.GetV2()
            hitB = [bkg.GetBkgIdx(dsTmp,hitR[i]) for i in range(n)]
            hitCh = [hitC[i] for i in range(n)]       # raw channels, used in weighting function
            hitC = [chMap[hitC[i]] for i in range(n)] # mapped channels, used in plotting

            hitW = [ 1 / expo["%d/%d" % (hitCh[i],hitB[i])] / (rateWin1[1]-rateWin1[0]) for i in range(n)]

            # h0,_,_,im0 = p0.hist2d(hitB, hitC, bins=[nbx, nby], range=[[xLo,xHi],[yLo,yHi]], cmap='jet')
            h0,_,_,im0 = p0.hist2d(hitB, hitC, weights=hitW, bins=[nbx, nby], range=[[xLo,xHi],[yLo,yHi]], cmap='jet')

            # p0.set_xlabel("bkgIdx, DS-%d" % ds, ha='right', x=1.)
            p0.set_xticks(np.arange(bLo, bHi, 5)+0.5)
            p0.set_xticklabels(np.arange(bLo,bHi, 5))
            p0.set_ylabel("cpd/channel", ha='right', y=1.)
            p0.set_yticks(np.arange(0, len(chNat))+0.5)
            p0.set_yticklabels(chLabels, fontsize=8)


            p1.cla()
            tCut = "isEnr && trapENFCal >= %.1f && trapENFCal <= %.1f" % (rateWin1[0], rateWin1[1])

            chEnr = [c for c in chList if det.getDetIDChan(dsTmp,c) > 100000] # enriched
            chLabels = ["%s/%s" % (det.getChanCPD(dsTmp,ch),ch) for ch in chEnr]
            chMap = {chEnr[i]:i for i in range(len(chEnr))}

            xLo, xHi, xpb = bLo, bHi+1, 1
            yLo, yHi = 0, len(chEnr)
            nbx, nby = int((xHi-xLo)/xpb), len(chEnr)

            n = tt.Draw("run:channel",tCut,"goff")
            hitR, hitC = tt.GetV1(), tt.GetV2()
            hitB = [bkg.GetBkgIdx(dsTmp,hitR[i]) for i in range(n)]
            hitCh = [hitC[i] for i in range(n)]       # raw channels, used in weighting function
            hitC = [chMap[hitC[i]] for i in range(n)] # mapped channels, used in plotting

            hitW = [ 1 / expo["%d/%d" % (hitCh[i],hitB[i])] / (rateWin1[1]-rateWin1[0]) for i in range(n)]

            # h1,_,_,im1 = p1.hist2d(hitB, hitC, bins=[nbx, nby], range=[[xLo,xHi],[yLo,yHi]], cmap='jet')
            h1,_,_,im1 = p1.hist2d(hitB, hitC, weights=hitW, bins=[nbx, nby], range=[[xLo,xHi],[yLo,yHi]], cmap='jet')

            # DEBUG: print the rate for a channel
            # print(h1.shape)
            # idx, ch = 11, 690
            # for i,v in enumerate(h1[:,idx]):
                # bIdx = i
                # rate = v
                # exp = expo["%d/%d" % (ch,bIdx)]
                # print("%d  %.2f  %.2f" % (bIdx, rate, exp))
            # return

            p1.set_xlabel("bkgIdx, DS-%s" % ds, ha='right', x=1.)
            p1.set_xticks(np.arange(bLo, bHi, 5)+0.5)
            p1.set_xticklabels(np.arange(bLo,bHi, 5))
            # p1.set_ylabel("cpd/channel", ha='right', y=1.)
            p1.set_yticks(np.arange(0, len(chEnr))+0.5)
            p1.set_yticklabels(chLabels, fontsize=8)

            cb0 = fig.colorbar(im0, ax=p0)
            cb1 = fig.colorbar(im1, ax=p1)
            cb0.set_label('cts/kg/d/keV', ha='right', rotation=270, labelpad=20)

            plt.tight_layout()
            # plt.show()
            plt.savefig(plotName3)


        tt.Reset()


def combineSpectra():

    from ROOT import TChain, TFile, TTree

    # this makes a combined enr/nat spectrum for all these combinations and energy ranges
    specList = [
        [[0,1,2,3,4,"5A","5B","5C"],[0,20,0.1]],
        [[0,1,2,3,4,"5A","5B","5C"],[0,50,0.1]],
        [[1,2,3,4,"5A","5B","5C"],[0,20,0.1]],
        [[1,2,3,4,"5A","5B","5C"],[0,50,0.1]],
        [[1,2,3,4,"5B","5C"],[0,20,0.1]],
        [[1,2,3,4,"5B","5C"],[0,50,0.1]]
    ]

    for spec in specList:

        dsList = spec[0]
        xLo, xHi, xpb = spec[1]

        plotName = "./plots/lat3-ds-e%d-" % pctTot
        for ds in dsList: plotName += str(ds)
        plotName += "-%.1fkev.pdf" % xHi

        # cType = "fr"   # before burst cut
        cType = "frb%d" % pctTot  # after burst cut <---- use this one
        # cType = "frb2" # after pass 2 burst cut

        tt = TChain("skimTree")

        for ds in dsList:
            dsNum = int(ds[0]) if isinstance(ds,str) else ds
            nBkg = bkg.dsMap()[dsNum]
            bLo, bHi = 0, nBkg
            if ds=="5A": bLo, bHi = 0, 79
            if ds=="5B": bLo, bHi = 80, 112
            if ds=="5C": bLo, bHi = 113, 121
            fList = ["%s/bkg/cut/%s/%s_ds%d_%d_*.root" % (dsi.dataDir, cType, cType, dsNum, bIdx) for bIdx in range(bLo, bHi+1)]
            for f in fList: tt.Add(f)

        fig = plt.figure(figsize=(8,7))
        p0 = plt.subplot(211)
        p1 = plt.subplot(212)

        tCut = "!isEnr"

        n = tt.Draw("trapENFCal",tCut,"goff")
        hitE = tt.GetV1()
        hitE = [hitE[i] for i in range(n)]
        x, hNat = wl.GetHisto(hitE, xLo, xHi, xpb, shift=False)

        p0.plot(x, hNat, lw=2, ls='steps', c='b', label='Natural')
        p0.axvline(1., c='g', lw=2, alpha=0.5, label="1.0 keV")
        p0.axvline(2., c='m', lw=2, alpha=0.5, label="2.0 keV")
        p0.set_ylabel("Counts/%.1f keV" % xpb, ha='right', y=1)
        p0.legend()

        tCut = "isEnr"

        n = tt.Draw("trapENFCal",tCut,"goff")
        hitE = tt.GetV1()
        hitE = [hitE[i] for i in range(n)]
        x, hEnr = wl.GetHisto(hitE, xLo, xHi, xpb, shift=False)

        n = tt.Draw("trapENFCal",tCut,"goff")
        hitE = tt.GetV1()
        hitE = [hitE[i] for i in range(n)]
        x, hA = wl.GetHisto(hitE, xLo, xHi, xpb, shift=False)

        p1.plot(x, hEnr, lw=2, ls='steps', c='b', label='Enriched')
        p1.axvline(1., c='g', lw=2, alpha=0.5, label="1.0 keV")
        p1.axvline(2., c='m', lw=2, alpha=0.5, label="2.0 keV")
        p1.set_xlabel("Energy (keV)", ha='right', x=1)
        p1.legend()

        plt.tight_layout()
        # plt.show()
        plt.savefig(plotName)

        tt.Reset()


if __name__=="__main__":
    main(sys.argv[1:])
