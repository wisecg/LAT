#!/usr/bin/env python3
import sys, os, math, glob
import numpy as np
import matplotlib.pyplot as plt
plt.style.use('../pltReports.mplstyle')

import waveLibs as wl
import dsi
bkg = dsi.BkgInfo()
det = dsi.DetInfo()

rateWin1 = [0, 5]
rateWin2 = [5, 20]
# kList = [5, 1.5]
kList = [5, 3, 1.5]

def main(argv):
    """
    1. Calculate typical enr/nat rates in a DS
    2. Find outliers in specific cpd/bIdx combinations
    3. Recalculate typical rates after rejecting outliers
    4? Find runs causing the outliers
    """
    # global rateWin1, rateWin2, kList


    # getRates()
    # getOutliers(True)
    # plotRates()
    # makeCutFiles()
    # plotSpecBeforeAfter()
    plotSpectraAfter()


def getRates():

    from ROOT import TFile, TTree, gROOT
    gROOT.ProcessLine("gErrorIgnoreLevel = 3001;") # suppress ROOT error messages

    cutType = "fr"

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
                fName = "%s/bkg/cut/%s/%s_ds%d_%d_ch%d.root" % (dsi.dataDir, cutType, cutType, dsNum, bIdx, ch)
                if not os.path.isfile(fName): continue

                tf = TFile(fName)
                tt = tf.Get("skimTree")
                nEvt = tt.GetEntries()
                live = livetime[ch]/86400
                expo = livetime[ch]*aMass/86400/1000
                if live == 0:
                    # if opt == "zombie":
                    # TODO: need to make sure we remove these hits!
                    # print("'Zombie' detector: DS-%s  bIdx %d  cpd %s  ch %d exposure 0, hits %d" % (ds, bIdx, cpd, ch, nEvt))
                    continue

                n1 = tt.Draw("Entry$:Iteration$","trapENFCal >= %.1f && trapENFCal <= %.1f" % (rateWin1[0], rateWin1[1]), "goff")
                n2 = tt.Draw("Entry$:Iteration$","trapENFCal >= %.1f && trapENFCal <= %.1f" % (rateWin2[0], rateWin2[1]), "goff")
                r1 = n1/expo/(rateWin1[1]-rateWin1[0]) # cts/kg-d-kev
                r2 = n2/expo/(rateWin2[1]-rateWin2[0])

                if opt == "verbose":
                    print("DS %s  bIdx %d  ch %d  cpd %s  exp %.1f - %d  %d  %.1f  %.1f" % (ds, bIdx, ch, cpd, exposure[ch], n1,n2,r1,r2))

                rateData[cpd].append([r1,r2,expo,bIdx,int(cpd)])

        np.savez('../data/cs4-rates-ds%s.npz' % ds, rateData, ds)


def getMuStd(vals, wts=None):
    """ Do a weighted average, and a weighted standard deviation. """
    mu = np.average(vals, weights=wts)
    std = math.sqrt( np.average((vals-mu)**2, weights=wts) )
    return mu, std


def getOutliers(verbose=False):
    """ Apply the "closeFence" method to enriched and natural detectors in each data set separately.
    Return a list of excluded [ds,cpd,bIdx]'s.
    https://math.stackexchange.com/questions/966331/why-john-tukey-set-1-5-iqr-to-detect-outliers-instead-of-1-or-2
    """
    if verbose:
        print("rateWin1: %.1f - %.1f keV" % (rateWin1[0], rateWin1[1]))
        print("rateWin2: %.1f - %.1f keV" % (rateWin2[0], rateWin2[1]))
        print("fence values:", kList)

    enrExc, natExc = [], []  # return [ds,cpd,bIdx] to be excluded
    enrRates, natRates = [], [] # return 'rates' objects: [:,0]=rate1, [:,1]=rate2, [:,2]=expo, [:,3]=bkgIdx, [:,4]=cpd, [:,5]=ds

    for ds in [0,1,2,3,4,"5A","5B","5C"]:
    # for ds in ["5A"]:

        f = np.load('../data/cs4-rates-ds%s.npz' % ds)
        rateData = f['arr_0'].item()

        dsNum = ds # this is a hack to keep the dataset in the numpy array
        if ds=="5A": dsNum=50
        if ds=="5B": dsNum=51
        if ds=="5C": dsNum=52

        # split rate data into enr/nat groups
        enr, nat = [], []
        for i, cpd in enumerate(sorted(rateData)):
            if len(rateData[cpd])==0:
                continue
            isEnr = True if det.allDetIDs[cpd] > 100000 else False
            for v in np.asarray(rateData[cpd]):
                t = np.append(v,dsNum)
                if isEnr: enr.append(t)
                else: nat.append(t)
        enr, nat = np.asarray(enr), np.asarray(nat) # [:,0]=rate1, [:,1]=rate2, [:,2]=expo, [:,3]=bkgIdx, [:,4]=cpd, [:,5]=ds

        # this is where the magic happens
        tmpEnrExc, tmpEnr = closeFence(ds, enr, kList, "Enr", False)
        tmpNatExc, tmpNat = closeFence(ds, nat, kList, "Nat", False)

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


def closeFence(ds, rates, kList, name="", verbose=False):
    """ Recursively applies the IQR / Tukey fence method to reject outliers, for the values in kList.
    Returns a list [[ds,cpd1,bIdx1],[ds,cpd2,bIdx2],...] and the updated averages.
    """
    excList = []

    dsNum = ds # this is a hack to keep the dataset in the numpy array
    if ds=="5A": dsNum=50
    if ds=="5B": dsNum=51
    if ds=="5C": dsNum=52

    for k in kList:

        # tag outliers from both rate windows
        iE1 = outliersIQR(rates[:,0], k, "hi", ignoreZeros=True)[0]
        iE2 = outliersIQR(rates[:,1], k, "hi", ignoreZeros=True)[0]
        iE = sorted(list(set(np.append(iE1, iE2))))

        # get the average before excluding
        avg1, std1 = getMuStd(rates[:,0], None)
        avg2, std2 = getMuStd(rates[:,1], None)
        if verbose:
            print("DS-%s %s rates:" % (ds, name))
            print("    rLo %.3f ± %-10.3f  %d/%d outliers" % (avg1, std1, len(iE1), len(rates)))
            print("    rHi %.3f ± %-10.3f  %d/%d outliers" % (avg2, std2, len(iE2), len(rates)))
            print("    Applying outlier cut, k=%.1f, %d/%d total outliers" % (k, len(iE), len(rates)))

        for i in iE:
            excList.append( np.asarray([dsNum, int(rates[:,4][i]), int(rates[:,3][i]) ])) # dsNum, cpd, bkgIdx

            if verbose:
                msg = ""
                if i in iE1: msg += "lo "
                if i in iE2: msg += "hi "
                print("    iE %-4d  %-5d  bIdx %-3d  rLo %-8.2f  rHi %-8.2f  %s" % (i, rates[:,4][i], rates[:,3][i], rates[:,0][i], rates[:,1][i], msg))

        # super verbose output
        # print("k value, ",k)
        # for idx in range(len(rates[:,0])):
        #     isOut = "out" if idx in iE else ""
        #     print("idx %-4d  %-5d  bIdx %-3d  rate %-8.2f  %s" % (idx, rates[:,4][idx], rates[:,3][idx], rates[:,0][idx], isOut))

        # exclude iE from rates, like the opposite of np.where
        rates = rates[~np.in1d(range(len(rates)),iE)]

    # if verbose:
        # print("DS-%-3s %s After rate %.3f ± %.3f  nIdx %d" % (ds, name, avg1, std1, len(rates[:,0])))

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
        plt.plot(xVals[i], dataBox[i], ".b", ms=3)

    plt.gca().set_xticklabels(dsList)

    plt.xlabel("Data Set", ha='right', x=1)
    plt.ylabel("Counts/kg-d-keV", ha='right', y=1)

    plt.show()

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


def getReduction():
    """ Calculate how much exposure we lose for a given outliers cut. """

    enrExc, natExc, enrRates, natRates = getOutliers(False)

    # print("enr")
    # for ds,cpd,bIdx in enrExc:
    #     print(ds,cpd,bIdx)

    # print("nat")
    # for ds,cpd,bIdx in natExc:
    #     print(ds,cpd,bIdx)


def makeCutFiles():

    from ROOT import TFile, TTree, MGTWaveform

    enrExc, natExc, enrRates, natRates = getOutliers(True)
    # enrExc, natExc: [:,0]=dsNum, [:,1]=cpd, [:,2]=bkgIdx
    # enrRates, natRates: [:,0]=rate1, [:,1]=rate2, [:,2]=expo, [:,3]=bkgIdx, [:,4]=cpd, [:,5]=ds

    cutType = "fr"
    outType = "frc"

    for ds in [0,1,2,3,4,"5A","5B","5C"]:
    # for ds in ["5B","5C"]:

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

        for bIdx in range(bLo, bHi+1):
            print("DS-%s  bIdx %d" % (ds, bIdx))

            # load the cut files
            for ch in sorted(chList):

                cpd = det.getChanCPD(dsNum,ch)

                if len(np.where((skipList == (dsTmp, int(cpd), bIdx)).all(axis=1))[0]) > 0:
                    # print("skipping det %s in bkgIdx %d" % (cpd, bIdx))
                    continue

                # skip nonexistent files.  other parts of chsel should tell us why these aren't here.
                fName = "%s/bkg/cut/%s/%s_ds%d_%d_ch%d.root" % (dsi.dataDir, cutType, cutType, dsNum, bIdx, ch)
                if not os.path.isfile(fName):
                    # print("no file for det %s in bkgIdx %d" % (cpd, bIdx))
                    continue

                tf = TFile(fName)
                tt = tf.Get("skimTree")
                nEvt = tt.GetEntries()

                # TODO: zombie file check (need to know livetime)

                outName = "%s/bkg/cut/%s/%s_ds%d_%d_ch%d.root" % (dsi.dataDir, outType, outType, dsNum, bIdx, ch)
                outFile = TFile(outName, "RECREATE")
                outTree = TTree()
                outTree = tt.CopyTree("")
                # print("Wrote %d entries." % outTree.GetEntries())

                if nEvt != outTree.GetEntries():
                    print("ERROR, number of entries don't match: input %d  output %d" % (nEvt, outTree.GetEntries()))

                outTree.Write()
                outFile.Close()
                tf.Close()


def plotSpecBeforeAfter():

    from ROOT import TChain

    dsList = [0,1,2,3,4,"5A","5B","5C"]
    # dsList = ["5C"]

    for ds in dsList:

        tB = TChain("skimTree") # before
        tA = TChain("skimTree") # after

        # for ds in dsList:
        dsNum = int(ds[0]) if isinstance(ds,str) else ds
        nBkg = bkg.dsMap()[dsNum]
        bLo, bHi = 0, nBkg
        if ds=="5A": bLo, bHi = 0, 79
        if ds=="5B": bLo, bHi = 80, 112
        if ds=="5C": bLo, bHi = 113, 121
        fB = ["%s/bkg/cut/fr/fr_ds%d_%d_*.root" % (dsi.dataDir, dsNum, bIdx) for bIdx in range(bLo, bHi+1)]
        for f in fB: tB.Add(f)
        fA = ["%s/bkg/cut/frc/frc_ds%d_%d_*.root" % (dsi.dataDir, dsNum, bIdx) for bIdx in range(bLo, bHi+1)]
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

        plt.savefig("../plots/cs4-comp-DS%s.pdf" % ds)

        tB.Reset()
        tA.Reset()


def plotSpectraAfter():

    from ROOT import TChain

    # dsList = [1,2,3,4,"5A","5B","5C"]
    dsList = [1,2,3,4,"5B","5C"]
    # dsList = [0]
    plotName = "../plots/cs4-spec-DS-14-5BC.pdf"

    xLo, xHi, xpb = 0, 20, 0.2
    # xLo, xHi, xpb = 0, 50, 0.1

    # for ds in dsList:

    tt = TChain("skimTree")

    for ds in dsList: # indent this block
        dsNum = int(ds[0]) if isinstance(ds,str) else ds
        nBkg = bkg.dsMap()[dsNum]
        bLo, bHi = 0, nBkg
        if ds=="5A": bLo, bHi = 0, 79
        if ds=="5B": bLo, bHi = 80, 112
        if ds=="5C": bLo, bHi = 113, 121
        fList = ["%s/bkg/cut/frc/frc_ds%d_%d_*.root" % (dsi.dataDir, dsNum, bIdx) for bIdx in range(bLo, bHi+1)]
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
    plt.savefig(plotName)

    tt.Reset()



if __name__=="__main__":
    main(sys.argv[1:])
