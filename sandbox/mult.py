#!/usr/bin/env python3
import sys, os, imp, glob
import numpy as np
import subprocess as sp
from scipy.optimize import curve_fit
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
sys.argv.append("-b")
ds = imp.load_source('DataSetInfo',os.environ['LATDIR']+'/DataSetInfo.py')
wl = imp.load_source('waveLibs',os.environ['LATDIR']+'/waveLibs.py')
calInfo = ds.CalInfo()

def main():

    # genSpectrum()
    # plotSpectrum()
    # genPeakSpec()
    # plotPeakSpec()
    # selectEvents()
    plotSelected()


def getSumEne(tree, theCut):
    """ get sum energy for all events passing cuts """
    sumArr = []
    tNames = ["Entry$","mH","channel","trapENFCal","sumEH","gain","isGood"]
    tVals = wl.GetVX(tree, tNames, theCut, False)
    nPass = len(tVals["Entry$"])
    prevEnt, sumE = -1, 0.
    for idx in range(nPass):
        ent = tVals["Entry$"][idx]
        mH = tVals["mH"][idx]
        chan = tVals["channel"][idx]
        tmp = tVals["sumEH"][idx]
        enf = tVals["trapENFCal"][idx]
        if enf > 99999: enf = 0
        if ent!=prevEnt and idx!=0:
            # print("%.2f" % sumE)
            sumArr.append(sumE)
            sumE = 0
        sumE += enf
        prevEnt = ent
        # print("%d  %d  %d  %.2f  %2f" % (ent,chan,mH,enf,sumE))
    # print("%.2f" % sumE)
    sumArr.append(sumE)
    sumArr = np.asarray(sumArr)
    return sumArr


def genSpectrum():
    from ROOT import TFile, TTree

    # 5 hr M1 calibration: https://majorana.npl.washington.edu/elog/Run+Elog/1703
    runList = calInfo.GetSpecialRuns("longCal",5)
    fileList = []
    for run in runList:
        fileList.extend(ds.getLATRunList([run],"%s/lat" % (ds.specialDir)))
    nFiles = len(fileList)

    xLo, xHi, bpx = 0., 4000., 2.
    nbx = int((xHi-xLo)/bpx)
    sum1 = np.zeros(nbx)
    sum2 = np.zeros(nbx)
    sum3 = np.zeros(nbx)
    sum4 = np.zeros(nbx)

    for iFile, f in enumerate(fileList):

        print("%d/%d %s" % (iFile,nFiles,f))
        sp.call("""free -g | grep "cache:" | awk '{print $4}'""", shell=True) # prints free RAM

        tf = TFile("%s/lat/%s" % (ds.specialDir,f))
        latTree = tf.Get("skimTree")

        sumArr = getSumEne(latTree, "mH==1 && gain==0 && isGood")
        y, x = np.histogram(sumArr, bins=nbx, range=(xLo,xHi))
        sum1 = np.add(sum1, y)

        sumArr = getSumEne(latTree, "mH==2 && gain==0 && isGood")
        y, x = np.histogram(sumArr, bins=nbx, range=(xLo,xHi))
        sum2 = np.add(sum2, y)

        sumArr = getSumEne(latTree, "mH==3 && gain==0 && isGood")
        y, x = np.histogram(sumArr, bins=nbx, range=(xLo,xHi))
        sum3 = np.add(sum3, y)

        sumArr = getSumEne(latTree, "mH==4 && gain==0 && isGood")
        y, x = np.histogram(sumArr, bins=nbx, range=(xLo,xHi))
        sum4 = np.add(sum4, y)

        tf.Close()

        if iFile > 200:
            break

    np.savez("../plots/longCalSumSpec.npz",x,sum1,sum2,sum3,sum4)


def plotSpectrum():

    f = np.load("../plots/longCalSumSpec.npz")
    x, sum1, sum2, sum3, sum4 = f['arr_0'], f['arr_1'], f['arr_2'], f['arr_3'], f['arr_4']

    fig = plt.figure(figsize=(9,6), facecolor='w')
    plt.semilogy(x[:-1], sum1, linewidth=0.7, alpha=1.0, ls='steps', color='r', label='mH=1')
    plt.semilogy(x[:-1], sum2, linewidth=0.8, alpha=0.8, ls='steps', color='b', label='mH=2')
    plt.semilogy(x[:-1], sum3, linewidth=0.9, alpha=0.6, ls='steps', color='m', label='mH=3')
    plt.semilogy(x[:-1], sum4, linewidth=1.0, alpha=0.4, ls='steps', color='c', label='mH=4')

    plt.xlabel("SumE (keV)", horizontalalignment='right', x=1.)
    plt.ylabel("Counts", horizontalalignment='right', y=1.)
    plt.legend(loc='best')

    plt.tight_layout()
    plt.savefig("../plots/mult-test.pdf")


def genPeakSpec():

    from ROOT import TFile, TTree

    # 5 hr M1 calibration: https://majorana.npl.washington.edu/elog/Run+Elog/1703
    runList = calInfo.GetSpecialRuns("longCal",5)
    runList = runList[:10] # truncate
    fileList = []
    for run in runList:
        fileList.extend(ds.getLATRunList([run],"%s/lat" % (ds.specialDir)))
    nFiles = len(fileList)

    bpx = 0.1

    e1Lo, e1Hi = 235, 242
    nb1 = int((e1Hi-e1Lo)/bpx)
    pk238 = np.zeros(nb1)

    e2Lo, e2Hi = 580, 586
    nb2 = int((e2Hi-e2Lo)/bpx)
    pk583 = np.zeros(nb2)

    for iFile, f in enumerate(fileList):
        print("%d/%d %s" % (iFile,nFiles,f))

        tf = TFile("%s/lat/%s" % (ds.specialDir,f))
        tree = tf.Get("skimTree")

        theCut = "mH==2 && gain==0 && sumEH > %f && sumEH < %f && isGood" % (e1Lo, e1Hi)
        sumArr = getSumEne(tree, theCut)
        y, x1 = np.histogram(sumArr, bins=nb1, range=(e1Lo,e1Hi))
        pk238 = np.add(pk238, y)

        theCut = "mH==3 && gain==0 && sumEH > %f && sumEH < %f && isGood" % (e2Lo, e2Hi)
        sumArr = getSumEne(tree, theCut)
        y, x2 = np.histogram(sumArr, bins=nb2, range=(e2Lo,e2Hi))
        pk583 = np.add(pk583, y)

        tf.Close()
        # if iFile > 10:
            # break

    np.savez("../plots/longCalPeaks.npz",x1,x2,pk238,pk583)


def roughSigma(ene):
    """ from gpxFitter, just used to set initial guesses """
    p0, p1, p2 = 0.2, 0.02, 0.0003
    return np.sqrt(p0**2. + p1**2. * ene + p2**2. * ene**2.)


def gaus(x, b, a, mu, sig):
    """ gaussian + flat bg """
    return b + a * np.exp(-(x-mu)**2. / (2. * sig**2.))


def plotPeakSpec():
    """ Plot sum peaks and fit them to gaussians.
    Print the peak widths (so we know how tight to make the cut)
    and the peak/background ratio (so we can estimate # accidentals)
    """
    f = np.load("../plots/longCalPeaks.npz")
    x1, x2, pk238, pk583 = f['arr_0'], f['arr_1'], f['arr_2'], f['arr_3']
    x1, x2 = x1[:-1], x2[:-1]

    fig = plt.figure(figsize=(9,6), facecolor='w')

    plt.plot(x1, pk238, linewidth=0.7, alpha=1.0, ls='steps', color='b', label='mH==2')

    p0 = (np.mean(pk238[:5]), max(pk238), 238., roughSigma(238.))
    popt,_ = curve_fit(gaus, x1, pk238, p0=p0)
    bpx = x1[1] - x1[0]
    bgRate = popt[0] * bpx
    mu, sig = popt[2], popt[3]
    idx = np.where((x1 > mu-3*sig) & (x1 < mu+3*sig))
    totCts = np.sum(pk238[idx])
    bgCts = bgRate * len(pk238[idx])
    pkCts = totCts - bgCts
    pbr = pkCts / bgCts
    print("%d  %.2f  %.2f  %.2f  %.2f  %d  %d  %d  P/B: %.2f" % (238, bpx, bgRate, mu, sig, totCts, bgCts, pkCts, pbr))
    print("sumLo, sumHi = %.2f, %.2f" % (mu-3*sig, mu+3*sig))

    plt.plot(x1, gaus(x1, *popt), 'r-', label='fit.  mu %.2f  sig %.2f\nP/B %.2f' % (popt[2], popt[3], pbr))

    plt.title("Peak-to-bkg ratio: %.3f" % pbr)
    plt.xlabel("SumE (keV)", horizontalalignment='right', x=1.)
    plt.ylabel("Counts", horizontalalignment='right', y=1.)
    plt.legend(loc='best')
    plt.tight_layout()
    plt.savefig("../plots/mult-m2-238.png")

    plt.cla()

    plt.plot(x2, pk583, linewidth=0.7, alpha=1.0, ls='steps', color='b', label='mH==3')

    p0 = (np.mean(pk583[:5]), max(pk583), 583., roughSigma(583.))
    popt,_ = curve_fit(gaus, x2, pk583, p0=p0)
    bpx = x2[1] - x2[0]
    bgRate = popt[0] * bpx
    mu, sig = popt[2], popt[3]
    idx = np.where((x2 > mu-3*sig) & (x2 < mu+3*sig))
    totCts = np.sum(pk583[idx])
    bgCts = bgRate * len(pk583[idx])
    pkCts = totCts - bgCts
    pbr = pkCts / bgCts
    print("%d  %.2f  %.2f  %.2f  %.2f  %d  %d  %d  P/B: %.2f" % (583, bpx, bgRate, mu, sig, totCts, bgCts, pkCts, pbr))
    print("sumLo, sumHi = %.2f, %.2f" % (mu-3*sig, mu+3*sig))

    plt.plot(x2, gaus(x2, *popt), 'r-', label='fit.  mu %.2f  sig %.2f\nP/B %.2f' % (popt[2], popt[3], pbr))

    plt.title("Peak-to-bkg ratio: %.3f" % pbr)
    plt.xlabel("SumE (keV)", horizontalalignment='right', x=1.)
    plt.ylabel("Counts", horizontalalignment='right', y=1.)
    plt.legend(loc='best')
    plt.tight_layout()
    plt.savefig("../plots/mult-m3-583.png")


def selectEvents():

    from ROOT import TFile, TTree

    # 5 hr M1 calibration: https://majorana.npl.washington.edu/elog/Run+Elog/1703
    runList = calInfo.GetSpecialRuns("longCal",5)
    # runList = runList[:10] # truncate
    fileList = []
    for run in runList:
        fileList.extend(ds.getLATRunList([run],"%s/lat" % (ds.specialDir)))
    nFiles = len(fileList)

    e1Lo, e1Hi = 235, 242
    sumLo1, sumHi1 = 237.25, 239.37

    e2Lo, e2Hi = 580, 586
    sumLo2, sumHi2 = 580.98, 584.29

    lo238, hi238 = [], []
    lo583, mid583, hi583 = [], [], []

    for iFile, f in enumerate(fileList):
        print("%d/%d %s" % (iFile,nFiles,f))
        tf = TFile("%s/lat/%s" % (ds.specialDir,f))
        tree = tf.Get("skimTree")

        theCut = "mH==2 && gain==0 && sumEH > %f && sumEH < %f && isGood" % (sumLo1, sumHi1)
        tNames = ["Entry$","mH","channel","trapENFCal","sumEH","gain","isGood"]
        tVals = wl.GetVX(tree, tNames, theCut)
        eList = sorted(set(tVals["Entry$"])) # we can do this because the hit lists aren't huge
        for iEnt in eList:
            idx = np.where(tVals["Entry$"]==iEnt)
            if len(idx[0]) != 2:
                continue # sometimes stragglers get through
            hitE = tVals["trapENFCal"][idx]
            chan = tVals["channel"][idx]
            hitLo, hitHi = min(hitE), max(hitE)
            chLo, chHi = chan[np.argmin(hitE)], chan[np.argmax(hitE)]
            # print("%d  %.2f  %d  %.2f  %d" % (iEnt, hitLo, chLo, hitHi, chHi))
            lo238.append(hitLo)
            hi238.append(hitHi)

        theCut = "mH==3 && gain==0 && sumEH > %f && sumEH < %f && isGood" % (sumLo2, sumHi2)
        tNames = ["Entry$","mH","channel","trapENFCal","sumEH","gain","isGood"]
        tVals = wl.GetVX(tree, tNames, theCut)
        eList = sorted(set(tVals["Entry$"]))
        for iEnt in eList:
            idx = np.where(tVals["Entry$"]==iEnt)
            if len(idx[0]) != 3: continue
            hitE = tVals["trapENFCal"][idx]
            chan = tVals["channel"][idx]

            iLo = np.argmin(hitE)
            iHi = np.argmax(hitE)
            iMid = list(set([0,1,2]) - set([iLo,iHi]))[0]
            # print(iEnt,hitE,chan, iLo, iHi, iMid)

            lo583.append(hitE[iLo])
            mid583.append(hitE[iMid])
            hi583.append(hitE[iHi])

    lo238, hi238 = np.asarray(lo238), np.asarray(hi238)
    lo583, hi583 = np.asarray(lo583), np.asarray(hi583)
    np.savez("../plots/longCal_hits.npz",lo238, hi238, lo583, mid583, hi583)


def plotSelected():

    f = np.load("../plots/longCal_hits.npz")
    lo238, hi238, lo583, mid583, hi583 = f['arr_0'], f['arr_1'], f['arr_2'], f['arr_3'], f['arr_4']
    critE238 = 114.8  # 238 peak
    critE583 = 405.36 # 583 peak

    fig = plt.figure(figsize=(9,6), facecolor='w')

    # 1
    xLo, xHi, bpx = 0., 242., 2.
    nb = int((xHi-xLo)/bpx)

    y, x = np.histogram(lo238, bins=nb, range=(xLo,xHi))
    plt.plot(x[:-1], y, linewidth=0.7, alpha=1.0, ls='steps', color='r', label='lo hits')

    y, x = np.histogram(hi238, bins=nb, range=(xLo,xHi))
    plt.plot(x[:-1], y, linewidth=0.7, alpha=1.0, ls='steps', color='b', label='hi hits')

    plt.axvline(critE238, color='green', linewidth=2., alpha=0.3, label='Crit. E: %.2f keV' % critE238)

    plt.xlabel("trapENFCal (keV)", horizontalalignment='right', x=1.0)
    plt.ylabel("Counts", horizontalalignment='right', y=1.0)
    plt.legend(loc='best')
    plt.tight_layout()
    plt.savefig("../plots/mult-m2.png")

    # 2
    plt.cla()

    xLo, xHi, bpx = 0., 20., 0.1
    nb = int((xHi-xLo)/bpx)

    y, x = np.histogram(lo238, bins=nb, range=(xLo,xHi))
    plt.plot(x[:-1], y, linewidth=0.7, alpha=1.0, ls='steps', color='r', label='lo hits')

    plt.xlabel("trapENFCal (keV)", horizontalalignment='right', x=1.0)
    plt.ylabel("Counts", horizontalalignment='right', y=1.0)
    plt.legend(loc='best')
    plt.tight_layout()
    plt.savefig("../plots/mult-m2-low.png")

    # 3
    plt.cla()

    xLo, xHi, bpx = 0, 585, 2.
    nb = int((xHi-xLo)/bpx)

    y, x = np.histogram(lo583, bins=nb, range=(xLo,xHi))
    plt.plot(x[:-1], y, linewidth=0.7, alpha=1.0, ls='steps', color='r', label='lo hits')

    y, x = np.histogram(mid583+hi583, bins=nb, range=(xLo,xHi))
    plt.plot(x[:-1], y, linewidth=0.7, alpha=1.0, ls='steps', color='b', label='mid+hi hits')

    y, x = np.histogram(hi583, bins=nb, range=(xLo,xHi))
    plt.plot(x[:-1], y, linewidth=0.7, alpha=1.0, ls='steps', color='g', label='hi hits')

    plt.axvline(critE583, color='green', linewidth=2., alpha=0.3, label='Crit. E: %.2f keV' % critE583)
    plt.axvline(177.64, color='green', linewidth=2., alpha=0.3, label='Crit. E: %.2f keV' % 177.64)

    plt.xlabel("trapENFCal (keV)", horizontalalignment='right', x=1.0)
    plt.ylabel("Counts", horizontalalignment='right', y=1.0)
    plt.legend(loc='best')
    plt.tight_layout()
    plt.savefig("../plots/mult-m3.png")

    # 4
    plt.cla()

    xLo, xHi, bpx = 0., 20., 0.1
    nb = int((xHi-xLo)/bpx)

    y, x = np.histogram(lo583, bins=nb, range=(xLo,xHi))
    plt.plot(x[:-1], y, linewidth=0.7, alpha=1.0, ls='steps', color='r', label='lo hits')

    plt.xlabel("trapENFCal (keV)", horizontalalignment='right', x=1.0)
    plt.ylabel("Counts", horizontalalignment='right', y=1.0)
    plt.legend(loc='best')
    plt.tight_layout()
    plt.savefig("../plots/mult-m3-low.png")





if __name__=="__main__":
    main()