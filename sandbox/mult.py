#!/usr/bin/env python3
import sys, os, imp, glob
import numpy as np
import subprocess as sp
import tinydb as db
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
    # plotSelected()
    # selectWideEvents()
    # plotWideEvents()
    # plotParams()
    # plotEfficiency()
    # selectWaveforms()
    # eventMovie()
    # plotSimTest()
    # generateSimSpec()
    # compareDataSimSpec()
    # getLowSimSpec()
    # plotLowSimSpec()
    # getMultip()
    # plotMultip()
    # getChannelRates()
    # getPhysProbability()
    # get2615()
    # plot2615()
    plotCrosstalk()


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
    plt.semilogy(x[1:], sum1, linewidth=0.7, alpha=1.0, ls='steps', color='r', label='mH=1')
    plt.semilogy(x[1:], sum2, linewidth=0.8, alpha=0.8, ls='steps', color='b', label='mH=2')
    plt.semilogy(x[1:], sum3, linewidth=0.9, alpha=0.6, ls='steps', color='m', label='mH=3')
    plt.semilogy(x[1:], sum4, linewidth=1.0, alpha=0.4, ls='steps', color='c', label='mH=4')

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
    x1, x2 = x1[1:], x2[1:]

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
    fs238, rn238 = [], []
    fs583, rn583 = [], []
    lo238chan, lo583chan = [], []

    for iFile, f in enumerate(fileList):
        print("%d/%d %s" % (iFile,nFiles,f))
        tf = TFile("%s/lat/%s" % (ds.specialDir,f))
        tree = tf.Get("skimTree")

        theCut = "mH==2 && gain==0 && sumEH > %f && sumEH < %f && isGood" % (sumLo1, sumHi1)
        tNames = ["Entry$","mH","channel","trapENFCal","sumEH","gain","isGood","fitSlo","riseNoise"]
        tVals = wl.GetVX(tree, tNames, theCut)
        eList = sorted(set(tVals["Entry$"])) # we can do this because the hit lists aren't huge
        for iEnt in eList:
            idx = np.where(tVals["Entry$"]==iEnt)
            if len(idx[0]) != 2:
                continue # sometimes stragglers get through
            hitE = tVals["trapENFCal"][idx]
            chan = tVals["channel"][idx]
            fitSlo = tVals["fitSlo"][idx]
            riseNoise = tVals["riseNoise"][idx]

            iLo = np.argmin(hitE)
            iHi = np.argmax(hitE)
            lo238.append(hitE[iLo])
            hi238.append(hitE[iHi])
            fs238.append(fitSlo[iLo])
            rn238.append(riseNoise[iLo])
            lo238chan.append(chan[iLo])

        theCut = "mH==3 && gain==0 && sumEH > %f && sumEH < %f && isGood" % (sumLo2, sumHi2)
        tNames = ["Entry$","mH","channel","trapENFCal","sumEH","gain","isGood","fitSlo","riseNoise"]
        tVals = wl.GetVX(tree, tNames, theCut)
        eList = sorted(set(tVals["Entry$"]))
        for iEnt in eList:
            idx = np.where(tVals["Entry$"]==iEnt)
            if len(idx[0]) != 3: continue
            hitE = tVals["trapENFCal"][idx]
            chan = tVals["channel"][idx]
            fitSlo = tVals["fitSlo"][idx]
            riseNoise = tVals["riseNoise"][idx]

            iLo = np.argmin(hitE)
            iHi = np.argmax(hitE)
            iMid = list(set([0,1,2]) - set([iLo,iHi]))[0]
            # print(iEnt,hitE,chan, iLo, iHi, iMid)
            lo583.append(hitE[iLo])
            mid583.append(hitE[iMid])
            hi583.append(hitE[iHi])
            fs583.append(fitSlo[iLo])
            rn583.append(riseNoise[iLo])
            lo583chan.append(chan[iLo])

    lo238, hi238 = np.asarray(lo238), np.asarray(hi238)
    lo583, hi583 = np.asarray(lo583), np.asarray(hi583)
    fs238, rn238 = np.asarray(fs238), np.asarray(rn238)
    fs583, rn583 = np.asarray(fs583), np.asarray(rn583)
    lo238chan, lo583chan = np.asarray(lo238chan), np.asarray(lo583chan)
    np.savez("../plots/longCal_hits.npz",lo238, hi238, lo583, mid583, hi583, fs238, rn238, fs583, rn583, lo238chan, lo583chan)


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
    plt.plot(x[1:], y, linewidth=0.7, alpha=1.0, ls='steps', color='r', label='lo hits')

    y, x = np.histogram(hi238, bins=nb, range=(xLo,xHi))
    plt.plot(x[1:], y, linewidth=0.7, alpha=1.0, ls='steps', color='b', label='hi hits')

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
    plt.plot(x[1:], y, linewidth=0.7, alpha=1.0, ls='steps', color='r', label='lo hits')

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
    plt.plot(x[1:], y, linewidth=0.7, alpha=1.0, ls='steps', color='r', label='lo hits')

    y, x = np.histogram(mid583+hi583, bins=nb, range=(xLo,xHi))
    plt.plot(x[1:], y, linewidth=0.7, alpha=1.0, ls='steps', color='b', label='mid+hi hits')

    y, x = np.histogram(hi583, bins=nb, range=(xLo,xHi))
    plt.plot(x[1:], y, linewidth=0.7, alpha=1.0, ls='steps', color='g', label='hi hits')

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
    plt.plot(x[1:], y, linewidth=0.7, alpha=1.0, ls='steps', color='r', label='lo hits')

    plt.xlabel("trapENFCal (keV)", horizontalalignment='right', x=1.0)
    plt.ylabel("Counts", horizontalalignment='right', y=1.0)
    plt.legend(loc='best')
    plt.tight_layout()
    plt.savefig("../plots/mult-m3-low.png")


def selectWideEvents():
    from ROOT import TFile, TTree

    # 5 hr M1 calibration: https://majorana.npl.washington.edu/elog/Run+Elog/1703
    runList = calInfo.GetSpecialRuns("longCal",5)
    # runList = runList[:10] # truncate
    fileList = []
    for run in runList:
        fileList.extend(ds.getLATRunList([run],"%s/lat" % (ds.specialDir)))
    nFiles = len(fileList)

    sumLo2, sumHi2 = 200, 500
    eLo2, fsLo2, chLo2 = [], [], []

    sumLo3, sumHi3 = 200, 500
    eLo3, fsLo3, chLo3 = [], [], []

    for iFile, f in enumerate(fileList):
        print("%d/%d %s" % (iFile,nFiles,f))
        tf = TFile("%s/lat/%s" % (ds.specialDir,f))
        tree = tf.Get("skimTree")

        theCut = "mH==2 && gain==0 && sumEH > %f && sumEH < %f && isGood" % (sumLo2, sumHi2)
        tNames = ["Entry$","mH","channel","trapENFCal","sumEH","gain","isGood","fitSlo","riseNoise"]
        tVals = wl.GetVX(tree, tNames, theCut)
        eList = sorted(set(tVals["Entry$"])) # we can do this because the hit lists aren't huge
        for iEnt in eList:
            idx = np.where(tVals["Entry$"]==iEnt)
            if len(idx[0]) != 2:
                continue # sometimes stragglers get through
            hitE = tVals["trapENFCal"][idx]
            chan = tVals["channel"][idx]
            fitSlo = tVals["fitSlo"][idx]
            riseNoise = tVals["riseNoise"][idx]
            iLo = np.argmin(hitE)
            eLo2.append(hitE[iLo])
            fsLo2.append(fitSlo[iLo])
            chLo2.append(chan[iLo])

        theCut = "mH==3 && gain==0 && sumEH > %f && sumEH < %f && isGood" % (sumLo3, sumHi3)
        tVals = wl.GetVX(tree, tNames, theCut)
        eList = sorted(set(tVals["Entry$"]))
        for iEnt in eList:
            idx = np.where(tVals["Entry$"]==iEnt)
            if len(idx[0]) != 3: continue
            hitE = tVals["trapENFCal"][idx]
            chan = tVals["channel"][idx]
            fitSlo = tVals["fitSlo"][idx]
            riseNoise = tVals["riseNoise"][idx]

            iLo = np.argmin(hitE)
            iHi = np.argmax(hitE)
            iMid = list(set([0,1,2]) - set([iLo,iHi]))[0]

            iLo = np.argmin(hitE)
            eLo3.append(hitE[iLo])
            fsLo3.append(fitSlo[iLo])
            chLo3.append(chan[iLo])

    eLo2, fsLo2, chLo2 = np.asarray(eLo2), np.asarray(fsLo2), np.asarray(chLo2)
    eLo3, fsLo3, chLo3 = np.asarray(eLo3), np.asarray(fsLo3), np.asarray(chLo3)
    np.savez("../plots/longCal_wide.npz", eLo2, fsLo2, chLo2, eLo3, fsLo3, chLo3)


def plotWideEvents():
    f = np.load("../plots/longCal_wide.npz")

    eLo2, fsLo2, chLo2 = f['arr_0'], f['arr_1'], f['arr_2']
    eLo3, fsLo3, chLo3 = f['arr_3'], f['arr_4'], f['arr_5']

    fig = plt.figure(figsize=(9,6))

    xLo, xHi, bpx = 0, 50, 0.3
    nb = int((xHi-xLo)/bpx)

    y, x = np.histogram(eLo2, bins=nb, range=(xLo,xHi))
    plt.plot(x[1:], y, color='r', ls='steps', label='m2')

    y, x = np.histogram(eLo3, bins=nb, range=(xLo,xHi))
    plt.plot(x[1:], y, color='b', ls='steps', label='m3')

    plt.legend(loc='best')
    plt.tight_layout()
    plt.savefig("../plots/wide-mult.png")


    plt.cla()
    plt.plot(eLo2, fsLo2, '.', c='black', markersize=0.5, label="wide-lo-m2", alpha=0.5)
    plt.ylim(0,300)
    plt.xlabel("trapENFCal (keV)", horizontalalignment='right', x=1.0)
    plt.ylabel("fitSlo", horizontalalignment='right', y=1.0)
    plt.legend(loc='best')
    plt.tight_layout()
    plt.savefig("../plots/mult-wide-m2.png")


    plt.cla()
    plt.plot(eLo3, fsLo3, '.', c='black', markersize=0.5, label="wide-lo-m3", alpha=0.5)
    plt.ylim(0,300)
    plt.xlabel("trapENFCal (keV)", horizontalalignment='right', x=1.0)
    plt.ylabel("fitSlo", horizontalalignment='right', y=1.0)
    plt.legend(loc='best')
    plt.tight_layout()
    plt.savefig("../plots/mult-wide-m3.png")


def plotParams():

    f = np.load("../plots/longCal_hits.npz")
    lo238, hi238, lo583, mid583, hi583 = f['arr_0'], f['arr_1'], f['arr_2'], f['arr_3'], f['arr_4']
    fs238, rn238, fs583, rn583 = f['arr_5'], f['arr_6'], f['arr_7'], f['arr_8']

    # 1
    fig = plt.figure(figsize=(9,6), facecolor='w')
    p1 = plt.subplot(111)

    xLo, xHi, bpX = 0, 50, 0.3
    yLo, yHi, bpY = 0, 300, 1.5
    nbX, nbY = int((xHi-xLo)/bpX), int((yHi-yLo)/bpY)
    _,_,_,im = p1.hist2d(lo238, fs238, bins=[nbX, nbY], range=[[xLo,xHi],[yLo,yHi]], norm=LogNorm(), cmap='jet')
    fig.colorbar(im, ax=p1)
    # p1.axhline(fsCut, color='black', linewidth=3)

    p1.set_xlabel("trapENFCal (keV)", horizontalalignment='right',x=1.0)
    p1.set_ylabel("fitSlo", horizontalalignment='right',y=1.0)
    plt.tight_layout()
    plt.savefig("../plots/mult-238-fitSlo.png")

    # 2
    fig.clf()
    p1 = plt.subplot(111)

    xLo, xHi, bpX = 0, 50, 0.3
    yLo, yHi, bpY = 0, 5, 0.01
    nbX, nbY = int((xHi-xLo)/bpX), int((yHi-yLo)/bpY)
    _,_,_,im = p1.hist2d(lo238, rn238, bins=[nbX, nbY], range=[[xLo,xHi],[yLo,yHi]], norm=LogNorm(), cmap='jet')
    fig.colorbar(im, ax=p1)

    p1.set_xlabel("trapENFCal (keV)", horizontalalignment='right',x=1.0)
    p1.set_ylabel("riseNoise", horizontalalignment='right',y=1.0)
    plt.tight_layout()
    plt.savefig("../plots/mult-238-riseNoise.png")

    # 3
    fig.clf()
    p1 = plt.subplot(111)

    xLo, xHi, bpX = 0, 50, 0.3
    yLo, yHi, bpY = 0, 300, 1.5
    nbX, nbY = int((xHi-xLo)/bpX), int((yHi-yLo)/bpY)
    _,_,_,im = p1.hist2d(lo583, fs583, bins=[nbX, nbY], range=[[xLo,xHi],[yLo,yHi]], norm=LogNorm(), cmap='jet')
    fig.colorbar(im, ax=p1)
    # p1.axhline(fsCut, color='black', linewidth=3)

    p1.set_xlabel("trapENFCal (keV)", horizontalalignment='right',x=1.0)
    p1.set_ylabel("fitSlo", horizontalalignment='right',y=1.0)
    plt.tight_layout()
    plt.savefig("../plots/mult-583-fitSlo.png")


    # 4
    fig.clf()
    p1 = plt.subplot(111)

    xLo, xHi, bpX = 0, 50, 0.3
    yLo, yHi, bpY = 0, 5, 0.01
    nbX, nbY = int((xHi-xLo)/bpX), int((yHi-yLo)/bpY)
    _,_,_,im = p1.hist2d(lo583, rn583, bins=[nbX, nbY], range=[[xLo,xHi],[yLo,yHi]], norm=LogNorm(), cmap='jet')
    fig.colorbar(im, ax=p1)

    p1.set_xlabel("trapENFCal (keV)", horizontalalignment='right',x=1.0)
    p1.set_ylabel("riseNoise", horizontalalignment='right',y=1.0)
    plt.tight_layout()
    plt.savefig("../plots/mult-583-riseNoise.png")


    # 5 - scatter, compare fitslo evts from 238 and 583
    fig.clf()
    plt.xlim(0,100)
    plt.ylim(0,300)
    plt.plot(lo238, fs238, '.', c='black', markersize=2., label="lo238", alpha=0.5)
    plt.plot(lo583, fs583, '.', c='red', markersize=3., label="lo583")
    plt.xlabel("trapENFCal (keV)", horizontalalignment='right', x=1.0)
    plt.ylabel("fitSlo", horizontalalignment='right', y=1.0)
    plt.legend(loc='best')
    plt.tight_layout()
    plt.savefig("../plots/mult-238-583-fscompare.png")


def plotEfficiency():

    f = np.load("../plots/longCal_hits.npz")
    lo238, fs238, lo238chan = f['arr_0'], f['arr_5'], f['arr_9']
    lo583, fs583, lo583chan = f['arr_2'], f['arr_7'], f['arr_10']

    fig = plt.figure(figsize=(9,6), facecolor='w')

    # load fitSlo constants for cal run range closest to the run range [22513, 22566]
    dsNum, modNum, calIdx = 5, 1, 11  # calIdx 11: [[22568,22635],22568,22841],
    calDB = db.TinyDB('../calDB.json')
    pars = db.Query()
    fsDict = ds.getDBRecord("fitSlo_ds%d_idx%d_m%d_Peak" % (dsNum, calIdx, modNum), False, calDB, pars)

    # 1. 238
    ePass, fsPass = [], []
    eFail, fsFail = [], []
    eAll, fsAll = [], []
    for idx in range(len(lo238chan)):
        ch = lo238chan[idx]
        if ch not in fsDict.keys(): # e.g. 692 is not in the fsDict
            continue
        fsCut = fsDict[ch][2]
        fsVal = fs238[idx]
        ene = lo238[idx]
        if fsVal <= fsCut:
            ePass.append(ene)
            fsPass.append(fsVal)
        else:
            eFail.append(ene)
            fsFail.append(fsVal)
        eAll.append(ene)
        fsAll.append(fsVal)

    plt.plot(ePass, fsPass, 'o', c='black', fillstyle='full', markersize=1., alpha=0.5, label="pass")
    plt.plot(eFail, fsFail, 'o', c='red', fillstyle='full', markersize=2., label="fail")
    plt.xlim(0,100)
    plt.ylim(0,300)
    plt.xlabel("trapENFCal (keV)", horizontalalignment='right', x=1.0)
    plt.ylabel("fitSlo", horizontalalignment='right', y=1.0)
    plt.legend(loc='best')
    plt.tight_layout()
    plt.savefig("../plots/mult-238-passFail.png")


    fig.clf()
    p1 = plt.subplot(111)

    xLo, xHi, bpx = 0., 50., 0.5
    nbx = int((xHi-xLo)/bpx)

    yAll, x = np.histogram(eAll, bins=nbx, range=(xLo,xHi))
    p1.plot(x[1:], yAll, linewidth=0.7, alpha=1.0, ls='steps', color='r', label='all mH==2')

    yPass, x = np.histogram(ePass, bins=nbx, range=(xLo,xHi))
    p1.plot(x[1:], yPass, linewidth=0.7, alpha=1.0, ls='steps', color='b', label='hits passing')

    p1.set_xlabel("trapENFCal (keV)", horizontalalignment='right', x=1.)
    p1.set_ylabel("Counts", horizontalalignment='right', y=1.)
    p1.legend(loc='best')

    p2 = p1.twinx()

    p2.plot(x[1:], 100.*yPass/yAll, '.r', markersize=5, label='efficiency')

    p2.set_ylim(0,110)
    p2.set_ylabel('% Efficiency', color='r', horizontalalignment='right', y=1.0)
    p2.tick_params('y', colors='black')

    plt.tight_layout()
    plt.savefig("../plots/mult-cutEff.png")


def selectWaveforms():

    from ROOT import TFile, TTree, MGTWaveform

    # 5 hr M1 calibration: https://majorana.npl.washington.edu/elog/Run+Elog/1703
    runList = calInfo.GetSpecialRuns("longCal",5)
    runList = runList[:10] # truncate
    fileList = []
    for run in runList:
        fileList.extend(ds.getLATRunList([run],"%s/lat" % (ds.specialDir)))
    nFiles = len(fileList)

    e1Lo, e1Hi = 235, 242
    sumLo1, sumHi1 = 237.25, 239.37
    # e2Lo, e2Hi = 580, 586
    # sumLo2, sumHi2 = 580.98, 584.29

    e238, fs238, wf238, ch238 = [], [], [], []

    for iFile, f in enumerate(fileList):
        print("%d/%d %s" % (iFile,nFiles,f))
        tf = TFile("%s/lat/%s" % (ds.specialDir,f))
        tree = tf.Get("skimTree")

        theCut = "mH==2 && gain==0 && sumEH > %f && sumEH < %f && isGood" % (sumLo1, sumHi1)
        tNames = ["Entry$","mH","channel","sumEH","gain","isGood","MGTWaveforms","trapENFCal","fitSlo"]
        tVals = wl.GetVX(tree, tNames, theCut)

        eList = sorted(set(tVals["Entry$"]))
        for iEnt in eList:
            idx = np.where(tVals["Entry$"]==iEnt)
            if len(idx[0]) != 2: continue # sometimes stragglers get through

            hitE = tVals["trapENFCal"][idx]
            chan = tVals["channel"][idx]
            fitSlo = tVals["fitSlo"][idx]
            wfs = tVals["MGTWaveforms"][idx]

            iLo = np.argmin(hitE)
            if hitE[iLo] > 50: continue

            e238.append(hitE[iLo])
            fs238.append(fitSlo[iLo])
            ch238.append(chan[iLo])
            wf238.append(wfs[iLo])

    np.savez("../plots/longCal-m2.npz",e238, fs238, ch238, wf238)


def eventMovie():
    """ adapted from LAT/sandbox/wave-movie.py """
    from matplotlib import animation

    f = np.load("../plots/longCal-m2.npz")
    e238, fs238, ch238, wf238 = f['arr_0'], f['arr_1'], f['arr_2'], f['arr_3']
    # nList = 100
    nList = len(wf238)

    outFile = "../plots/movie-longCal-lo238.mp4"

    # load fitSlo constants for cal run range closest to the run range [22513, 22566]
    dsNum, modNum, calIdx = 5, 1, 11  # calIdx 11: [[22568,22635],22568,22841],
    calDB = db.TinyDB('../calDB.json')
    pars = db.Query()
    fsDict = ds.getDBRecord("fitSlo_ds%d_idx%d_m%d_Peak" % (dsNum, calIdx, modNum), False, calDB, pars)

    fig = plt.figure(figsize=(11,5), facecolor='w')
    a1 = plt.subplot(111)
    a1.set_xlabel("time (ns)")
    a1.set_ylabel("ADC")
    p1, = a1.plot(np.ones(1), np.ones(1), color='blue')

    def init():
        p1.set_data([],[])
        return p1,

    def animate(iList):
        energy = e238[iList]
        chan = ch238[iList]
        fitSlo = fs238[iList]
        signal = wf238[iList]

        if chan not in fsDict.keys(): return p1, # e.g. 692 is not in the fsDict
        fsCut = fsDict[chan][2]
        p1.set_color('blue') if fitSlo < fsCut else p1.set_color('red')

        waveRaw = signal.GetWaveRaw()
        waveTS = signal.GetTS()
        p1.set_ydata(waveRaw)
        p1.set_xdata(waveTS)
        titleStr = "Entry %d  Channel %d  trapENFCal %.2f  fitSlo %.2f  fsCut %.2f" % (iList,chan,energy,fitSlo,fsCut)
        plt.title(titleStr)
        # print(titleStr)

        xmin, xmax = np.amin(waveTS), np.amax(waveTS)
        a1.set_xlim([xmin,xmax])
        if energy < 400:
            a1.set_ylim(0,400)
        else:
            ymin, ymax = np.amin(waveRaw), np.amax(waveRaw)
            a1.set_ylim([ymin-abs(0.1*ymin),ymax+abs(0.1*ymax)])

        if iList%500 == 0 and iList!=0:
            print("%d / %d entries saved (%.2f %% done)." % (iList,nList,100*(float(iList)/nList)))
        return p1,

    anim = animation.FuncAnimation(fig, animate, init_func=init, frames=nList, interval=0, blit=True)
    anim.save(outFile, fps=20, extra_args=['-vcodec', 'libx264'])


def plotSimTest():
    """ Using code from sims people, create a simulated spectrum. """
    from ROOT import TChain

    config = "DS5"
    utils = ds.simUtils(config)
    module = "M1"
    basePath = "/global/projecta/projectdirs/majorana/sim/MJDG41003GAT/"
    sourceType = "linesource"
    partClass = "%sCalSource" % module
    segment = "A224_Z88"
    simPath = "%s/MJDemonstrator/%s/%s/%s" % (basePath,sourceType,partClass,segment)
    specPath = "../plots/"

    simList = sorted(glob.glob("%s/processed_*.root" % simPath))
    simChain = TChain("simTree")
    for f in simList[:2]: simChain.Add(f)

    auxList = sorted(glob.glob("%s/aux_processed_*.root" % simPath))
    auxChain = TChain("auxTree_%s" % config)
    for f in auxList[:2]: auxChain.Add(f)

    simChain.AddFriend(auxChain)

    # Count the number of primaries to get the correct factor for finding the efficiency.
    nPrimaries = sum([int(fName.split('_')[-2]) for fName in simList[:2]])

    # Get detectors active in this module & dataset
    detList = []
    for iD, det in enumerate(utils.GetDetectorList(module)):
        if utils.activeDets[module][config][iD] == 1:
            detList.append(det)

    # Set up 3 histograms: All events, granularity cut, multisite cut
    xLo, xHi, bpx = 0, 3000, 1.
    nb = int((xHi-xLo)/bpx)
    xRaw = np.zeros(nb)
    xGran = np.zeros(nb)
    xPSA = np.zeros(nb)

    fig = plt.figure(figsize=(9,6))

    for iD, det in enumerate(detList):

        theCut = "isGoodDet_%s*(fWaveformID==%s)/%d" % (config,det,nPrimaries)
        nRaw = simChain.Draw("fEnergy*1000.0", theCut, "GOFF")
        xArr = simChain.GetV1()
        xArr = np.asarray([xArr[i] for i in range(nRaw)])
        y, x = np.histogram(xArr, bins=nb, range=(xLo, xHi))
        xRaw = np.add(xRaw, y)

        theCut += "*(mH_%s==1)" % config
        nGran = simChain.Draw("fEnergy*1000.0", theCut, "GOFF")
        xArr = simChain.GetV1()
        xArr = np.asarray([xArr[i] for i in range(nGran)])
        y, x = np.histogram(xArr, bins=nb, range=(xLo, xHi))
        xGran = np.add(xGran, y)

        theCut += "*(fDtHeuristic < %s)" % utils.GetDTCutoff(module, det)
        nPSA = simChain.Draw("fEnergy*1000.0", theCut, "GOFF")
        xArr = simChain.GetV1()
        xArr = np.asarray([xArr[i] for i in range(nPSA)])
        y, x = np.histogram(xArr, bins=nb, range=(xLo, xHi))
        xPSA = np.add(xPSA, y)

        print("%s  raw %d  gran %d  psa %d" % (det, nRaw, nGran, nPSA))

        if iD > 5:
            break

    plt.semilogy(x[1:], xRaw, linewidth=0.7, alpha=1.0, ls='steps', color='r', label='raw')
    plt.semilogy(x[1:], xGran, linewidth=0.7, alpha=1.0, ls='steps', color='b', label='+mH==1')
    plt.semilogy(x[1:], xPSA, linewidth=0.7, alpha=1.0, ls='steps', color='g', label='+PSA')
    plt.legend(loc='best')
    plt.xlabel("fEnergy (keV)", horizontalalignment='right', x=1.0)
    plt.ylabel("Counts (arb)", horizontalalignment='right', y=1.0)
    plt.savefig("../plots/mult-sim-test.png")


def generateSimSpec():
    # FIXME!  HIT energy is being plotted instead of SUM.

    from ROOT import TChain

    config, module = "DS5", "M1"
    basePath = "/global/projecta/projectdirs/majorana/sim/MJDG41003GAT/"
    sourceType, partClass, segment = "linesource", "%sCalSource" % module, "A224_Z88"
    simPath = "%s/MJDemonstrator/%s/%s/%s" % (basePath,sourceType,partClass,segment)

    simInfo = ds.SimInfo(config)
    detList = simInfo.GetActiveDets(config, module)

    nLimit = 20
    simList = sorted(glob.glob("%s/processed_*.root" % simPath))
    simChain = TChain("simTree")
    for f in simList[:nLimit]: simChain.Add(f) # limit for now

    nPrimaries = sum([int(fName.split('_')[-2]) for fName in simList[:nLimit]])

    auxList = sorted(glob.glob("%s/aux_processed_*.root" % simPath))
    auxChain = TChain("auxTree_%s" % config)
    for f in auxList[:nLimit]: auxChain.Add(f)
    simChain.AddFriend(auxChain)

    xLo, xHi, bpx = 0, 4000, 2.
    nb = int((xHi-xLo)/bpx)
    sm1, sm2, sm3, sm4 = np.zeros(nb), np.zeros(nb), np.zeros(nb), np.zeros(nb)

    # fill sim histograms
    for iD, det in enumerate(detList):

        theCut = "isGoodDet_%s * (fWaveformID==%s) * (mH_%s==1) / %d" % (config, det, config, nPrimaries)
        np1 = simChain.Draw("fEnergy*1000.0", theCut, "GOFF")
        xArr = simChain.GetV1()
        xArr = np.asarray([xArr[i] for i in range(np1)])
        y,x = np.histogram(xArr, bins=nb, range=(xLo, xHi))
        sm1 = np.add(sm1, y)

        theCut = "isGoodDet_%s * (fWaveformID==%s) * (mH_%s==2) / %d" % (config, det, config, nPrimaries)
        np2 = simChain.Draw("fEnergy*1000.0", theCut, "GOFF")
        xArr = simChain.GetV1()
        xArr = np.asarray([xArr[i] for i in range(np2)])
        y,x = np.histogram(xArr, bins=nb, range=(xLo, xHi))
        sm2 = np.add(sm2, y)

        theCut = "isGoodDet_%s * (fWaveformID==%s) * (mH_%s==3) / %d" % (config, det, config, nPrimaries)
        np3 = simChain.Draw("fEnergy*1000.0", theCut, "GOFF")
        xArr = simChain.GetV1()
        xArr = np.asarray([xArr[i] for i in range(np3)])
        y,x = np.histogram(xArr, bins=nb, range=(xLo, xHi))
        sm3 = np.add(sm3, y)

        theCut = "isGoodDet_%s * (fWaveformID==%s) * (mH_%s==4) / %d" % (config, det, config, nPrimaries)
        np4 = simChain.Draw("fEnergy*1000.0", theCut, "GOFF")
        xArr = simChain.GetV1()
        xArr = np.asarray([xArr[i] for i in range(np4)])
        y,x = np.histogram(xArr, bins=nb, range=(xLo, xHi))
        sm4 = np.add(sm4, y)

        print(det,np1,np2,np3,np4)

    np.savez("../plots/longCal-sim-mH.npz",x, sm1, sm2, sm3, sm4)


def compareDataSimSpec():
    # both files were generated with this binning:
    # xLo, xHi, bpx = 0, 4000, 2.

    f1 = np.load("../plots/longCal-sim-mH.npz")
    x, sm1, sm2, sm3, sm4 = f1['arr_0'], f1['arr_1'], f1['arr_2'], f1['arr_3'], f1['arr_4']

    f2 = np.load("../plots/longCalSumSpec.npz")
    x, dm1, dm2, dm3, dm4 = f2['arr_0'], f2['arr_1'], f2['arr_2'], f2['arr_3'], f2['arr_4']

    x = x[1:]

    # scale the sims by the 2615 peak
    idx = np.where((x > 2610) & (x < 2620))
    sim2615 = np.amax(sm1[idx])
    data2615 = np.amax(dm1[idx])
    dsr = data2615/sim2615

    fig = plt.figure(figsize=(10,7))
    p0 = plt.subplot2grid((3,1), (0,0), rowspan=2)
    p1 = plt.subplot2grid((3,1), (2,0), sharex=p0)

    p0.semilogy(x, dm1, ls='steps', color='r', label='data m=1')
    p0.semilogy(x, sm1*dsr, ls='steps', color='b', alpha=0.6, label='sim m=1')
    p0.set_xlabel("Energy (keV)", horizontalalignment='right',x=1.0)
    p0.set_ylabel("Counts (arb)")
    p0.legend()
    p1.plot(x, sm1-dm1)
    p1.set_ylabel("Residual (sim-data)")
    plt.tight_layout()
    plt.savefig("../plots/mult-sim-m1.png")

    p0.cla()
    p1.cla()
    p0.semilogy(x, dm2, ls='steps', color='r', label='data m=2')
    p0.semilogy(x, sm2*dsr, ls='steps', color='b', alpha=0.6, label='sim m=2')
    p0.set_xlabel("Energy (keV)", horizontalalignment='right',x=1.0)
    p0.set_ylabel("Counts (arb)")
    p0.legend()
    p1.plot(x, sm2-dm2)
    p1.set_ylabel("Residual (sim-data)")
    plt.tight_layout()
    plt.savefig("../plots/mult-sim-m2.png")

    p0.cla()
    p1.cla()
    p0.semilogy(x, dm3, ls='steps', color='r', label='data m=3')
    p0.semilogy(x, sm3*dsr, ls='steps', color='b', alpha=0.6, label='sim m=3')
    p0.set_xlabel("Energy (keV)", horizontalalignment='right',x=1.0)
    p0.set_ylabel("Counts (arb)")
    p0.legend()
    p1.plot(x, sm3-dm3)
    p1.set_ylabel("Residual (sim-data)")
    plt.tight_layout()
    plt.savefig("../plots/mult-sim-m3.png")

    p0.cla()
    p1.cla()
    p0.semilogy(x, dm4, ls='steps', color='r', label='data m=4')
    p0.semilogy(x, sm4*dsr, ls='steps', color='b', alpha=0.6, label='sim m=4')
    p0.set_xlabel("Energy (keV)", horizontalalignment='right',x=1.0)
    p0.set_ylabel("Counts (arb)")
    p0.legend()
    p1.plot(x, sm4-dm4)
    p1.set_ylabel("Residual (sim-data)")
    plt.tight_layout()
    plt.savefig("../plots/mult-sim-m4.png")


def getLowSimSpec():

    from ROOT import TChain

    config, module = "DS5", "M1"
    basePath = "/global/projecta/projectdirs/majorana/sim/MJDG41003GAT/"
    sourceType, partClass, segment = "linesource", "%sCalSource" % module, "A224_Z88"
    simPath = "%s/MJDemonstrator/%s/%s/%s" % (basePath,sourceType,partClass,segment)

    simInfo = ds.SimInfo(config)
    detList = simInfo.GetActiveDets(config, module)

    nLimit = 2
    simList = sorted(glob.glob("%s/processed_*.root" % simPath))
    simChain = TChain("simTree")
    for f in simList[:nLimit]: simChain.Add(f) # limit for now

    nPrimaries = sum([int(fName.split('_')[-2]) for fName in simList[:nLimit]])

    auxList = sorted(glob.glob("%s/aux_processed_*.root" % simPath))
    auxChain = TChain("auxTree_%s" % config)
    for f in auxList[:nLimit]: auxChain.Add(f)
    simChain.AddFriend(auxChain)

    xLo, xHi, bpx = 0, 50, 0.2
    nb = int((xHi-xLo)/bpx)
    sm1 = np.zeros(nb)

    # fill sim histograms
    for iD, det in enumerate(detList):
        print(det)

        theCut = "isGoodDet_%s * (fWaveformID==%s) * (mH_%s==1) / %d" % (config, det, config, nPrimaries)
        np1 = simChain.Draw("fEnergy*1000.0", theCut, "GOFF")
        xArr = simChain.GetV1()
        xArr = np.asarray([xArr[i] for i in range(np1)])
        y,x = np.histogram(xArr, bins=nb, range=(xLo, xHi))
        sm1 = np.add(sm1, y)

    np.savez("../plots/sim-lowe-spec.npz",x,sm1)


def plotLowSimSpec():

    f = np.load("../plots/sim-lowe-spec.npz")
    x, sm1 = f['arr_0'], f['arr_1']

    fig = plt.figure(figsize=(9,6))

    plt.plot(x[1:], sm1, ls='steps', color='red', label='sim spec, mH==1, all dets')
    plt.xlabel("Energy (keV)", horizontalalignment='right', x=1.)
    plt.legend(loc='best')
    plt.tight_layout()
    plt.savefig("../plots/sim-lowEspec.png")


def getMultip():

    from ROOT import TFile, TTree, GATDataSet

    # 5 hr M1 calibration: https://majorana.npl.washington.edu/elog/Run+Elog/1703
    runList = calInfo.GetSpecialRuns("longCal",5)
    runList = runList[:15] # this gets 10k mH=4 events
    fileList = []
    for run in runList:
        fileList.extend(ds.getLATRunList([run],"%s/lat" % (ds.specialDir)))
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

    xLo, xHi, bpx = 0, 39, 1
    nbx = int((xHi-xLo)/bpx)
    mult = np.zeros(nbx)

    for iFile, f in enumerate(fileList):

        tf = TFile("%s/lat/%s" % (ds.specialDir,f))
        latTree = tf.Get("skimTree")

        tNames = ["Entry$","mH","channel","trapENFCal","gain","isGood"]
        theCut = ""
        nPass = latTree.Draw("mH",theCut,"GOFF")
        mH = latTree.GetV1()
        mH = [mH[i] for i in range(nPass)]
        if len(mH) != latTree.GetEntries():
            # make sure # in draw matches # entries (don't double count hits)
            print("ERROR.  len(mH) %d  tree entries %d" % (len(mH), latTree.GetEntries()))
            return

        y, x = np.histogram(mH, bins=nbx, range=(xLo, xHi))
        mult = np.add(mult, y)

        m4 = 100 * mult[4] / 10000.
        print("%d/%d %s  m4: %d  %.2f%%" % (iFile,nFiles,f,mult[4],m4))

        tf.Close()

    np.savez("../plots/longCal-mH.npz",x,mult,runTime)


def plotMultip():
    f = np.load("../plots/longCal-mH.npz")
    x, mult, runTime = f['arr_0'], f['arr_1'], f['arr_2']

    rates = mult/runTime

    for i in range(1,6):
        cts = mult[i]
        rate = cts/runTime
        err = np.sqrt(cts)/runTime
        pctErr = 100/np.sqrt(cts)
        print("mH=%d  cts %-8d  rate %.2f pm %-8.5f  (%.5f%%)" % (i, cts, rate, err, pctErr))

    fig = plt.figure(figsize=(9,6))
    # plt.semilogy(x[:-1], rates, ls='steps', color='black')

    plt.bar(x[:-1], rates)

    ax = plt.gca()
    ax.set_yscale('log')
    plt.xticks(x[:12])
    plt.xlim(0,12)
    plt.ylim(0.01)

    plt.xlabel("mH", horizontalalignment='right', x=1.0)
    plt.ylabel("Rate (Hz)")

    plt.title("mH. runTime %.0f sec.  mH=4 rate: %.2f pm %.4f Hz" % (runTime, rates[4], np.sqrt(mult[4])/runTime))

    plt.savefig("../plots/longCal-mH.png")


def getChannelRates():
    # from plotChanNoiseRate2
    chRates = {
        584: (0.02504, 0.65272), 592: (0.04016, 1.12538), 598: (0.04214, 7.43857), 608: (0.04133, 0.89446),
        610: (0.05058, 1.16650), 614: (0.02029, 1.03391), 624: (0.02900, 0.93005), 626: (0.03535, 0.96318),
        628: (0.04316, 1.48850), 632: (0.02492, 0.57832), 640: (0.03153, 0.73312), 648: (0.03063, 0.57880),
        658: (0.01900, 0.52327), 660: (0.02465, 0.46342), 662: (0.02693, 0.60194), 672: (0.01873, 0.50940),
        674: (0.02377, 0.62129), 678: (0.02555, 0.58299), 680: (0.01812, 0.55507), 688: (0.04286, 1.03254),
        690: (0.01906, 0.62837), 692: (0.02236, 0.77410), 694: (0.03285, 0.83018),
        1106: (0.00210, 0.26005), 1120: (0.01229, 0.46665), 1124: (0.00514, 0.24493), 1128: (0.00198, 0.23542),
        1170: (0.01055, 0.43628), 1172: (0.00535, 0.45415), 1174: (0.00204, 0.12108), 1176: (0.00595, 0.20878),
        1204: (0.00472, 0.21207), 1208: (0.01247, 0.37574), 1232: (0.00201, 0.08488), 1236: (0.00234, 0.07277),
        1298: (0.00078, 0.09937), 1302: (0.01199, 0.32350), 1330: (0.01064, 0.40703), 1332: (0, 0)
        }
    rNoise = sum([chRates[key][1] + chRates[key][0] for key in chRates]) # un-subtract rate from 40-50 kev

    xLabels = sorted(chRates.keys())
    y = [chRates[ch][0] + chRates[ch][1] for ch in xLabels]
    x = np.arange(0, len(chRates), 1)

    noiseAvg = np.mean(y)
    noiseTot = np.sum(y)

    fig = plt.figure(figsize=(9,6))
    plt.bar(x, y)
    plt.xticks(x, xLabels, rotation=85.)
    plt.axhline(noiseAvg, color='green', label='Avg. Noise Rate: %.3f Hz' % noiseAvg)
    plt.plot(np.nan, np.nan, color='w', label='Total Noise Rate: %.3f Hz' % noiseTot)
    plt.legend(loc='best')
    plt.savefig("../plots/mult-channel.png")


def getPhysProbability():
    from scipy.stats import poisson

    rMult = [0, 344.73, 69.60, 12.14, 1.96, 0.58] # mult. rates

    # from plotChanNoiseRate2
    chRates = {
        584: (0.02504, 0.65272), 592: (0.04016, 1.12538), 598: (0.04214, 7.43857), 608: (0.04133, 0.89446),
        610: (0.05058, 1.16650), 614: (0.02029, 1.03391), 624: (0.02900, 0.93005), 626: (0.03535, 0.96318),
        628: (0.04316, 1.48850), 632: (0.02492, 0.57832), 640: (0.03153, 0.73312), 648: (0.03063, 0.57880),
        658: (0.01900, 0.52327), 660: (0.02465, 0.46342), 662: (0.02693, 0.60194), 672: (0.01873, 0.50940),
        674: (0.02377, 0.62129), 678: (0.02555, 0.58299), 680: (0.01812, 0.55507), 688: (0.04286, 1.03254),
        690: (0.01906, 0.62837), 692: (0.02236, 0.77410), 694: (0.03285, 0.83018),
        1106: (0.00210, 0.26005), 1120: (0.01229, 0.46665), 1124: (0.00514, 0.24493), 1128: (0.00198, 0.23542),
        1170: (0.01055, 0.43628), 1172: (0.00535, 0.45415), 1174: (0.00204, 0.12108), 1176: (0.00595, 0.20878),
        1204: (0.00472, 0.21207), 1208: (0.01247, 0.37574), 1232: (0.00201, 0.08488), 1236: (0.00234, 0.07277),
        1298: (0.00078, 0.09937), 1302: (0.01199, 0.32350), 1330: (0.01064, 0.40703), 1332: (0, 0)
        }
    rNoise = sum([chRates[key][1] + chRates[key][0] for key in chRates]) # un-subtract rate from 40-50 kev
    print("raw noise rate:  %.3f" % rNoise) # 28.709

    pNHits = [0]

    # calculate probability of pure physics event vs. physics + noise
    for n in range(1, 6):
        print("N ==",n)

        dt = 4e-6 if n <= 2 else 4e-6 * (n-1)

        pPhys = poisson.pmf(1, rMult[n] * dt) # remember it's the probability of getting ONE mH==n event

        pNoise = poisson.pmf(n, rNoise * dt) # this is prob. of getting n noise events

        pairs = [ (j, n-j) for j in range(1, n) ]
        pMix = 0
        for k, l in pairs:
            phys = poisson.pmf(1, rMult[k] * dt) # one mH==k event in dt
            noise = poisson.pmf(l, rNoise * dt) # l noise events in dt
            pMix += phys*noise
            print("(%d,%d)   phys %.4e  noise %.4e  prod %.4e  mix %.4e" % (k, l, phys, noise, phys*noise, pMix))

        print("Result: phys %.4e  noise %.4e   mix %.4e" % (pPhys, pNoise, pMix))

        pTot = pPhys + pNoise + pMix
        fPhys, fNoise, fMix = 100*pPhys/pTot, 100*pNoise/pTot, 100*pMix/pTot
        print("------> phys %.4f%%    noise %.4f%%   mix %.4f%%" % (fPhys, fNoise, fMix))

        print(" ")

        pNHits.append(pNoise)

    pNTot = sum(pNHits)
    print("expected mH fraction, all noise evts (< 10 kev)")
    for i in range(1,6):
        print("mH==%d, %.2e" % (i, pNHits[i]/pNTot))


def get2615():
    from ROOT import TFile, TTree

    # 5 hr M1 calibration: https://majorana.npl.washington.edu/elog/Run+Elog/1703
    runList = calInfo.GetSpecialRuns("longCal",5)
    # runList = runList[:10] # truncate
    fileList = []
    for run in runList:
        fileList.extend(ds.getLATRunList([run],"%s/lat" % (ds.specialDir)))
    nFiles = len(fileList)

    chList = ds.GetGoodChanListNew(5)

    xLo, xHi, xpb = 2605.0, 2625.0, 0.2
    nb = int((xHi-xLo)/xpb)
    h2615 = []

    xTalkEvts = []

    for iFile, f in enumerate(fileList):
        print("%d/%d %s" % (iFile,nFiles,f))
        tf = TFile("%s/lat/%s" % (ds.specialDir,f))
        tree = tf.Get("skimTree")

        tNames = ["channel","mH","trapENFCal"]
        tree.SetBranchStatus('*',0)
        for name in tNames: tree.SetBranchStatus(name,1)

        for iEnt in range(tree.GetEntries()):
            tree.GetEntry(iEnt)

            nHits = tree.channel.size()
            hitList = [i for i in range(nHits) if int(tree.channel.at(i)) in chList]
            hitE = np.asarray([tree.trapENFCal.at(i) for i in hitList])

            idx2615 = np.where((hitE > xLo) & (hitE < xHi))
            if len(idx2615[0]) == 0: continue
            hit2615 = hitE[idx2615][0]
            h2615.append(hit2615)

            # now save events where the hit is within peak mean + 2 sigma (only 0.38% bkg contamination)
            if not 2612.051 < hit2615 < 2616.971: continue
            hitCh = np.asarray([tree.channel.at(i) for i in hitList])
            xTalkEvts.append([hitE, hitCh])

        tf.Close()

    y, x = np.histogram(h2615, bins=nb, range=(xLo, xHi))
    np.savez("../plots/longCal-2615.npz", x, y, xTalkEvts)


def plot2615():
    f = np.load("../plots/longCal-2615.npz")
    x, h2615 = f['arr_0'], f['arr_1']
    x = x[1:]

    fig = plt.figure(figsize=(9,6))

    plt.plot(x, h2615, ls='steps', color='blue', label='pk2615')

    p0 = (np.mean(h2615[:5]), max(h2615), 2614.511, roughSigma(2614.511))
    popt,_ = curve_fit(gaus, x, h2615, p0=p0)
    bpx = x[1] - x[0]
    bgRate = popt[0] * bpx
    mu, sig = popt[2], popt[3]
    idx = np.where((x > mu-2*sig) & (x < mu+2*sig))  # yes, 2 sigma
    totCts = np.sum(h2615[idx])
    bgCts = bgRate * len(h2615[idx])
    pkCts = totCts - bgCts
    pbr = pkCts / bgCts
    plt.plot(x, gaus(x, *popt), 'r-', label='fit.  mu %.2f  sig %.2f\nP/B %.2f' % (mu,sig,pbr))

    plt.axvline(mu+2*sig, color='green', label='2 sigma')
    plt.axvline(mu-2*sig, color='green')
    plt.plot(np.nan, np.nan, 'w', label='bkg in peak: %.2f%%' % (100 * bgCts / totCts))

    plt.xlabel("trapENFCal (keV)", horizontalalignment='right', x=1.0)
    plt.legend()
    plt.savefig("../plots/longCal-2615-peak.png")


def plotCrosstalk():

    f = np.load("../plots/longCal-2615.npz")
    evts, nTot = f['arr_2'], f['arr_3']
    nTot = len(evts)

    eAny, eLoOnly = [], []

    nHaveNoise = len(evts)
    nOnlyNoise = 0
    for evt in evts:
        hitE, hitCh = evt[0], evt[1]
        mGood = len(hitE)

        idxLo = np.where(hitE <= 10.)
        mNoise = len(idxLo[0])
        for hit in hitE[idxLo]:
            eLoOnly.append(hit)
            # eAny.append(hit)

        idxMid = np.where((hitE > 10.) & (hitE < 2612.051))
        mMid = len(idxMid[0])
        for hit in hitE[idxMid]:
            eAny.append(hit)

        idxHi = np.where(hitE > 2612.051)
        mHi = len(idxHi[0])

        # print("mG %d  mN %d  mM %d  mHi %d" % (mGood, mNoise, mMid, mHi))

        if mNoise > 0 and mMid == 0:
            nOnlyNoise += 1

        # return

    print("%d total evts w/ hit=2615, %d have any hits < 10 keV, %d have only hits < 10 keV" % (nTot, nHaveNoise, nOnlyNoise))

    xLo, xHi, xpb = -5., 30., 0.2
    nb = int((xHi-xLo)/xpb)

    fig = plt.figure(figsize=(9,6))

    yAny, x = np.histogram(eAny, bins=nb, range=(xLo, xHi))
    plt.plot(x[1:], yAny, ls='steps', color='green', linewidth=2, label='eAny')

    yOnly, x = np.histogram(eLoOnly, bins=nb, range=(xLo, xHi))
    plt.plot(x[1:], yOnly, ls='steps', color='red', label='eOnly')

    plt.legend()
    plt.tight_layout()
    plt.savefig("../plots/longCal-2615-eLoHits.png")

    x, y = wl.GetHisto(mAll, nb, xLo, xHi, xpb)
    # plt.semilogy(x, y, ls='steps', color='r', label='mAll')
    plt.bar(x-xpb/2., y, 0.95, color='r', log=True, label='All Hits E < 2614: %d' % len(mAll))

    x, y = wl.GetHisto(mMid, nb, xLo, xHi, xpb)
    plt.bar(x-xpb/2., y, 0.9, color='b', log=True, label='Hits 10 < E < 2614: %d' % len(mMid))

    x, y = wl.GetHisto(mNoise, nb, xLo, xHi, xpb)
    plt.bar(x-xpb/2., y, 0.85, color='g', log=True, label='Hits E < 10: %d' % len(mNoise))

    plt.xlabel("multiplicity", horizontalalignment='right', x=1.)
    plt.xlim(xLo, xHi)
    plt.legend()
    plt.tight_layout()
    # plt.savefig("../plots/longCal-2615-mult.png")
    plt.show()




if __name__=="__main__":
    main()
