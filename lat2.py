#!/usr/local/bin/python
"""
===================== LAT2.py =====================

Calibrate energy parameters in LAT data.

Modes:
    -cal  : Scans a calibration range and updates parameters in a pandas file.
    -upd  : Updates an input file with data from the calibration database file.
    -test : IDEA: Quick check if a given range has enough events in each channel to calibrate.

v1: 07 Aug 2017

========= C. Wiseman (USC), B. Zhu (LANL) =========
"""
import sys, time, os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import DataSetInfo as ds
import waveLibs as wl
from ROOT import gROOT, TFile, TTree
from scipy.optimize import curve_fit

homePath = os.path.expanduser('~')
bgDir = homePath + "/project/lat"
calDir = homePath + "/project/cal-lat"

def main(argv):

    print "========================================"
    print "LAT2 started:",time.strftime('%X %x %Z')
    startT = time.clock()
    gROOT.ProcessLine("gErrorIgnoreLevel = 3001;") # suppress ROOT error messages

    dsNum, subNum, runNum = -1, -1, -1
    fCal, fUpd, fSub, fRun = 0,0,0,0
    fPaths = [".","."]
    if len(argv)==0: return
    for i,opt in enumerate(argv):
        if opt == "-cal":
            fCal = True
            print "Calibration mode."
        if opt == "-sub":
            fSub = True
            print "File update mode."
        if opt == "-f":
            fRun = True
            print "Scanning DS-%d, run %d" % (dsNum, runNum)
        if opt == "-s":
            fSub, dsNum, subNum = True, int(argv[i+1]), int(argv[i+2])
            print "Scanning DS-%d sub-range %d" % (dsNum, subNum)
        if opt == "-p":
            fPaths = [argv[i+1], argv[i+2]]
            print "Manual I/O paths set."

    if fCal: calibrateRuns(dsNum,subNum,runNum,fPaths)
    if fUpd: updateFile(dsNum,subNum, runNum,fPaths)

    stopT = time.clock()
    print "Stopped:",time.strftime('%X %x %Z'),"\nProcess time (min):",(stopT - startT)/60


def calibrateRuns(dsNum,subNum,runNum,fPaths):
    """ ./lat2.py -cal [options] """

    # File I/O
    inPath = "waveSkimDS5_run21975_NLC2.root"  # placeholder, latskim file goes here
    inFile = TFile(inPath)

    latTree = inFile.Get("skimTree")
    theCut = inFile.Get("theCut").GetTitle()

    # parameters to calibrate:
    iPar = ["trapENM"]
    # iPar = ["trapENM","fitAmp","lat","latFC","latAFC"]
    oPar = ["trapENMCal","fitE","latE","latFCE","latAFCE"]

    # keep a figure
    fig = plt.figure(figsize=(9,6), facecolor='w')
    p0 = plt.subplot(111)
    plt.show(block=False)

    # loop over channels
    for ch in ds.DetID[dsNum]:

        if ch%2==1: continue

        start = time.clock()

        # can do up to 4 vars in one draw
        n = latTree.Draw("trapENM", theCut + " && trapENM > 100 && channel==%d" % ch, "GOFF")
        if n==0:
            print "ch %d, no events found." % ch
            continue
        vec = latTree.GetV1()
        lst = []
        for i in range(n): lst.append(vec[i])
        arr = np.asarray(lst)
        print "channel",ch,": ",len(arr),"entries."

        # peak finder - only assumes the 238 peak is the highest peak in the spectrum.

        # coarse step
        h1, x1 = np.histogram(arr,bins=1000)
        x1 = x1[:-1]
        pk1, ct1 = wl.GetPeaks(h1, x1, 30)

        if len(pk1)<1:
            print "Coarse step couldn't find peaks, ch %d.  Continuing ..." % ch
            continue
        bigPk = pk1[ np.argmax(ct1) ]

        # fine step
        lo, hi = bigPk - 0.01*bigPk, bigPk + 0.02*bigPk
        h2, x2 = np.histogram(arr, bins=150, range=(lo, hi))
        x2 = x2[:-1]
        pk2, ct2 = wl.GetPeaks(h2, x2, 15)

        if len(pk2) < 1:
            print "Fine step couldn't find peaks, ch %d.  Continuing ..." % ch
            continue
        bigPk = pk2[ np.argmax(ct2) ]
        bigCt = ct2[ np.argmax(ct2) ]

        # fit step
        # from http://nucleardata.nuclear.lu.se/toi/radSearch.asp
        # Pb212 : 238.6322, Ra224 : 240.9866

        pars, guess = [0.,0.,0.,0.], [bigCt, bigPk, 0.8, 2.]
        try:
            pars,_ = curve_fit(wl.peakModel238, x2, h2, p0=guess)
        except ValueError:
            print "Channel %d, ValueError. ydata or xdata contain nan's." % ch
        except RuntimeError:
            print "Channel %d, RuntimeError.  Leastsq minimization failed." % ch

        # calculate stuff
        mu, sig = pars[1], pars[2]
        fwhm = pars[2] * 2.3548
        const = 238.6322 / mu
        fitSpeed = time.clock()-start
        idx = np.where((x2 >= mu-3.*sig) & (x2 <= mu+3.*sig))
        nCts = np.sum(h2[idx])
        fit = wl.peakModel238(x2, *pars)
        chi2NDF = np.sum( np.square(h2[idx] - fit[idx]) / fit[idx] ) / (len(fit[idx]) - 4.)

        title = "Chan %d  Uncal.Peak %.3f  Cal.Const %.3f  Cal.FWHM %.3f\n  Speed %.3f  PeakCts %d  Chi2NDF %.2f" % (ch,mu,const,fwhm*const,fitSpeed,nCts,chi2NDF)

        print title

        value = raw_input()
        if value=='q': break
        p0.cla()
        p0.margins(x=0,y=0)
        p0.plot(x2, h2, ls='steps-post')
        for pk in pk2: p0.axvline(pk, color='green')
        p0.plot(x2, fit, color='red')
        p0.set_title(title)
        plt.tight_layout()
        plt.pause(0.0000001)


def updateFile(dsNum,subNum,runNum,fPaths):
    """ ./lat2.py -upd [options] """

    print "hi"


if __name__ == "__main__":
    main(sys.argv[1:])