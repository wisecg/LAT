#!/usr/local/bin/python
"""
===================== LAT2.py =====================

Calibrate energy parameters in LAT data.

Two modes:
    -cal : Scans a calibration range and updates parameters in a pandas file.
    -upd : Updates an input file with data from the calibration database file.

v1: 07 Aug 2017

========= C. Wiseman (USC), B. Zhu (LANL) =========
"""
import sys, time, os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import DataSetInfo as ds
import waveLibs as wl
from ROOT import TFile, TTree
from scipy.optimize import curve_fit

homePath = os.path.expanduser('~')
bgDir = homePath + "/project/lat"
calDir = homePath + "/project/cal-lat"

def main(argv):

    print "========================================"
    print "LAT2 started:",time.strftime('%X %x %Z')
    startT = time.clock()

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
        print len(arr),"entries."

        # peak finder - assume the 238 peak is the highest peak in the spectrum.

        # coarse
        h1, x1 = np.histogram(arr,bins=1000)
        x1 = x1[:-1]

        fig = plt.figure(figsize=(9,6), facecolor='w')
        plt.plot(x1, h1, ls='steps-post')
        plt.show()

        pk1, ct1 = wl.GetPeaks(h1, x1, 200)
        if len(pk1)<1:
            print "Coarse couldn't find peaks, ch %d.  Continuing ..." % ch
            continue
        bigPk = pk1[ np.argmax(ct1) ]

        plt.axvline(bigPk, color='green')
        plt.show()

        continue

        # fine
        lo, hi = bigPk - 0.05*bigPk, bigPk + 0.05*bigPk
        h2, x2 = np.histogram(arr, bins=100, range=(lo, hi))
        x2 = x2[:-1]
        pk2, ct2 = wl.GetPeaks(h2, x2, 10)
        if len(pk2) < 1:
            print "Fine couldn't find peaks, ch %d.  Continuing ..." % ch
            continue
        bigPk = pk2[ np.argmax(ct2) ]

        if len(pk2)==2:
            print "Ch. %d  Ratio of counts in big peak (%d) to little peak (%d) = %.2f  (Scan = %.2f sec)" % (ch,ct2[0], ct2[1], float(ct2[0])/ct2[1], time.clock()-start)
        else:
            print "Ch %d found %d peaks.  (Scan = %.2f sec)" % (ch, len(pk2), time.clock()-start)

        fig = plt.figure(figsize=(9,6), facecolor='w')
        plt.plot(x2, h2, ls='steps-post')
        for pk in pk2: plt.axvline(pk, color='green')
        plt.show()


        # fit




        # return


def updateFile(dsNum,subNum,runNum,fPaths):
    """ ./lat2.py -upd [options] """

    print "hi"


if __name__ == "__main__":
    main(sys.argv[1:])