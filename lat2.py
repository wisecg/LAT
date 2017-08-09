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
    """ ./lat2.py -cal [options]
    Source: http://nucleardata.nuclear.lu.se/toi/radSearch.asp
    Pb212 : 238.6322, Ra224 : 240.9866
    """

    # File I/O
    inPath = "waveSkimDS5_run21975_NLC2.root"  # placeholder, latskim file goes here
    inFile = TFile(inPath)

    latTree = inFile.Get("skimTree")
    theCut = inFile.Get("theCut").GetTitle()

    # parameters to calibrate:
    iPar = ["trapENM","trapENF"]
    # iPar = ["trapENM","fitAmp","lat","latFC","latAFC"]
    oPar = ["trapENMCal","fitE","latE","latFCE","latAFCE"]

    # keep a figure
    fig = plt.figure(figsize=(9,6), facecolor='w')
    p0 = plt.subplot(111)
    plt.show(block=False)

    # loop over channels
    for ch in ds.DetID[dsNum]:

        if ch%2==1: continue
        value = raw_input()
        if value=='q': break

        # can do up to 4 vars in one draw.  TODO: add multisite cut?
        n = latTree.Draw("trapENM:trapENF", theCut + " && trapENM > 100 && channel==%d" % ch, "GOFF")
        if n==0:
            print "Chan %d - no events found." % ch
            continue

        pk1 = peakInfo(latTree.GetV1(), n, ch)

        print ch, pk1.fail()


        # title = "Chan %d  Pk %.3f  Const %.3f  CalFWHM %.3f  PeakCts %d  Chi2 %.2f" % (ch,mu,const,fwhm*const,nCts,chi2NDF)
        # print title
        #
        # p0.cla()
        # p0.margins(x=0,y=0)
        # p0.plot(xvals, hist, color='blue', ls='steps-post', label='data')
        # p0.axvline(bigPk, color='orange', label='pkfinder')
        # p0.axvline(mu, color='magenta', label='centroid')
        # p0.plot(xvals, fit, color='red', label='bestfit')
        # p0.set_title(title)
        # p0.legend()
        # plt.tight_layout()
        plt.pause(0.00001)

class peakInfo:
    """ Stores the results from 'find238Peak' """
    def __init__(self, vec, n, ch):

        r = find238Peak(vec, n, ch)

        self.Fail = True
        self.Ch, self.Mu, self.Fwhm = 0,0,0
        self.Spd, self.Cts, self.Chi2 = 0,0,0
        self.Fit, self.Hist, self.Xvals = 0,0,0
        self.BigPk = 0
        if len(r)!=0:
            self.Fail = False
            self.Ch, self.Mu, self.Fwhm = r[0], r[1], r[2]
            self.Spd, self.Cts, self.Chi2 = r[3], r[4], r[5]
            self.Fit, self.Hist, self.Xvals = r[6], r[7], r[8]
            self.BigPk = r[9]

    def fail(self): return self.Fail
    def ch(self): return self.Ch
    def mu(self): return self.Mu
    def fwhm(self): return self.Fwhm
    def spd(self): return self.Spd
    def cts(self): return self.Cts
    def chi2(self): return self.Chi2
    def fit(self): return self.Fit
    def hist(self): return self.Hist
    def xvals(self): return self.Xvals
    def bigPk(self): return self.BigPk


def find238Peak(vec, nEnt, ch):
    """ Assumes that the 238 peak is the highest peak in the spectrum.
        This means you have to feed it a cal file truncated at 250 kev.
    """
    start = time.clock()

    # read in raw data from a tree.GetVX() object
    lst = []
    for i in range(nEnt): lst.append(vec[i])
    arr = np.asarray(lst)

    # coarse histogram
    h1, x1 = np.histogram(arr,bins=1500)
    x1 = x1[:-1]
    pks, cts = wl.GetPeaks(h1, x1, 30)

    if len(pks) < 1:
        print "Chan. %d - Couldn't find peaks." % ch
        return []
    bigPk = pks[ np.argmax(cts) ]
    bigCt = cts[ np.argmax(cts) ]

    # fine histogram
    lo, hi = bigPk - 0.01*bigPk, bigPk + 0.02*bigPk
    hist, xvals = np.histogram(arr, bins=150, range=(lo,hi))
    xvals = xvals[:-1]

    # fit step
    pars, guess = [0.,0.,0.,0.], [bigCt, bigPk, 0.8, 2.]
    try:
        pars,_ = curve_fit(wl.peakModel238, xvals, hist, p0=guess)
    except ValueError:
        print "Chan. %d - ValueError. ydata or xdata contain nan's." % ch
        return []
    except RuntimeError:
        print "Chan. %d - RuntimeError.  Leastsq minimization failed." % ch
        return []

    # calculate stuff
    mu, sig = pars[1], pars[2]
    fwhm = pars[2] * 2.3548
    const = 238.6322 / mu
    fitSpeed = time.clock()-start
    idx = np.where((xvals >= mu-3.*sig) & (xvals <= mu+3.*sig))
    nCts = np.sum(hist[idx])
    fit = wl.peakModel238(xvals, *pars)
    chi2NDF = np.sum( np.square(hist[idx] - fit[idx]) / fit[idx] ) / (len(fit[idx]) - 4.)

    return [ch,mu,const,fwhm,fitSpeed,nCts,chi2NDF,fit,hist,xvals,bigPk]

def updateFile(dsNum,subNum,runNum,fPaths):
    """ ./lat2.py -upd [options] """

    print "hi"


if __name__ == "__main__":
    main(sys.argv[1:])