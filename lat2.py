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
import tinydb as db
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
    fCal, fUpd, fSub, fRun, fTest, fBat = 0,0,0,0,0,0
    fPaths = [".","."]
    if len(argv)==0: return
    for i,opt in enumerate(argv):
        if opt == "-test":
            fTest = True
            print "Test mode."
            dbStuff()
        if opt == "-cal":
            fCal = True
            print "Calibration mode."
        if opt == "-sub":
            fSub = True
            print "File update mode."
        if opt == "-b":
            fBat = True
            print "Batch mode."
        if opt == "-f":
            fRun = True
            print "Scanning DS-%d, run %d" % (dsNum, runNum)
        if opt == "-s":
            fSub, dsNum, subNum = True, int(argv[i+1]), int(argv[i+2])
            print "Scanning DS-%d sub-range %d" % (dsNum, subNum)
        if opt == "-p":
            fPaths = [argv[i+1], argv[i+2]]
            print "Manual I/O paths set."

    if fCal: calibrateRuns(dsNum,subNum,runNum,fPaths,fBat)
    if fUpd: updateFile(dsNum,subNum, runNum,fPaths)

    stopT = time.clock()
    print "Stopped:",time.strftime('%X %x %Z'),"\nProcess time (min):",(stopT - startT)/60


def calibrateRuns(dsNum,subNum,runNum,fPaths,batMode):
    """ ./lat2.py -cal [options]
    Source: http://nucleardata.nuclear.lu.se/toi/radSearch.asp
    Pb212 : 238.6322, Ra224 : 240.9866
    """

    inPath = "latSkimDS5_run21975.root"
    inFile = TFile(inPath)
    latTree = inFile.Get("skimTree")
    theCut = inFile.Get("theCut").GetTitle()

    # this is the record we'll return
    calRecord = {}

    fig = plt.figure(figsize=(9,6), facecolor='w')
    p0 = plt.subplot(111)
    if not batMode: plt.show(block=False)

    # loop over channels
    for ch in sorted(ds.DetID[dsNum]):

        # skip dumb stuff
        if ch%2==1: continue
        if ch in ds.PMon[dsNum]: continue
        if ds.CPD[dsNum][ch] == 1 or ds.CPD[dsNum][ch] == 2: continue

        # set const default
        calRecord[ch] = [-999,-999,-999,-999]

        # TODO: change to multisite cut 'nMS'
        n = latTree.Draw("trapENF:fitAmp:latAF:latAFC", theCut + " && trapENF > 100 && channel==%d && avse >-1" % ch, "GOFF")
        if n==0:
            print "Chan %-4d - TotCts 0" % ch
            continue

        # do fitting and pull parameters
        pk1 = peakInfo(latTree.GetV1(), n, ch)
        if pk1.fail(): continue
        pk2 = peakInfo(latTree.GetV2(), n, ch)
        if pk2.fail(): continue
        pk3 = peakInfo(latTree.GetV3(), n, ch)
        if pk3.fail(): continue
        pk4 = peakInfo(latTree.GetV4(), n, ch)
        if pk4.fail(): continue
        par1 = pk1.GetResults()
        par2 = pk2.GetResults()
        par3 = pk3.GetResults()
        par4 = pk4.GetResults()

        # write to the record
        calRecord[ch] = [par1[1],par2[1],par3[1],par4[1]]

        # print status message
        totCts, pkCts, fitChi2 = par1[3], par1[4], par1[5]
        print "Chan %-4d - TotCts %d  PeakCts %d  Chi2 %.2f  calRecord[%d]:[%.3f,%.3f,%.3f,%.3f]" % (ch,totCts,pkCts,fitChi2,ch,calRecord[ch][0],calRecord[ch][1],calRecord[ch][2],calRecord[ch][3])

        # plotting
        if not batMode:
            mu, const, fwhm, totCts, cts, chi2, bigPk = pk1.GetResults()
            xvals, hist, fit = pk1.GetArrays()

            title = "Chan %-4d - Pk %.3f  Const %.3f  CalFWHM %.3f  TotCts %d  PeakCts %d  Chi2 %.2f" % (ch,mu,const,fwhm*const,totCts,cts,chi2)
            print title

            value = raw_input()
            if value=='q': break
            p0.cla()
            p0.margins(x=0,y=0)
            p0.plot(xvals, hist, color='blue', ls='steps-post', label='data')
            p0.axvline(bigPk, color='orange', label='pkfinder')
            p0.axvline(mu, color='magenta', label='centroid')
            p0.plot(xvals, fit, color='red', label='bestfit')
            p0.set_title(title)
            p0.legend()
            plt.tight_layout()
            plt.pause(0.00001)

    return calRecord


class peakInfo:
    """ Stores the results from find238Peak.
        Mainly this is useful for error handling.
    """
    def __init__(self, vec, n, ch):

        result = find238Peak(vec, n, ch)

        self.Fail = True
        self.Ch, self.Mu, self.Fwhm, self.Const = 0,0,0,0
        self.Spd, self.totCts, self.Cts, self.Chi2 = 0,0,0,0
        self.Fit, self.Hist, self.Xvals = 0,0,0
        self.BigPk = 0
        if len(result)!=0:
            self.Fail = False
            self.Ch, self.Mu, self.Fwhm, self.Const, self.Spd, self.totCts, self.Cts, self.Chi2, self.Fit, self.Hist, self.Xvals, self.BigPk = result

    def fail(self): return self.Fail
    def GetResults(self): return [self.Mu, self.Const, self.Fwhm, self.totCts, self.Cts, self.Chi2, self.BigPk]
    def GetArrays(self): return [self.Xvals, self.Hist, self.Fit]


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
        print "Chan %-4d - TotCts %d, peakdet fail." % (ch, nEnt)
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
        print "Chan. %-4d - ValueError. ydata or xdata contain nan's." % ch
        return []
    except RuntimeError:
        print "Chan. %-4d - RuntimeError.  Leastsq minimization failed." % ch
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

    return [ch,mu,fwhm,const,fitSpeed,nEnt,nCts,chi2NDF,fit,hist,xvals,bigPk]


def dbStuff():
    """ ./lat2.py -test """

    calDB = db.TinyDB('calDB.json') # set DB file
    pars = db.Query()

    key1 = 'ds5_idx999'

    if len( calDB.search(pars.key == 'ds5_idx999') )==0:
        print "record doesn't exist, adding it ...."
        vals = calibrateRuns(5,-1,21975,[".","."],True)
        calDB.insert({'key':key1,"vals":vals})

    # returns a list - multiple keys are possible in principle
    search = calDB.search(pars.key == 'ds5_idx999')
    record = search[0]
    print record['key']
    # print record['vals']
    print record['vals']['640']



def updateFile(dsNum,subNum,runNum,fPaths):
    """ ./lat2.py -upd [options] """

    print "hi"


if __name__ == "__main__":
    main(sys.argv[1:])