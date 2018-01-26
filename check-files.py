#!/usr/bin/env python
import sys, itertools
sys.argv.append("-b")
import DataSetInfo as ds
import tinydb as db
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import scipy.special as sp

""" TODO:
Make this check integrity of all files:
TTrees work, have nonzero entries, branches all work, etc.
"""

def main():
    checkFiles()

def checkFiles():
    """ ./job-panda.py -checkFiles """
    from ROOT import TFile, TTree
    import os.path, imp
    import DataSetInfo as ds

    # Check BG skim and waveskim files
    # dsMap = {0:75,1:51,2:7,3:24,4:18,5:112}
    dsMap = {2:7}
    for dsNum in dsMap:
        for sub in range(dsMap[dsNum]+1):

            # check skims
            fileName = "%s/skimDS%d_%d_low.root" % (ds.skimDir,dsNum,sub)
            if not os.path.isfile(fileName):
                print("file not found! name:", fileName)
                continue
            f1 = TFile(fileName)
            t1 = f1.Get("skimTree")
            n1 = t1.GetEntries()
            print("DS %d  sub %d  skim entries %d" % (dsNum, sub, n1))
            if n1==0:
                print("no skim entries found! file:", fileName)
                continue

            # check waveskims
            fileName = "%s/waveSkimDS%d_%d.root" % (ds.skimDir,dsNum,sub)
            if not os.path.isfile(fileName):
                print("file not found! name:", fileName)
                continue
            f2 = TFile(fileName)
            t2 = f2.Get("skimTree")
            n2 = t2.GetEntries()
            print("DS %d  sub %d  wave entries %d" % (dsNum, sub, n2))
            if n2==0:
                print("no waveskim entries found! file:", fileName)
                continue

    # Check CAL skim and waveskim files
    calList = getCalRunList(dsNum=2) # None checks all ds's
    for run in calList:
        dsNum=-1
        for key in ds.dsRanges:
            if ds.dsRanges[key][0] <= run <= ds.dsRanges[key][1]:
                dsNum=key

        # check skims
        fileName = "%s/skimDS%d_run%d_low.root" % (ds.calSkimDir,dsNum,run)
        if not os.path.isfile(fileName):
            print("file not found! name:", fileName)
            continue
        f1 = TFile(fileName)
        t1 = f1.Get("skimTree")
        n1 = t1.GetEntries()
        print("DS %d  run %d  skim entries %d" % (dsNum, run, n1))
        if n1==0:
            print("no skim entries found! file:", fileName)
            continue

        # check waveskims
        fileName = "%s/waveSkimDS%d_run%d.root" % (ds.calWaveDir,dsNum,run)
        if not os.path.isfile(fileName):
            print("file not found! name:", fileName)
            continue
        f2 = TFile(fileName)
        t2 = f2.Get("skimTree")
        n2 = t2.GetEntries()
        print("DS %d  run %d  wave entries %d" % (dsNum, run, n2))
        if n2==0:
            print("no waveskim entries found! file:", fileName)
            continue

if __name__=="__main__":
    main()