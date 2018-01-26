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

Search log files for common errors.
"""

def main():
    checkFiles()
    # checkLogErrors()
    # checkLogErrors2()

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


def checkLogErrors():
    """ ./job-panda.py -checkLogs
        This isn't really complete. but you get the idea.  Error checking via bash inside python is kind of a PITA.
        Maybe it would be better to just have python look at the files directly.
    """

    # Shell commands
    c1 = "ls -F ./logs/ | grep -v / | wc -l" # count total files
    c2 = "grep -rIl \"Done! Job Panda\" ./logs/ | wc -l" # count completed files
    c3 = "grep -rL \"Done! Job Panda\" ./logs/"  # negative file matching, can also count # fails.  gives a file list
    c4 = "grep -rIl \"Segmentation\" ./logs/" # segfaults
    c5 = "grep -rIl \"bad_alloc\" ./logs/" # memory errors

    # using sp to deal with a pipe is kind of annoying
    p1 = sp.Popen('ls -F ./logs/'.split(), stdout=sp.PIPE)
    p2 = sp.Popen('grep -v /'.split(), stdin=p1.stdout, stdout=sp.PIPE)
    p3 = sp.Popen('wc -l'.split(), stdin=p2.stdout,stdout=sp.PIPE)
    output = p3.communicate()[0]
    num = int(output.strip('\n'))
    print(num)

    # make a dummy bash script that runs all the shell commands.  who knows if this is smart or not
    outFile = open('logCheck.sh','w+')
    dummyScript = "#!/bin/bash \n %s \n %s \n %s \n %s \n %s \n" % (c1,c2,c3,c4,c5)
    outFile.write(dummyScript)
    outFile.close()
    sh('chmod a+x logCheck.sh')
    sh('./logCheck.sh')
    os.remove('logCheck.sh')


def checkLogErrors2():
    """ Usage: ./job-panda -checkLogs2
        Globs together log files and then searches for "Error", returning the failed ROOT files.
    """
    print("Checking log errors ...")

    ErrList = []
    for fl in glob.glob("./logs/*"):
        fErr = open(fl,'r').read()
        if 'Error' in open(fl, 'r').read():
            print(ErrList.append(fl))

    for errFile in ErrList:
        fErr = open(errFile,'r')
        for lineErr in fErr:
            if '/lat.py -b' in lineErr:
                print('Error from: ', lineErr)
            if 'Error' in lineErr:
                print(lineErr)



if __name__=="__main__":
    main()