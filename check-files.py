#!/usr/bin/env python
import sys, os
sys.argv.append("-b")
import tinydb as db

import dsi
bkg = dsi.BkgInfo()
cal = dsi.CalInfo()

def main(argv):
    """ NOTE: don't use globs when getting files.
    Manually make sure everything is here.
    """
    global checkCal
    checkCal = False

    for i,opt in enumerate(argv):
        if opt == "-c": checkCal = True
        if opt == "-s": checkSkim()
        if opt == "-w": checkWave()
        if opt == "-p": checkSplit()
        if opt == "-l": checkLAT()
        if opt == "-all":
            checkSkim()
            checkWave()
            checkSplit()
            checkLAT()

def checkSkim():

    # make file list
    fileList = []

    # bkg
    if not checkCal:
        dsMap = bkg.dsMap()
        for ds in dsMap:
            for sub in range(dsMap[ds]+1):
                skimFile = "%s/skimDS%d_%d_low.root" % (dsi.skimDir, ds, sub)
                if os.path.isfile(skimFile):
                    fileList.append(skimFile)
                else:
                    print("File not found:",skimFile)
                    continue
    # cal
    else:
        for key in cal.GetKeys():
            ds = int(key[2])
            for cIdx in range(cal.GetIdxs(key)): # 0-indexed
                cRuns = cal.GetCalList(key, cIdx)
                for run in cRuns:
                    skimFile = "%s/skimDS%d_run%d_low.root" % (dsi.calSkimDir, ds, run)
                    if os.path.isfile(skimFile):
                        fileList.append(skimFile)
                    else:
                        print("File not found:",skimFile)
                        continue

    for skimFile in fileList:
        print(skimFile)



def checkWave():
    print("hi")


def checkSplit():
    print("hi")


def checkLAT():
    print("hi")


if __name__=="__main__":
    main(sys.argv[1:])