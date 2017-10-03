#!/usr/bin/env python
from ROOT import TFile, TTree
import os.path, imp
dsi = imp.load_source('DataSetInfo', '../DataSetInfo.py')

def main():
    """ Two ideas here:
    1) Scan over all lat-related files, both BG and Cal:
        - skims
        - wave-skims
        - split-skims
        - split-lats
        - merged lats
    and make sure they all exist, all have tree entries, etc.

    2) Scan over calibration run lists from calDB, and see how many
    entries we have per channel.  Could help to find the sweet spot
    of how many cal runs we need to process per range.

    also, should it check if we have events down to 0.7 kev or whatever?

    also, this should be moved to job-panda once it's working well.
    """
    # checkBkgFiles()
    checkCalFiles()


def checkBkgFiles():

    dsMap = {0:75,1:51,2:7,3:24,4:18,5:112}
    for ds in dsMap:
        for sub in range(dsMap[ds]+1):

            # check skims
            fileName = "/global/homes/w/wisecg/project/bg-skim/skimDS%d_%d_low.root" % (ds,sub)
            if not os.path.isfile(fileName):
                print "file not found! name:", fileName
                continue
            f1 = TFile(fileName)
            t1 = f1.Get("skimTree")
            n1 = t1.GetEntriesFast()
            print "DS %d  sub %d  skim entries %d" % (ds, sub, n1)
            if n1==0:
                print "no skim entries found! file:", fileName
                continue

            # check waveskims
            fileName = "/global/homes/w/wisecg/project/bg-waves/waveSkimDS%d_%d.root" % (ds,sub)
            if not os.path.isfile(fileName):
                print "file not found! name:", fileName
                continue
            f2 = TFile(fileName)
            t2 = f2.Get("skimTree")
            n2 = t2.GetEntriesFast()
            print "DS %d  sub %d  wave entries %d" % (ds, sub, n2)
            if n2==0:
                print "no waveskim entries found! file:", fileName
                continue


def getCalRunList(dsNum=None,subNum=None):
    """ ./job-panda.py -cal (-ds [dsNum] -sub [dsNum] [calIdx])
        Create a calibration run list, using the CalInfo object in DataSetInfo.py .
        Note that the -sub option is re-defined here to mean a calibration range idx.
        Note that running with -cal alone will create a list for all datasets (mega mode).
    """
    runLimit = 10 # yeah I'm hardcoding this, sue me.
    calList = []
    calInfo = dsi.CalInfo()
    calKeys = calInfo.GetKeys(dsNum)

    for key in calKeys:
        print "key:",key

        # -cal (mega mode)
        if dsNum==None:
            for idx in range(calInfo.GetIdxs(key)):
                lst = calInfo.GetCalList(key,idx,runLimit)
                print lst
                calList += lst
        # -ds
        elif subNum==None:
            for idx in range(calInfo.GetIdxs(key)):
                lst = calInfo.GetCalList(key,idx,runLimit)
                print lst
                calList += lst
        # -sub
        else:
            lst = calInfo.GetCalList(key,subNum,runLimit)
            if lst==None: continue
            print lst
            calList += lst

    # remove any duplicates, but there probably aren't any
    calList = sorted(list(set(calList)))

    return calList


def checkCalFiles():

    calList = getCalRunList(dsNum=3)

    for run in calList:
        dsNum=-1
        for key in dsi.dsRanges:
            if dsi.dsRanges[key][0] <= run <= dsi.dsRanges[key][1]:
                dsNum=key

        # check skims
        fileName = "/global/homes/w/wisecg/project/cal-skim/skimDS%d_run%d_low.root" % (dsNum,run)
        if not os.path.isfile(fileName):
            print "file not found! name:", fileName
            continue
        f1 = TFile(fileName)
        t1 = f1.Get("skimTree")
        n1 = t1.GetEntriesFast()
        print "DS %d  run %d  skim entries %d" % (dsNum, run, n1)
        if n1==0:
            print "no skim entries found! file:", fileName
            continue

        # check waveskims
        fileName = "/global/homes/w/wisecg/project/cal-waves/waveSkimDS%d_run%d.root" % (dsNum,run)
        if not os.path.isfile(fileName):
            print "file not found! name:", fileName
            continue
        f2 = TFile(fileName)
        t2 = f2.Get("skimTree")
        n2 = t2.GetEntriesFast()
        print "DS %d  run %d  wave entries %d" % (dsNum, run, n2)
        if n2==0:
            print "no waveskim entries found! file:", fileName
            continue


if __name__ == "__main__":
    main()