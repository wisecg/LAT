#!/usr/common/usg/software/python/2.7.9/bin/python
#!/usr/local/bin/python
#!/bin/sh
""":" `which python` "$0" "$@"; exit 1 ":""" # dammit, this works on PDSF but not my laptop

import os
from ROOT import TFile, TTree, TEntryList, gDirectory, TNamed, std, TObject

def main():

    # splitTree()
    # printFileSizes()
    # printCalFileSizes()


def printCalFileSizes():

    calMap = {
    0:[2931,2940, 6854,6863],
    1:[9497,9503, 14149,14155],
    2:[14568,14574, 15789,15794],
    3:[16911,16920, 17950,17959],
    4:[60001014,60001023, 60001855,60001864],
    5:[19055,19064, 22513,22523]
    }

    for dsNum in calMap:
        totalSize = 0
        runList = []

        for i in xrange(0,len(calMap[dsNum]),2):
            lower = calMap[dsNum][i]
            upper = calMap[dsNum][i+1]
            for j in xrange(lower,upper+1):
                runList.append(j)

        for run in runList:

            # print dsNum, run

            inPath = "/global/homes/w/wisecg/project/v1-cal-waveskim/waveSkimDS%d_run%d.root" % (dsNum, run)

            inFile = TFile()
            try:
                infile = TFile(inPath)
                inFile = TFile(inPath)
                gatTree = inFile.Get("skimTree")
                nEnt = gatTree.GetEntriesFast()

                fileSize = os.path.getsize(inPath)/1e6 # mb
                totalSize += fileSize

                print "DS%d, run %d: entries %d  size %.0f" % (dsNum,run,nEnt,fileSize)
            except:
                AttributeError
                print "file ain't exist"

        print "Total Size: %.0f\n" % totalSize


def printFileSizes():
    dsMap = {0:76, 1:51, 2:7, 3:24, 4:22, 5:112}

    for dsNum, subDS in dsMap.iteritems():

        totalSize = 0

        for subNum in range(subDS+1):

            waveSkimPath = "/global/homes/w/wisecg/project/waveskim/waveSkimDS%d_%d.root" % (dsNum, subNum)
            skimPath = "/global/homes/w/wisecg/project/skim/skimDS%d_%d_low.root" % (dsNum, subNum)

            wFile = TFile(waveSkimPath)
            theCut = wFile.Get("theCut").GetTitle()
            wTree = wFile.Get("skimTree")
            wEnt = wTree.GetEntriesFast()

            sFile = TFile(skimPath)
            sTree = sFile.Get("skimTree")
            # sTree.Draw(">>elist",theCut,"entrylist")
            # elist = gDirectory.Get("elist")
            # sTree.SetEntryList(elist)
            # nList = elist.GetN()
            sEnt = sTree.GetEntriesFast()

            sSize = os.path.getsize(skimPath)/1e6 # mb
            wSize = os.path.getsize(waveSkimPath)/1e6

            if wEnt > 100000:
                print "DS%d-%d:  sEnt %-8d  wEnt %-8d  s/w %-4.1f  sSize %-4d  wSize %-4d" % (dsNum,subNum,sEnt,wEnt,sEnt/wEnt,sSize,wSize)

            # totalSize += fileSize

            # if nEnt > 100000:
            # print "DS%d-%d: entries %d  size %.0f" % (dsNum,subNum,nEnt,fileSize)
                # print "%d %.1f" % (nEnt,fileSize)
                # print "%.1f," % fileSize

        # print "Total Size: %.0f\n" % totalSize


def splitTree():

    dsNum, subNum = 0, 25

    inPath = "/global/homes/w/wisecg/project/waveskim/waveSkimDS%d_%d.root" % (dsNum, subNum)
    outPath = "./cutSkimDS%d_%d.root" % (dsNum, subNum)

    inFile = TFile(inPath)
    bigTree = inFile.Get("skimTree")

    theCut = inFile.Get("theCut").GetTitle()

    bigTree.Draw(">>elist",theCut,"entrylist")
    elist = gDirectory.Get("elist")
    bigTree.SetEntryList(elist)
    nList = elist.GetN()

    outFile = TFile(outPath,"RECREATE")

    lilTree = TTree()
    lilTree.SetMaxTreeSize(150000000)
    lilTree = bigTree.CopyTree("")

    lilTree.Write("",TObject.kOverwrite)

    thisCut = TNamed("theCut",theCut)
    thisCut.Write()


if __name__ == "__main__":
    main()

