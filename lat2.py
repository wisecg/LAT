#!/usr/bin/env python3
"""
============== lat2.py: (L)ANL (A)nalysis (T)oolkit ==============

Tunes cut values for LAT files.
Usage: ./lat2.py [options]

================ C. Wiseman (USC), B. Zhu (LANL) ================
"""
import sys
import tinydb as db

import dsi
bkg = dsi.BkgInfo()
cal = dsi.CalInfo()

skipDS6Cal = True # ignore DS6 cal runs until they're processed

def main(argv):

    getCalStats()


def getCalStats():
    """ Load HV changes to see which cal ranges we can chain together.
        Get an idea of statistics in each file.
        Need 100 counts / keV / bin - gives sqrt(100) error.

        If this counts as channel selection, should i move it to chan-sel.py ?
    """
    ds = "5A"

    # load dataset background run list
    # set up a list of HV-change run boundaries

    runList = bkg.getRunList(ds)
    print("DS:",ds,"  first run: %d  last run: %d" % (runList[0], runList[-1]))

    subRanges = bkg.GetSubRanges(ds,opt='hv')
    for sub in subRanges:
        print(sub)

    hvBounds = []
    if len(subRanges) == 0:
        print("No HV changes in this DS.")
        hvBounds.append((runList[0], runList[-1]))
    else:
        runLo = runList[0]

        # damn.  i really need to know the exact run number
        # the hv changed, not just as it applies to the bkg idx.
        # ---> chan-sel.py

        # for sub in subRanges:
        #     runChg = sub[2]
        #
        #     runBeforeChg = -1
        #     for run in runList:
        #         if run < runChg:
        #             runBeforeChg = run
        #     if runBeforeChg == -1:
        #         print("No runs before change")
        #     runHi = runBeforeChg
        #
        #     hvBounds.append((runLo, runHi))
        #     runLo = runChg

    # for runLo, runHi in hvBounds:
        # print(runLo, runHi)

    # load calibration sub-ranges
    # dsList = [0,1,2,3,4,"5A","5B","5C"]
    # for ds in dsList:
    #     subr = bkg.GetSubRanges(ds,opt='hv')
    #     for s in subr:
    #         print(s)
        # [1 51 14342 14372 31 56]
        # [1 51 14386 14387 2 56] <---- huh, probably just remove this one?
        # [4 1 60000970 60000991 22 0]
        # [4 1 60000992 60001000 9 0]
        # ['5A' 16 19733 19747 15 (4, 3, True, True)]
        # ['5A' 16 19771 19773 3 (4, 3, True, True)]
        # ['5A' 61 21565 21567 3 (9, 8, True, False)]
        # ['5A' 61 21568 21587 19 (9, 8, True, False)]
        # ['5A' 78 22248 22304 19 (12, 10, True, False)]
        # ['5A' 78 22316 22333 18 (12, 10, True, False)]
        # ['5B' 91 22944 22946 3 (14, 12, True, False)]
        # ['5B' 91 22952 22986 29 (14, 12, True, False)]


    # load calibration files
    # fileList = []
    # dsMap = cal.GetKeys()
    # for key in dsMap:
    #     ds = int(key[2])
    #     if skipDS6Cal and ds==6: continue
    #     for sub in range(cal.GetIdxs(key)):
    #         cRuns = cal.GetCalList(key, sub)
    #         covLo, covHi = GetCalRunCoverage(self,key,idx)
    #         for run in cRuns:
    #             latList = dsi.getSplitList("%s/latSkimDS%d_run%d*" % (dsi.calLatDir, ds, run), run)
    #             tmpList = [f for idx, f in sorted(latList.items())]
    #             fileList.extend(tmpList)






if __name__=="__main__":
    main(sys.argv[1:])