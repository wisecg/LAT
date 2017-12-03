#!/usr/bin/env python
import sys, imp
sys.argv.append("-b")
ds = imp.load_source('DataSetInfo','../DataSetInfo.py')
wl = imp.load_source('waveLibs','../waveLibs.py')
import tinydb as db
"""
Remember, ApplyChannelCuts generates a file for each channel in each bkgIdx.

What did we cut out when we ran LAT3::ApplyChannelCuts ?
    - From outright channel selection (2 DS5 beges)
    - From no data in cut tuning
    - Files w/ no statistics after cuts (not actually a cut ...)

If the bkgIdx goes over multiple calIdx's, and we don't have data for all of them,
we don't lose the whole bkgIdx, just the ranges that we don't have valid cuts for.

fitSlo - uses calIdx
    keys: fitSlo_ds[i]_idx[j]_m[k]_Peak
    vals: {[chan]:[1%, 5%, 90%, 95%, 99%]}

riseNoise - uses calIdx
    keys: riseNoise_ds[i]_idx[j]_m[k]_SoftPlus
    vals: {[chan]:[1%, 5%, 90%, 95%, 99%]}

wfStd - uses M1 calIdx (including DS5-M2)
    keys: wfstd_ds[i]_idx[j]_mod[k]
    vals: {[chan]:[y/n, expo, aThresh, a, b, c, d, e, base, n, m, chi2]}

threshold - uses bkgIdx
    keys: thresh_ds[i]_bkgidx[j]
    vals: {[chan]:[50 pct mean, sigma]}

FUTURE: write a routine that identifies all active channels in a dataset
and makes a plot of active channels vs. run number.
So we're super sure that ds.GetGoodChanList is complete.
"""

def main():

    # quickSearch()
    channelSelection()


def quickSearch():
    """ Use a regexp to quickly search the DB ... very handy. """
    calDB = db.TinyDB('../calDB.json')
    pars = db.Query()
    recList = calDB.search(pars.key.matches("wf"))
    print len(recList)
    for idx in range(len(recList)):
        key = recList[idx]['key']
        vals = recList[idx]['vals']
        print key
        for ch in vals:
            print ch, vals[ch]
            return


def channelSelection():

    calDB = db.TinyDB('../calDB.json')
    pars = db.Query()

    for dSet in [(0,1),(1,1),(2,1),(3,1),(4,2),(5,1),(5,2)]:
        dsNum, modNum = dSet[0], dSet[1]
        print "DS-%d, M-%d" % (dsNum, modNum)

        chList = ds.GetGoodChanList(dsNum)
        if dsNum==5 and modNum==1: chList = [ch for ch in chList if ch < 1000 and ch!=692]
        if dsNum==5 and modNum==2: chList = [ch for ch in chList if ch > 1000 and ch!=1232]

        nBkg = ds.dsMap[dsNum]
        bkgRuns = ds.bkgRunsDS[dsNum]
        nCal = wl.getNCalIdxs(dsNum, modNum)
        calInfo = ds.CalInfo()

        for bkgIdx in range(nBkg+1):

            runLo = bkgRuns[bkgIdx][0]
            runHi = bkgRuns[bkgIdx][-1]

            calIdxLo = calInfo.GetCalIdx("ds%d_m%d" % (dsNum, modNum), runLo)
            calIdxHi = calInfo.GetCalIdx("ds%d_m%d" % (dsNum, modNum), runHi)

            for calIdx in range(calIdxLo, calIdxHi+1):

                calRunLo = calInfo.master["ds%d_m%d" % (dsNum, modNum)][calIdx][1]
                calRunHi = calInfo.master["ds%d_m%d" % (dsNum, modNum)][calIdx][2]

                fsD = wl.getDBCalRecord("fitSlo_ds%d_idx%d_m%d_Peak" % (dsNum, calIdx, modNum), False, calDB, pars)
                rnD = wl.getDBCalRecord("riseNoise_ds%d_idx%d_m%d_SoftPlus" % (dsNum, calIdx, modNum), False, calDB, pars)
                wfD = wl.getDBCalRecord("wfstd_ds%d_idx%d_mod%d" % (dsNum, calIdx, modNum), False, calDB, pars)

                goodFS = True if fsD!=0 else False
                goodRN = True if rnD!=0 else False
                goodWF = True if wfD!=0 else False

                print "DS%d  bkgIdx %d  calIdx %d  FS %d  RN %d  WF %d" % (dsNum, bkgIdx, calIdx, int(goodFS), int(goodRN), int(goodWF))

                for ch in chList:

                    goodChanFS = True if goodFS and fsD[ch][2] > 0 else False
                    goodChanRN = True if goodRN and rnD[ch][3] > 0 else False
                    goodChanWF = True if goodWF and ch in wfD.keys() and wfD[ch][0]==u'y' else False

                # now think about how to print out this information.
                # maybe print bad keys and bad channels, and keep a counter?



if __name__=="__main__":
    main()