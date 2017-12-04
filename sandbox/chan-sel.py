#!/usr/bin/env python
import sys, imp
sys.argv.append("-b")
ds = imp.load_source('DataSetInfo','../DataSetInfo.py')
wl = imp.load_source('waveLibs','../waveLibs.py')
import tinydb as db
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
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


def getErrorCode(arr):
    binString = ''.join(['1' if x else '0' for x in arr])
    code = int(binString, 2)
    return code


def unpackErrorCode(code, padLength=3):
    tmp = "{:b}".format(code)
    while len(tmp) != padLength:
        tmp = '0'+tmp
    arr = [bool(int(i)) for i in tmp]
    return arr


def channelSelection():

    calDB = db.TinyDB('../calDB.json')
    pars = db.Query()

    np.set_printoptions(threshold=np.inf) # print full numpy array

    for dSet in [(0,1),(1,1),(2,1),(3,1),(4,2),(5,1),(5,2)]:
    # for dSet in [(4,2)]:
        dsNum, modNum = dSet[0], dSet[1]
        print "DS-%d, M-%d" % (dsNum, modNum)

        chList = ds.GetGoodChanList(dsNum)
        if dsNum==5 and modNum==1: chList = [ch for ch in chList if ch < 1000 and ch!=692]
        if dsNum==5 and modNum==2: chList = [ch for ch in chList if ch > 1000 and ch!=1232]

        chStats = {ch:[] for ch in chList}
        xTicks = []

        nBkg = ds.dsMap[dsNum]
        bkgRuns = ds.bkgRunsDS[dsNum]
        nCal = wl.getNCalIdxs(dsNum, modNum)
        calInfo = ds.CalInfo()

        for bkgIdx in range(nBkg+1):

            runLo = bkgRuns[bkgIdx][0]
            runHi = bkgRuns[bkgIdx][-1]

            calIdxLo = calInfo.GetCalIdx("ds%d_m%d" % (dsNum, modNum), runLo)
            calIdxHi = calInfo.GetCalIdx("ds%d_m%d" % (dsNum, modNum), runHi)

            # if calIdxHi != calIdxLo:
                # print "DS-%d  bkgIdx %d  Multiple calIdxs: %d - %d" % (dsNum, bkgIdx, calIdxLo, calIdxHi)

            for calIdx in range(calIdxLo, calIdxHi+1):

                xTicks.append(bkgIdx + 0.1*(calIdx - calIdxLo))

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

                    goodList = [goodChanFS, goodChanRN, goodChanWF]
                    code = getErrorCode(goodList)
                    chStats[ch].append(code)


        data = np.asarray([chStats[ch] for ch in chStats])
        xTicks = [ int(t) if float(t).is_integer() else t for t in xTicks ]
        # yTicks = [ch for ch in chStats]
        yTicks = ["P%sD%s" % (str(ds.CPD[dsNum][ch])[1], str(ds.CPD[dsNum][ch])[2]) for ch in chStats]

        # print "x:",xTicks
        # print "y:",yTicks
        # print "data:"
        # print data

        a, b = int(len(xTicks)/3.), int(len(yTicks)/2.)
        if a < 12: a = 12
        if b < 8:  b = 8
        print "a",a,"b",b

        fig = plt.figure(figsize=(a,b), facecolor='w')
        ax = plt.subplot(111)

        im = plt.imshow(np.asarray(data),interpolation='nearest')
        ax.set_xlabel('bkgIdx',x=1., ha='right')
        ax.set_ylabel('Detector',rotation=0)
        ax.yaxis.set_label_coords(-0.01,1.01)
        ax.set_title("DS-%d, Module %d" % (dsNum, modNum))

        # Set the major ticks at the centers and minor tick at the edges
        xlocs = np.arange(0, len(xTicks), 1)
        ax.xaxis.set_ticks(xlocs + 0.5, minor=True)
        ax.xaxis.set(ticks=xlocs, ticklabels=xTicks)
        ax.tick_params(axis='x',labelsize=8)
        ylocs = np.arange(0, len(yTicks), 1)
        ax.yaxis.set_ticks(ylocs + 0.5, minor=True)
        ax.yaxis.set(ticks=ylocs, ticklabels=yTicks)
        ax.grid(True, which='minor')

        values = np.unique(data.ravel())
        colors = [im.cmap(im.norm(value)) for value in values]

        patches = []
        for idx, val in enumerate(values):
            tmp = unpackErrorCode(val)
            lab = "FS %d RN %d WF %d" % (int(tmp[0]), int(tmp[1]), int(tmp[2]))
            patches.append(mpatches.Patch(color=colors[idx], label=lab))

        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.85, box.height])
        ax.legend(handles=patches, bbox_to_anchor=(1, 1), loc=2, borderaxespad=0. )

        plt.savefig("../plots/dc-results-DS%d-M%d.pdf" % (dsNum,modNum))



if __name__=="__main__":
    main()