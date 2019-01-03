#!/usr/bin/env python3
import numpy as np
import tinydb as db

import matplotlib.pyplot as plt
plt.style.use('../clint.mpl')

import dsi
bkg = dsi.BkgInfo()
cal = dsi.CalInfo()

def main():

    makeSubDSPlot()


def makeSubDSPlot():
    """
    lat-expo uses a 5-layer loop to subdivide data runs.
    Let's try to graphically illustrate that.
    """

    ds = "5A"

    rBkg, rCalM1, rCalM2 = [], [], []

    # set DS stuff
    dsNum = int(ds[0]) if isinstance(ds,str) else ds
    nBkg = bkg.dsMap()[dsNum]
    bLo, bHi = 0, nBkg
    if ds=="5A": bLo, bHi = 0, 79
    if ds=="5B": bLo, bHi = 80, 112
    if ds=="5C": bLo, bHi = 113, 121
    bkgRanges = bkg.getRanges(ds)

    calDB = db.TinyDB('%s/calDB-v2.json' % (dsi.latSWDir))
    pars = db.Query()

    # 1. loop over modules
    mods = [1]
    if dsNum == 4: mods = [2]
    if dsNum == 5: mods = [1,2]
    for mod in mods:

        calKey = "ds%d_m%d" % (dsNum, mod)
        if ds == "5C": calKey = "ds5c"
        if calKey not in cal.GetKeys(dsNum):
            print("Error: Unknown cal key:",calKey)
            return

        rFirst, rLast = bkgRanges[bLo][0], bkgRanges[bHi][-1]

        nCal = cal.GetNCalIdxs(dsNum, mod)

        for i in range(nCal):
            if mod==1: rCalM1.append(cal.GetCalList(calKey,i)[0])
            if mod==2: rCalM2.append(cal.GetCalList(calKey,i)[0])

        # 2. loop over bkgIdx
        for bIdx in range(bLo, bHi+1):

            rLo, rHi = bkgRanges[bIdx][0], bkgRanges[bIdx][-1]
            rBkg.append(rLo)

            # # 3. loop over sub-bIdx
            subRanges = bkg.GetSubRanges(ds, bIdx)
            if len(subRanges) == 0: subRanges.append((rLo, rHi))
            for sbIdx, (subLo, subHi) in enumerate(subRanges):

                # 4. loop over cIdx's in this sub-bIdx
                cIdxLo, cIdxHi = cal.GetCalIdx(calKey, subLo), cal.GetCalIdx(calKey, subHi)
                for i, cIdx in enumerate(range(cIdxLo, cIdxHi+1)):

                    # get the run coverage of this sub-sub-bIdx
                    if cIdxLo==cIdxHi:
                        covLo, covHi = subLo, subHi
                    else:
                        runList = bkg.getRunList(ds, bIdx)  # change 1
                        subList = [r for r in runList if subLo <= r <= subHi and cal.GetCalIdx(calKey,r) == cIdx]
                        if len(subList)==0: continue
                        covLo, covHi = subList[0], subList[-1]

                    # if bIdx==78:
                        # print(ds, smod, bIdx, "bkg", sbIdx, subLo, subHi, "cal", cIdx, covLo, covHi)


    fig = plt.figure()
    p0 = plt.subplot(211)
    p0.axes.get_xaxis().set_visible(False)
    p0.axes.get_yaxis().set_visible(False)
    p0.set_frame_on(False)

    p1 = plt.subplot(212)
    p1.axes.get_xaxis().set_visible(False)
    p1.axes.get_yaxis().set_visible(False)
    p1.set_frame_on(False)

    # 1. whole-ds plot, with bkgIdx positions

    p0.plot((rFirst, rLast), (1, 1), '-k')

    p0.plot((rFirst, rFirst), (0.9, 1.1), '-k')
    p0.plot((rLast, rLast), (0.9, 1.1), '-k')

    for r in rBkg:
        p0.plot((r,r),(0.95,1.05),'-r', lw=1)

    p0.set_ylim(0.8, 1.2)
    p0.set_xlim(rFirst-100, rLast+100)

    # p0.annotate("%d" % rFirst, xy=(rFirst, 0.9), xytext=(rFirst-200, 0.85), fontsize=15)
    # p0.annotate("%d" % rLast, xy=(rLast, 0.9), xytext=(rLast-200, 0.85), fontsize=15)
    # p0.annotate("DS-%s" % ds, xy=(rLast, 0.9), xytext=(rLast-200, 0.85), fontsize=15)
    p0.annotate("Bkg    Cal", xy=(rLast, 0.9), xytext=(rLast-200, 0.85), fontsize=15)

    # p0.annotate('bIdx %d' % 78, xy=(rBkg[78],0.95), xytext=(rBkg[78]-800, 0.85), fontsize=15, \
        # arrowprops=dict(arrowstyle="->"))


    # 2. whole-ds plot, with calIdx positions

    p1.plot((rFirst, rLast), (1, 1), '-k')

    p1.plot((rFirst, rFirst), (0.9, 1.1), '-k')
    p1.plot((rLast, rLast), (0.9, 1.1), '-k')

    for r in rCalM1:
        p1.plot((r,r),(0.95,1.05),'-g', lw=2)

    for r in rCalM2:
        p1.plot((r,r),(0.95,1.05),'-m', lw=2)

    p1.set_ylim(0.8, 1.2)
    p1.set_xlim(rFirst-100, rLast+100)

    # plt.show()
    plt.savefig("../plots/bkgCalIdxs2.pdf")




    return

    # 2. bIdx 78, which has multiple sbIdx and scIdx's
    p1.set_ylim(0.6, 1.35)

    # Module 1                         |  these  | are the fine-grained run limits
    # 5A 1 78 bkg 0 22248 22250 cal 11 22248 22250
    # 5A 1 78 bkg 1 22266 22304 cal 11 22266 22280 * note: no data between 22280 - 22304!
    # 5A 1 78 bkg 1 22266 22304 cal 12 22304 22304
    # 5A 1 78 bkg 2 22316 22333 cal 12 22316 22333
    # # Module 2
    # 5A 2 78 bkg 0 22248 22250 cal 9 22248 22250
    # 5A 2 78 bkg 1 22266 22304 cal 9 22266 22280
    # 5A 2 78 bkg 1 22266 22304 cal 10 22304 22304
    # 5A 2 78 bkg 2 22316 22333 cal 10 22316 22333

    # whole bkgIdx (vert)
    p1.plot((22248, 22248), (0.75, 1.25), '-k')
    p1.plot((22333, 22333), (0.75, 1.25), '-k')
    # p1.annotate("bIdx 78", xy=(22248, 0.9), xytext=(22248-5, 1.25), fontsize=15)

    # m1/m2 bkg run ranges (horiz)
    p1.plot((22248, 22250), (1.1, 1.1), '-k')
    p1.plot((22266, 22280), (1.1, 1.1), '-k')
    p1.plot((22304-0.5, 22304), (1.1, 1.1), '-k')
    p1.plot((22316, 22333), (1.1, 1.1), '-k')
    p1.plot((22248, 22250), (0.9, 0.9), '-k')
    p1.plot((22266, 22280), (0.9, 0.9), '-k')
    p1.plot((22304-0.5, 22304), (0.9, 0.9), '-k')
    p1.plot((22316, 22333), (0.9, 0.9), '-k')
    p1.annotate("M1", xy=(22248, 1.1), xytext=(22248-6, 1.08), fontsize=15)
    p1.annotate("M2", xy=(22248, 0.9), xytext=(22248-6, 0.88), fontsize=15)

    # m1/m2 bkg run boundaries (vert)
    p1.plot((22250, 22250), (0.75, 1.25), '-r', alpha=0.6)
    p1.plot((22266, 22266), (0.75, 1.25), '-r', alpha=0.6)
    p1.plot((22304, 22304), (0.75, 1.25), '-r', alpha=0.6)
    p1.plot((22316, 22316), (0.75, 1.25), '-r', alpha=0.6)

    p1.annotate("cIdx", xy=(22248, 0.9), xytext=(22248-8, 1.3), fontsize=15)
    p1.annotate("11", xy=(22280-3, 1.32), xytext=(22248-1, 1.3), fontsize=15, arrowprops=dict(arrowstyle="->"))
    p1.annotate("12", xy=(22333, 1.32), xytext=(22282-1, 1.3), fontsize=15, arrowprops=dict(arrowstyle="->"))

    p1.annotate("cIdx", xy=(22248, 0.9), xytext=(22248-8, 0.65), fontsize=15)
    p1.annotate("9", xy=(22280-3, 0.672), xytext=(22248-1, 0.65), fontsize=15, arrowprops=dict(arrowstyle="->"))
    p1.annotate("10", xy=(22333, 0.672), xytext=(22282-1, 0.65), fontsize=15, arrowprops=dict(arrowstyle="->"))

    # m1 cal coverage (vert)
    p1.plot((22248, 22248), (1.02, 1.2), '-g')
    p1.plot((22280, 22280), (1.02, 1.2), '-g')
    p1.plot((22316, 22316), (1.02, 1.2), '-g')
    p1.plot((22333, 22333), (1.02, 1.2), '-g')

    # m2 cal coverage (vert)
    p1.plot((22250, 22250), (0.8, 0.98), '-m')
    p1.plot((22266, 22266), (0.8, 0.98), '-m')
    p1.plot((22280, 22280), (0.8, 0.98), '-m')
    p1.plot((22316, 22316), (0.8, 0.98), '-m')


    # p1.legend(loc=1)
    plt.show()


if __name__=="__main__":
    main()