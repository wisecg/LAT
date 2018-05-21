#!/usr/bin/env python3
import sys, os
import numpy as np

import matplotlib as mpl
mpl.use('Agg')
sys.argv.append("-b")
import matplotlib.pyplot as plt
plt.style.use('../pltReports.mplstyle')
from matplotlib.colors import LogNorm, Normalize
from matplotlib import gridspec

import waveLibs as wl
import dsi
det = dsi.DetInfo()
bkg = dsi.BkgInfo()

def main(argv):

    # p1()
    # p2()
    p3()


def loadCutFiles(dsList,cutType):
    from ROOT import TChain
    fList = []
    for ds in dsList:
        print("Loading DS%s cut files ..." % str(ds))
        dsNum = int(ds[0]) if isinstance(ds, str) else int(ds)
        chList = det.getGoodChanList(dsNum)
        bkgRanges = bkg.getRanges(ds)
        for bIdx in bkgRanges:
            for ch in chList:
                cutFile = "%s/%s/%s_ds%d_%d_ch%d.root" % (dsi.cutDir,cutType,cutType,dsNum,bIdx,ch)
                if not os.path.isfile(cutFile):
                    # print("Couldn't find file:",cutFile)
                    continue
                fList.append(cutFile)
    tt = TChain("skimTree")
    for f in fList: tt.Add(f)
    return tt


def p1():
    tfr = loadCutFiles([1,2,3,4,"5B","5C"],"fr")

    n = tfr.Draw("trapENFCal:fitSlo:riseNoise","trapENFCal < 50","goff")
    hitE, fSlo, rise = tfr.GetV1(), tfr.GetV2(), tfr.GetV3()
    hitE = np.asarray([hitE[i] for i in range(n)])
    fSlo = np.asarray([fSlo[i] for i in range(n)])
    rise = np.asarray([rise[i] for i in range(n)])

    fig = plt.figure()

    plt.cla()
    plt.plot(hitE, fSlo, ".b", ms=2)
    plt.xlabel("Energy (keV)", ha='right', x=1)
    plt.ylabel("fitSlo", ha='right',y=1)
    plt.savefig("../plots/p1-fSlo.png")

    plt.cla()
    plt.plot(hitE, rise, ".b", ms=2)
    plt.xlabel("Energy (keV)", ha='right', x=1)
    plt.ylabel("riseNoise", ha='right',y=1)
    plt.savefig("../plots/p1-rise.png")

    plt.cla()
    x, y = wl.GetHisto(hitE, 0, 50, 0.1)
    plt.plot(x, y, ls='steps', lw=1., c='b')
    plt.xlabel("Energy (keV)", ha='right', x=1)
    plt.ylabel("Counts", ha='right',y=1)
    plt.savefig("../plots/p1-espec.png")


def p2():

    tfs = loadCutFiles([1,2,3,4,"5B","5C"],"fs")
    n = tfs.Draw("trapENFCal:fitSlo:riseNoise","trapENFCal < 50","goff")
    hitE_fs, fSlo_fs, rise_fs = tfs.GetV1(), tfs.GetV2(), tfs.GetV3()
    hitE_fs = np.asarray([hitE_fs[i] for i in range(n)])
    fSlo_fs = np.asarray([fSlo_fs[i] for i in range(n)])
    rise_fs = np.asarray([rise_fs[i] for i in range(n)])

    tfr = loadCutFiles([1,2,3,4,"5B","5C"],"fr")
    n = tfr.Draw("trapENFCal:fitSlo:riseNoise","trapENFCal < 50","goff")
    hitE_fr, fSlo_fr, rise_fr = tfr.GetV1(), tfr.GetV2(), tfr.GetV3()
    hitE_fr = np.asarray([hitE_fr[i] for i in range(n)])
    fSlo_fr = np.asarray([fSlo_fr[i] for i in range(n)])
    rise_fr = np.asarray([rise_fr[i] for i in range(n)])

    fig = plt.figure()

    plt.cla()
    plt.plot(hitE_fs, fSlo_fs, ".r", ms=3, label='fs')
    plt.plot(hitE_fr, fSlo_fr, ".b", ms=2, label='fr')
    plt.legend()
    plt.xlabel("Energy (keV)", ha='right', x=1)
    plt.ylabel("fitSlo", ha='right',y=1)
    plt.savefig("../plots/p1-fSlo.png")

    plt.cla()
    plt.plot(hitE_fs, rise_fs, ".r", ms=3, label='fs')
    plt.plot(hitE_fr, rise_fr, ".b", ms=2, label='fr')
    plt.legend()
    plt.xlabel("Energy (keV)", ha='right', x=1)
    plt.ylabel("fitSlo", ha='right',y=1)
    plt.savefig("../plots/p1-rise.png")

    plt.cla()
    plt.plot(*wl.GetHisto(hitE_fs, 0, 50, 0.1), ls='steps', lw=1., c='r', label='fs')
    plt.plot(*wl.GetHisto(hitE_fr, 0, 50, 0.1), ls='steps', lw=1., c='b', label='fr')
    plt.legend()
    plt.xlabel("Energy (keV)", ha='right', x=1)
    plt.ylabel("Counts", ha='right',y=1)
    plt.savefig("../plots/p1-espec.png")


def p3():

    ds = "5B"

    tfr = loadCutFiles([ds],"fr")
    n = tfr.Draw("trapENFCal:run:channel", "trapENFCal < 50", "goff")
    hitE, run, chan = tfr.GetV1(), tfr.GetV2(), tfr.GetV3()
    hitE = np.asarray([hitE[i] for i in range(n)])
    runs = np.asarray([run[i] for i in range(n)])
    chan = np.asarray([chan[i] for i in range(n)])

    # rate under 10 kev, by channel

    # for every run
        # for every channel
            # count the number of hits under 10 kev

    runList = list(sorted(set(runs)))
    chMap = list(sorted(set(chan)))
    for run in runList:
        for ch in chMap:
            print(ch)
        return


    runLo, runHi = np.amin(run), np.amax(run)
    nRuns = runHi-runLo+1

    chDict = {chMap[i]:i for i in range(len(chMap))}
    chan = [chDict[chan] for chan in chan]
    yLo, yHi = 0, len(chMap)
    nbx, nby = int((xHi-xLo)/xpb), len(chMap)


    # plt.hist2d(run)


    # this is really too fine grained
    # runLo, runHi = np.amin(run), np.amax(run)
    # nRuns = runHi-runLo+1
    # eLo, eHi, epb = 0, 50, 0.5
    # nbe = int((eHi-eLo)/epb)
    # plt.hist2d(run, hitE, bins=[nRuns,nbe], range=[[runLo,runHi],[eLo,eHi]], norm=LogNorm(), cmap='jet')
    # plt.hist2d(run, hitE, bins=[nRuns,nbe], range=[[runLo,runHi],[eLo,eHi]], cmap='jet')
    # plt.colorbar()
    # plt.xlabel("Run number", ha='right', x=1)
    # plt.ylabel("Energy (keV)", ha='right', y=1)
    # plt.savefig("../plots/p1-runvsch.png")

    # goodChans = det.getGoodChanList(ds)

    # this doesn't show bad runs
    # fig = plt.figure()
    # p1 = plt.subplot(111)
    # chMap = list(sorted(set(chan)))
    # chDict = {chMap[i]:i for i in range(len(chMap))}
    # chan = [chDict[chan] for chan in chan]
    # xLo, xHi, xpb = 0, 50, 0.5
    # yLo, yHi = 0, len(chMap)
    # nbx, nby = int((xHi-xLo)/xpb), len(chMap)
    # _,_,_,im1 = p1.hist2d(hitE, chan, bins=[nbx,nby], range=[[xLo,xHi],[yLo,yHi]], cmap='jet')
    # p1.set_xlabel("Energy (keV)", ha='right', x=1.)
    # # xticklabels = ["%.0f" % t for t in np.arange(xLo, xHi, xpb)]
    # # p1.set_xticklabels(xticklabels)
    # yticks = np.arange(0, len(chMap))
    # p1.draw(fig.canvas.get_renderer()) # workaround to print tick labels
    # print([tl.get_text() for tl in p1.get_xticklabels()])
    # p1.set_ylabel("channel", ha='right', y=1.)
    # p1.set_yticks(yticks)
    # p1.set_yticklabels(chMap, fontsize=12)
    # fig.colorbar(im1, ax=p1, fraction=0.037, pad=0.04)
    # plt.savefig("../plots/p1-evsch.png")





if __name__=="__main__":
    main(sys.argv[1:])