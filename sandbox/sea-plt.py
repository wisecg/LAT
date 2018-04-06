#!/usr/bin/env python3
import sys, os, imp, glob
import numpy as np
import tinydb as db

import matplotlib as mpl
mpl.use('Agg')
sys.argv.append("-b")
import matplotlib.pyplot as plt
plt.style.use('../pltReports.mplstyle')
from matplotlib.colors import LogNorm, Normalize

dsi = imp.load_source('dsi',os.environ['LATDIR']+'/dsi.py')
wl = imp.load_source('waveLibs',os.environ['LATDIR']+'/waveLibs.py')

def main():

    # scrapeData()
    # plotData()
    plotData2()
    # getData3()
    # plotData3()


def scrapeData():
    from ROOT import TChain

    ds, mod = 1, 1
    print("Loading cut data ...")
    bkg = dsi.BkgInfo()
    cal = dsi.CalInfo()
    nSub = bkg.dsMap()[ds]
    chList = dsi.GetGoodChanList(ds)
    calDB = db.TinyDB('../calDB.json')
    pars = db.Query()

    fsD, rnSD, rnCD = {}, {}, {}
    nCal = cal.GetNCalIdxs(ds,mod)
    for iC in range(nCal+1):
        fsD[iC] = dsi.getDBRecord("fitSlo_ds%d_idx%d_m%d_Peak" % (ds, iC, mod), False, calDB, pars)
        rnSD[iC] = dsi.getDBRecord("riseNoise_ds%d_idx%d_m%d_SoftPlus" % (ds, iC, mod), False, calDB, pars)
        rnCD[iC] = dsi.getDBRecord("riseNoise_ds%d_idx%d_m%d_Continuum" % (ds, iC, mod), False, calDB, pars)
    thD = {}
    for iB in range(nSub+1):
        thD[iB] = dsi.getDBRecord("thresh_ds%d_bkgidx%d" % (ds, iB), False, calDB, pars)

    print("Looping over sub-ranges ...")
    hitE, fSlo, chans, runs, fsCut, thMu, thSig = [], [], [], [], [], [], []
    rise, riseCut = [], []
    for sub in range(nSub+1):
        tc = TChain("skimTree")
        fRaw = sorted(glob.glob("%s/latSkimDS%d_%d_*.root" % (dsi.latDir, ds, sub)))
        for f in fRaw: tc.Add(f)

        print("%d/%d %d" % (sub, nSub, tc.GetEntries()))

        # for some reason range 44 in ds1 is corrupted?
        n = tc.Draw("trapENFCal:fitSlo:channel:run","","goff")
        if n==0:
            print("skipped",sub)
            continue

        # n = tc.Draw("run:channel:trapENFCal:fitSlo","","goff")
        # t1, t2, t3, t4 = tc.GetV1(), tc.GetV2(), tc.GetV3(), tc.GetV4()

        tn = ["run","channel","trapENFCal","fitSlo","riseNoise"]
        vals = wl.GetVX(tc,tn,"")
        t1, t2, t3, t4, t5 = vals["run"], vals["channel"], vals["trapENFCal"], vals["fitSlo"], vals["riseNoise"]
        n = len(t1)

        pRun = -1
        for i in range(n):
            run = int(t1[i])
            if run != pRun:
                cIdx = cal.GetCalIdx("ds%d_m%d" % (ds,mod), run)
                bIdx = bkg.GetBkgIdx(ds, run)
                tmpFS = fsD[cIdx]
                tmpTH = thD[bIdx]
                tmpRNC = rnCD[cIdx]
                tmpRNS = rnSD[cIdx]
                fsCutChan = list(tmpFS.keys())
                thCutChan = list(tmpTH.keys())
                rnCutChan = list(tmpRNC.keys())
            pRun = run

            chan = int(t2[i])
            if chan not in chList: continue
            if chan not in fsCutChan: continue
            if chan not in thCutChan: continue
            if chan not in rnCutChan: continue

            fsVal = tmpFS[chan][2]
            thVal = tmpTH[chan][0] + 3*tmpTH[chan][1]

            a = max(tmpRNS[chan][0], tmpRNC[chan][4])
            b = tmpRNS[chan][1]
            c = tmpRNS[chan][2]
            d = tmpRNS[chan][3]
            if d == 0: continue
            rnVal = a + b * np.log(1 + np.exp(t3[i] - c/d))

            hitE.append(t3[i])
            fSlo.append(t4[i])
            runs.append(run)
            chans.append(chan)
            fsCut.append(tmpFS[chan][2])
            thMu.append(tmpTH[chan][0])
            thSig.append(tmpTH[chan][1])
            rise.append(t5[i])
            riseCut.append(rnVal)

    # # fCut = sorted(glob.glob("%s/results_v1/fs/fitSlo-DS%d-*.root" % (dsi.dataDir, dsNum)))

    print(len(hitE))
    np.savez("../plots/sea-plt-2.npz", hitE, fSlo, chans, runs, fsCut, thMu, thSig, rise, riseCut)


def plotData():

    f = np.load("../plots/sea-plt.npz")

    hitE, fSlo, chans, runs = f['arr_0'], f['arr_1'], f['arr_2'], f['arr_3']
    fsCut, thMu, thSig = f['arr_4'], f['arr_5'], f['arr_6']

    hitPass, fsPass, hitFail, fsFail = [], [], [], []
    for i in range(len(hitE)):
        # if chans[i] in [578, 580, 582, 592, 598, 600, 608, 610, 626, 632, 640, 648, 664, 672, 690, 692]:
        if chans[i] not in [592]:
            # bad? 578
            continue
        if fSlo[i] < fsCut[i] and hitE[i] > (thMu[i]+3*thSig[i]):
            hitPass.append(hitE[i])
            fsPass.append(fSlo[i])
        else:
            hitFail.append(hitE[i])
            fsFail.append(fSlo[i])

    fig = plt.figure()
    xLo, xHi, xpb = 0.5, 20, 0.1

    # x, yAll = wl.GetHisto(hitE, xLo, xHi, xpb)
    # x, yPass = wl.GetHisto(hitPass, xLo, xHi, xpb)
    # plt.semilogy(x, yAll, lw=2, ls='steps', c='k', label='all')
    # plt.semilogy(x, yPass, lw=2, ls='steps', c='b', label='pass')
    # plt.axvline(1.0, color='g', label='1 keV')
    # plt.axvline(1.5, color='m', label='1.5 keV')
    # plt.legend()
    # plt.savefig("../plots/sea-test-2.png")

    # chMap = list(sorted(set(chans)))
    # # print(chMap)
    # x, yCts = wl.GetHisto(chans, chMap[0], chMap[-1], 1)
    # idx = np.where(yCts!=0)
    # x, yCts = x[idx]-0.5, yCts[idx]
    # x = [int(ch) for ch in x]
    # xb = np.arange(0,len(x),1)
    # plt.bar(xb, yCts, 0.95, color='blue', label='chans')
    # plt.xticks(xb)
    # plt.gca().set_xticklabels(x, fontsize=12, rotation='vertical')
    # plt.xlabel("channel", ha='right', x=1.)
    # plt.legend(fontsize=14, ncol=2)
    # plt.savefig("../plots/sea-chans.png")

    # plt.cla()
    # plt.plot(hitPass, fsPass, '.', c='k', ms=5., label='pass')
    # plt.plot(hitFail, fsFail, '.', c='r', ms=5., label='fail')
    # plt.xlabel("Energy (keV)", ha='right', x=1.)
    # plt.ylabel("fitSlo", ha='right', y=1.)
    # plt.ylim(0, 500)
    # plt.xlim(0, 20)
    # plt.legend(loc=1)
    # plt.savefig("../plots/sea-fsPass.png")

    # for c in [578, 580, 582, 592, 598, 600, 608, 610, 626, 632, 640, 648, 664, 672, 690, 692]:
    #     print(c)
    #     hitPass, fsPass, hitFail, fsFail = [], [], [], []
    #     for i in range(len(hitE)):
    #         if chans[i] != c: continue
    #         if fSlo[i] < fsCut[i] and hitE[i] > (thMu[i]+3*thSig[i]):
    #             hitPass.append(hitE[i])
    #             fsPass.append(fSlo[i])
    #         else:
    #             hitFail.append(hitE[i])
    #             fsFail.append(fSlo[i])
    #     plt.cla()
    #     plt.plot(hitPass, fsPass, '.', c='k', ms=5., label='pass')
    #     plt.plot(hitFail, fsFail, '.', c='r', ms=5., label='fail')
    #     plt.xlabel("Energy (keV)", ha='right', x=1.)
    #     plt.ylabel("fitSlo", ha='right', y=1.)
    #     plt.axvline(1.0, color='b', label='1 keV')
    #     plt.axvline(1.5, color='m', label='1.5 keV')
    #     plt.ylim(0, 500)
    #     plt.xlim(0, 20)
    #     plt.legend(loc=1)
    #     plt.savefig("../plots/sea-fsPass-ch%d.png" % c)

    hitPass, fsPass, hitFail, fsFail = [], [], [], []
    for i in range(len(hitE)):
        if chans[i] in [610, 690, 692, 600, 598, 578]:
            continue
        if hitE[i] < 1.0: continue
        if fSlo[i] < fsCut[i] and hitE[i] > (thMu[i]+3*thSig[i]):
            hitPass.append(hitE[i])
            fsPass.append(fSlo[i])
        else:
            hitFail.append(hitE[i])
            fsFail.append(fSlo[i])
    fig = plt.figure()
    xLo, xHi, xpb = 0.5, 20, 0.1
    plt.cla()
    plt.plot(hitFail, fsFail, '.', c='r', ms=10., label='fail')
    plt.plot(hitPass, fsPass, '.', c='k', ms=10., label='pass')
    plt.axvline(1.0, c='b', label='1 keV')
    plt.xlabel("Energy (keV)", ha='right', x=1.)
    plt.ylabel("fitSlo", ha='right', y=1.)
    plt.ylim(0, 500)
    plt.xlim(0, 20)
    plt.legend(loc=1)
    plt.savefig("../plots/sea-fsPass-2.png")


def plotData2():

    f = np.load("../plots/sea-plt-2.npz")

    hitE, fSlo, chans, runs = f['arr_0'], f['arr_1'], f['arr_2'], f['arr_3']
    fsCut, thMu, thSig = f['arr_4'], f['arr_5'], f['arr_6']
    rise, riseCut = f['arr_7'], f['arr_8']

    # hitPass, rnPass, hitFail, rnFail = [], [], [], []
    # for i in range(len(hitE)):
    #     if chans[i] in [610, 690, 692, 600, 598, 578, 592]:
    #         continue
    #     # if rise[i] < riseCut[i] and fSlo[i] < fsCut[i] and hitE[i] > (thMu[i]+3*thSig[i]):
    #     # if fSlo[i] < fsCut[i] and hitE[i] > (thMu[i]+3*thSig[i]):
    #     # eThresh = thMu[i]+3*thSig[i]
    #     # if hitE[i] > (eThresh):
    #         # print("%d  e %.2f  t %.2f" % (chans[i], hitE[i],eThresh))
    #     if rise[i] < riseCut[i]:
    #         hitPass.append(hitE[i])
    #         rnPass.append(rise[i])
    #     else:
    #         hitFail.append(hitE[i])
    #         rnFail.append(rise[i])

    # hitPass, rnPass = [], []
    # for i in range(len(hitE)):
    #     if chans[i] in [610, 690, 692, 600, 598, 578, 592]:
    #         continue
    #     if rise[i] < riseCut[i] and hitE[i] > (thMu[i]+3*thSig[i]):
    #         rnPass.append(rise[i])
    #         hitPass.append(hitE[i])


    fig = plt.figure()
    xLo, xHi, xpb = 0.5, 20, 0.1
    # plt.semilogy(hitE,rise,".b",ms=5,label='all') # this shows the untagged pulser evts.  whew.

    # for the purposes of the talk, why not just draw a line at 10?
    n = len(hitE)
    hitPass = [hitE[i] for i in range(n) if rise[i] < 10]
    rnPass = [rise[i] for i in range(n) if rise[i] < 10]
    hitFail = [hitE[i] for i in range(n) if rise[i] > 10]
    rnFail = [rise[i] for i in range(n) if rise[i] > 10]

    plt.semilogy(hitPass,rnPass,".b", ms=5, label='pass')
    plt.semilogy(hitFail,rnFail,".r", ms=5, label='fail')
    plt.axhline(10,color='red')

    plt.ylim(0.5, 3e2)
    plt.xlim(0, 30)
    plt.xlabel("Energy (keV)",ha='right',x=1.)
    plt.ylabel("riseNoise",ha='right',y=1.)
    plt.legend(loc=1)
    plt.savefig("../plots/sea-rn-simpleCut.png")


def getData3():
    from ROOT import TChain

    # fCut = sorted(glob.glob("%s/results_v1/fs/fitSlo-DS%d-*.root" % (dsi.dataDir, dsNum)))
    # fCut = sorted(glob.glob("%s/results_v1/fs/fitSlo-DS%d-*.root" % (dsi.dataDir, dsNum)))

    print("drawing LAT files ...")
    ch = TChain("skimTree")
    ds = 1
    fileList = []
    bkg = dsi.BkgInfo()
    det = dsi.DetInfo()
    goodChans = det.getGoodChanList(ds)
    dsMap = bkg.dsMap()
    for sub in range(dsMap[ds]+1):
        # if 38 < sub < 45: continue
        # if 33 < sub < 51: continue
        print("Sub:",sub)
        latList = dsi.getSplitList("%s/latSkimDS%d_%d*" % (dsi.latDir, ds, sub), sub)
        tmpList = [f for idx, f in sorted(latList.items())]
        for f in tmpList: ch.Add(f)
    print(ch.GetEntries())
    n = ch.Draw("trapENFCal:riseNoise:channel","trapENFCal < 250","goff")
    t1, t2, t3 = ch.GetV1(), ch.GetV2(), ch.GetV3()
    t1 = np.asarray([t1[i] for i in range(n) if t3[i] in goodChans])
    t2 = np.asarray([t2[i] for i in range(n) if t3[i] in goodChans])

    # passing v1 of riseNoise cut
    print("drawing v1 files ...")
    ch2 = TChain("skimTree")
    ch2.Add("~/project/results_v1/rn/riseNoise-DS1-*.root")
    n = ch2.Draw("trapENFCal:riseNoise:channel","trapENFCal < 250","goff")
    t4, t5, t6 = ch2.GetV1(), ch2.GetV2(), ch2.GetV3()
    t4 = np.asarray([t4[i] for i in range(n) if t6[i] in goodChans])
    t5 = np.asarray([t5[i] for i in range(n) if t6[i] in goodChans])

    np.savez("../plots/sea-data3.npz",t1,t2,t4,t5)

def plotData3():

    f = np.load("../plots/sea-data3.npz")
    t1 = f['arr_0']
    t2 = f['arr_1']
    t4 = f['arr_2']
    t5 = f['arr_3']

    fig = plt.figure()
    plt.semilogy(t1, t2, ".r", ms=2)
    plt.semilogy(t4, t5, ".b", ms=3)
    plt.ylim(0.5, 3e2)
    plt.xlim(0, 50)

    plt.savefig("../plots/sea-riseNoise.png")


if __name__=="__main__":
    main()