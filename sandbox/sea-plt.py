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

# dsi = imp.load_source('dsi',os.environ['LATDIR']+'/dsi.py')
dsi = imp.load_source('dsi',os.environ['LATDIR']+'/sandbox/DataSetInfo.py')
wl = imp.load_source('waveLibs',os.environ['LATDIR']+'/waveLibs.py')

# load threshold data
import tinydb as db
dsNum, bkgIdx = 5, 83
calDB = db.TinyDB('../calDB.json')
pars = db.Query()
thD = dsi.getDBRecord("thresh_ds%d_bkgidx%d" % (dsNum, bkgIdx), False, calDB, pars)

# load threshold data
import tinydb as db
dsNum, bkgIdx = 5, 83
calDB = db.TinyDB('../calDB.json')
pars = db.Query()
thD = dsi.getDBRecord("thresh_ds%d_bkgidx%d" % (dsNum, bkgIdx), False, calDB, pars)

# load fitSlo vals for cal run range closest to the run range [22513, 22566]
dsNum, modNum, calIdx = 5, 1, 11  # calIdx 11: [[22568,22635],22568,22841],
fsD = dsi.getDBRecord("fitSlo_ds%d_idx%d_m%d_Peak" % (dsNum, calIdx, modNum), False, calDB, pars)

# load riseNoise vals
rnSD = dsi.getDBRecord("riseNoise_ds%d_idx%d_m%d_SoftPlus" % (dsNum, calIdx, modNum), False, calDB, pars)
rnCD = dsi.getDBRecord("riseNoise_ds%d_idx%d_m%d_Continuum" % (dsNum, calIdx, modNum), False, calDB, pars)

def main():

    # scrapeData()
    # plotData()
    # plotData2()
    # getData3()
    # plotData3()

    # plotThreshCut()

    # getCalData()
    # plotCalSlowness()
    plotSlowness2()

    # plotSloHits()


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

    plt.semilogy(hitPass,rnPass,".k", ms=5, label='pass')
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


def plotThreshCut():
    """ Adapted from LAT/sandbox/mult2.py, and dependent on input:
    "../plots/mult2-dtVals-ene.npz"
    """
    f = np.load("../plots/mult2-dtVals-ene.npz")
    dtVals = f['arr_0'].item()
    mH = 1
    n = len(dtVals[mH])

    # plot 1: all hits under 10 keV
    # chan = [dtVals[mH][i][1][0] for i in range(n) if dtVals[mH][i][2][0]<10]
    # hitE = [dtVals[mH][i][2][0] for i in range(n) if dtVals[mH][i][2][0]<10]
    chan = [dtVals[mH][i][1][0] for i in range(n) if dtVals[mH][i][1][0] not in [598,1172] and dtVals[mH][i][2][0]<10]
    hitE = [dtVals[mH][i][2][0] for i in range(n) if dtVals[mH][i][1][0] not in [598,1172] and dtVals[mH][i][2][0]<10]

    # plot 2: only take hits above the detector threshold.
    n = len(hitE)
    chanTh, hitETh = [], []
    for i in range(n):
        if chan[i] in thD.keys() and hitE[i] > thD[chan[i]][0] + 3*thD[chan[i]][1]:
            chanTh.append(chan[i])
            hitETh.append(hitE[i])

    # make the figure

    f = plt.figure(figsize=(20,8))
    p1 = plt.subplot(121)
    p2 = plt.subplot(122)

    chMap = list(sorted(set(chan)))
    chDict = {chMap[i]:i for i in range(len(chMap))}
    chan = [chDict[chan] for chan in chan]
    chanTh = [chDict[chan] for chan in chanTh]

    xLo, xHi, xpb = 0.5, 5., 0.1
    yLo, yHi = 0, len(chMap)
    nbx, nby = int((xHi-xLo)/xpb), len(chMap)


    h1,xedg1,yedg1 = np.histogram2d(hitE,chan,bins=[nbx,nby], range=[[xLo,xHi],[yLo,yHi]])
    h2,xedg2,yedg2 = np.histogram2d(hitETh,chanTh,bins=[nbx,nby], range=[[xLo,xHi],[yLo,yHi]])
    h1, h2 = h1.T, h2.T
    hMin, hMax = np.amin(h1), np.amax(h2)
    hMin, hMax = 2, hMax*1.3

    # arr[arr > 255] = x
    h1[h1 < 0.1] = -1 # avoid LogNorm? set zeros to white?
    h2[h2 < 0.1] = -1

    im1 = p1.imshow(h1,cmap='jet',vmin=hMin,vmax=hMax)#,norm=LogNorm())
    xticklabels = ["%.1f" % t for t in np.arange(0, 5.5, 0.5)]
    yticks = np.arange(0, len(chMap))

    p1.set_xlabel("Energy (keV)", ha='right', x=1.)
    p1.set_xticklabels(xticklabels)
    p1.set_ylabel("channel", ha='right', y=1.)
    p1.set_yticks(yticks)
    p1.set_yticklabels(chMap, fontsize=12)
    f.colorbar(im1, ax=p1, fraction=0.037, pad=0.04)

    im2 = p2.imshow(h2,cmap='jet',vmin=hMin,vmax=hMax)#,norm=LogNorm())
    p2.set_xlabel("Energy (keV)", ha='right', x=1.)
    p2.set_xticklabels(xticklabels)
    p2.set_ylabel("channel", ha='right', y=1.)
    p2.set_yticks(yticks)
    p2.set_yticklabels(chMap, fontsize=12)
    f.colorbar(im2, ax=p2, fraction=0.037, pad=0.04)

    plt.tight_layout()

    plt.savefig("../plots/sea-dtm1-thresh.png")


def getCalData():
    """ load a small subset of the long calibration files
        and apply v1 of the fitSlo and threshold cuts to the hits.
        also uses the old DataSetInfo (see top of this file)
    """
    from ROOT import TFile, TTree
    fileList = []
    calInfo = dsi.CalInfo()
    runList = calInfo.GetSpecialRuns("longCal",5) # 54 total
    runList = runList[:2]
    # runList = runList[:] # histats file, see filename below
    for run in runList: fileList.extend(dsi.getLATRunList([run],"%s/lat" % (dsi.specialDir)))
    nFiles = len(fileList)
    chList = dsi.GetGoodChanListNew(dsNum=5)

    chan, hitE, fSlo, rise = [], [], [], []

    for f in fileList:
        print(f)
        tf = TFile("%s/longCal/%s" % (dsi.specialDir,f))
        tt = tf.Get("skimTree")
        n = tt.Draw("channel:trapENFCal:fitSlo:riseNoise","trapENFCal < 250","goff")
        ch, he, fs, rn = tt.GetV1(), tt.GetV2(), tt.GetV3(), tt.GetV4()
        ch = np.asarray([ch[i] for i in range(n)])
        he = np.asarray([he[i] for i in range(n)])
        fs = np.asarray([fs[i] for i in range(n)])
        rn = np.asarray([rn[i] for i in range(n)])
        chan.extend(ch)
        hitE.extend(he)
        fSlo.extend(fs)
        rise.extend(rn)
        tf.Close()

    np.savez("../plots/sea-cal.npz",chan,hitE,fSlo,rise)


def plotCalSlowness():

    f = np.load('../plots/sea-cal.npz')
    chanIn, hitEIn, fSloIn, riseIn = f['arr_0'], f['arr_1'], f['arr_2'], f['arr_3']

    # average slowness values
    fsVals = {
        584: 102.5, 592: 75.5, 608: 73.5, 610: 76.5, 614: 94.5, 624: 69.5,
        626: 81.5, 628: 102.5, 632: 81.5, 640: 73.5, 648: 74.5, 658: 75.5,
        660: 127.5, 662: 84.5, 672: 80.5, 678: 82.5, 680: 86.5, 688: 77.5,
        690: 80.5, 694: 80.5
        }

    # only plot channels we have fitSlo and threshold cuts for
    chList = []
    for ch in (list(fsD.keys())+list(thD.keys())+list(fsVals.keys())):
        if ch in fsD.keys() and ch in thD.keys() and ch in fsVals.keys():
            chList.append(ch)
    chList = sorted(list(set(chList)))

    # fill all hits and hits passing fitSlo and 3sig threshold cuts
    hitE, hitEPass, hitEFail = [], [], []
    fSlo, fSloPass, fSloFail = [], [], []
    for i in range(len(chanIn)):
        if chanIn[i] not in chList: continue
        hitE.append(hitEIn[i])
        fSlo.append(fSloIn[i])

        fsCut = fsD[chanIn[i]][2] # [1%, 5%, 90%, 95%, 99%]
        thCut = thD[chanIn[i]][0] + 3*thD[chanIn[i]][1]

        # calculate the shifted fitSlo value to make the plot more clear
        fSloShift = fSloIn[i] - fsVals[chanIn[i]]

        if fSloIn[i] < fsCut and hitEIn[i] > thCut:
            hitEPass.append(hitEIn[i])
            fSloPass.append(fSloShift)
        else:
            hitEFail.append(hitEIn[i])
            fSloFail.append(fSloShift)

    plt.cla()
    fig = plt.figure()
    xLo, xHi, xpb = 0, 250, 1

    x, y1 = wl.GetHisto(hitE,xLo,xHi,xpb)
    plt.semilogy(x, y1, 'r',lw=2.,ls='steps',label='all hits')

    x, y2 = wl.GetHisto(hitEPass,xLo,xHi,xpb)
    plt.semilogy(x, y2, 'b',lw=2.,ls='steps',label='hits passing cuts')

    plt.xlabel("Energy (keV)", ha='right', x=1)
    plt.ylabel("Counts", ha='right', y=1)
    plt.legend(loc=1)
    plt.savefig("../plots/sea-cal-fsCut.png")

    plt.cla()
    plt.plot(hitEFail,fSloFail,".r",ms=5,label='fail')
    plt.plot(hitEPass,fSloPass,".k",ms=5,label='pass')
    plt.axvline(1.0,color='b',label='1 keV')
    plt.xlabel("Energy (keV)", ha='right', x=1)
    plt.ylabel("fitSlo (channel-shifted)", ha='right', y=1)
    plt.legend(loc=1)
    plt.xlim(0,20)
    plt.ylim(-50,500)
    plt.savefig("../plots/sea-cal-fsThrScatter.png")


def plotSlowness2():

    # fitSlo results from tuneFitSlo.
    fsVals = {
        584: 102.5, 592: 75.5, 608: 73.5, 610: 76.5, 614: 94.5, 624: 69.5,
        626: 81.5, 628: 102.5, 632: 81.5, 640: 73.5, 648: 74.5, 658: 75.5,
        660: 127.5, 662: 84.5, 672: 80.5, 678: 82.5, 680: 86.5, 688: 77.5,
        690: 80.5, 694: 80.5
        }

    chList = list(fsVals.keys())

    # peak fit results from mult2.  could retune, but it shouldn't matter much
    mHT = 2
    sumPks = [
        (2,238,237.28,239.46), (2,583,581.26,584.46), (2,2615,2610.57,2618.01),
        (3,238,237.13,239.43), (3,583,581.04,584.36), (3,2615,2610.10,2617.92)
        ]

    # f1 = np.load("../plots/mult4-hitE.npz")
    f1 = np.load("../plots/mult4-hitE-histats.npz")
    runTime, hitList, hitData, eCut = f1['arr_0'], f1['arr_1'], f1['arr_2'], f1['arr_3']

    hitE, chan, fSlo = [], [], []
    for i in range(len(hitData)):
        if hitData[i][0]==mHT and 237.28 < hitData[i][1] < 239.46:
            hitE.extend(hitList[i][0])
            chan.extend(hitList[i][1])
            fSlo.extend(hitList[i][2])
    n = len(hitE)
    hitE = [hitE[i] for i in range(n) if chan[i] in fsVals.keys()]
    fSloShift = [fSlo[i]-fsVals[chan[i]] for i in range(n) if chan[i] in chList]
    print("peak evts:",n)

    # plt.cla()
    # plt.plot(hitE, fSloShift, '.', c='k', ms=1.)
    # plt.ylim(-50, 400)
    # plt.savefig("../plots/mult4-fitSlo-shift.png")

    cLo, cHi = 0, 230

    hitE2, chan2, fSlo2 = [], [], []
    for i in range(len(hitData)):
        # if hitData[i][0]==mHT and cLo < hitData[i][1] < cHi:
        if hitData[i][1] < 230:
            hitE2.extend(hitList[i][0])
            chan2.extend(hitList[i][1])
            fSlo2.extend(hitList[i][2])
    n = len(hitE2)
    hitE2 = [hitE2[i] for i in range(n) if chan2[i] in fsVals.keys()]
    fSloShift2 = [fSlo2[i]-fsVals[chan2[i]] for i in range(n) if chan2[i] in chList]
    print("cont evts:",n)

    # plt.cla()
    # plt.plot(hitE2, fSloShift2, '.', c='k', ms=1.)
    # plt.ylim(-50, 400)
    # plt.savefig("../plots/mult4-fitSlo-shift2.png")

    plt.cla()
    xLo, xHi, xpb = -50, 300, 1

    x, y238 = wl.GetHisto(fSloShift,xLo,xHi,xpb)
    x, yCon = wl.GetHisto(fSloShift2,xLo,xHi,xpb)

    plt.semilogy(x, y238/np.sum(y238),'b',lw=2.,ls='steps',label='hits m=2,s=238 pk')
    plt.semilogy(x, yCon/np.sum(yCon),'r',lw=2.,ls='steps',label='all hits %d - %d keV' % (cLo, cHi))
    plt.xlabel("fitSlo (shifted)",ha='right',x=1)
    plt.ylabel("Counts (arb)",ha='right',y=1)
    plt.legend()
    plt.savefig("../plots/sea-fitSlo-hist.png")


def plotSloHits():
    """ adapted from mult2::plotSloHits
    Fixes the bug I found 14 Apr 18, where mult2::plotSloHits was looking
    at the wrong values for fSlo and erroneously seeing a flat distribution.
    """
    # average slowness values
    fsVals = {
        584: 102.5, 592: 75.5, 608: 73.5, 610: 76.5, 614: 94.5, 624: 69.5,
        626: 81.5, 628: 102.5, 632: 81.5, 640: 73.5, 648: 74.5, 658: 75.5,
        660: 127.5, 662: 84.5, 672: 80.5, 678: 82.5, 680: 86.5, 688: 77.5,
        690: 80.5, 694: 80.5
        }
    # load a unique channel list
    chList = []
    for ch in list(fsVals.keys())+list(thD.keys())+list(fsD.keys()):
        if ch in fsVals.keys() and ch in fsD.keys() and ch in thD.keys():
            chList.append(ch)
    chList = sorted(list(set(chList)))

    # load data (from mult2::getSumEvents)
    # f = np.load("../plots/mult2-sumLoEvts.npz")
    f = np.load("../plots/mult2-sumLoEvts-histats.npz")
    runTime, sumEvts = f['arr_0'], f['arr_1'].item()
    mH, pk = 2, 238
    evts = sumEvts[mH][pk]

    # fill arrays
    hitE, fSlo, fSloShift = [], [], []
    ePass, sPass = [], []
    eCut, sCut = [], []

    for i in range(len(evts)):
        chanIn = evts[i][0]
        if chanIn not in chList: continue
        hitEIn = evts[i][2]
        if hitEIn < 0.7 or hitEIn < thD[chanIn][0] + 3*thD[chanIn][1]: continue

        fSloIn = evts[i][4]
        fShift = fSloIn - fsVals[chanIn]

        fCut = 1.5 * fsD[chanIn][2] # v1 fitSlo vals look too tight.

        hitE.append(hitEIn)
        fSlo.append(fSloIn)
        fSloShift.append(fShift)

        if fSloIn < fCut:
            ePass.append(hitEIn)
            # sPass.append(fSloIn)
            sPass.append(fShift)
        else:
            eCut.append(hitEIn)
            # sCut.append(fSloIn)
            sCut.append(fShift)

    # ========================================================================================
    # -- plot 1: 2d scatter plot, fSlo vs hitE, hits passing/failing v1 of fitSlo cut.
    fig = plt.figure()

    plt.cla()
    # plt.plot(hitE, fSlo, ".k", ms=1., label='all')
    # plt.plot(hitE, fSloShift, ".k", ms=1., label='all')
    plt.plot(ePass, sPass, ".k", ms=0.5, label='pass')
    plt.plot(eCut, sCut, '.r', ms=1.0, label='cut')

    plt.ylim(-100, 300)
    plt.xlim(0, 30)

    plt.xlabel("hitE (keV)", ha='right', x=1.)
    plt.ylabel("fitSlo", ha='right', y=1.)
    plt.legend(loc=1)
    plt.savefig("../plots/sea-fitSloScatter.png")

    # -- plot 2: 1-d energy histo, hits passing/failing
    plt.cla()

    xLo, xHi, xpb = 0, 15, 0.2
    x, hTot = wl.GetHisto(hitE, xLo, xHi, xpb)
    x, hPass = wl.GetHisto(ePass, xLo, xHi, xpb)
    x, hCut = wl.GetHisto(eCut, xLo, xHi, xpb)

    plt.plot(x, hTot, c='b', lw=1.5, ls='steps', label='m=%d, pk %d' % (mH, pk))
    plt.plot(x, hPass, c='g', lw=1.5, ls='steps', label='fitSlo pass')
    plt.plot(x, hCut, c='r', lw=1.5, ls='steps', label='fitSlo cut')

    plt.xlabel("hitE (keV)", ha='right', x=1.)
    plt.ylabel("Counts / %.1f kev" % xpb, ha='right', y=1.)
    plt.legend(loc=1)
    plt.savefig("../plots/sea-m2s238-hitE.png")

    # -- plot 3: efficiency of v1 slowness cut, assuming all m2s238 hits are fast.
    plt.cla()

    from statsmodels.stats import proportion
    idx = np.where((hTot > 0) & (hPass > 0))
    ci_low, ci_upp = proportion.proportion_confint(hPass[idx], hTot[idx], alpha=0.1, method='beta')
    sloEff = hPass[idx] / hTot[idx]

    nPad = len(hPass)-len(hPass[idx])
    sloEff = np.pad(sloEff, (nPad,0), 'constant', constant_values=0)
    ci_low = np.pad(ci_low, (nPad,0), 'constant', constant_values=0)
    ci_upp = np.pad(ci_upp, (nPad,0), 'constant', constant_values=0)

    plt.plot(x, sloEff, '.b', ms=10., label='efficiency')
    plt.errorbar(x, sloEff, yerr=[sloEff - ci_low, ci_upp - sloEff], color='k', linewidth=0.8, fmt='none')

    idx = np.where(np.abs(ci_upp) > 0)
    xCut = x[6]
    plt.axvline(xCut,color='g',label='cutoff: %.2f keV' % xCut)

    plt.xlabel("hitE (keV)", ha='right', x=1.)
    plt.ylabel("Efficiency", ha='right', y=1.)
    plt.legend(loc=4)

    plt.savefig("../plots/sea-m2s238-eff.png")


if __name__=="__main__":
    main()