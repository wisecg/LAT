#!/usr/bin/env python3
import sys, os, imp, glob
import numpy as np
import subprocess as sp
import tinydb as db
from scipy.optimize import curve_fit
import matplotlib
matplotlib.use('Agg')
sys.argv.append("-b")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
ds = imp.load_source('DataSetInfo',os.environ['LATDIR']+'/DataSetInfo.py')
wl = imp.load_source('waveLibs',os.environ['LATDIR']+'/waveLibs.py')
calInfo = ds.CalInfo()

def main():

    # getDTMult()
    # plotDTMult()
    # chanNoiseRate()
    # plotNoiseRate()
    # chanNoiseRate2()
    # plotChanNoiseRate2()
    m4Prob()


def getDTMult():
    from ROOT import TFile, TTree

    # 5 hr M1 calibration: https://majorana.npl.washington.edu/elog/Run+Elog/1703
    runList = calInfo.GetSpecialRuns("longCal",5)
    fileList = []
    for run in runList:
        fileList.extend(ds.getLATRunList([run],"%s/lat" % (ds.specialDir)))

    fileList = fileList[:100]
    nFiles = len(fileList)

    dtVals = {mH:[] for mH in range(50)}

    for iFile, f in enumerate(fileList):

        print("%d/%d %s" % (iFile,nFiles,f))

        tf = TFile("%s/lat/%s" % (ds.specialDir,f))
        latTree = tf.Get("skimTree")

        bNames = ["startClockTime_s","clockTime_s","tOffset","trapENFCal","mH","fitSlo","channel","gain","isGood","sumEH"]
        latTree.SetBranchStatus('*',0)
        for name in bNames:
            latTree.SetBranchStatus(name,1)
        latTree.GetEntry(0)
        startTime = latTree.startClockTime_s

        dt, prev = {}, {}

        for iEnt in range(latTree.GetEntries()):
            latTree.GetEntry(iEnt)
            clockTime = latTree.clockTime_s
            nHit = latTree.channel.size()
            mH = latTree.mH
            evtTime = clockTime-startTime

            if mH not in prev:
                prev[mH] = 0.

            if prev[mH] != 0:
                dt[mH] = evtTime - prev[mH]
                # print("%-4d  %-3d  evt %-10.9f  dt[%d] %-10.9f  pr[%d] %-10.9f" % (iEnt, mH, evtTime, mH, dt[mH], mH, prev[mH]))
                dtVals[mH].append(dt[mH])
            # else:
                # print("%-4d  %-3d  evt %-10.9f  pr[%d] %-10.9f" % (iEnt, mH, evtTime, mH, prev[mH]))

            prev[mH] = evtTime


    np.savez("../plots/dtVals.npz",dtVals)


def expoDist(x, beta, amp):
    """ The waiting time between events in a poisson process is an exponential distribution.
    https://nicolewhite.github.io/2015/05/23/understanding-waiting-times.html """
    return amp/beta * np.exp(-x / beta)


def plotDTMult():

    f = np.load("../plots/dtVals.npz")
    dtVals = f['arr_0'].item() # trick to recover a dict

    fig = plt.figure(figsize=(9,6))

    xLo, xHi, xpb = 0., 2., 0.02
    nb = int((xHi-xLo)/xpb)

    # plt.hist(dtVals[1], bins=nb, range=(xLo,xHi), histtype='step', linewidth=0.5, color='r', log=True)
    # plt.plot(np.nan, np.nan, c='r', label='hist m=1')

    cmap = plt.cm.get_cmap('hsv',5)
    for idx in range(1,3):

        y,x = np.histogram(dtVals[idx], bins=nb, range=(xLo, xHi))

        popt,_ = curve_fit(expoDist, x[1:], y)
        plt.semilogy(x, expoDist(x, *popt), 'r-')

        y = np.insert(y,0,0)
        plt.semilogy(x, y, ls='steps', color=cmap(idx), label="mH==%d" % idx)

    plt.ylim(2,1e7)
    plt.legend(loc='best')
    plt.xlabel("delta-t since last m=N event (sec)", horizontalalignment='right', x=1.0)
    plt.ylabel("Counts (arb)", horizontalalignment='right', y=1.0)
    plt.tight_layout()
    plt.savefig("../plots/dt-expoFit.png")


def chanNoiseRate():
    from ROOT import TFile, TTree

    runList = calInfo.GetSpecialRuns("longCal",5)
    fileList = []
    for run in runList:
        fileList.extend(ds.getLATRunList([run],"%s/lat" % (ds.specialDir)))
    fileList = fileList[:200]
    nFiles = len(fileList)

    xLo, xHi, xpb = 0., 50., 0.2
    nb = int((xHi-xLo)/xpb)

    chList = ds.GetGoodChanList(5)
    chSpec = {ch:np.zeros(nb) for ch in chList}

    for iFile, f in enumerate(fileList):
        print("%d/%d %s" % (iFile,nFiles,f))
        tf = TFile("%s/lat/%s" % (ds.specialDir,f))
        latTree = tf.Get("skimTree")

        theCut = "gain==0 && isGood && trapENFCal < 50"
        tNames = ["channel","trapENFCal","isGood","gain"]
        tVals = wl.GetVX(latTree, tNames, theCut)
        chData = tVals["channel"]
        nPass = len(chData)

        for ch in chList:
            idx = np.where(chData == ch)
            y,x = np.histogram(tVals["trapENFCal"][idx], bins=nb, range=(xLo, xHi) )
            chSpec[ch] += y

        tf.Close()

    np.savez("../plots/cal-loChanSpec.npz",x,chSpec)


def plotNoiseRate():
    f = np.load("../plots/cal-loChanSpec.npz")
    x = f['arr_0']
    x = x[1:]
    chSpec = f['arr_1'].item() # trick to recover a dict

    chList = ds.GetGoodChanList(5)

    fig = plt.figure(figsize=(9,6))

    xLo, xHi, xpb = 0., 50., 0.2
    nb = int((xHi-xLo)/xpb)

    cmap = plt.cm.get_cmap('hsv',len(chList)+1)
    for idx, ch in enumerate(chSpec):
        if np.sum(chSpec[ch])==0:
            print("no hits in ch",ch)
            continue

        avg50 = np.mean(chSpec[ch][np.where((x>=40) & (x<=50))])
        nPhys = avg50 / xpb
        nNoise = np.sum(chSpec[ch][np.where((x>=0) & (x<=10))]) - nPhys
        print("ch %d  phys %d  noise %d  N/P %.2f" % (ch, nPhys, nNoise, nNoise/nPhys))

        plt.semilogy(x, chSpec[ch], ls='steps', color=cmap(idx), label='ch %d' % ch)

    plt.legend(loc='best')
    plt.tight_layout()
    plt.savefig("../plots/dt-test.png")


    plt.cla()

    ch = chList[7]

    plt.plot(x, chSpec[ch], ls='steps', color='red', label='ch %d' % ch)

    avg50 = np.mean(chSpec[ch][np.where((x>=40) & (x<=50))])
    nPhys = avg50 / xpb
    nNoise = np.sum(chSpec[ch][np.where((x>=0) & (x<=10))]) - nPhys
    print("ch %d  phys %d  noise %d  N/P %.2f" % (ch, nPhys, nNoise, nNoise/nPhys))

    plt.axhline(avg50,color='g',label='avg 40-50')

    plt.legend(loc='best')
    plt.tight_layout()
    plt.savefig("../plots/noise-test.png")


def chanNoiseRate2():
    """ save energy & timestamp for channels s/t we can look at
    the dt distributions for two different energy regions. """
    from ROOT import TFile, TTree, TChain, GATDataSet

    runList = calInfo.GetSpecialRuns("longCal",5)
    runList = runList[:10] # limit
    fileList = []
    for run in runList:
        fileList.extend(ds.getLATRunList([run],"%s/lat" % (ds.specialDir)))
    nFiles = len(fileList)

    chList = ds.GetGoodChanListNew(5)
    chData = {ch:[] for ch in chList}

    runTime = 0
    for run in runList:
        gds = GATDataSet()
        gatPath = gds.GetPathToRun(run, GATDataSet.kGatified)
        tf = TFile(gatPath)
        gatTree = tf.Get("mjdTree")
        gatTree.GetEntry(0)
        runTime += gatTree.timeinfo.stopTime - gatTree.timeinfo.startTime
        tf.Close()
        print(run, runTime)

    for iFile, f in enumerate(fileList):
        print("%d/%d %s" % (iFile,nFiles,f))
        tf = TFile("%s/lat/%s" % (ds.specialDir,f))
        latTree = tf.Get("skimTree")

        theCut = "gain==0 && isGood && trapENFCal < 50"
        tNames = ["channel","trapENFCal","isGood","gain","startClockTime_s","clockTime_s","tOffset","mH","Entry$"]
        tVals = wl.GetVX(latTree, tNames, theCut)

        for idx in range(len(tVals["channel"])):
            ent = tVals["Entry$"][idx]
            ch = tVals["channel"][idx]
            ene = tVals["trapENFCal"][idx]
            toff = tVals["tOffset"][idx]/1e9
            ts = tVals["clockTime_s"][idx] - tVals["startClockTime_s"][idx] + toff
            mH = tVals["mH"][idx]
            chData[ch].append((ene,ts,mH,ent))

        tf.Close()

    np.savez("../plots/cal-loChan-E-Ts.npz",runTime,chData)


def plotChanNoiseRate2():
    f = np.load("../plots/cal-loChan-E-Ts.npz")
    runTime = f['arr_0']
    chData = f['arr_1'].item()

    chList = ds.GetGoodChanListNew(5)

    fig = plt.figure(figsize=(9,6))

    for ch in chList:

        plt.cla()

        if len(chData[ch]) < 10:
            print("no data for ch %d, continuing." % ch)
            continue

        eVals = [chData[ch][idx][0] for idx in range(1,len(chData[ch]))]
        dtVals = [chData[ch][idx][1] - chData[ch][idx-1][1] for idx in range(1,len(chData[ch]))]

        eVals, dtVals = np.asarray(eVals), np.asarray(dtVals)

        e10, dt10, e20, dt20, e30, dt30, e40, dt40, e50, dt50 = [], [], [], [], [], [], [], [], [], []
        p10, p20, p30, p40, p50 = 0., 0., 0., 0., 0.

        for idx in range(len(chData[ch])):
            ene = chData[ch][idx][0]
            ts = chData[ch][idx][1]
            if 0 < ene < 10:
                dt10.append(ts - p10)
                e10.append(ene)
                p10 = ts
            if 10 < ene < 20:
                dt20.append(ts - p20)
                e20.append(ene)
                p20 = ts
            if 20 < ene < 30:
                dt30.append(ts - p30)
                e30.append(ene)
                p30 = ts
            if 30 < ene < 40:
                dt40.append(ts - p40)
                e40.append(ene)
                p40 = ts
            if 40 < ene < 50:
                dt50.append(ts - p50)
                e50.append(ene)
                p50 = ts

        # compute physics and noise rate
        xLo, xHi, xpb = 0, 50, 0.2
        nb = int((xHi-xLo)/xpb)
        y, x = np.histogram(eVals, bins=nb, range=(xLo, xHi))
        x = x[1:]
        avg50 = np.mean(y[np.where((x>40) & (x<50))])
        nPhys = avg50 / xpb
        nNoise = np.sum(y[np.where((x>0) & (x<10))]) - nPhys

        rPhys, rNoise = nPhys/runTime, nNoise/runTime

        titleStr = "Ch %d  phys: %.5f Hz, noise: %.5f Hz  N/P %.2f" % (ch, rPhys, rNoise, nNoise/nPhys)
        plt.title(titleStr)
        print(titleStr)

        # plot dt since last event in each energy range
        cmap = plt.cm.get_cmap('hsv',5)
        plt.plot(e10,dt10,'.',color=cmap(0), label='10')
        plt.plot(e20,dt20,'.',color=cmap(1), label='20')
        plt.plot(e30,dt30,'.',color=cmap(2), label='30')
        plt.plot(e40,dt40,'.',color=cmap(3), label='40')
        plt.plot(e50,dt50,'.',color=cmap(4), label='50')
        plt.xlim(0, 50)
        plt.ylim(-0.1, 10.)
        plt.xlabel("trapENFCal", horizontalalignment='right', x=1.)
        plt.ylabel("dt (sec) since last evt in e-bin", horizontalalignment='right', y=1.)

        plt.legend(loc='best')
        plt.savefig("../plots/dt-chan%d-ebins.png" % ch)


def m4Prob():
    from ROOT import TFile, TTree

    # from plotChanNoiseRate2
    chRates = {
        584: (0.02504, 0.65272), 592: (0.04016, 1.12538), 598: (0.04214, 7.43857), 608: (0.04133, 0.89446),
        610: (0.05058, 1.16650), 614: (0.02029, 1.03391), 624: (0.02900, 0.93005), 626: (0.03535, 0.96318),
        628: (0.04316, 1.48850), 632: (0.02492, 0.57832), 640: (0.03153, 0.73312), 648: (0.03063, 0.57880),
        658: (0.01900, 0.52327), 660: (0.02465, 0.46342), 662: (0.02693, 0.60194), 672: (0.01873, 0.50940),
        674: (0.02377, 0.62129), 678: (0.02555, 0.58299), 680: (0.01812, 0.55507), 688: (0.04286, 1.03254),
        690: (0.01906, 0.62837), 692: (0.02236, 0.77410), 694: (0.03285, 0.83018),
        1106: (0.00210, 0.26005), 1120: (0.01229, 0.46665), 1124: (0.00514, 0.24493), 1128: (0.00198, 0.23542),
        1170: (0.01055, 0.43628), 1172: (0.00535, 0.45415), 1174: (0.00204, 0.12108), 1176: (0.00595, 0.20878),
        1204: (0.00472, 0.21207), 1208: (0.01247, 0.37574), 1232: (0.00201, 0.08488), 1236: (0.00234, 0.07277),
        1298: (0.00078, 0.09937), 1302: (0.01199, 0.32350), 1330: (0.01064, 0.40703), 1332: (0, 0)
        }

    runList = calInfo.GetSpecialRuns("longCal",5)
    runList = runList[:1] # limit
    fileList = []
    for run in runList:
        fileList.extend(ds.getLATRunList([run],"%s/lat" % (ds.specialDir)))
    nFiles = len(fileList)

    pMon = ds.PMon[5] # ds5 pulser monitor list

    for iFile, f in enumerate(fileList):
        print("%d/%d %s" % (iFile,nFiles,f))
        tf = TFile("%s/lat/%s" % (ds.specialDir,f))
        latTree = tf.Get("skimTree")

        theCut = "gain==0 && isGood && mH==4 && sumEH > 200 && sumEH < 1000"
        tNames = ["Entry$","channel","trapENFCal","mH","tOffset","clockTime_s","startClockTime_s","sumEH"]
        tVals = wl.GetVX(latTree, tNames, theCut)

        prevEnt = 0
        for idx in range(len(tVals["channel"])):
            ent = tVals["Entry$"][idx]
            ch = tVals["channel"][idx]
            ene = tVals["trapENFCal"][idx]
            mH = tVals["mH"][idx]
            toff = tVals["tOffset"][idx]/1e9
            ts = tVals["clockTime_s"][idx] - tVals["startClockTime_s"][idx] + toff

            rPhys, rNoise = 0., 0.
            if ch in chRates:
                rPhys, rNoise = chRates[ch][0], chRates[ch][1]

            ignore = "x" if ch in pMon else ""

            if ent != prevEnt:
                print(" ")

            print("%d  %d  m%d  e %-8.2f  t %-8.5f  rp %.4f  rn %.4f  %s" % (ent,ch,mH,ene,ts,rPhys,rNoise,ignore))

            prevEnt = ent

        tf.Close()


if __name__=="__main__":
    main()