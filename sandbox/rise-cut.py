#!/usr/bin/env python3
"""
Should be able to quickly merge this code w/ lat2.py
once we know what we need.

RiseNoise should be set w/ ~1 cal run at each calIdx.
    maybe re-use the softplus function fit
    adapt code from LAT3
    verify the cut here,

then move the db-writing part to LAT2,
as well as the cut file generation.
"""
import sys, os, time
import numpy as np
from scipy.optimize import curve_fit

import matplotlib as mpl
mpl.use('Agg')
sys.argv.append("-b")
import matplotlib.pyplot as plt
plt.style.use('../pltReports.mplstyle')
from matplotlib.colors import LogNorm, Normalize
# from matplotlib import gridspec

import dsi
bkg = dsi.BkgInfo()
cal = dsi.CalInfo()
det = dsi.DetInfo()
skipDS6Cal=True
import waveLibs as wl


def main(argv):

    ds, cIdx, mod = None, None, None
    for i, opt in enumerate(argv):

        # call scanRunsRise directly (used by job-panda)
        if opt=="-scan":
            ds, key, mod, cIdx = int(argv[i+1]), argv[i+2], int(argv[i+3]), int(argv[i+4])
            scanRunsRise(ds,key,mod,cIdx)

        if opt=="-set":
            setRiseCut()


def scanRunsRise(ds, key, mod, cIdx):
    from ROOT import TFile, TTree

    rLim, eLim = 4, 250
    print("Limiting to",rLim,"runs and a",eLim,"keV hit upper limit.")

    # load file and channel list
    fileList = []
    calRuns = cal.GetCalList(key,cIdx,runLimit=rLim) # should not need much for riseNoise
    for run in calRuns:
        latList = dsi.getSplitList("%s/latSkimDS%d_run%d*" % (dsi.calLatDir, ds, run), run)
        tmpList = [f for idx, f in sorted(latList.items())]
        fileList.extend(tmpList)
    chList = det.getGoodChanList(ds)

    print("Scanning DS:%d  calIdx %d  mod %d  key %s  nFiles:%d" % (ds, cIdx, mod, key, len(fileList)), time.strftime('%X %x %Z'))
    outFile = "%s/rise_%s_c%d.npz" % (dsi.effDir, key, cIdx)
    print("Saving output in:",outFile)

    # this is what we'll output for every calIdx
    hitE, chan, rise = [], [], []

    # loop over LAT cal files
    scanStart = time.time()
    prevRun = 0
    evtCtr, totCtr, totRunTime = 0, 0, 0
    for iF, f in enumerate(fileList):

        print("%d/%d %s" % (iF, len(fileList), f))
        tf = TFile(f)
        tt = tf.Get("skimTree")

        tt.GetEntry(0)
        run = tt.run
        if run!=prevRun:
            calIdx = cal.GetCalIdx(key,run)
            start = tt.startTime_s
            stop = tt.stopTime_s
            runTime = stop-start
            if runTime < 0 or runTime > 9999:
                print("run time error, run",run,"start",start,"stop")
            else:
                totRunTime += runTime

            # find thresholds for this run,
            # to calculate sumET and mHT in the loop.

            n = tt.Draw("channel:threshKeV:threshSigma","","goff")
            thrC, thrM, thrS = tt.GetV1(), tt.GetV2(), tt.GetV3()
            tmpThresh = {}
            for i in range(n):
                if thrC[i] not in chList:
                    continue
                if thrC[i] in tmpThresh.keys():
                    continue
                if thrM[i] < 9999:
                    thrK = thrM[i] + 3*thrS[i]
                    tmpThresh[thrC[i]] = [run,thrM[i],thrS[i],thrK]
            for ch in chList:
                if ch not in tmpThresh.keys():
                    tmpThresh[ch] = [-1,-1,-1,-1]

        prevRun = run
        # continue

        # loop over tree
        for iE in range(tt.GetEntries()):
            tt.GetEntry(iE)
            if tt.EventDC1Bits != 0: continue
            # totCtr += 1

            n = tt.channel.size()
            chTmp = np.asarray([tt.channel.at(i) for i in range(n)])
            idxRaw = [i for i in range(tt.channel.size()) if tt.channel.at(i) in chList]
            hitERaw = np.asarray([tt.trapENFCal.at(i) for i in idxRaw])

            # get indexes of hits above threshold (use thresholds from THIS CAL RUN)
            idxList = [i for i in range(tt.channel.size())
                if tt.channel.at(i) in chList
                and tt.trapENFCal.at(i) > tmpThresh[tt.channel.at(i)][3]
                and 0.7 < tt.trapENFCal.at(i) < eLim
                ]

            # save riseNoise data
            for i in idxList:
                hitE.append(tt.trapENFCal.at(i))
                chan.append(tt.channel.at(i))
                rise.append(tt.riseNoise.at(i))

    print("done.")
    hitE, chan, rise = np.asarray(hitE), np.asarray(chan), np.asarray(rise)
    print(len(hitE),'total entries')

    for ch in chList:
        idx = np.where(chan==ch)
        idx2 = np.where(hitE[idx] < 10)
        print(ch, "nTot",len(hitE[idx]), "nCts under 10 keV:",len(hitE[idx2]), "nCts<10/0.5 keV: ",len(hitE[idx2])/20)

    np.savez(outFile,hitE,chan,rise)
    print("Done:",time.strftime('%X %x %Z'),", %.2f sec/file." % ((time.time()-scanStart)/len(fileList)))


def loadRiseData(key):
    """ Load files generated by scanRunsRise, return data in a dict.
    To avoid confusion, must specify a key from runsCal.json .
    """
    if key not in cal.GetKeys():
        print("Unknown key!")
        return None
    else:
        print("Loading eff data for key:",key)

    # output dict
    eff = {}
    eff["hitE"] = {}  # {ci: [hitE1, hitE2 , ...] }
    eff["chan"] = {}  # {ci: [chan1, chan2 , ...] }
    eff["rise"] = {}  # {ci: [rise1, rise2, ...] }
    for ci in range(cal.GetIdxs(key)):
        eFile = "%s/rise_%s_c%d.npz" % (dsi.effDir, key, ci)
        if not os.path.isfile(eFile):
            print("File not found:",eFile)
            continue
        f = np.load(eFile)
        eff["hitE"][ci] = np.asarray(f['arr_0'])
        eff["chan"][ci] = np.asarray(f['arr_1'])
        eff["rise"][ci] = np.asarray(f['arr_2'])

    return eff


def setRiseCut():

    makePlots = False

    dsList = [1] # still need to generate the other ds's

    # loop over ds's, separated by cal key
    for ds in dsList:
        for calKey in cal.GetKeys(ds):

            chList = det.getGoodChanList(ds)
            mod = -1
            if "m1" in calKey:
                mod = 1
                chList = [ch for ch in chList if ch < 1000]
            if "m2" in calKey:
                mod = 2
                chList = [ch for ch in chList if ch > 1000]

            eff = loadRiseData(calKey)
            nCal = cal.GetNCalIdxs(ds,mod)

            # loop over calIdx's
            for ci in range(nCal):

                dbKey = "riseNoise_%s_ci%d_pol" % (calKey,ci)
                dbVals = {ch : None for ch in chList}

                for ch in chList:

                    cTmp = eff["chan"][ci]
                    idx = np.where(cTmp==ch)

                    hitE = eff["hitE"][ci][idx]
                    rise = eff["rise"][ci][idx]

                    # make sure we have hits, and that riseNoise vals are good
                    if len(hitE)==0 or len(rise[np.where(rise > 0)])==0:
                        print("No data, ch",ch)
                        continue

                    # fit the data to a pol1
                    popt, pcov = curve_fit(wl.pol1, hitE, rise)

                    # move the y-intercept up (c), and calculate how many events you keep
                    start = time.time()
                    nTot = len(hitE)
                    a, b, c = popt
                    fitPass = False
                    for i in range(1000):
                        c99 = c + 0.05*i
                        evtPass = np.asarray([[hitE[i],rise[i]] for i in range(nTot) if rise[i] <= wl.pol1(hitE[i],a,b,c99)])
                        evtFail = np.asarray([[hitE[i],rise[i]] for i in range(nTot) if rise[i] > wl.pol1(hitE[i],a,b,c99)])
                        nPass = len(evtPass)
                        if nPass/nTot > 0.995: # keep 99.5
                            fitPass = True
                            break

                    # print("DS%d ci%d ch%d  99pct fit: %.2f sec  a %-9.2e  b %-5.3f  c99 %.3f  Pass:" % (ds,ci,ch,time.time()-start,a,b,c99),fitPass)

                    dbVals[ch] = [a,b,c99,fitPass]

                    if makePlots:
                        plt.cla()
                        xLo, xHi, xpb = 0, 250, 1
                        # nbx = int((xHi-xLo)/xpb)
                        # yLo, yHi, ypb = 0, 5, 0.05
                        # nby = int((yHi-yLo)/ypb)
                        # plt.hist2d(hitE, rise, bins=[nbx, nby], range=[[xLo,xHi],[yLo,yHi]], cmap='jet',norm=LogNorm(),alpha=0.5)
                        cpd = det.getChanCPD(ds,ch)
                        plt.plot(np.nan, np.nan, ".w", label="ch %d, C%sP%sD%s" % (ch, cpd[0],cpd[1],cpd[2]))
                        xFit = np.arange(xLo, xHi, 0.1)
                        plt.plot(xFit, wl.pol1(xFit, *popt), 'r-', label="a %.4f b %.4f c %.4f" % tuple(popt))
                        plt.plot(xFit, wl.pol1(xFit, a,b,c99), 'g-', label="a %-9.4f b %-9.4f c99 %.4f" % (a,b,c99))
                        plt.plot(evtPass[:,0], evtPass[:,1], ".b", ms=1, label="pass")
                        plt.plot(evtFail[:,0], evtFail[:,1], ".r", ms=1, label="fail")
                        plt.xlabel("Energy (keV)", ha='right', x=1)
                        plt.ylabel("riseNoise", ha='right', y=1)
                        leg = plt.legend(loc='best', fontsize=12)
                        leg.get_frame().set_alpha(0.5)
                        plt.savefig("../plots/rise-ds%d-ci%d-ch%d.png" % (ds, ci, ch))
                        # return

                # final db check
                print(dbKey)
                for ch in sorted(dbVals):
                    print(ch, dbVals[ch])

                return


def plotStability():
    print("you need to plot each channel, each calIdx")


if __name__=="__main__":
    main(sys.argv[1:])

