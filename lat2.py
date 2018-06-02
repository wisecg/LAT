#!/usr/bin/env python3
"""
====== LAT2.py =======
Tunes cut parameters,
calculates cut efficiencies.
C. Wiseman, USC
v1. 17 Apr 2018
======================
"""
import sys, os, time
import numpy as np
import tinydb as db
from scipy.optimize import curve_fit

import matplotlib as mpl
mpl.use('Agg')
sys.argv.append("-b")
import matplotlib.pyplot as plt
plt.style.use('pltReports.mplstyle')
from matplotlib.colors import LogNorm, Normalize

import waveLibs as wl
import dsi
bkg = dsi.BkgInfo()
cal = dsi.CalInfo()
det = dsi.DetInfo()
skipDS6Cal = True # ignore DS6 cal runs until they're processed

def main(argv):

    global writeDB
    writeDB = False
    ds, cIdx, mod = None, None, None

    for i, opt in enumerate(argv):

        # set dataset and options
        if opt=="-ds":
            ds = argv[i+1]
        if opt=="-ci":
            ds = int(argv[i+1])
            cIdx = int(argv[i+2])
        if opt=="-m":
            mod = int(argv[i+1])
        if opt=="-db":
            writeDB = True

        # wrapper function for scanRunsSlo
        if opt=="-load":
            loadRuns(ds,cIdx,mod)

        # call scanRunsSlo directly (used by lat-jobs)
        if opt=="-scan":
            ds, key, mod, cIdx = int(argv[i+1]), argv[i+2], int(argv[i+3]), int(argv[i+4])
            scanRunsSlo(ds,key,mod,cIdx)

        # call scanRunsRise directly (used by lat-jobs)
        if opt == "-rscan":
            ds, key, mod, cIdx = int(argv[i+1]), argv[i+2], int(argv[i+3]), int(argv[i+4])
            scanRunsRise(ds,key,mod,cIdx)

        # set fitSlo cut
        if opt=="-fs":
            setSloCut()

        # set riseNoise cut
        if opt=="-rn":
            setRiseCut()

        # check riseNoise stability
        if opt=="-rs":
            riseStability()

        # update riseNoise DB entries for problem channels
        if opt=="-rc":
            badRiseChans()

        # generate cut files (specify cut type: (blank) fs rn fr)
        if opt=="-cut":
            if ds is None:
                for d in [0,1,2,3,4,"5A","5B","5C"]: applyCuts(d, argv[i+1])
            else:
                # cutType = "-th" if argv[i+1]!="-b" else argv[i+1]
                cutType = argv[i+1]
                applyCuts(ds,cutType)


def loadRuns(dsIn=None,subIn=None,modIn=None):

    # loop over datasets, skipping DS6 cal runs till they're processed
    for ds in [0,1,2,3,4,5,6]:
        if skipDS6Cal is True and ds==6:
            continue

        if dsIn is not None and ds!=dsIn:
            continue

        # loop over keys in this DS
        for key in cal.GetKeys(ds):

            mod = -1
            if "m1" in key: mod = 1
            if "m2" in key: mod = 2

            # loop over cIdx's for this key
            for cIdx in range(cal.GetIdxs(key)):
                if subIn is not None and cIdx!=subIn:
                    continue
                if modIn is not None and mod!=modIn:
                    continue

                # now that ds, cIdx, and module are determined, we can call scanRuns
                scanRunsSlo(ds, key, mod, cIdx)


def scanRunsSlo(ds, key, mod, cIdx):
    from ROOT import TFile, TTree

    # load file and channel list
    fileList = []
    calRuns = cal.GetCalList(key,cIdx)
    for run in calRuns:
        latList = dsi.getSplitList("%s/latSkimDS%d_run%d*" % (dsi.calLatDir, ds, run), run)
        tmpList = [f for idx, f in sorted(latList.items())]
        fileList.extend(tmpList)
    chList = det.getGoodChanList(ds)

    print("Scanning DS:%d  calIdx %d  mod %d  key %s  nFiles:%d" % (ds, cIdx, mod, key, len(fileList)), time.strftime('%X %x %Z'))
    outFile = "%s/eff_%s_c%d.npz" % (dsi.effDir, key, cIdx)
    print("Saving output in:",outFile)

    # declare the output stuff
    evtIdx, evtSumET, evtHitE, evtChans, evtSlo, evtRise, evtToE = [], [], [], [], [], [], []
    thrCal = {ch:[] for ch in chList}
    fLo, fHi, fpb = -200, 400, 1
    nbf = int((fHi-fLo)/fpb)+1
    fSloSpec = {ch:[np.zeros(nbf) for i in range(3)] for ch in chList} # 0-10, 10-200, 236-240

    # loop over LAT cal files
    scanStart = time.time()
    prevRun = 0
    evtCtr, totCtr, totRunTime = 0, 0, 0
    for iF, f in enumerate(fileList):

        print("%d/%d %s" % (iF, len(fileList), f))
        tf = TFile(f)
        tt = tf.Get("skimTree")

        # histogram these for each file (and then add to total)
        fs10 = {ch:[] for ch in chList}
        fs200 = {ch:[] for ch in chList}
        fs238 = {ch:[] for ch in chList}

        # increment the run time and fill the output dict of thresholds
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
            # save them into the output dict (so we can compare w/ DB later).

            n = tt.Draw("channel:threshKeV:threshSigma","","goff")
            chan, thrM, thrS = tt.GetV1(), tt.GetV2(), tt.GetV3()
            tmpThresh = {}
            for i in range(n):
                if chan[i] not in chList:
                    continue
                if chan[i] in tmpThresh.keys():
                    continue
                if thrM[i] < 9999:
                    thrK = thrM[i] + 3*thrS[i]
                    tmpThresh[chan[i]] = [run,thrM[i],thrS[i],thrK]
            for ch in chList:
                if ch not in tmpThresh.keys():
                    tmpThresh[ch] = [-1,-1,-1,-1]

            # fill the output dict
            for ch in tmpThresh:
                thrCal[ch].append(tmpThresh[ch]) # [run, thrM, thrS, thrK]

        prevRun = run
        # continue

        # loop over tree
        for iE in range(tt.GetEntries()):
            tt.GetEntry(iE)
            if tt.EventDC1Bits != 0: continue
            totCtr += 1

            n = tt.channel.size()
            chTmp = np.asarray([tt.channel.at(i) for i in range(n)])
            idxRaw = [i for i in range(tt.channel.size()) if tt.channel.at(i) in chList]
            hitERaw = np.asarray([tt.trapENFCal.at(i) for i in idxRaw])

            # get indexes of hits above threshold (use thresholds from THIS CAL RUN)
            idxList = [i for i in range(tt.channel.size())
                if tt.channel.at(i) in chList
                and tt.trapENFCal.at(i) > tmpThresh[tt.channel.at(i)][3]
                and 0.7 < tt.trapENFCal.at(i) < 9999
                ]
            hitE = np.asarray([tt.trapENFCal.at(i) for i in idxList])

            # calculate mHT and sumET
            mHT, sumET = len(hitE), sum(hitE)

            # save fitSlo data for 0-10 and 10-200 kev ranges for each channel
            for i in idxList:
                en = tt.trapENFCal.at(i)
                ch = tt.channel.at(i)
                fs = tt.fitSlo.at(i)
                if fLo < fs < fHi:
                    if 0 < en < 10:
                        fs10[ch].append(fs)
                    if 10 < en < 200:
                        fs200[ch].append(fs)
                    if 236 < en < 240:
                        fs238[ch].append(fs)

            # Save m2s238 events to output, skip everything else
            if mHT!=2: continue
            if not 237.28 < sumET < 239.46: continue
            hitChans = np.asarray([tt.channel.at(i) for i in idxList])
            hitSlo = np.asarray([tt.fitSlo.at(i) for i in idxList])
            hitRise = np.asarray([tt.riseNoise.at(i) for i in idxList])
            hitkvorrT = np.asarray([tt.kvorrT.at(i) for i in idxList])
            evtIdx.append([run,iE,calIdx])
            evtSumET.append(sumET)
            evtHitE.append(hitE)
            evtChans.append(hitChans)
            evtSlo.append(hitSlo)
            evtRise.append(hitRise)
            evtToE.append(hitkvorrT/hitE)

            evtCtr += 1

        # fill the fitSlo histograms w/ the events from this file
        for ch in chList:
            x, h1 = wl.GetHisto(fs10[ch],fLo,fHi,fpb)
            x, h2 = wl.GetHisto(fs200[ch],fLo,fHi,fpb)
            x, h3 = wl.GetHisto(fs238[ch],fLo,fHi,fpb)
            fSloSpec[ch][0] = np.sum([fSloSpec[ch][0], h1], axis=0)
            fSloSpec[ch][1] = np.sum([fSloSpec[ch][1], h2], axis=0)
            fSloSpec[ch][2] = np.sum([fSloSpec[ch][2], h3], axis=0)
            # n1 = np.sum(fSloSpec[ch][0])
            # n2 = np.sum(fSloSpec[ch][1])
            # n3 = np.sum(fSloSpec[ch][2])
            # print("ch:%d  n10 %d  n200 %d  n238 %d" % (ch, n1, n2, n3))

    # get average threshold for each channel in this file list
    thrFinal = {chan:[] for chan in thrCal}
    for chan in thrCal:
        thrVals = []
        for iT in range(len(thrCal[chan])):
            run, thrM, thrS, thrK = thrCal[chan][iT]
            # print("%d  %d  %.3f  %.3f  %.3f" % (chan,run,thrM,thrS,thrK))
            if thrK > -1:
                thrVals.append(thrK)
        thrVals = np.asarray(thrVals)
        thrAvg = np.mean(thrVals)
        thrDev = np.std(thrVals)
        # print("%d  %.3f  %.3f" % (chan, thrAvg, thrDev))
        thrFinal[chan] = [thrAvg,thrDev]

    # print to screen the final thresholds, stdev, and an error message if necessary
    print("Detector Thresholds:")
    for chan in sorted(thrFinal):
        thKeV = thrFinal[chan][0]
        thE = thrFinal[chan][1]
        errString = ""
        if thE/thKeV > 0.5:
            errString = ">50pct error:  thE/thKeV=%.2f" % (thE/thKeV)
        if thKeV > 2:
            errString = ">2kev"
        print("%d  %.3f  %.3f  %s" % (chan,thKeV,thE,errString))

    # save output
    np.savez(outFile, evtIdx, evtSumET, evtHitE, evtChans, thrCal, thrFinal, evtCtr, totCtr, totRunTime, fSloSpec, x, evtSlo, evtRise, evtToE)

    # output stats
    print("Done:",time.strftime('%X %x %Z'),", %.2f sec/file." % ((time.time()-scanStart)/len(fileList)))
    print("  m2s238 evts:",evtCtr, "total evts:",totCtr, "runTime:",totRunTime)


def loadSloData(key):
    """ Load files generated by scanRuns, return data in a dict.
    To avoid confusion, must specify a key from runsCal.json .
    """
    if key not in cal.GetKeys():
        print("Unknown key!")
        return None
    else:
        print("Loading eff data for key:",key)

    # output dict
    eff = {}
    eff["hitE"] = []  # [hitE1, hitE2 , ...] (remove sub-list of input format)
    eff["chan"] = []  # [chan1, chan2 , ...]
    eff["fSlo"] = []  # [fSlo1, fSlo2, ...]
    eff["rise"] = []  # [rise1, rise2, ...]
    eff["run"]  = []  # [run1, run2, ...]
    eff["cIdx"] = []  # [cIdx1, cIdx2, ...]
    eff["spec"] = []  # array of fitSlo histo dicts (i should have used pandas probably)
    eff["specX"] = [] # x values for "spec" histos (all the same)
    for ci in range(cal.GetIdxs(key)):
        eFile = "%s/eff_%s_c%d.npz" % (dsi.effDir, key, ci)
        if not os.path.isfile(eFile):
            print("File not found:",eFile)
            continue
        f = np.load(eFile)
        evtIdx = f['arr_0']          # m2s238 event [[run,iE,cIdx] , ...]
        evtSumET = f['arr_1']        # m2s238 event [sumET , ...]
        evtHitE = f['arr_2']         # m2s238 event [[hitE1, hitE2] , ...]
        evtChans = f['arr_3']        # m2s238 event [[chan1, chan2] , ...]
        thrCal = f['arr_4'].item()   # {ch : [run,thrM,thrS,thrK] for ch in goodList(ds)}
        thrFinal = f['arr_5'].item() # {ch : [thrAvg, thrDev] for ch in goodList(ds)}
        evtCtr = f['arr_6']          # num m2s238 evts
        totCtr = f['arr_7']          # num total evts
        runTime = f['arr_8']         # cal run time
        fSloSpec = f['arr_9'].item() # fitSlo histos (all hits) {ch:[h10, h200, h238] for ch in chList}
        fSloX = f['arr_10']          # xVals for fitSlo histos
        evtSlo = f['arr_11']         # m2s238 event [[fSlo1, fSlo2], ...]
        evtRise = f['arr_12']        # m2s238 event [[rise1, rise2], ...]

        # remove the hit pair
        for i in range(len(evtHitE)):
            eff["hitE"].extend(evtHitE[i])
            eff["chan"].extend(evtChans[i])
            eff["fSlo"].extend(evtSlo[i])
            eff["rise"].extend(evtRise[i])
            eff["run"].extend([evtIdx[i][0], evtIdx[i][0]])
            eff["cIdx"].extend([evtIdx[i][2], evtIdx[i][2]])
        eff["spec"].append(fSloSpec)
        eff["specX"] = fSloX # this doesn't change

    # convert to numpy arrays and return
    for key in eff:
        if key=="spec": continue
        eff[key] = np.asarray(eff[key])
    return eff


def setSloCut():
    """ ./lat2.py -fs
    Set slow pulse cut for each detector in each DS.
    Uses the combined m2s238 events from ALL datasets, in each detector.
    If makePlots is true, output plots w/ the efficiency functions for each detector.
    Sandbox version: LAT/sandbox/slo-cut.py :: combineDSEff

    NOTE: uses a 90% cut on m2s238 events.  Someone may ask for this to change someday,
          it might be worth making it a parameter at the top of the code.

    DB entries:
    "fitSlo_[calKey]_idx[ci]_m2s238" : {ch : [fsCut, fs200] for ch in chList}}
        and
    "fitSlo_cpd_eff" : {cpd:[fsShiftCut, nBin, amp, sig, mu, ampE, sigE, muE] for cpd in detList}
    """
    from statsmodels.stats import proportion
    if writeDB:
        dbFile = '%s/calDB-v2.json' % (dsi.latSWDir)
        print("Writing results to DB :",dbFile)
        calDB = db.TinyDB(dbFile)
        pars = db.Query()

    makePlots = False
    printTable = False
    writeEffData = False

    dsList = [0,1,2,3,4,5]
    # dsList = [1]
    detList = det.allDets
    detIDs = det.allDetIDs

    # parameter limits and binning
    xAllLo, xAllHi, xpbAll = 0, 250, 0.5        # all hits < 250
    xPassLo, xPassHi, xpbPass = 0, 50, 1        # "low energy" region
    fAllLo, fAllHi, fpbAll = -200, 400, 1       # raw fitSlo vals
    fShiftLo, fShiftHi, fpbShift = -50, 50, 1   # shifted fitSlo

    nbxAll = int((xAllHi-xAllLo)/xpbAll)
    nbxPass = int((xPassHi-xPassLo)/xpbPass)
    nbyShift = int((fShiftHi-fShiftLo)/fpbShift)
    nbyAll = int((fAllHi-fAllLo)/fpbAll) + 1

    # overall hit and efficiency plots (low-E region)
    xTot, hPassAll = wl.GetHisto([], xPassLo, xPassHi, xpbPass, shift=False)
    xTot, hFailAll = wl.GetHisto([], xPassLo, xPassHi, xpbPass, shift=False)
    xTot, hTotAll = wl.GetHisto([], xPassLo, xPassHi, xpbPass, shift=False)

    # ** these two dicts are what we use to fill the DB for calIdx/channel pairs **
    shiftVals = {} # store the UNshifted 10-200 fitSlo peak vals
    shiftCut = {}  # store the 90% vals of the SHIFTED fitSlo peaks

    # shiftSpec = {cpd:np.zeros(nbyShift+1) for cpd in detList} # these are the 10-200 fitSlo hits (unused)

    # save the hits in the "low energy" region
    hitE = {cpd:[] for cpd in detList}
    fSlo = {cpd:[] for cpd in detList} # SHIFTED

    # save the efficiency plot data
    effData = {}

    # load the slowness data we got from scanRunsSlow()
    for ds in dsList:
        for key in cal.GetKeys(ds):

            # get channels in this DS and map back to CPD
            chList = det.getGoodChanList(ds)
            mod = -1
            if "m1" in key:
                mod = 1
                chList = [ch for ch in chList if ch < 1000]
            if "m2" in key:
                mod = 2
                chList = [ch for ch in chList if ch > 1000]
            cpdList = [det.getChanCPD(ds,ch) for ch in chList]
            chMap = {det.getChanCPD(ds,ch):ch for ch in chList}

            eff = loadSloData(key)
            nCal = cal.GetNCalIdxs(ds,mod)
            shiftVals[key] = {ci:None for ci in range(nCal)}

            # loop over calIdx's
            for ci in range(nCal):

                # save the fs shift value for each cpd/ch in each calIdx, and fill the hit lists for plotting
                shiftVals[key][ci] = {cpd:None for cpd in cpdList}

                for cpd in cpdList:

                    # load fitSlo hist of all cal hits in this channel 10-200 keV
                    ch = chMap[cpd]
                    h1 = eff["spec"][ci][ch][1]
                    x1 = eff["specX"]
                    # shiftSpec[cpd] = np.add(shiftSpec[cpd], h1) # unused
                    if np.sum(h1)==0:
                        # print("ci %d  cpd %d  no counts" % (ci, cpd))
                        shiftVals[key][ci][cpd] = -1
                        continue

                    # get mode (maximum) of the 10-200 hits and save it
                    fMax = x1[np.argmax(h1)]
                    shiftVals[key][ci][cpd] = fMax

                    # fill the hit lists
                    idx = np.where((eff["cIdx"]==ci) & (eff["chan"]==ch))
                    fSlo[cpd].extend([f - fMax for f in eff["fSlo"][idx]])
                    hitE[cpd].extend([e for e in eff["hitE"][idx]])

    # values for a bar chart of all det counts vs cts under 10 keV
    ctsAll = {cpd:0 for cpd in detList}
    ctsU10 = {cpd:0 for cpd in detList}

    # it figures
    fig1 = plt.figure(1, figsize=(20,15)) # diagnostic m2s238 plot
    p1 = plt.subplot(221)
    p2 = plt.subplot(222)
    p3 = plt.subplot(223)
    p4 = plt.subplot(224)
    fig2 = plt.figure(2) # hit spectrum plot
    fig3 = plt.figure(3) # efficiency plot
    p31 = plt.subplot2grid((3,1), (0,0), rowspan=2)
    p32 = plt.subplot2grid((3,1), (2,0))

    if printTable: print("CPD  amp  c     loc      scale   e1keV  n10/bin")

    # loop over the combined DS m2s238 populations for each detector and fill shiftCut
    for cpd in detList:

        # if cpd!='114': continue

        hTmp = np.asarray(hitE[cpd])
        idx = np.where(hTmp < 10)
        ctsAll[cpd] = len(hitE[cpd])
        ctsU10[cpd] = len(hTmp[idx])

        # plot fitSlo 1D, calculate the 90% value
        xS, hSlo = wl.GetHisto(fSlo[cpd], fShiftLo, fShiftHi, fpbShift)
        if np.sum(hSlo)==0:
            if printTable: print(cpd)
            continue
        max, avg, std, pct, wid = wl.getHistInfo(xS,hSlo)

        # ***************************
        fsShiftCut = pct[2] # 90% val
        # ***************************

        # zoom on low-E region & fit erf
        hitPass, hitFail = [], []
        for i in range(len(hitE[cpd])):
            if fSlo[cpd][i] <= fsShiftCut: hitPass.append(hitE[cpd][i])
            else: hitFail.append(hitE[cpd][i])
        hitPass, hitFail = np.asarray(hitPass), np.asarray(hitFail)

        # get the average number of counts per bin under 10 kev
        nBin = len(hitPass[np.where(hitPass < 10)])/(10/xpbPass)

        # fill the pass/fail/tot histograms
        xELow, hPass = wl.GetHisto(hitPass, xPassLo, xPassHi, xpbPass, shift=False)
        xELow, hFail = wl.GetHisto(hitFail, xPassLo, xPassHi, xpbPass, shift=False)
        hTot = np.add(hPass, hFail)
        hTotAll = np.add(hTotAll, hTot)
        hPassAll = np.add(hPassAll, hPass)
        hFailAll = np.add(hFailAll, hFail)

        # efficiency for low-E region (xPassLo, xPassHi)

        # start by only looking at points where we have >0 hits passing
        idxP = np.where(hPass > 0)
        sloEff = hPass[idxP] / hTot[idxP]
        xEff = xELow[idxP]

        # calculate confidence intervals for each point (error bars)
        ci_low, ci_upp = proportion.proportion_confint(hPass[idxP], hTot[idxP], alpha=0.1, method='beta')

        # limit the fit energy range
        idxE = np.where((xEff < xPassHi) & (xEff > 1))
        xEff, sloEff, ci_low, ci_upp = xEff[idxE], sloEff[idxE], ci_low[idxE], ci_upp[idxE]

        # ** this is essential, otherwise the blue dots show up on the right side of
        # where the bin center would be, and throw off the fit by half a bin width
        xEff -= xpbPass/2.

        # fit to constrained weibull (c, loc, scale, amp)
        b3 = ((0,-15,0,0),(np.inf,np.inf,np.inf,1.)) # this one is working best
        popt, pcov = curve_fit(wl.weibull, xEff, sloEff, bounds=b3)

        # save pars and errors
        perr = np.sqrt(np.diag(pcov))
        c, loc, sc, amp = popt
        cE, locE, scE, ampE = perr
        eff1 = wl.weibull(1.,*popt)

        if printTable:
            print("%s  %-3.1f  %-4.1f  %-7.2f  %-5.2f  %-3.2f  %d" % (cpd, amp, c, loc, sc, eff1, nBin))

        # fill output dict
        shiftCut[cpd] = [fsShiftCut, nBin, amp, c, loc, sc, ampE, cE, locE, scE]

        # save eff data
        effData[cpd] = [xEff, sloEff, ci_low, ci_upp, hPass, hFail, hTot, xELow]

        # make plots for each detector
        if makePlots:

            # if cpd!='114': continue
            # continue

            # NOTE: we could save the full hitE spec into the efficiency data object if needed.

            plt.figure(1)
            plt.cla()

            # plot fs vs hitE (2D) (xAllLo, xAllHi)
            p1.cla()
            p1.hist2d(hitE[cpd], fSlo[cpd], bins=[nbxAll, nbyShift], range=[[xAllLo,xAllHi],[fShiftLo,fShiftHi]], cmap='jet',norm=LogNorm())
            p1.axhline(fsShiftCut, c='r', lw=3)
            p1.set_xlabel("Energy (keV)", ha='right', x=1)
            p1.set_ylabel("fitSlo", ha='right', y=1)
            p1.set_xlim(xAllLo, xAllHi)
            p1.set_ylim(fShiftLo, fShiftHi)

            # plot SHIFTED fs
            p2.cla()
            p2.plot(hSlo, xS, ls='steps', c='k', label='cpd %s' % cpd)
            p2.axhline(fsShiftCut, c='r', label="90%% value: %.0f" % fsShiftCut)
            p2.set_ylabel("fitSlo", ha='right', x=1)
            p2.set_xlabel("Counts", ha='right', y=1)
            p2.legend(loc=1)
            p2.set_ylim(fShiftLo, fShiftHi)

            # plot hitE (ALL)
            p3.cla()
            x, hHit = wl.GetHisto(hitE[cpd], xAllLo, xAllHi, xpbAll)
            p3.plot(x, hHit, ls='steps', c='b', label="cpd %s" % cpd)
            p3.set_xlabel("Energy (keV)", ha='right', x=1)
            p3.set_ylabel("Counts/%.1f keV" % xpbAll, ha='right', y=1)
            p3.set_xlim(xAllLo, xAllHi)
            p3.legend(loc=1)

            # zoom in on low-e region and plot pass/fail
            p4.cla()
            p4.plot(xELow, hTot, ls='steps', c='k', lw=2., label='all m2s238 hits')
            p4.plot(xELow, hPass, ls='steps', c='b', lw=2., label='cpd %s pass' % cpd)
            p4.plot(xELow, hFail, ls='steps', c='r', lw=2., label='fail')
            p4.axvline(1., c='g', lw=1, label="1 kev")
            p4.set_xlabel("Energy (keV)", ha='right', x=1)
            p4.set_ylabel("Counts/%.1f keV" % xpbPass, ha='right', y=1)
            p4.set_xlim(xPassLo, xPassHi)
            p4.legend(loc=1)

            plt.tight_layout()
            plt.savefig("./plots/lat2-%s.png" % cpd)

            # plot efficiency vs energy.
            plt.figure(3)
            p31.cla()

            p31.plot(xEff, sloEff, '.b', ms=10., label='C%sP%sD%s  nBin %.1f' % (cpd[0],cpd[1],cpd[2], nBin))
            p31.errorbar(xEff, sloEff, yerr=[sloEff - ci_low, ci_upp - sloEff], color='k', linewidth=0.8, fmt='none')

            xFunc = np.arange(xPassLo, xPassHi, 0.1)
            p31.plot(xFunc, wl.weibull(xFunc, *popt), 'g-', label='weibull: c %.1f  loc %.1f  sc %.1f  a %.3f ' % tuple(popt))
            p31.axvline(1.,color='b', lw=1., label='1keV eff: %.2f' % wl.weibull(1.,*popt))

            p31.set_xlim(xPassLo, xPassHi)
            p31.set_ylim(0,1)
            p31.set_xlabel("hitE (keV)", ha='right', x=1)
            p31.set_ylabel("Efficiency", ha='right', y=1)
            p31.legend(loc=4, fontsize=10)

            p32.cla()
            p32.set_xlim(xPassLo, xPassHi)

            hResid = 100*(wl.weibull(xEff, *popt) - sloEff)
            meanRes, stdRes = np.mean(hResid), np.std(hResid)
            p32.axhline(meanRes, c='b', alpha=0.3, label='mean:%.2f%%' % meanRes)
            p32.axhline(meanRes+stdRes, c='m', alpha=0.3, label='std:%.2f%%' % stdRes)
            p32.axhline(meanRes-stdRes, c='m', alpha=0.3)

            p32.plot(xEff, hResid, ".g")
            p32.errorbar(xEff, np.zeros(len(xEff)), yerr=100*np.asarray([sloEff - ci_low, ci_upp - sloEff]), \
                 color='k', linewidth=0.8, fmt='none')
            p32.axvline(1.,color='b', lw=1.)
            p32.set_ylabel("Resid(%)")

            plt.legend(loc=1,ncol=2,fontsize=8)

            plt.tight_layout()
            plt.savefig("./plots/lat2-eff-%s.png" % cpd)

    # make total (all-detector) plots
    if makePlots:
        # plot a bar of all det counts vs counts under 10 keV
        plt.figure(2)
        plt.cla()
        x = np.arange(0,len(detList),1)
        hAll = [ctsAll[cpd] for cpd in detList]
        plt.bar(x, hAll, 0.95, color='b', label='all m2s238 hits')
        hLow = [ctsU10[cpd] for cpd in detList]
        plt.bar(x, hLow, 0.95, color='r', label='m2s238 E<10 keV')
        plt.gca().set_ylim(1)
        plt.gca().set_yscale('log')
        # plt.xlabel("channel", ha='right', x=1.)
        xticks = np.arange(0, len(detList))
        plt.xticks(xticks)
        plt.gca().set_xticklabels(detList, fontsize=8, rotation=90)
        # plt.ylabel("Counts, mHT=2, sumET=238 hits", ha='right', x=1.)
        plt.legend(loc=1)
        plt.savefig("./plots/lat2-totCts.png")

        # plot overall hit spectrum
        plt.figure(2)
        plt.cla()
        plt.plot(xTot, hTotAll, ls='steps', c='k', label="Total Hits")
        plt.plot(xTot, hPassAll, ls='steps', c='b', label="Pass")
        plt.plot(xTot, hFailAll, ls='steps', c='r', label="Fail")
        plt.axvline(1., c='g', lw=1, label='1 kev')
        plt.xlabel("Energy (keV)", ha='right', x=1)
        plt.ylabel("Counts/%.1f keV" % xpbPass, ha='right', y=1)
        plt.xlim(xPassLo, xPassHi)
        plt.legend(loc=1)
        plt.tight_layout()
        plt.savefig("./plots/lat2-totHits.png")


        # plot overall efficiency (xPassLo, xPassHi)

        plt.figure(3)
        p31.cla()

        # start by only looking at points where we have >0 hits passing
        idxP = np.where(hPassAll > 0)
        sloEff = hPassAll[idxP] / hTotAll[idxP]
        xEff = xELow[idxP]

        # calculate confidence intervals for each point (error bars)
        ci_low, ci_upp = proportion.proportion_confint(hPassAll[idxP], hTotAll[idxP], alpha=0.1, method='beta')

        # limit the fit energy range
        idxE = np.where((xEff < xPassHi) & (xEff > 1))
        xEff, sloEff, ci_low, ci_upp = xEff[idxE], sloEff[idxE], ci_low[idxE], ci_upp[idxE]

        # ** this is essential, otherwise the blue dots show up on the right side of
        # where the bin center would be, and throw off the fit by half a bin width
        xEff -= xpbPass/2.

        # fit to constrained weibull (c, loc, scale, amp)
        b3 = ((0,-15,0,0),(np.inf,np.inf,np.inf,1.)) # this one is working best
        popt, pcov = curve_fit(wl.weibull, xEff, sloEff, bounds=b3)

        # save pars and errors
        perr = np.sqrt(np.diag(pcov))
        c, loc, sc, amp = popt
        cE, locE, scE, ampE = perr
        eff1 = wl.weibull(1.,*popt)

        print("eLo %.1f  eHi %.1f  c %.3f (%.3f)  loc %.3f (%.3f)  sc %.3f (%.3f)  amp %.3f (%.3f)  eff1 %.3f" % (xPassLo,xPassHi,c,cE,loc,locE,sc,scE,amp,ampE,eff1))

        # fill the plot
        p31.cla()

        nBin = np.sum(hPassAll[np.where(xELow <= 10)])/(10/xpbPass)

        p31.plot(xEff, sloEff, '.b', ms=10., label='all dets. nBin %.1f' % nBin)
        p31.errorbar(xEff, sloEff, yerr=[sloEff - ci_low, ci_upp - sloEff], color='k', linewidth=0.8, fmt='none')

        xFunc = np.arange(xPassLo, xPassHi, 0.1)
        p31.plot(xFunc, wl.weibull(xFunc, *popt), 'g-', label='weibull: c %.1f  loc %.1f  sc %.1f  a %.3f ' % tuple(popt))
        p31.axvline(1.,color='b', lw=1., label='1keV eff: %.2f' % wl.weibull(1.,*popt))

        p31.set_xlim(xPassLo, xPassHi)
        p31.set_ylim(0,1)
        p31.set_xlabel("hitE (keV)", ha='right', x=1)
        p31.set_ylabel("Efficiency", ha='right', y=1)
        p31.legend(loc=4, fontsize=10)

        p32.cla()
        p32.set_xlim(xPassLo, xPassHi)

        hResid = 100*(wl.weibull(xEff, *popt) - sloEff)
        meanRes, stdRes = np.mean(hResid), np.std(hResid)
        p32.axhline(meanRes, c='b', alpha=0.3, label='mean:%.2f%%' % meanRes)
        p32.axhline(meanRes+stdRes, c='m', alpha=0.3, label='std:%.2f%%' % stdRes)
        p32.axhline(meanRes-stdRes, c='m', alpha=0.3)

        p32.plot(xEff, hResid, ".g")
        p32.errorbar(xEff, np.zeros(len(xEff)), yerr=100*np.asarray([sloEff - ci_low, ci_upp - sloEff]), \
             color='k', linewidth=0.8, fmt='none')
        p32.axvline(1.,color='b', lw=1.)
        p32.set_ylabel("Resid(%)")

        plt.legend(loc=1,ncol=2,fontsize=8)
        plt.tight_layout()
        plt.savefig("./plots/lat2-effTot.png")

    # save efficiency for faster plotting
    if writeEffData:
        np.savez('./data/lat2-eff-data.npz',effData)

    # this is here to emphasize looking at the 'makePlots' output before writing to the DB.
    if not writeDB:
        return

    # --------------------------------------------------------------
    # Detectors to cut.  This is from inspecting the 'makePlots' output.
    # Criteria: nBin must be higher than 4 (50% statistical error in the bins)
    # NOTE: these detectors could be brought back if we included the DS6 cal runs
    # to bump up the stats in M2.
    cutDets = ['211','214','221','261','274'] # 2018/5/11 cgw verified not in DB

    # --------------------------------------------------------------
    # Now that we've filled shiftVals and shiftCut, translate back to datasets and channels

    # reminder:
    # shiftVals = {} # store the UNshifted 10-200 fitSlo peak vals
    # shiftCut = {}  # store the 90% vals of the SHIFTED fitSlo peaks, and the fit results
    #     {cpd:[fsShiftCut, nBin, amp, c, loc, sc, ampE, cE, locE, scE]}

    for key in shiftVals:

        # if "ds5_m2" not in key: continue

        print(key)
        ds = int(key[2])
        chList = det.getGoodChanList(ds)
        if "m1" in key:
            chList = [ch for ch in chList if ch < 1000]
        elif "m2" in key:
            chList = [ch for ch in chList if ch > 1000]
        print(chList)

        for ci in shiftVals[key]:
            dbKey = "fitSlo_%s_idx%d_m2s238" % (key, ci)
            dbVals = {}
            # print("pass 1:")
            # print(dbKey)
            for ch in chList:

                cpd = det.getChanCPD(ds, ch)

                fs200 = shiftVals[key][ci][cpd]  # watch out for "-1", it signifies a bad fit (or no data)
                fsShiftCut, nBin, amp, c, loc, sc, ampE, cE, locE, scE = shiftCut[cpd]

                if fs200 > 0:
                    fsCut = fs200 + fsShiftCut # ****** this is the 90% cut value for this channel ******
                else:
                    fsCut, nBin = -1, -1

                # apply the detector cut
                if cpd in cutDets:
                    dbVals[ch] = None
                    # print(ch, cpd, None)
                    continue

                dbVals[ch] = [fsCut, fs200]

                # print(ch, cpd, wl.niceList(dbVals[ch]))

            # final review
            # print("pass 2:")
            # print(dbKey)
            # for ch in dbVals:
            #     cpd = det.getChanCPD(ds, ch)
            #     if dbVals[ch] is not None:
            #         print(ch, cpd, wl.niceList(dbVals[ch]))
            #     else:
            #         print(ch, cpd, None)

            # fill the DB
            if writeDB:
                dsi.setDBRecord({"key":dbKey, "vals":dbVals}, forceUpdate=True, calDB=calDB, pars=pars)
                print("DB filled:",dbKey)


    # Finally, write the fit function parameters as a separate DB entry for each cpd (all DS's)
    if writeDB:
        dbKey = "fitSlo_cpd_eff"
        dbVals = shiftCut
        dsi.setDBRecord({"key":dbKey, "vals":dbVals}, forceUpdate=True, calDB=calDB, pars=pars)
        print("DB filled:",dbKey)


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
    """ ./lat2.py [-db] -rn

    riseNoise DB entries:
    {"key":"riseNoise_%s_ci%d_pol", "vals": {ch:[a,b,c99,c,fitPass] for ch in goodList} }
    """
    makePlots = False
    if writeDB:
        dbFile = '%s/calDB-v2.json' % (dsi.latSWDir)
        print("Writing results to DB :",dbFile)
        calDB = db.TinyDB(dbFile)
        pars = db.Query()

    dsList = [0,1,2,3,4,5]
    # dsList = [0]

    # loop over ds's, separated by cal key
    for ds in dsList:
        for calKey in cal.GetKeys(ds):
            print("Scanning",calKey,"...")

            # if calKey!="ds5c": continue

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

            # set up 2d histo object for diagnostics
            riseHist = {ci:{} for ci in range(nCal)}

            # loop over calIdx's
            for ci in range(nCal):

                # if ci!=4: continue

                dbKey = "riseNoise_%s_ci%d_pol" % (calKey,ci)
                dbVals = {ch : None for ch in chList}

                # set up 2d histo object for diagnostics
                riseHist[ci] = {ch:None for ch in chList}

                for ch in chList:

                    # print(ds, ci, ch)
                    # if ch!=1204: continue

                    cTmp = eff["chan"][ci]
                    idx = np.where(cTmp==ch)

                    hitE = eff["hitE"][ci][idx]
                    rise = eff["rise"][ci][idx]

                    # make sure we have at least as many hits as parameters,
                    # and that riseNoise vals are good
                    if len(hitE)<10 or len(rise[np.where(rise > 0)])==0:
                        # print("No data, ch",ch)
                        continue

                    # fit the data to a pol1
                    popt, pcov = curve_fit(wl.pol1, hitE, rise)

                    # move the y-intercept up (c), and calculate how many events you keep.
                    # The version in rise-cut.py will make a plot showing the raw fitting (so we don't do it here.)
                    # s1 = time.time()
                    nTot = len(hitE)
                    a, b, c = popt
                    fitPass = False
                    for i in range(1000):
                        c99 = c + 0.1*i
                        nPass = 0
                        for i in range(nTot):
                            if rise[i] <= wl.pol1(hitE[i],a,b,c99):
                                nPass += 1
                        if nPass/nTot > 0.995: # keep 99.5% of cal events
                            fitPass = True
                            break
                    # print("DS%d ci%d ch%d  99pct fit: %.2f sec  a %-9.2e  b %-5.3f  c99 %.3f  Pass:" % (ds,ci,ch,time.time()-s1,a,b,c99),fitPass)

                    dbVals[ch] = [a,b,c99,c,fitPass]

                    # save a 2d hist for this channel/calIdx
                    xLo, xHi, xpb = 0, 250, 2
                    yLo, yHi, ypb = -5, 10, 0.1
                    nbx = int((xHi-xLo)/xpb)
                    nby = int((yHi-yLo)/ypb)

                    # make sure we're not getting too many overflow counts
                    idx = np.where((rise < yLo) | (rise > yHi))
                    if len(idx[0])/len(rise) > 0.05:
                        print("Warning, getting overflow counts, chan %d  nTot %d  nOVF %d  maxRise %.1f  minRise %.1f" % (ch,len(idx[0]),len(rise),rise[idx].max(),rise[idx].min()))

                    hist,_,_ = np.histogram2d(hitE,rise,bins=[nbx,nby], range=[[xLo,xHi],[yLo,yHi]])
                    riseHist[ci][ch] = hist

                    if makePlots:
                        plt.cla()
                        xLo, xHi, xpb = 0, 250, 1
                        cpd = det.getChanCPD(ds,ch)
                        plt.plot(np.nan, np.nan, ".w", label="ch %d, C%sP%sD%s" % (ch, cpd[0],cpd[1],cpd[2]))
                        xFit = np.arange(xLo, xHi, 0.1)
                        plt.plot(xFit, wl.pol1(xFit, *popt), 'r-', label="a %.4f b %.4f c %.4f" % tuple(popt))
                        plt.plot(xFit, wl.pol1(xFit, a,b,c99), 'g-', label="a %.4f b %.4f c99 %.4f" % (a,b,c99))
                        plt.plot(evtPass[:,0], evtPass[:,1], ".b", ms=1, label="pass")
                        plt.plot(evtFail[:,0], evtFail[:,1], ".r", ms=1, label="fail")
                        plt.xlabel("Energy (keV)", ha='right', x=1)
                        plt.ylabel("riseNoise", ha='right', y=1)
                        leg = plt.legend(loc='best', fontsize=12)
                        leg.get_frame().set_alpha(0.5)
                        plt.savefig("./plots/rise-ds%d-ci%d-ch%d.png" % (ds, ci, ch))
                        # return

                # final db check
                if not writeDB: print(dbKey)
                # for ch in sorted(dbVals):
                #     if dbVals[ch] is not None:
                #         print(ch, "%-9.2e  %-9.2e  %.2f " % (*dbVals[ch])
                #     else:
                #         print(ch, None)

                # fill the DB
                if writeDB:
                    dsi.setDBRecord({"key":dbKey, "vals":dbVals}, forceUpdate=True, calDB=calDB, pars=pars)
                    print("DB filled:",dbKey)

            # save 2d histos for this cal key
            histOutput = "./data/lat2-rise-%s.npz" % calKey
            print("Saving hist output:",histOutput)
            np.savez(histOutput, riseHist, calKey)


def riseStability():
    """ ./lat2.py -rs
    Track problem channels in riseNoise and return a (suggested)
    list of channels to cut, along w/ a diagnostic plot. """

    dsList = [0,1,2,3,4,5]
    # dsList = [3]

    makeStabilityPlot = True
    makeChannelPlots = True

    calDB = db.TinyDB('%s/calDB-v2.json' % (dsi.latSWDir))
    pars = db.Query()

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
            nCal = cal.GetNCalIdxs(ds,mod)

            # load DB vals : {calIdx: {ch:[a,b,c99,c,fitPass] for ch in goodList} }}
            dbVals = {}
            for ci in range(nCal):
                dbVals[ci] = dsi.getDBRecord("riseNoise_%s_ci%d_pol" % (calKey,ci),False,calDB,pars)

            # average a,b,c for ALL detectors, all calIdx's
            allA, allB, allC = [], [], []
            for ci in range(nCal):
                for ch in dbVals[ci]:
                    if dbVals[ci][ch] is not None:
                        allA.append(dbVals[ci][ch][0])
                        allB.append(dbVals[ci][ch][1])
                        allC.append(dbVals[ci][ch][2])
            avgA, stdA = np.mean(allA), np.std(allA)
            avgB, stdB = np.mean(allB), np.std(allB)
            avgC, stdC = np.mean(allC), np.std(allC)

            fig = plt.figure(figsize=(18,6))
            p1 = plt.subplot(131)
            p2 = plt.subplot(132)
            p3 = plt.subplot(133)
            cmap = plt.cm.get_cmap('tab20',len(chList)+1)

            # these are the threshold values for a candidate for a bad riseNoise channel
            chk = {'a':[-500,2000], 'b':[0,200], 'c':[50,200]}

            checkList = [] # this is what we fill with bad stuff

            for i, ch in enumerate(chList):
                x = [ci for ci in range(nCal) if dbVals[ci][ch] is not None]
                yA, yB, yC = [], [], []
                for ci in range(nCal):
                    if dbVals[ci][ch] is not None:
                        valA = dbVals[ci][ch][0]/avgA
                        valB = dbVals[ci][ch][1]/avgB
                        valC = dbVals[ci][ch][2]/avgC
                        yA.append(valA)
                        yB.append(valB)
                        yC.append(valC)

                        # save bad vals
                        if not (chk['a'][0] < valA*100 < chk['a'][1] \
                            and chk['b'][0] < valB*100 < chk['b'][1] \
                            and chk['c'][0] < valC*100 < chk['c'][1]):
                            checkList.append([ci,ch])

                        # if fit is bad, always check it
                        if not dbVals[ci][ch][4]:
                            checkList.append([ci,ch])

                if makeStabilityPlot:
                    yA, yB, yC = np.asarray(yA), np.asarray(yB), np.asarray(yC)
                    p1.plot(x, yA*100, ".", ms=10, c=cmap(i), label=ch)
                    p2.plot(x, yB*100, ".", ms=10, c=cmap(i))
                    p3.plot(x, yC*100, ".", ms=10, c=cmap(i))

            if makeStabilityPlot:
                p1.axhline(100, c='g', alpha=0.5, label='avgA %.2e' % avgA)
                p1.axhline(chk['a'][0], c='r', alpha=0.5, label='bad:%d' % chk['a'][0])
                p1.axhline(chk['a'][1], c='r', alpha=0.5, label='bad:%d' % chk['a'][1])
                p1.set_xlabel("calIdx", ha='right', x=1)
                p1.set_ylabel("Pct.Deviation from Avg.", ha='right', y=1)
                p1.legend(loc='best', ncol=4, fontsize=10)

                p2.axhline(100, c='g', alpha=0.5, label='avgB %.2e' % avgB)
                p2.axhline(chk['b'][0], c='r', alpha=0.5, label='bad:%d' % chk['b'][0])
                p2.axhline(chk['b'][1], c='r', alpha=0.5, label='bad:%d' % chk['b'][1])
                p2.set_xlabel("calIdx", ha='right', x=1)
                p2.legend(loc='best', fontsize=10)

                p3.axhline(100, c='g', alpha=0.5, label='avgC %.2f' % avgC)
                p3.axhline(chk['c'][0], c='r', alpha=0.5, label='bad:%d' % chk['c'][0])
                p3.axhline(chk['c'][1], c='r', alpha=0.5, label='bad:%d' % chk['c'][1])
                p3.set_xlabel("calIdx", ha='right', x=1)
                p3.legend(loc='best', fontsize=10)

                plt.tight_layout()
                plt.savefig("./plots/rise-stability-%s.png" % calKey)


            # ============================================================
            # if we have candidates for removal, check this diagnostic plot
            # before you make a final decision.

            if len(checkList)==0:
                print("No bad channels found!")
                return

            if not makeChannelPlots:
                print(checkList)
                continue

            f = np.load("./data/lat2-rise-%s.npz" % calKey)
            riseHist = f['arr_0'].item()

            # match the 2d hists in setRiseCut
            xLo, xHi, xpb = 0, 250, 2
            yLo, yHi, ypb = -5, 10, 0.1
            nbx = int((xHi-xLo)/xpb)
            nby = int((yHi-yLo)/ypb)
            _, xe, ye = np.histogram2d([],[],bins=[nbx,nby], range=[[xLo,xHi],[yLo,yHi]])

            fig = plt.figure(figsize=(18,6))
            p1 = plt.subplot(131)
            p2 = plt.subplot(132)
            p3 = plt.subplot(133)

            # print a good channel, just for comparison
            # checkList = [[0,578]]

            print("Removal candidates:")
            for ci, ch in checkList:
                cpd = det.getChanCPD(ds,ch)
                # print("ci",ci,"ch",ch,"cpd",cpd)

                if riseHist[ci][ch] is None:
                    print("No data found! ci, ch",ci,ch)
                    continue

                hRise = riseHist[ci][ch]

                # -- 1. 2d hist, riseNoise vs hitE, with fit and 99.9% cut shown
                p1.cla()
                p1.set_aspect('auto')
                x, y = np.meshgrid(xe, ye)
                p1.pcolormesh(x, y, hRise.T, norm=LogNorm())
                p1.plot(np.nan, np.nan, '.w', label='cIdx %d ch%d C%sP%sD%s' % (ci, ch, tuple(cpd)))

                a,b,c99,c,fitPass = dbVals[ci][ch]
                if not fitPass:
                    print("warning: this riseNoise fit failed!")
                xFit = np.arange(xLo, xHi, 0.1)
                p1.plot(xFit, wl.pol1(xFit, a, b, c), 'g-', label="a %.1e b %.1e c %.4f" % (a,b,c))
                p1.plot(xFit, wl.pol1(xFit, a, b, c99), 'r-', label="a %.1e b %.1e c99 %.4f" % (a,b,c99))

                p1.set_xlabel("Energy (keV)", ha='right', x=1)
                p1.set_ylabel("riseNoise", ha='right', y=1)

                p1.legend(loc=4)

                # -- 2. 1d projection, energy
                p2.cla()
                yE = np.sum(hRise, axis=1)
                xE = xe[:-1] + 0.5*(xe[1]-xe[0]) # center the bin edges
                p2.plot(xE, yE, "b", ls='steps')
                p2.set_xlabel("Energy (keV)", ha='right', x=1)

                # -- 3. 1d projection, riseNoise. show the E=0 and E=238 cut vals
                p3.cla()
                yR = np.sum(hRise, axis=0)
                xR = ye[:-1] + 0.5*(ye[1]-ye[0])
                p3.plot(xR, yR, "b", ls='steps')

                p3.axvline(wl.pol1(0,a,b,c99), c='r', label='@E=0 : %.2f' % wl.pol1(0,a,b,c99))
                p3.axvline(wl.pol1(238,a,b,c99), c='r', label='@E=238 : %.2f' % wl.pol1(238,a,b,c99))
                p3.set_xlabel("riseNoise", ha='right', x=1)
                p3.legend(loc=1, fontsize=10)

                plt.tight_layout()
                plt.savefig("./plots/rise-%s-ci%d-ch%d.png" % (calKey,ci,ch))

            print("Cut candidates,",calKey)
            print(checkList)


def badRiseChans():
    """ ./lat2.py [-db] -rc
    Using the diagnostic plots from riseStability,
    manually identify channels which fail the fit.
    Update the riseNoise DB entries for these channels/calIdxs to be declared bad.
        (Just change the last parameter to False.)
    This could maybe go in lat-expo.py but we're just gonna turn around
    and use it in applyCuts anyway.
    """
    # the corresponding plots are saved in ./plots/rise/ for reference

    # NOTE: these are calIdx's.
    removeList = {}
    removeList["ds0_m1"] = {
        692:[26,27]   # HF burst.  this calIdx isn't used in DS0 anyway.  OK cgw 2018/06/01
        }
    removeList["ds1_m1"] = {
        594:list(range(29,56+1)),   # 2nd HF population starting @ 50 keV (C1P7D3)
        692:[56]                    # too much curvature
        }
    removeList["ds3_m1"] = {
        594:list(range(0,8+1))      # 2nd HF population starting @ 50 keV (C1P7D3)
        }
    removeList["ds4_m2"] = {
        1106:[1,4,7,8], # HF burst
        1136:[4,7,8],   # "
        1144:[7],       # too much curvature
        1296:[4,7,8],   # HF burst
        1298:[4]        # "
        }
    removeList["ds5_m1"] = {
        584:[7],    # threshold noise causes too much curvature
        608:[7,8],  # "
        632:[7],    # "
        662:[8],    # "
        692:[7,8]   # "
        }
    removeList["ds5_m2"] = {
        1232:[4,5,6,7],     # threshold noise causes too much curvature
        1236:[4,5,6,7,8],   # "
        1298:[4,5,6,7],     # "
        1330:[4,6,7,8],     # "
        1332:[4]            # "
        }

    # load DB vals : {calIdx: {ch:[a,b,c99,c,fitPass] for ch in goodList} }}
    dbFile = '%s/calDB-v2.json' % (dsi.latSWDir)
    calDB = db.TinyDB(dbFile)
    pars = db.Query()
    for calKey in removeList:
        for ch in removeList[calKey]:
            badIdx = removeList[calKey][ch]
            print("Updating DB entry for %s ch %d, badIdxs:" % (calKey,ch),badIdx)

            for bI in badIdx:
                dbKey = "riseNoise_%s_ci%d_pol" % (calKey,bI)
                dbVals = dsi.getDBRecord(dbKey,False,calDB,pars)

                # dbVals[ch][4] = False # this marks the entry bad in the DB

                print(dbVals[ch])
                # return

                if writeDB:
                    # dsi.setDBRecord({"key":dbKey, "vals":dbVals}, forceUpdate=True, dbFile=dbFile, calDB=calDB, pars=pars)
                    dsi.setDBRecord({"key":dbKey, "vals":dbVals}, True, dbFile, calDB, pars)


def applyCuts(ds, cutType):
    """ ./lat2.py -ds [N] -cut [cutType]
    DS values:
        0, 1, 2, 3, 4, 5A, 5B, 5C
    Cut types:
        -cut [blank]: threshold cut only
        -cut fs : threshold + fitSlo
        -cut rn : threshold + riseNoise
        -cut fr : threshold + fitSlo + riseNoise
    """
    from ROOT import gROOT, TFile, TChain, TTree, TNamed, MGTWaveform
    gROOT.ProcessLine("gErrorIgnoreLevel = 3001;")

    # if this is set, don't overwrite good files.
    cleanupMode = False
    print("CLEANUP MODE?",cleanupMode)

    # if this is set, check that we get all files we should. (instead of writing new ones)
    checkMode = False
    print("CHECK MODE?", checkMode)

    # NOTE: input for DS5 must be 5A, 5B, or 5C, not 5.
    dsNum = int(ds[0]) if isinstance(ds, str) else int(ds)
    print("Generating cut files for DS-%s (%d) ..." % (ds, dsNum))

    calDB = db.TinyDB('%s/calDB-v2.json' % (dsi.latSWDir))
    pars = db.Query()
    dsMap = bkg.dsMap() # number of bIdx's
    bkgRanges = bkg.getRanges(ds)
    calKeys = cal.GetKeys(dsNum)

    # have to treat modules separately
    mods = [1]
    if dsNum == 4: mods = [2]
    if dsNum == 5: mods = [1,2]
    for mod in mods:
        chList = det.getGoodChanList(dsNum, mod)

        for bIdx in bkgRanges:
            bkgDict, calDict, _,_ = dsi.GetDBCuts(ds,bIdx,mod,cutType,calDB,pars)

            # # debug block
            # if bIdx!=32:continue
            # if bIdx==32:
            #     bkgDict, calDict = dsi.GetDBCuts(ds,bIdx,mod,cutType,calDB,pars,True)
            #     for ch in sorted(bkgDict):
            #         print(ch, bkgDict[ch])
            #     for ch in calDict:
            #         print(ch, calDict[ch])

            # get the list of LAT files for this bkgIdx
            latList = dsi.getSplitList("%s/latSkimDS%d_%d*" % (dsi.latDir, dsNum, bIdx), bIdx)
            fileList = [f for idx, f in sorted(latList.items())]
            skimTree = TChain("skimTree")
            for f in fileList: skimTree.Add(f)
            # f = TFile(fileList[0])
            # theCut = f.Get("theCut").GetTitle() # don't need to reapply this

            for ch in chList:
                cpd = det.getChanCPD(dsNum,ch)

                # check the TCut string
                thData = True if ch in bkgDict.keys() else False
                fsData = True if ch in calDict.keys() and "fitSlo" in calDict[ch] else False
                rnData = True if ch in calDict.keys() and "riseNoise" in calDict[ch] else False

                # print("DS%d  bIdx %d  ch %d  cpd %s  th %d  fs %d  rn %d" % (dsNum,bIdx,ch,cpd,int(thData),int(fsData),int(rnData)))

                # ** thresholds are REQUIRED **
                if not thData: continue
                chanCut = "channel==%d && %s" % (ch, bkgDict[ch])
                outFile = "%s/th/th_ds%d_%d_ch%d.root" % (dsi.cutDir,dsNum,bIdx,ch)

                if cutType=="fs" and fsData:
                    outFile = "%s/fs/fs_ds%d_%d_ch%d.root" % (dsi.cutDir,dsNum,bIdx,ch)
                    chanCut += "&& fitSlo>0 && %s" % calDict[ch]

                elif cutType=="rn" and rnData:
                    outFile = "%s/rn/rn_ds%d_%d_ch%d.root" % (dsi.cutDir,dsNum,bIdx,ch)
                    chanCut += "&& %s" % calDict[ch]

                elif cutType=="fr" and fsData and rnData:
                    outFile = "%s/fr/fr_ds%d_%d_ch%d.root" % (dsi.cutDir,dsNum,bIdx,ch)
                    chanCut += "&& %s" % calDict[ch]

                else:
                    continue

                # debug: just print cuts
                # if ch in bkgDict.keys():
                #     print(ch, chanCut)
                # else:
                #     print(ch, None)
                # continue

                # do file integrity checks (instead of making new output)
                if checkMode:

                    if not os.path.isfile(outFile):
                        print("File not found:",outFile)
                        print("DS%d  bIdx %d  ch %d  cpd %s  th %d  fs %d  rn %d" % (dsNum,bIdx,ch,cpd,int(thData),int(fsData),int(rnData)))
                        exit(1)

                    cutFile = TFile(outFile)
                    if cutFile.IsZombie():
                        print("Zombieeeah, eeah, eeaah,",outFile)
                        exit(1)
                    elif cutFile.TestBit(TFile.kRecovered):
                        print("Recovered",outFile)
                        exit(1)
                    elif cutFile.GetListOfKeys().GetSize()==0:
                        print("No keys:",outFile)
                        exit(1)

                    tt = cutFile.Get("skimTree")
                    nEnt = tt.GetEntries()
                    n = tt.Draw("trapENFCal:riseNoise:fitSlo","trapENFCal < 250","goff")
                    t1, t2, t3 = tt.GetV1(), tt.GetV2(), tt.GetV3()
                    t1 = np.asarray([t1[i] for i in range(n)])
                    t2 = np.asarray([t2[i] for i in range(n)])
                    print("%s: nEnt %-8d nDraw %-8d" % (outFile.split("/")[-1], nEnt, n))

                    cutFile.Close()
                    continue

                # if a file exists, only recreate it if it's bad.
                makeNew = True
                if cleanupMode:
                    # let's try an idea from a Rene Brun forum post from 2002.  Ugggg
                    if os.path.isfile(outFile):
                        makeNew = False
                        fTmp = TFile(outFile)
                        if fTmp.IsZombie():
                            print("zombie:",outFile)
                            makeNew = True
                        elif fTmp.TestBit(TFile.kRecovered):
                            print("recovered:",outFile)
                            makeNew = True
                        elif fTmp.GetListOfKeys().GetSize()==0:
                            print("no keys:",outFile)
                            makeNew = True
                        fTmp.Close()

                # Wow, such a crazy amount of work to get to this block.
                if makeNew:
                    print("   Writing to:",outFile)
                    print("   Chan",ch,"Cut:",chanCut,"\n")
                    outFile = TFile(outFile,"RECREATE")
                    outTree = TTree()
                    outTree = skimTree.CopyTree(chanCut)
                    outTree.Write()
                    cutUsed = TNamed("chanCut",chanCut)
                    cutUsed.Write()
                    print("Wrote",outTree.GetEntries(),"entries.")
                    outFile.Close()




if __name__=="__main__":
    main(sys.argv[1:])
