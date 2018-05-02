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
from statsmodels.stats import proportion
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

        # wrapper function for scanRuns
        if opt=="-load":
            loadRuns(ds,cIdx,mod)

        # call scanRuns directly (used by job-panda)
        if opt=="-scan":
            ds, key, mod, cIdx = int(argv[i+1]), argv[i+2], int(argv[i+3]), int(argv[i+4])
            scanRuns(ds,key,mod,cIdx)

        # set fitSlo cut
        if opt=="-fs":
            setSloCut()

        # generate cut files (specify cut type)
        if opt=="-cut":
            applyCuts(ds,argv[i+1])


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
                scanRuns(ds, key, mod, cIdx)


def scanRuns(ds, key, mod, cIdx):
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


def loadScanData(key):
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
    If writeDB is set, fills the database with:
        "fitSlo_[calKey]_idx[ci]_m2s238" : {ch : [fsCut, fs200, nBin] for ch in chList}}
        and
        "fitSlo_cpd_eff" : {cpd:[fsShiftCut, nBin, amp, sig, mu, ampE, sigE, muE] for cpd in detList}
    Sandbox version: LAT/sandbox/slo-cut.py :: combineDSEff
    """
    makePlots = False
    printTable = True
    if writeDB:
        dbFile = '%s/calDB-v2.json' % (dsi.latSWDir)
        print("Writing results to DB :",dbFile)
        calDB = db.TinyDB(dbFile)
        pars = db.Query()

    dsList = [0,1,2,3,4,5]
    detList = det.allDets
    detIDs = det.allDetIDs

    yLo, yHi, ypb = -200, 400, 1
    nby = int((yHi-yLo)/ypb)
    shiftSpec = {cpd:np.zeros(nby+1) for cpd in detList} # these are the 10-200 hits
    shiftVals = {} # this stores the 10-200 fitSlo vals
    hitE = {cpd:[] for cpd in detList}
    fSlo = {cpd:[] for cpd in detList}

    # overall hit and efficiency plots
    nTot = 0
    xLo, xHi, xpbE = 0, 30, 0.5
    xE, hPassAll = wl.GetHisto([], xLo, xHi, xpbE)
    xE, hFailAll = wl.GetHisto([], xLo, xHi, xpbE)
    xE, hTotAll = wl.GetHisto([], xLo, xHi, xpbE)

    # loop over multiple ds's, separated by cal key
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

            eff = loadScanData(key)
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
                    shiftSpec[cpd] = np.add(shiftSpec[cpd], h1)
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

    if printTable: print("CPD  amp  sig   e50%  e1keV  n10/bin")

    shiftCut = {}

    # loop over the combined DS m2s238 populations for each detector and fill shiftCut
    for cpd in detList:
        # if cpd!='114': continue

        hTmp = np.asarray(hitE[cpd])
        idx = np.where(hTmp < 10)
        ctsAll[cpd] = len(hitE[cpd])
        ctsU10[cpd] = len(hTmp[idx])

        xLo, xHi, xpb = 0, 250, 1
        nbx = int((xHi-xLo)/xpb)
        fLo, fHi, fpb = -50, 50, 1
        nby = int((fHi-fLo)/fpb)

        # plot fitSlo 1D, calculate the 90% value
        xS, hSlo = wl.GetHisto(fSlo[cpd], fLo, fHi, fpb)
        if np.sum(hSlo)==0:
            if printTable: print(cpd)
            continue
        max, avg, std, pct, wid = wl.getHistInfo(xS,hSlo)

        # ****
        fsShiftCut = pct[2] # 90% val
        # ****

        # zoom on low-E region & fit erf
        hitPass, hitFail = [], []
        for i in range(len(hitE[cpd])):
            if fSlo[cpd][i] <= fsShiftCut: hitPass.append(hitE[cpd][i])
            else: hitFail.append(hitE[cpd][i])

        xLo, xHi, xpb = 0, 30, 0.5
        xE, hPass = wl.GetHisto(hitPass, xLo, xHi, xpb)
        xE, hFail = wl.GetHisto(hitFail, xLo, xHi, xpb)
        hTot = np.add(hPass, hFail)

        nTot += 1
        hTotAll = np.add(hTotAll, hTot)
        hPassAll = np.add(hPassAll, hPass)
        hFailAll = np.add(hFailAll, hFail)

        idx = np.where((hTot > 0) & (hPass > 0))
        sloEff = hPass[idx] / hTot[idx]
        nPad = len(hPass)-len(hPass[idx])
        sloEff = np.pad(sloEff, (nPad,0), 'constant', constant_values=0)
        ci_low, ci_upp = proportion.proportion_confint(hPass[idx], hTot[idx], alpha=0.1, method='beta')
        ci_low = np.pad(ci_low, (nPad,0), 'constant', constant_values=0)
        ci_upp = np.pad(ci_upp, (nPad,0), 'constant', constant_values=0)
        idx2 = np.where(xE > 1.)
        # erf params: mu,sig,amp
        bnd = (0,[np.inf,np.inf,1])
        popt,pcov = curve_fit(wl.logisticFunc, xE[idx2], sloEff[idx2], bounds=bnd)
        perr = np.sqrt(np.diag(pcov))
        mu, sig, amp = popt
        muE, sigE, ampE = perr
        xErf = np.arange(0, xE[-1], 0.1)
        yErf = wl.logisticFunc(xErf, *popt)

        hitPass = np.asarray(hitPass)
        nBin = len(hitPass[np.where(hitPass < 10)])/((10/xpb))
        eff1 = wl.logisticFunc(1.,*popt)

        if printTable:
            print("%s  %-3.1f  %-4.1f  %-4.2f  %-3.2f  %d" % (cpd, amp, sig, mu, eff1, nBin))

        # fill output dict
        shiftCut[cpd] = [fsShiftCut, nBin, amp, sig, mu, ampE, sigE, muE]

        if makePlots:
            # if cpd!='114': continue

            plt.figure(1)
            plt.cla()

            # plot fs vs hitE (2D)
            p1.cla()
            xLo, xHi, xpb = 0, 250, 1
            nbx = int((xHi-xLo)/xpb)
            fLo, fHi, fpb = -50, 50, 1
            nby = int((fHi-fLo)/fpb)
            p1.hist2d(hitE[cpd], fSlo[cpd], bins=[nbx, nby], range=[[xLo,xHi],[fLo,fHi]], cmap='jet',norm=LogNorm())
            p1.axhline(v90, c='r', lw=3)
            p1.set_xlabel("Energy (keV)", ha='right', x=1)
            p1.set_ylabel("fitSlo", ha='right', y=1)

            # plot fs
            p2.cla()
            p2.plot(hSlo, xS, ls='steps', c='k', label='cpd %s' % cpd)
            p2.axhline(v90, c='r', label="90%% value: %.0f" % v90)
            p2.set_xlabel("fitSlo", ha='right', x=1)
            p2.legend(loc=1)

            # plot hitE
            p3.cla()
            x, hHit = wl.GetHisto(hitE[cpd], xLo, xHi, xpb)
            p3.plot(x, hHit, ls='steps', c='b', label="cpd %s" % cpd)
            p3.set_xlabel("Energy (keV)", ha='right', x=1)
            p3.set_ylabel("Counts/%.1f keV" % xpb, ha='right', y=1)
            p3.legend(loc=1)

            # zoom in on low-e region and plot pass/fail
            p4.cla()
            p4.plot(xE, hTot, ls='steps', c='k', lw=2., label='all m2s238 hits')
            p4.plot(xE, hPass, ls='steps', c='b', lw=2., label='cpd %s pass' % cpd)
            p4.plot(xE, hFail, ls='steps', c='r', lw=2., label='fail')
            p4.set_xlabel("Energy (keV)", ha='right', x=1)
            p4.set_ylabel("Counts/%.1f keV" % xpb, ha='right', y=1)
            p4.legend(loc=1)

            # save figure 1
            plt.tight_layout()
            plt.savefig("./plots/lat2-%s.png" % cpd)

            # plot efficiency vs energy.
            plt.cla()
            plt.figure(3)
            p31.cla()
            p31.plot(xE, sloEff, '.b', ms=10., label='C%sP%sD%s' % (cpd[0],cpd[1],cpd[2]))
            p31.errorbar(xE, sloEff, yerr=[sloEff - ci_low, ci_upp - sloEff], color='k', linewidth=0.8, fmt='none')
            p31.plot(xErf, yErf, 'r-', label="m %.1f s %.2f a %.2f" % tuple(popt))
            p31.axvline(1.,color='g',label='1keV eff: %.2f' % wl.logisticFunc(1.,*popt))
            p31.plot(np.nan, np.nan, 'w', label='nBin %d' % nBin)
            p31.set_xlabel("hitE (keV)", ha='right', x=1)
            p31.set_ylabel("Efficiency", ha='right', y=1)
            p31.legend(loc=4)

            p32.cla()
            hResid = wl.logisticFunc(xE, *popt) - sloEff
            p32.plot(xE, hResid, ".b")
            p32.errorbar(xE, hResid, yerr=[sloEff - ci_low, ci_upp - sloEff], color='k', linewidth=0.8, fmt='none')

            plt.tight_layout()
            plt.savefig("./plots/lat2-eff-%s.png" % cpd)

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
        plt.plot(xE, hTotAll, ls='steps', c='k', label="Total Hits")
        plt.plot(xE, hPassAll, ls='steps', c='b', label="Pass")
        plt.plot(xE, hFailAll, ls='steps', c='r', label="Fail")
        plt.xlabel("Energy (keV)", ha='right', x=1)
        plt.ylabel("Counts/%.1f keV" % xpbE, ha='right', y=1)
        plt.legend(loc=1)
        plt.tight_layout()
        plt.savefig("./plots/lat2-totHits.png")

        # plot overall efficiency
        # NOTE: fitting one logistic to all detectors doesn't fit very well below 2 keV.
        # Better to use the individual detector efficiency functions.
        eLow = 2.

        plt.figure(3)
        p31.cla()

        idx = np.where((hTotAll > 0) & (hPassAll > 0))
        sloEff = hPassAll[idx] / hTotAll[idx]
        nPad = len(hPassAll)-len(hPassAll[idx])
        sloEff = np.pad(sloEff, (nPad,0), 'constant', constant_values=0)
        ci_low, ci_upp = proportion.proportion_confint(hPassAll[idx], hTotAll[idx], alpha=0.1, method='beta')
        ci_low = np.pad(ci_low, (nPad,0), 'constant', constant_values=0)
        ci_upp = np.pad(ci_upp, (nPad,0), 'constant', constant_values=0)
        idx = np.where(xE > eLow)
        bnd = (0,[np.inf,np.inf,1])
        popt,pcov = curve_fit(wl.logisticFunc, xE[idx], sloEff[idx], bounds=bnd)

        p31.plot(xE, sloEff, '.b', ms=10., label='Efficiency')
        p31.errorbar(xE, sloEff, yerr=[sloEff - ci_low, ci_upp - sloEff], color='k', linewidth=0.8, fmt='none')
        xErf = np.arange(0, xE[-1], 0.1)
        p31.plot(xErf, wl.logisticFunc(xErf, *popt), 'r-', label="m %.1f s %.2f a %.2f" % tuple(popt))
        p31.axvline(eLow,color='g',label='%.1fkeV eff: %.2f' % (eLow, wl.logisticFunc(1.,*popt)))
        p31.set_xlabel("Energy (keV)", ha='right', x=1)
        p31.set_ylabel("Efficiency", ha='right', y=1)
        p31.legend(loc=4)

        p32.cla()
        hResid = wl.logisticFunc(xE, *popt) - sloEff
        p32.plot(xE, hResid, ".b")
        p32.errorbar(xE, hResid, yerr=[sloEff - ci_low, ci_upp - sloEff], color='k', linewidth=0.8, fmt='none')

        plt.tight_layout()
        plt.savefig("./plots/lat2-effTot.png")

        # print error of overall logistic fit
        perr = np.sqrt(np.diag(pcov))
        mu, sig, amp = popt
        muE, sigE, ampE = perr
        print("eLow %.1f  mu %.3f pm %3f  sig %.3f pm %.3f  amp %.3f pm %.3f" % (eLow, mu,muE,sig,sigE,amp,ampE))

    # --------------------------------------------------------------
    # Now that we've filled shiftVals, translate it back to datasets and channels
    for key in shiftVals:
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
            # print(dbKey)
            for ch in chList:
                cpd = det.getChanCPD(ds, ch)
                fs200 = shiftVals[key][ci][cpd]  # watch out for "-1", it signifies a bad fit (or no data)
                if fs200 > 0:
                    fsCut = fs200 + shiftCut[cpd][0]
                    nBin = shiftCut[cpd][1]
                else:
                    fsCut, nBin = -1, -1
                dbVals[ch] = [fsCut, fs200, nBin]
                # print(cpd, ch, fs200, fsCut, nBin)

            # final review
            # print(dbKey)
            # for ch in dbVals:
                # print(ch, dbVals[ch])

            # fill the DB
            if writeDB:
                print(dbKey)
                dsi.setDBRecord({"key":dbKey, "vals":dbVals}, forceUpdate=True, calDB=calDB, pars=pars)
                print("DB filled.")

    # Finally, write the erf function parameters as a separate DB entry for each cpd (all DS's)
    if writeDB:
        dbKey = "fitSlo_cpd_eff"
        dbVals = shiftCut
        print(dbKey)
        dsi.setDBRecord({"key":dbKey, "vals":dbVals}, forceUpdate=True, calDB=calDB, pars=pars)
        print("DB filled.")


def applyCuts(ds, cutType):
    """ ./lat2.py -ds [N] -cut [cutType]"""

    # from ROOT import gROOT
    # gROOT.ProcessLine("gErrorIgnoreLevel = 3001;")

    # NOTE: input for DS5 must be 5A, 5B, or 5C, not 5.
    dsNum = int(ds[0]) if isinstance(ds, str) else int(ds)

    calDB = db.TinyDB('%s/calDB-v2.json' % (dsi.latSWDir))
    pars = db.Query()
    dsMap = bkg.dsMap() # number of sub-ranges
    bkgRanges = bkg.getRanges(ds)
    calKeys = cal.GetKeys(dsNum)

    mods = [1]
    if dsNum == 4: mods = [2]
    if dsNum == 5: mods = [1,2]

    for mod in mods:
        chList = det.getGoodChanList(dsNum, mod)
        # print("DS",ds,"Module",mod,"chans:",chList)

        fileList = []
        for sub in bkgRanges:

            # latList = dsi.getSplitList("%s/latSkimDS%d_%d*" % (dsi.latDir, dsNum, sub), sub)
            # tmpList = [f for idx, f in sorted(latList.items())]
            # fileList.extend(tmpList)

            firstRun, lastRun = bkgRanges[sub][0], bkgRanges[sub][-1]

            calKey = "ds%d_m%d" % (dsNum, mod)
            if ds == "5C": calKey = "ds5c"
            if calKey not in calKeys:
                print("Error: Unknown cal key:",calKey)
                return
            cIdxLo = cal.GetCalIdx(calKey, firstRun)
            cIdxHi = cal.GetCalIdx(calKey, lastRun)
            print(ds,mod,sub,firstRun,cIdxLo,lastRun,cIdxHi)

            # create a dict of cuts for each channel, covering each calIdx within each bkgIdx
            cutDict = {}
            for cIdx in range(cIdxLo, cIdxHi+1):

                runCovMin = cal.master[calKey][cIdx][1]
                runCovMax = cal.master[calKey][cIdx][2]
                runLo = firstRun if runCovMin < firstRun else runCovMin
                runHi = lastRun if lastRun < runCovMax else runCovMax
                runCut = "run>=%d && run<=%d" % (runLo, runHi)
                print("  ci",cIdx,runCut)

                fsD = dsi.getDBRecord("fitSlo_%s_idx%d_m2s238" % (calKey, cIdx), False, calDB, pars)
                fsCut = None

                for ch in sorted(fsD):

                    # TODO: thresholds

                    # fitSlo: check the 90% value is positive
                    if fsD[ch][2] > 0:
                        fsCut = "fitSlo<%.2f" % fsD[ch][2]

                    # TODO: riseNoise

                    # set the combination channel cut
                    if cutType == "fs" and fsCut!=None:
                        chanCut = "(%s && %s)" % (runCut, fsCut)

                    # create dict entry for this channel or append to existing, taking care of parentheses and OR's.
                    if ch in cutDict.keys() and chanCut!=None:
                        cutDict[ch] += " || %s" % chanCut
                    elif ch not in cutDict.keys() and chanCut!=None:
                        cutDict[ch] = "(%s" % chanCut

            for ch in cutDict:
                cutDict[ch] += ")" # close the parens for each channel entry

            # -- finally, loop over each channel we have an entry for, get its cut, and create an output file. --
            for ch in sorted(cutDict):

                theCut = ""
                chanCut = theCut + "&& channel==%d " % ch

                if cutType == "fs":
                    # outFile = "~/project/cuts/%sfs/%sfitSlo-DS%d-%d-ch%d.root" % (dString, dString, dsNum, bkgIdx, ch)
                    chanCut += "&& fitSlo>0 && %s" % cutDict[ch]

                print(ch, chanCut)

            # f = TFile(file0)
            # theCut = f.Get("theCut").GetTitle()






if __name__=="__main__":
    main(sys.argv[1:])
