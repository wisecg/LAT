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
            ds = int(argv[i+1])
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

        # set fitSlo cut for a ds
        if opt=="-fs":
            setSloCut(ds)

    getCombinedEff()


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
    evtIdx, evtSumET, evtHitE, evtChans, evtSlo, evtRise = [], [], [], [], [], []
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
            evtIdx.append([run,iE,calIdx])
            evtSumET.append(sumET)
            evtHitE.append(hitE)
            evtChans.append(hitChans)
            evtSlo.append(hitSlo)
            evtRise.append(hitRise)

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
    np.savez(outFile, evtIdx, evtSumET, evtHitE, evtChans, thrCal, thrFinal, evtCtr, totCtr, totRunTime, fSloSpec, x, evtSlo, evtRise)

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


def getCombinedEff():
    """ Identical to sandbox/slo-cut.py::getCombinedEff, except it doesn't make plots.
        Goes by CPD instead of channel number, since CPD never changes through the DS's.
    """
    detList = det.allDets
    detIDs = det.allDetIDs

    yLo, yHi, ypb = -200, 400, 1
    nby = int((yHi-yLo)/ypb)
    shiftSpec = {cpd:np.zeros(nby+1) for cpd in detList} # these are the 10-200 hits (unlike the function above)
    shiftVals = {} # this stores the 10-200 fitSlo vals
    hitE = {cpd:[] for cpd in detList}
    fSlo = {cpd:[] for cpd in detList}

    # loop over scanRuns output from ALL DS's
    for ds in [0,1,2,3,4,5]:
    # for ds in [1]:
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

    print("CPD  amp  sig   e50%  e1keV  n10/bin")

    for cpd in detList:
        # if cpd!='114': continue

        xLo, xHi, xpb = 0, 250, 1
        nbx = int((xHi-xLo)/xpb)
        fLo, fHi, fpb = -50, 50, 1
        nby = int((fHi-fLo)/fpb)

        # plot fitSlo 1D, calculate the 90% value
        x, hSlo = wl.GetHisto(fSlo[cpd], fLo, fHi, fpb)
        if np.sum(hSlo)==0:
            print(cpd)
            continue
        max, avg, std, pct, wid = wl.getHistInfo(x,hSlo)
        v90 = pct[2]

        # zoom on low-E region & fit erf
        hitPass, hitFail = [], []
        for i in range(len(hitE[cpd])):
            if fSlo[cpd][i] <= v90: hitPass.append(hitE[cpd][i])
            else: hitFail.append(hitE[cpd][i])

        xLo, xHi, xpb = 0, 30, 0.5
        x, hPass = wl.GetHisto(hitPass, xLo, xHi, xpb)
        x, hFail = wl.GetHisto(hitFail, xLo, xHi, xpb)
        hTot = np.add(hPass, hFail)

        idx = np.where((hTot > 0) & (hPass > 0))
        sloEff = hPass[idx] / hTot[idx]
        nPad = len(hPass)-len(hPass[idx])
        sloEff = np.pad(sloEff, (nPad,0), 'constant', constant_values=0)
        ci_low, ci_upp = proportion.proportion_confint(hPass[idx], hTot[idx], alpha=0.1, method='beta')
        ci_low = np.pad(ci_low, (nPad,0), 'constant', constant_values=0)
        ci_upp = np.pad(ci_upp, (nPad,0), 'constant', constant_values=0)

        idx = np.where(x > 0.5)
        bnd = (0,[np.inf,np.inf,1])
        popt,pcov = curve_fit(wl.threshFunc, x[idx], sloEff[idx], bounds=bnd)
        perr = np.sqrt(np.diag(pcov))
        mu, sig, amp = popt

        hitPass = np.asarray(hitPass)
        nBin = len(hitPass[np.where(hitPass < 10)])/((10/xpb))
        eff1 = wl.threshFunc(1.,*popt)

        # identify the cutoff point where the erf is outside the error bar

        print("%s  %-3.1f  %-4.1f  %-4.2f  %-3.2f  %d" % (cpd, amp, sig, mu, eff1, nBin))


def setSloCut(ds):
    """ Set slow pulse cut for each detector in each DS.
    Current method uses the combined m2s238 events from ALL datasets.
    {"key":"fitSlo_[calKey]_idx[ci]_m2s238","value":{ch:[] for ch in chList}}
    """

    if writeDB:
        print("Writing results to DB ...")
    calDB = db.TinyDB('%s/calDB.json' % (dsi.latSWDir))
    pars = db.Query()

    # treat each cal key separately
    for key in cal.GetKeys(ds):

        chList = det.getGoodChanList(ds)
        mod = -1
        if "m1" in key:
            mod = 1
            chList = [ch for ch in chList if ch < 1000]
        if "m2" in key:
            mod = 2
            chList = [ch for ch in chList if ch > 1000]

        eff = loadScanData(key)
        nCal = cal.GetNCalIdxs(ds,mod)

        shiftDict = {ci:None for ci in range(nCal)}
        yLo, yHi, ypb = -50, 50, 1
        nby = int((yHi-yLo)/ypb)
        shiftSpec = {ch:np.zeros(nby+1) for ch in chList}

        # find the fs shift value for each ch in each calIdx
        for ci in range(nCal):
            shiftDict[ci] = {ch:[] for ch in chList}

            for ch in chList:

                # load fitSlo hist of ALL cal hits in this channel 10-200 keV
                h1 = eff["spec"][ci][ch][1]
                x1 = eff["specX"]

                if np.sum(h1)==0:
                    print("ci %d  ch %d  no counts" % (ci, ch))
                    shiftDict[ci][ch].extend([-1,-1,-1])
                    continue

                # get mode (maximum) of the 10-200 hits
                b = np.argmax(h1)
                fMax = x1[b] # maximum x (fitSlo) value

                # get the 50% width of the 10-200 hits
                bLo, bHi = -1, -1
                for j in range(len(h1)):
                    if bLo==-1 and h1[j] >= h1[b]/2.:
                        bLo = j
                    if j > b and h1[j] <= h1[b]/2.:
                        bHi = j
                        break
                fLo, fHi = x1[bLo], x1[bHi]
                fWid = fHi - fLo

                # save the max and width of h1
                shiftDict[ci][ch].extend([fMax, fWid])
                # print("ci %-2d  ch %d  h10-200 nCts %-7d  max %-4d  fLo %-4d  fHi %-4d  wid %-4d" % (ci,ch,np.sum(h1),fMax,fLo,fHi,fWid))

                # shift the m2s238 events
                idx = np.where((eff["cIdx"]==ci) & (eff["chan"]==ch))
                thisFS = [f - fMax for f in eff["fSlo"][idx]] #
                if len(thisFS) == 0: continue

                # fill the shifted m2s238 fitSlo histograms
                x2, hSlo = wl.GetHisto(thisFS, yLo, yHi, ypb)
                shiftSpec[ch] = np.add(shiftSpec[ch], hSlo)
                # print("ci %-2d  ch %d  m2s238  nCts %-7d" % (ci,ch,np.sum(hSlo)))

        # now find the 90% value from the shifted m2s238 fs histograms
        for ch in chList:
            max, avg, std, pct, wid = wl.getHistInfo(x2,shiftSpec[ch])
            pct90 = pct[2]
            for ci in shiftDict:
                shiftDict[ci][ch].extend([pct90])

        # find the unshifted value for each ch, each calIdx and fill the DB
        for ci in shiftDict:

            dbKey = "fitSlo_%s_idx%d_m2s238" % (key, ci)
            dbVals = {}

            for ch in chList:
                max = shiftDict[ci][ch][0] # fs max
                v90 = shiftDict[ci][ch][2] # m2s238 90% val
                cut90 = shiftDict[ci][ch][0] + shiftDict[ci][ch][2]
                # print("ds %d  cIdx %d  ch%d  max %-4d  v90 %-3d  cut90 %d" % (ds, ci, ch, max, v90, cut90))
                dbVals[ch] = [cut90, max, v90] # watch out for bad values < 0

            # final review
            print(dbKey)
            # for key in dbVals:
                # print(key, dbVals[key])

            # fill the DB
            if writeDB:
                dsi.setDBRecord({"key":dbKey, "vals":dbVals}, forceUpdate=True, calDB=calDB, pars=pars)
                print("DB filled.")


if __name__=="__main__":
    main(sys.argv[1:])