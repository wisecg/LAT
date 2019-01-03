#!/usr/bin/env python3
"""
===================== lat-expo.py ======================
Apply PSA and burst cuts from LAT2 and LAT3, to generate
final spectrum files.
Using output of ds_livetime, calculate the final exposure
and generate combined spectrum efficiency curves for
each DS.
=================== C. Wiseman (USC) ===================
"""
import sys, os, time
import tinydb as db
import numpy as np
import pandas as pd

# LAT libraries
import dsi
bkg = dsi.BkgInfo()
cal = dsi.CalInfo()
det = dsi.DetInfo()
import waveLibs as wl


def main(argv):

    # a very important parameter, use 90 or 95
    global pctTot
    # pctTot = 90
    pctTot = 95
    print("Using pctTot ==",pctTot)


    # NOTE: the options outline a rough 'procedure' to use this code.
    for i, opt in enumerate(argv):

        # check DB cuts
        if opt=="-cov":
            # manual
            getPSACutRuns(argv[i+1],argv[i+2], verbose=True) # ds, cutType

            # batch <-- use this one
            # for ds in [0,1,2,3,4,"5A","5B","5C",6]:
                # getPSACutRuns(ds,"fr")

        # finally, calculate exposure!
        if opt=="-xp":
            getExposure()

        # make final output files
        if opt=="-f":
            makeFinalFiles()

        # make waveform movies
        if opt=="-wfm":
            makeMovies()

        # get efficiency functions
        if opt == "-eff":
            getEfficiency()

        # get efficiency functions in ROOT histogram
        if opt == "-root":
            getEfficiencyROOT()


def getPSACutRuns(ds, cutType, verbose=False):
    """ ./lat-expo.py -cov [ds] [cutType]
    Check which bkgIdx's we have good cut values for,
    accounting for multiple bkg and cal sub-Idx's.
    This is a 5-layered loop, which is pretty badass.

    ** For BOTH of these, the DB has already been updated to reflect these changes **

    ========= fitSlo detectors cut:  (lat2::setSloCut) =========
    # Detectors to cut.  This is from inspecting the 'makePlots' output.
    # Criteria: nBin must be higher than 4 (50% statistical error in the bins)
    # NOTE: these detectors could probably be brought back if we included the
    # DS6 cal runs to bump up the stats in M2.
    cutDets = {
        90: ['211','214','221','254','261','274'],
        95: ['211','214','221','261','274']
    }
    fitSlo can ALSO be bad in a calIdx when fsCut and nBin are == -1.
    dbVals[ch] = [fsCut, fs200, nBin]

    ========= riseNoise chan/calIdx cut: (lat2::badRiseChans) =========
    # the corresponding plots are saved in ./plots/rise/ for reference
    removeList = {}
    removeList["ds0_m1"] = {
        692:[26,27]   # HF burst
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
    """
    # NOTE: input for DS5 must be 5A, 5B, or 5C, not 5.
    dsNum = int(ds[0]) if isinstance(ds, str) else int(ds)
    print("Getting PSA cut run/ch vals for DS-%s (%d) ..." % (ds, dsNum))

    debugMode = False
    if debugMode:
        print("DEBUG MODE.  Not writing output file!")

    if cutType not in ["fr","fs","rn","thr"]:
        print("Unknown cut type:",cutType,"... exiting ...")
        return

    # set up output file (deprecated)
    # writeFile = True
    # dsLabel = str(dsNum)
    # if isinstance(ds, str):
    #     if ds=="5A": dsLabel = "5a"
    #     if ds=="5B": dsLabel = "5b"
    #     if ds=="5C": dsLabel = "5c"
    # outFile = "./data/dbCut_%s_%s.txt" % (cutType,dsLabel)

    dRanges = {} # this is the list of ch/runs to exclude

    calDB = db.TinyDB('%s/calDB-v2.json' % (dsi.latSWDir))
    pars = db.Query()
    dsMap = bkg.dsMap() # number of bIdx's
    bkgRanges = bkg.getRanges(ds)

    # 1. loop over modules
    mods = [1]
    if dsNum == 4: mods = [2]
    if dsNum >= 5: mods = [1,2]
    for mod in mods:

        calKey = "ds%d_m%d" % (dsNum, mod)
        if ds == "5C": calKey = "ds5c"
        if calKey not in cal.GetKeys(dsNum):
            print("Error: Unknown cal key:",calKey)
            return

        chList = det.getGoodChanList(dsNum, mod)
        for ch in chList:
            dRanges[ch] = []

        # 2. loop over bkgIdx
        for i, bIdx in enumerate(bkgRanges):

            # bkgDict, calDict are only filled when we have good entries, bkgCov, calCov are always filled.
            bkgDict, calDict, bkgCov, calCov = dsi.GetDBCuts(ds,bIdx,mod,cutType,calDB,pars,pctTot,False)

            rFirst, rLast = bkgRanges[bIdx][0], bkgRanges[bIdx][-1]
            dsSub = ds if ds in ["5A","5B","5C"] else int(ds)
            subRanges = bkg.GetSubRanges(dsSub, bIdx)
            if len(subRanges) == 0: subRanges.append((rFirst, rLast))

            # 3. loop over sub-bkgIdx
            for sbIdx, (runLo, runHi) in enumerate(subRanges):

                # 4. loop over cIdx's in this sub-bkgIdx
                cIdxLo, cIdxHi = cal.GetCalIdx(calKey, runLo), cal.GetCalIdx(calKey, runHi)
                for i, cIdx in enumerate(range(cIdxLo, cIdxHi+1)):

                    # get the run coverage of this sub-sub-bkgIdx
                    if cIdxLo==cIdxHi:
                        covLo, covHi = runLo, runHi
                    else:
                        runList = bkg.getRunList(ds, bIdx)
                        subList = [r for r in runList if runLo <= r <= runHi and cal.GetCalIdx(calKey,r) == cIdx]
                        if len(subList)==0:
                            if verbose:
                                print("No good runs in this sub-sub-bkgIdx")
                            continue
                        covLo, covHi = subList[0], subList[-1]

                    # 5. loop over channels
                    for ch in chList:

                        goodThr = True if bkgCov[ch][sbIdx] else False
                        goodSlo = True if calCov[ch][0][i+1] else False
                        goodRise = True if calCov[ch][1][i+1] else False

                        exclude = False

                        if cutType == "fr" and not (goodThr and goodSlo and goodRise): exclude = True
                        elif cutType == "fs" and not (goodThr and goodSlo): exclude = True
                        elif cutType == "rn" and not (goodThr and goodRise): exclude = True
                        elif cutType == "thr" and not (goodThr): exclude = True
                        excludeMsg = ""
                        if exclude is True:
                            excludeMsg = " exclude, runs %d - %d" % (covLo, covHi)
                            dRanges[ch].extend([covLo, covHi])

                        if verbose:
                            if ch==692: # it seems that this channel always failed auto-thresh in bkg runs
                                print("%s  bIdx %d  sbIdx %d  cIdx %d  ch %d  cpd %d  th %d  fs %d  rn %d  %s" % (calKey,bIdx,sbIdx,cIdx,ch,int(det.getChanCPD(dsNum,ch)), int(goodThr),int(goodSlo),int(goodRise),excludeMsg))

    # Create output suitable for ds_livetime (deprecated)
    # if writeFile:
    #     with open(outFile, 'w') as f:
    #         for ch in sorted(dRanges):
    #             if len(dRanges[ch]) > 0:
    #                 outStr = "%d" % ch
    #                 for r in wl.niceList(dRanges[ch], "%d", "i"): outStr += " %d" % r
    #                 outStr += "\n"
    #                 f.write(outStr)
    #     print("Wrote cut file:",outFile)

    # Create numpy output for getExposure
    if not debugMode:
        np.savez('./data/lat-psa%dRunCut-ds%s.npz' % (pctTot, ds), dRanges)


def getExposure():
    """./lat-expo.py -xp
    NOTE: pctTotal has to be changed in lat3 to match its setting here.
    TODO: Add uncertainty from active mass?
    """
    import lat3
    from ROOT import TFile, TTree

    # be careful, you have to change pctTotal in lat3 to the right value
    enrExc, natExc, _, _ = lat3.getOutliers(False, usePass2=False)
    # format of enrExc, natExc: [:,0]=dsNum, [:,1]=cpd, [:,2]=bkgIdx

    cutType = "fr"
    burstType = "frb"

    # dsList = [0,1,2,3,4,"5A","5B","5C",6]
    # dsList = [1,2,3,4,"5A","5B","5C"]
    # dsList = [1,2,3,4,"5B"]
    dsList = [0]

    # output
    dsExpo = {} # {ds: [dsEnrExp, dsNatExp]}
    dsUnc = {} # {ds: [dsEnrUnc, dsNatUnc]}
    detExpo = {ds:{cpd:0 for cpd in det.allDets} for ds in dsList}

    writeRuns = False # Flag for writting run-lists that detectors are active
    bExclude253 = True # Flag for excluding C2P5D3 from exposure calculations cuz it sucks
    if bExclude253:
        print("Excluding C2P5D3 because it sucks!")

    # store all vals so we can add the uncertainties in quadrature
    rawTot, psaTot, burstTot, expTot = {}, {}, {}, {}
    grandTotEnrVals, grandTotNatVals = [], []

    for ds in dsList:

        f = np.load('./data/lat-psa%dRunCut-ds%s.npz' % (pctTot, ds))
        psaRuns = f['arr_0'].item() # {ch: [runLo1, runHi1, runLo2, runHi2, ...]}

        dsNum = int(ds[0]) if isinstance(ds,str) else ds
        nBkg = bkg.dsMap()[dsNum]
        bLo, bHi = 0, nBkg
        if ds=="5A": bLo, bHi = 0, 79
        if ds=="5B": bLo, bHi = 80, 112
        if ds=="5C": bLo, bHi = 113, 121
        runRanges = bkg.getRanges(dsNum)
        chList = det.getGoodChanList(dsNum)

        # build burst cut
        dsTmp = ds
        if ds=="5A": dsTmp=50
        if ds=="5B": dsTmp=51
        if ds=="5C": dsTmp=52
        iE = np.where(enrExc[:,0]==dsTmp)
        iN = np.where(natExc[:,0]==dsTmp)
        skipList = np.vstack((enrExc[iE], natExc[iN]))

        # load ds_livetime output
        tl = TFile("./data/ds_%s_livetime.root" % str(ds))
        lt = tl.Get("dsTree")

        # totals for this DS
        rawTot[ds] = {ch:[] for ch in chList}
        psaTot[ds] = {ch:[] for ch in chList}
        burstTot[ds] = {ch:[] for ch in chList}
        expTot[ds] = {ch:[] for ch in chList} # after cuts
        runTot = {det.getChanCPD(dsNum,ch):[] for ch in chList} # for run lists

        # loop over bIdx's
        for bIdx in range(bLo, bHi+1):

            rLo, rHi = runRanges[bIdx][0], runRanges[bIdx][-1]

            psaCutRuns = {ch:[] for ch in chList}
            for ch in chList:
                if len(psaRuns[ch]) > 0:
                    for i in range(0,len(psaRuns[ch]),2):
                        psaCutRuns[ch].extend([r for r in range(psaRuns[ch][i],psaRuns[ch][i+1]+1) if rLo <= r <= rHi])

            burstCutRuns = {ch:False for ch in chList}
            for ch in chList:
                cpd = det.getChanCPD(dsNum,ch)
                iSkip = np.where((skipList == (dsTmp, int(cpd), bIdx)).all(axis=1))
                if len(iSkip[0]) > 0:
                    burstCutRuns[ch] = True

            n = lt.Draw("run:channel:livetime","run>=%d && run<=%d" % (rLo, rHi), 'goff')
            ltRun, ltChan, ltLive = lt.GetV1(), lt.GetV2(), lt.GetV3()
            for i in range(n):
                ch = ltChan[i]
                cpd = det.getChanCPD(dsNum,ch)
                detID = det.getDetIDChan(dsNum,ch)

                aMass = det.allActiveMasses[detID]
                aMassUnc = det.allActiveMassesUnc[detID]

                expo = ltLive[i]*aMass/86400/1000

                rawTot[ds][ch].append(expo)

                if ltRun[i] in psaCutRuns[ch]:
                    psaTot[ds][ch].append(expo)
                    continue

                if burstCutRuns[ltChan[i]] is True:
                    burstTot[ds][ch].append(expo)
                    continue

                # Flag for skipping C2P5D3 here
                if bExclude253 and int(cpd) == 253: continue
                expTot[ds][ch].append(expo)
                runTot[cpd].append((int(ltRun[i]),expo))

        # get values for this DS

        rawEnrExp, rawEnrUnc = 0, 0
        psaEnrExp, psaEnrUnc = 0, 0
        burstEnrExp, burstEnrUnc = 0, 0
        finalEnrExp, finalEnrUnc = 0, 0
        rawNatExp, rawNatUnc = 0, 0
        psaNatExp, psaNatUnc = 0, 0
        burstNatExp, burstNatUnc = 0, 0
        finalNatExp, finalNatUnc = 0, 0

        for ch in chList:
            isEnr = True if det.getDetIDChan(dsNum,ch) > 100000 else False
            detID = det.getDetIDChan(dsNum,ch)
            aMass = det.allActiveMasses[detID]
            aMassUnc = det.allActiveMassesUnc[detID]


            re = np.sum(rawTot[ds][ch])
            pe = np.sum(psaTot[ds][ch])
            be = np.sum(burstTot[ds][ch])
            ee = np.sum(expTot[ds][ch])
            if isEnr:
                rawEnrExp += re
                rawEnrUnc += re * (aMassUnc/aMass)
                psaEnrExp += pe
                psaEnrUnc += pe * (aMassUnc/aMass)
                burstEnrExp += be
                burstEnrUnc += be * (aMassUnc/aMass)

                # This block is only for rejecting C2P5D3
                # Consider 253 to be rejected by PSA cut -- Set the burst cut for this detector to 0
                if bExclude253 and int(cpd)== 253:
                    psaEnrExp += expTot[ds][ch]
                    psaTot[ds][ch] += expTot[ds][ch]
                    burstEnrExp -= burstTot[ds][ch]
                    burstTot[ds][ch] -= burstTot[ds][ch]

                finalEnrExp += ee
                finalEnrUnc += ee * (aMassUnc/aMass)

            else:
                rawNatExp += re
                rawNatUnc += re * (aMassUnc/aMass)
                psaNatExp += pe
                psaNatUnc += pe * (aMassUnc/aMass)
                burstNatExp += be
                burstNatUnc += be * (aMassUnc/aMass)
                finalNatExp += ee
                finalNatUnc += ee * (aMassUnc/aMass)

            detExpo[ds][cpd] += np.sum(expTot[ds][ch]) - np.sum(psaTot[ds][ch]) - np.sum(burstTot[ds][ch])
            print('DS{} Ch{} C{}P{}D{} - Exposure: {exp}'.format(ds, ch, *str(cpd), exp=detExpo[ds][cpd]))

        print("DS-%s" % ds)
        print("Enriched (kg-d)")

        # ugly latex version
        print("DS%s & %.2f $\pm$ %.2f & %.2f $\pm$ %.2f & %.2f $\pm$ %.2f & %.2f $\pm$ %.2f \\\\" % (str(ds), rawEnrExp,rawEnrUnc, psaEnrExp,psaEnrUnc, burstEnrExp,burstEnrUnc, finalEnrExp,finalEnrUnc))

        print("Natural (kg-d)")
        print("DS%s & %.2f $\pm$ %.2f & %.2f $\pm$ %.2f & %.2f $\pm$ %.2f & %.2f $\pm$ %.2f \\\\" % (str(ds), rawNatExp,rawNatUnc, psaNatExp,psaNatUnc, burstNatExp,burstNatUnc, finalNatExp,finalNatUnc))

        # pretty terminal version
        # print("Raw: %.3f ± %.3f" % (rawEnrExp, rawEnrUnc))
        # print("PSA: %.3f ± %.3f" % (psaEnrExp, psaEnrUnc))
        # print("Burst: %.3f ± %.3f" % (burstEnrExp, burstEnrUnc))
        # print("Final: %.3f ± %.3f" % (finalEnrExp, finalEnrUnc))
        # print("Natural (kg-d)")
        # print("Raw: %.3f ± %.3f" % (rawNatExp, rawNatUnc))
        # print("PSA: %.3f ± %.3f" % (psaNatExp, psaNatUnc))
        # print("Burst: %.3f ± %.3f" % (burstNatExp, burstNatUnc))
        # print("Final: %.3f ± %.3f" % (finalNatExp, finalNatUnc))

        grandTotEnrVals.append([finalEnrExp, finalEnrUnc])
        grandTotNatVals.append([finalNatExp, finalNatUnc])

        dsExpo[ds] = [finalEnrExp, finalNatExp]
        dsUnc[ds] = [finalEnrUnc, finalNatUnc]

        # Write exposure lists
        if writeRuns:
            for k in runTot:
                with open('./data/expLists/DS{}_C{}P{}D{}.txt'.format(ds, *str(k)), 'w') as f:
                    for (r,j) in runTot[k]:
                        f.write('{},{}\n'.format(r,j))

    grandTotEnr = np.sum([v[0] for v in grandTotEnrVals]) / 365.25
    grandTotEnrUnc = np.sqrt(np.sum([v[1]**2 for v in grandTotEnrVals])) / 365.25

    grandTotNat = np.sum([v[0] for v in grandTotNatVals]) / 365.25
    grandTotNatUnc = np.sqrt(np.sum([v[1]**2 for v in grandTotNatVals])) / 365.25

    print("\nTotals for DS:",dsList)
    print("Enriched (kg-y): %.3f ± %.3f" % (grandTotEnr, grandTotEnrUnc))
    print("Natural  (kg-y): %.3f ± %.3f" % (grandTotNat, grandTotNatUnc))

    # save output
    np.savez("./data/expo-totals-e%d.npz" % (pctTot), dsExpo, dsUnc, detExpo)


def makeFinalFiles():
    """ ./lat-expo.py -f
    TChain the 'frb(pctTot)' files for each dataset and make output files for spec-fit and others.
    Save the enriched & natural exposure into the files too
    """
    from ROOT import TChain, TFile, TTree, TNamed, MGTWaveform

    # here is where you might put some last minute cuts
    tCut = "tOffset < 4000"

    cutType = "frb%d" % pctTot
    # outType = "final%d" % pctTot
    outType = "final%dt" % pctTot

    dsList = [0,1,2,3,4,"5A","5B","5C",6]
    # dsList = [1,2,3,4,"5A","5B","5C"]
    # dsList = [1,2,3,4,"5B"]
    # dsList = [0]

    # maybe some detectors have noise features that the PSA and burst cuts aren't able to remove.
    # finalDetCut = [
    #     [0,'152'],    # enr, noise wall at ~2.8 keV
    #     [0,'141'],    # nat, noise wall at ~2.1 keV
    #     ['5A','253'], # enr, noise wall at ~5.2 keV
    # ]

    for ds in dsList:

        outName = "%s/bkg/cut/%s/%s_DS%s.root" % (dsi.dataDir, outType, outType, ds)
        print("Writing final LAT output:",outName)

        dummyTree = TChain("skimTree")

        dsNum = int(ds[0]) if isinstance(ds,str) else ds
        nBkg = bkg.dsMap()[dsNum]
        bLo, bHi = 0, nBkg
        if ds=="5A": bLo, bHi = 0, 79
        if ds=="5B": bLo, bHi = 80, 112
        if ds=="5C": bLo, bHi = 113, 121
        chList = det.getGoodChanList(dsNum)

        for bIdx in range(bLo, bHi+1):
            for ch in sorted(chList):

                fName = "%s/bkg/cut/%s/%s_ds%d_%d_ch%d.root" % (dsi.dataDir, cutType, cutType, dsNum, bIdx, ch)
                if not os.path.isfile(fName): continue

                # cpd = det.getChanCPD(dsNum,ch)
                # if [ds,cpd] in finalDetCut: continue

                dummyTree.Add(fName)


        outFile = TFile(outName, "RECREATE")
        outTree = TTree()
        outTree = dummyTree.CopyTree(tCut)
        outTree.Write()

        f = np.load("./data/expo-totals-e%d.npz" % pctTot)
        dsExpo = f['arr_0'].item()
        enrExpTot = dsExpo[ds][0]
        natExpTot = dsExpo[ds][1]
        enrExp = TNamed("enrExp (kg-d)","%.4f" % enrExpTot)
        natExp = TNamed("natExp (kg-d)","%.4f" % natExpTot)
        enrExp.Write()
        natExp.Write()

        outFile.Close()


def makeMovies():
    """ ./lat-expo.py -wfm
    Makes .mp4 movies of waveforms in our final spectra.
    Default is only to save under 10 keV.
    Requires ffmpeg.  (brew install ffmpeg)
    """
    from ROOT import TFile, TTree, MGTWaveform
    import matplotlib.pyplot as plt
    plt.style.use('./clint.mpl')
    from matplotlib import animation
    import pywt

    dsList = [0,1,2,3,4,"5A","5B","5C",6]
    wfLimit = None

    # dsList = [0]
    # wfLimit = 500

    for ds in dsList:

        outFile = "./plots/final%d-wfs-DS%s.mp4" % (pctTot, ds)
        print("Saving WFs for DS-%s: %s" % (ds, outFile))

        inFile = "%s/bkg/cut/final%dt/final%dt_DS%s.root" % (dsi.dataDir, pctTot, pctTot, ds)
        tf = TFile(inFile)
        tt = tf.Get("skimTree")

        # make an entry list
        tCut = "trapENFCal >= 1 && trapENFCal <= 10"

        # simple entry list
        # n = tt.Draw("Entry$:Iteration$",tCut,"goff")
        # evt, itr = tt.GetV1(), tt.GetV2()
        # evtList = [[int(evt[i]),int(itr[i])] for i in range(n)]

        # slightly fancy entry list, sorted by energy
        n = tt.Draw("Entry$:Iteration$:trapENFCal",tCut,"goff")
        evt, itr, ene = tt.GetV1(), tt.GetV2(), tt.GetV3()
        ene = [ene[i] for i in range(n)]
        eHit = np.argsort(ene) # sort by ascending energy
        # eHit = np.flip(np.argsort(ene), axis=0) # sort by descending energy
        evtList = [[int(evt[i]),int(itr[i])] for i in eHit]

        # set up the figure, the axis, and the plot element we want to animate
        fig = plt.figure(figsize=(10,6))
        nFont = 14
        a1 = plt.subplot(111)
        a1.set_xlabel("Time (ns)", ha='right', x=1, fontsize=nFont)
        a1.set_ylabel("Voltage (ADC)", ha='right', y=1, fontsize=nFont)

        p1, = a1.plot(np.ones(1), np.ones(1), c='b')
        p2, = a1.plot(np.ones(1), np.ones(1), c='k',alpha=0.6)

        # initialization: plot the background of each frame
        def init():
            p1.set_data([],[])
            return p1,

        # this is the loop over entries
        def animate(i):

            iE = evtList[i][0]
            iH = evtList[i][1]
            tt.GetEntry(iE)

            run = tt.run
            chan = tt.channel.at(iH)
            hitE = tt.trapENFCal.at(iH)
            fSlo = tt.fitSlo.at(iH)
            rise = tt.riseNoise.at(iH)

            wf = tt.MGTWaveforms.at(iH)
            truncLo, truncHi = 0, 2
            if ds==6 or ds==2: truncLo = 4
            signal = wl.processWaveform(wf,truncLo,truncHi)

            waveBLSub = signal.GetWaveBLSub()
            waveTS = signal.GetTS()
            p1.set_ydata(waveBLSub)
            p1.set_xdata(waveTS)

            # wavelet denoised
            wp = pywt.WaveletPacket(waveBLSub, 'db2', 'symmetric', maxlevel=4)
            new_wp = pywt.WaveletPacket(data=None, wavelet='db2', mode='symmetric')
            new_wp['aaa'] = wp['aaa'].data
            waveDenoised = new_wp.reconstruct(update=False)
            diff = len(waveDenoised) - len(waveBLSub)
            if diff > 0: waveDenoised = waveDenoised[diff:]
            p2.set_ydata(waveDenoised)
            p2.set_xdata(waveTS)

            plt.title("Run %d  Ch %d  iE %d/%d  trapENFCal %.1f  fitSlo %.1f  riseNoise %.1f" % (run,chan,i,len(evtList),hitE,fSlo,rise), fontsize=nFont)

            # dynamically scale the axes
            xmin, xmax = np.amin(waveTS), np.amax(waveTS)
            a1.set_xlim([xmin,xmax])

            ymin, ymax = np.amin(waveDenoised), np.amax(waveDenoised)
            yLo = -7 if ymin-abs(0.1 * ymin) < 7 else ymin-abs(0.1 * ymin)
            yHi = ymax + abs(0.5+ymax)
            a1.set_ylim([yLo,yHi])

            # print("%d / %d  Run %d  nCh %d  chan %d  trapE %.1f  samp %d" % (i,nList,run,nChans,chan,energy,wf.GetLength()))
            if i % 500 == 0 and i != 0:
                print("%d / %d entries saved (%.2f %% done)." % (i,len(evtList),100*(float(i)/len(evtList))))
            return p1,

        # make the movie
        nWF = wfLimit if wfLimit is not None else len(evtList)
        anim = animation.FuncAnimation(fig, animate, init_func=init, frames=nWF, interval=0, blit=False)
        anim.save(outFile, fps=20)#, extra_args=['-vcodec', 'libx264'])


def getEfficiency():
    """ ./lat-expo.py -eff
    Same structure as getPSACutRuns, looping over sub-sub-bkgIdx's.
    Wow, it's a 6-layer loop.  Can I get a degree now?
    """
    import lat3
    from ROOT import TFile, TTree
    import matplotlib.pyplot as plt
    plt.style.use('./clint.mpl')

    dsList = [0,1,2,3,4,"5A","5B","5C",6]  # default
    # dsList = ["5A"]
    # dsList = [3]

    # mode = "trig"  # trigger efficiency only
    # mode = "slo"   # slowness efficiency only
    mode = "all"   # does all PS <-- use this one
    npzOut = "./data/lat-expo-efficiency-%s-e%d.npz" % (mode, pctTot)

    # pandas output for cosmo calculation
    pndOut = "./data/lat-runTimes.h5"
    pndDict = {"ds":[], "run":[], "det":[], "start":[], "stop":[], "lt":[], "rt":[], "expo":[]}

    calDB = db.TinyDB('%s/calDB-v2.json' % (dsi.latSWDir))
    pars = db.Query()
    enrExc, natExc, _, _ = lat3.getOutliers(verbose=False, usePass2=False)

    bExclude253 = True # Flag for excluding C2P5D3 from exposure calculations cuz it sucks
    if bExclude253:
        print("Excluding C2P5D3 because it sucks!")

    debugMode = False
    if debugMode:
        print("WARNING: debug mode.  Not saving output. dsList:",dsList)

    # efficiency output
    # xLo, xHi = 0, 50
    xLo, xHi = 0, 200 # Higher range
    xEff = np.arange(xLo, xHi, 0.01)
    totEnrEff = {ds:np.zeros(len(xEff)) for ds in dsList}
    totEnrEffLo = {ds:np.zeros(len(xEff)) for ds in dsList}
    totEnrEffHi = {ds:np.zeros(len(xEff)) for ds in dsList}
    totNatEff = {ds:np.zeros(len(xEff)) for ds in dsList}
    totNatEffLo = {ds:np.zeros(len(xEff)) for ds in dsList}
    totNatEffHi = {ds:np.zeros(len(xEff)) for ds in dsList}
    detEff = {ds:{} for ds in dsList}
    detEffLo = {ds:{} for ds in dsList}
    detEffHi = {ds:{} for ds in dsList}
    detExpo = {ds:{} for ds in dsList}
    for ds in dsList:
        detEff[ds] = {cpd:np.zeros(len(xEff)) for cpd in det.allDets}
        detEffLo[ds] = {cpd:np.zeros(len(xEff)) for cpd in det.allDets}
        detEffHi[ds] = {cpd:np.zeros(len(xEff)) for cpd in det.allDets}
        detExpo[ds] = {cpd:0 for cpd in det.allDets}

    # recalculate these, make sure they match getExposure
    enrExp = {ds:0 for ds in dsList}
    natExp = {ds:0 for ds in dsList}

    # 1. loop over dataset
    for ds in dsList:

        # set DS stuff
        dsNum = int(ds[0]) if isinstance(ds,str) else ds
        nBkg = bkg.dsMap()[dsNum]
        bLo, bHi = 0, nBkg
        if ds=="5A": bLo, bHi = 0, 79
        if ds=="5B": bLo, bHi = 80, 112
        if ds=="5C": bLo, bHi = 113, 121
        bkgRanges = bkg.getRanges(ds)

        # get psa cut runs and detector fitSlo efficiencies
        f = np.load('./data/lat-psa%dRunCut-ds%s.npz' % (pctTot,ds))
        psaRuns = f['arr_0'].item() # {ch: [runLo1, runHi1, runLo2, runHi2, ...]}

        fsD = dsi.getDBRecord("fitSlo_cpd_eff%d" % pctTot, False, calDB, pars)
        fsU = dsi.getDBRecord("fitSlo_cpd_effHi%d" % pctTot, False, calDB, pars) # upper
        fsL = dsi.getDBRecord("fitSlo_cpd_effLo%d" % pctTot, False, calDB, pars) # lower

        # get burst cut
        dsTmp = ds
        if ds=="5A": dsTmp=50
        if ds=="5B": dsTmp=51
        if ds=="5C": dsTmp=52
        iE = np.where(enrExc[:,0]==dsTmp)
        iN = np.where(natExc[:,0]==dsTmp)
        skipList = np.vstack((enrExc[iE], natExc[iN]))
        # print(skipList)

        # load ds_livetime output
        # tl = TFile("./data/ds_%s_livetime.root" % str(ds))
        tl = TFile("./data/ds_%s_output.root" % str(ds)) # these have extra info
        lt = tl.Get("dsTree")

        # 2. loop over modules
        mods = [1]
        if dsNum == 4: mods = [2]
        if dsNum >= 5: mods = [1,2]
        for mod in mods:

            calKey = "ds%d_m%d" % (dsNum, mod)
            if ds == "5C": calKey = "ds5c"
            if calKey not in cal.GetKeys(dsNum):
                print("Error: Unknown cal key:",calKey)
                return

            print("Scanning DS-%s, m%d ..." % (ds, mod))

            chList = det.getGoodChanList(dsNum, mod)

            # save total efficiency for each channel in this DS
            totEff = {ch:np.zeros(len(xEff)) for ch in chList}
            totEffLo = {ch:np.zeros(len(xEff)) for ch in chList}
            totEffHi = {ch:np.zeros(len(xEff)) for ch in chList}
            trigEff = {ch:np.zeros(len(xEff)) for ch in chList}
            fSloEff = {ch:np.zeros(len(xEff)) for ch in chList}

            # 3. loop over bkgIdx
            for i, bIdx in enumerate(bkgRanges):

                # load bkg (trigger) and cal (PSA) cut coverage
                _,_, bkgCov, calCov = dsi.GetDBCuts(ds,bIdx,mod,"fr",calDB,pars,pctTot,False)

                rLo, rHi = bkgRanges[bIdx][0], bkgRanges[bIdx][-1]

                # get psa cut runs
                psaCutRuns = {ch:[] for ch in chList}
                for ch in chList:
                    if len(psaRuns[ch]) > 0:
                        for i in range(0,len(psaRuns[ch]),2):
                            psaCutRuns[ch].extend([r for r in range(psaRuns[ch][i],psaRuns[ch][i+1]+1) if rLo <= r <= rHi])

                # get burst cut runs
                burstCutRuns = {ch:False for ch in chList}
                for ch in chList:
                    cpd = det.getChanCPD(dsNum,ch)
                    iSkip = np.where((skipList == (dsTmp, int(cpd), bIdx)).all(axis=1))
                    if len(iSkip[0]) > 0:
                        burstCutRuns[ch] = True

                # 4. loop over sub-bIdx
                subRanges = bkg.GetSubRanges(ds, bIdx)
                if len(subRanges) == 0: subRanges.append((rLo, rHi))
                for sbIdx, (subLo, subHi) in enumerate(subRanges):

                    # load trigger efficiencies
                    key = "thresh_ds%d_bkg%d_sub%d" % (dsNum, bIdx, sbIdx)
                    thD = dsi.getDBRecord(key, False, calDB, pars)

                    # 5. loop over cIdx's in this sub-bIdx
                    cIdxLo, cIdxHi = cal.GetCalIdx(calKey, subLo), cal.GetCalIdx(calKey, subHi)
                    for i, cIdx in enumerate(range(cIdxLo, cIdxHi+1)):

                        # get the run coverage of this sub-sub-bIdx
                        if cIdxLo==cIdxHi:
                            covLo, covHi = subLo, subHi
                        else:
                            runList = bkg.getRunList(ds, bIdx)
                            subList = [r for r in runList if subLo <= r <= subHi and cal.GetCalIdx(calKey,r) == cIdx]
                            if len(subList)==0: continue
                            covLo, covHi = subList[0], subList[-1]

                        # calculate exposure for this sub-sub-bIdx
                        subExpo = {ch:0 for ch in chList}
                        subExpoLo = {ch:0 for ch in chList} # save upper and lower limits due to AM uncertainty
                        subExpoHi = {ch:0 for ch in chList}
                        n = lt.Draw("run:channel:livetime","run>=%d && run<=%d" % (covLo, covHi), 'goff')
                        ltRun, ltChan, ltLive = lt.GetV1(), lt.GetV2(), lt.GetV3()
                        ltRun = [ltRun[j] for j in range(n)]
                        ltChan = [ltChan[j] for j in range(n)]
                        ltLive = [ltLive[j] for j in range(n)]

                        nf = lt.Draw("runtime:unixStart:unixStop","run>=%d && run<=%d" % (covLo, covHi), 'goff')
                        ltRT, ltStart, ltStop = lt.GetV1(), lt.GetV2(), lt.GetV3()
                        ltRT = [ltRT[j] for j in range(n)]
                        ltStart = [ltStart[j] for j in range(n)]
                        ltStop = [ltStop[j] for j in range(n)]

                        for j in range(n):
                            ch = ltChan[j]
                            detID = det.getDetIDChan(dsNum,ch)
                            aMass = det.allActiveMasses[detID]
                            expo = ltLive[j]*aMass/86400/1000

                            aMassUnc = det.allActiveMassesUnc[detID]
                            expoUnc = ltLive[j]*aMassUnc/86400/1000
                            expoLo = expo - expoUnc
                            expoHi = expo + expoUnc

                            # print("chan %d  expo %.3f  expoUnc %.3f  (%.3f pct)" % (ch, expo, expoUnc, 100*expoUnc/expo))

                            # since we're splitting by module, ignore channels in the other module
                            # print(ds, mod, int(cpd[0]))
                            # exit()
                            if ch < 1000 and mod!=1: continue
                            if ch > 1000 and mod!=2: continue

                            if ltRun[j] in psaCutRuns[ch]:
                                continue
                            if burstCutRuns[ltChan[j]] is True:
                                continue

                            # skip all runs removed by low-e analysis
                            subExpo[ch] += expo
                            subExpoLo[ch] += expoLo
                            subExpoHi[ch] += expoHi

                            # skip C2P5D3
                            cpd = det.getChanCPD(dsNum,ch)
                            if bExclude253 and int(cpd) == 253: continue
                            detExpo[ds][cpd] += expo

                            if detID > 100000:
                                enrExp[ds] += expo
                            else:
                                natExp[ds] += expo

                            # fill pandas output
                            pndDict["ds"].append(dsTmp)
                            pndDict["run"].append(ltRun[j])
                            pndDict["det"].append(int(cpd))
                            pndDict["start"].append(ltStart[j])
                            pndDict["stop"].append(ltStop[j])
                            pndDict["lt"].append(ltLive[j])
                            pndDict["rt"].append(ltRT[j])
                            pndDict["expo"].append(expo)


                        # 6. loop over channels
                        for ch in chList:

                            goodThr = True if bkgCov[ch][sbIdx] else False
                            goodSlo = True if calCov[ch][0][i+1] else False
                            goodRise = True if calCov[ch][1][i+1] else False
                            if not (goodThr and goodSlo and goodRise):
                                continue

                            # finally, get trigger, fitSlo, and riseNoise efficiencies and scale by exposure

                            # trigger
                            mu, sig, isGood = thD[ch]
                            if isGood != 0:
                                print("error, bad threshold, ch",ch)
                                exit(1)
                            effThresh = mu + 3*sig
                            idx = np.where(xEff >= effThresh)
                            nPad = len(xEff) - len(xEff[idx])
                            tEff = wl.erFunc(xEff[idx],mu,sig,1)
                            tEff = np.pad(tEff, (nPad,0), 'constant')
                            tEffLo, tEffHi = tEff, tEff

                            tEff = np.multiply(tEff, subExpo[ch])
                            tEffLo = np.multiply(tEffLo, subExpoLo[ch])
                            tEffHi = np.multiply(tEffHi, subExpoHi[ch])

                            trigEff[ch] += tEff

                            # fitSlo & riseNoise
                            cpd = int(det.getChanCPD(dsNum,ch))
                            c, loc, scale, amp = fsD[cpd][3], fsD[cpd][4], fsD[cpd][5], fsD[cpd][2]
                            fEff = wl.weibull(xEff,c,loc,scale,amp)

                            # lower and upper bounds for fitSlo efficiency, multiplied by riseNoise efficiency
                            fEffHi = wl.weibull(xEff,fsU[cpd][3], fsU[cpd][4], fsU[cpd][5], fsU[cpd][2]) * 0.995
                            fEffLo = wl.weibull(xEff,fsL[cpd][3], fsL[cpd][4], fsL[cpd][5], fsL[cpd][2]) * 0.995

                            riseEff = 0.995 # riseNoise is defined to be 99.5% efficient, no energy dependence
                            fEff = np.multiply(fEff, riseEff)

                            fSloEff[ch] += fEff

                            # total efficiency
                            if mode == "trig":
                                totEff[ch] += tEff
                            elif mode == "slo":
                                totEff[ch] += fEff
                                totEffHi[ch] += fEffHi
                                totEffLo[ch] += fEffLo
                            elif mode == "all":
                                totEff[ch] += np.multiply(tEff, fEff) # this is what we want

                                # these curves are just the fitSlo uncertainty
                                # totEffHi[ch] += np.multiply(tEff, fEffHi)
                                # totEffLo[ch] += np.multiply(tEff, fEffLo)

                                # these curves are the outer limits of fitSlo + active mass uncertainty << use these
                                totEffHi[ch] += np.multiply(tEffHi, fEffHi)
                                totEffLo[ch] += np.multiply(tEffLo, fEffLo)
                            else:
                                print("Unknown mode! exiting ...")
                                exit()

            # if debugMode:
            #     for ch in chList[1:2]:
            #         # plt.plot(xEff, trigEff[ch], '-')
            #         # plt.plot(xEff, fSloEff[ch], '-')
            #         plt.plot(xEff, totEff[ch], '-r',label='Best, %s' % det.getChanCPD(ds,ch))
            #         plt.plot(xEff, totEffLo[ch], '-g',label="Lo")
            #         plt.plot(xEff, totEffHi[ch], '-b',label='Hi')
            #     plt.legend()
            #     plt.xlim(0,10)
            #     plt.show()
            #     exit()

            expTot = 0
            for ch in chList:
                cpd = det.getChanCPD(dsNum,ch)
                detEff[ds][cpd] = totEff[ch]
                detEffLo[ds][cpd] = totEffLo[ch]
                detEffHi[ds][cpd] = totEffHi[ch]
                # print(ch, cpd, detEff[ds][cpd][500:520], totEff[ch][500:520])

                expTot += detExpo[ds][cpd]
                # print(ds, cpd, detExpo[ds][cpd])
            # print(expTot)

            # plt.plot(xEff,detEff[ds]['164'])
            # plt.show()
            # exit()

        # ========= Done w/ 5 layer loop.  whew! =======

        # get total enr/nat efficiency for this DS
        for cpd in detEff[ds]:

            if det.isEnr(cpd):
                totEnrEff[ds] += detEff[ds][cpd]
                totEnrEffLo[ds] += detEffLo[ds][cpd]
                totEnrEffHi[ds] += detEffHi[ds][cpd]
            else:
                totNatEff[ds] += detEff[ds][cpd]
                totNatEffLo[ds] += detEffLo[ds][cpd]
                totNatEffHi[ds] += detEffHi[ds][cpd]

    # done w/ loop over datasets.
    finalEnrEff = np.zeros(len(xEff))
    finalEnrEffLo = np.zeros(len(xEff))
    finalEnrEffHi = np.zeros(len(xEff))
    finalNatEff = np.zeros(len(xEff))
    finalNatEffLo = np.zeros(len(xEff))
    finalNatEffHi = np.zeros(len(xEff))

    # plot overall enriched/natural efficiency
    finalEnrExp, finalNatExp = 0, 0
    finalEnrEffLo
    for ds in dsList:
        finalEnrEff += totEnrEff[ds]
        finalEnrEffLo += totEnrEffLo[ds]
        finalEnrEffHi += totEnrEffHi[ds]

        finalNatEff += totNatEff[ds]
        finalNatEffHi += totNatEffHi[ds]
        finalNatEffLo += totNatEffLo[ds]

        finalEnrExp += enrExp[ds]/365.25
        finalNatExp += natExp[ds]/365.25

    enr1 = finalEnrEff[np.where(xEff > 1.)][0]/365.25
    enr1p5 = finalEnrEff[np.where(xEff > 1.5)][0]/365.25

    plt.plot(xEff, finalEnrEff/365.25, '-r', lw=2, ls='steps', label="Enriched: %.2f kg-y" % finalEnrExp)
    plt.plot(xEff, finalEnrEffLo/365.25, '-r', lw=1, ls='steps')
    plt.plot(xEff, finalEnrEffHi/365.25, '-r', lw=1, ls='steps')

    plt.plot(xEff, finalNatEff/365.25, '-b', lw=2, ls='steps', label="Natural: %.2f kg-y" % finalNatExp)
    plt.plot(xEff, finalNatEffLo/365.25, '-b', lw=1, ls='steps')
    plt.plot(xEff, finalNatEffHi/365.25, '-b', lw=1, ls='steps')

    plt.axvline(1.0, c='g', alpha=0.5, label="1.0 keV enr: %.2f kg-y" % enr1 )
    plt.axvline(1.5, c='m', alpha=0.5, label="1.5 keV enr: %.2f kg-y" % enr1p5)

    plt.xlabel("Energy (keV)", ha='right', x=1)
    plt.ylabel("Exposure (kg-y)", ha='right', y=1)
    plt.legend(loc=1, bbox_to_anchor=(0., 0.5, 1, 0.2))
    plt.xlim(0,30)
    plt.tight_layout()
    # plt.show()
    plt.savefig('./plots/lat-eff%d-finalEff.pdf'% pctTot)

    # for ds in dsList:
    #     plt.cla()
    #     plt.plot(xEff, enrEff, '-r', label="Enriched, %.3f kg-y" % enrExp[ds])
    #     plt.plot(xEff, natEff, '-b', label="Natural,  %.3f kg-y" % natExp[ds])
    #     plt.axvline(1.0, c='g', alpha=0.5, label="1.0 keV")
    #     plt.axvline(1.5, c='m', alpha=0.5, label="1.5 keV")
    #
    #     plt.xlabel("Energy (keV)", ha='right', x=1.)
    #     plt.ylabel("Exposure (kg-d)", ha='right', y=1.)
    #     # plt.ylabel("Acceptance", ha='right', y=1.)
    #     plt.legend()
    #     plt.tight_layout()
    #     plt.show()
    #     plt.savefig("./plots/lat-expo-eff-%s-ds%s.pdf" % (mode, ds))

    # Finally, save output.
    if debugMode:
        print("I'm not saving output, I'm in debug mode, dsList:",dsList)
    if not debugMode:
        print("Saving output ...")
        np.savez(npzOut, xEff, totEnrEff, totNatEff, enrExp, natExp, finalEnrEff, finalNatEff, finalEnrExp, finalNatExp, totEnrEffLo, totEnrEffHi, totNatEffLo, totNatEffHi, detEff, detExpo)
        df = pd.DataFrame.from_dict(pndDict)
        try: os.remove(pndOut)
        except OSError: pass
        df.to_hdf(pndOut, key='runTimeInfo')


def getEfficiencyROOT():
    """
    Uses the trigger efficiency only file to build total efficiencies
    and uncertainties. Total efficiencies and uncertainties are calculated
    by a toy MC in lat-eff.py
    Note: because of the amount of data the Toy MC takes up, the binning is in
    0.1 keV bins rather than 0.01 keV bins, must rebin the efficiencies manually
    """
    from ROOT import TFile, TH1D
    import matplotlib.pyplot as plt
    # load trigger efficiency npz file and convert to TH1D for RooFit.
    f = np.load(dsi.latSWDir+'/data/lat-expo-efficiency-trig-e95.npz')

    xEff = f['arr_0']
    eMin, eMax, nBins = 0, 200, len(xEff)
    totEnrEff = f['arr_1'].tolist()
    totNatEff = f['arr_2'].tolist()
    enrExp = f['arr_3'].tolist()
    natExp = f['arr_4'].tolist()
    detEff = f['arr_13'].item()
    detExpo = f['arr_14'].item()
    # Have to manually set this because depending on the version of python,
    # sometimes the dictionary keys don't automatically get read
    dsList = [0, 1, 2, 3, 4, '5A', '5B', '5C', 6]

    # Grab total fitSlo efficiency from DB
    calDB = db.TinyDB('%s/calDB-v2.json' % (dsi.latSWDir))
    pars = db.Query()
    dbKey = "fitSlo_Combined_m2s238_eff95"
    fsN = dsi.getDBRecord(dbKey, False, calDB, pars)
    enrpars = fsN[0]
    natpars = fsN[1]
    fEffEnr = wl.weibull(xEff, *enrpars[:4])
    fEffNat = wl.weibull(xEff, *natpars[:4])

    # Correct efficiency by riseNoise efficiency
    riseEff = 0.995 # riseNoise is defined to be 99.5% efficient, no energy dependence
    fEffEnr = np.multiply(fEffEnr, riseEff)
    fEffNat = np.multiply(fEffNat, riseEff)

    totalEnrEff, totalNatEff = {}, {}
    totalEnrExp, totalNatExp = {}, {}
    for ds in dsList:
        dsNum = int(ds[0]) if isinstance(ds,str) else ds
        mods = [1]
        if dsNum == 4: mods = [2]
        if dsNum >= 5: mods = [1,2]
        totalEnrEff.setdefault(ds, np.zeros(len(xEff)))
        totalNatEff.setdefault(ds, np.zeros(len(xEff)))
        totalEnrExp.setdefault(ds, 0.)
        totalNatExp.setdefault(ds, 0.)

        for mod in mods:
            expTot = 0
            chList = det.getGoodChanList(dsNum, mod)

            for ch in chList:
                cpd = det.getChanCPD(dsNum,ch)
                intcpd = int(cpd)
                if intcpd == 253:
                    continue
                if det.allDetIDs[cpd] > 100000:
                    totalEnrEff[ds] += detEff[ds][cpd]
                    totalEnrExp[ds] += detExpo[ds][cpd]
                    # print(ds, cpd, detExpo[ds][cpd], totalEnrExp[ds])
                else:
                    totalNatEff[ds] += detEff[ds][cpd]
                    totalNatExp[ds] += detExpo[ds][cpd]
                    # print(ds, cpd, detExpo[ds][cpd], totalNatExp[ds])

                expTot += detExpo[ds][cpd]
                # print(ds, cpd, detExpo[ds][cpd])
            # print(expTot)

    # For debugging purposes
    # print(totalEnrEff)
    # print(totalNatEff)
    # print(totalEnrExp)
    # print(totalNatExp)
    # print(len(totalEnrEff[1][::10]))

    # Load Toy MC data now
    dfEnr = pd.concat(
            [pd.read_hdf(os.environ['LATDIR']+'/data/ToyMCEff_{}.h5'.format(i), 'Enr')
            for i in range(1,11)])
    dfNat = pd.concat(
            [pd.read_hdf(os.environ['LATDIR']+'/data/ToyMCEff_{}.h5'.format(i), 'Nat')
            for i in range(1,11)])

    zVal = 1.645 # this is for 90% CI for Gaussian (we checked this is valid)
    toyEnrEff = dfEnr.mean(axis=0).values
    toyNatEff = dfNat.mean(axis=0).values
    toyEnrStd = dfEnr.std(axis=0).values
    toyNatStd = dfNat.std(axis=0).values

    eMin, eMax = 0, 200
    xVals = np.arange(eMin, eMax, 0.1)
    nBins = len(xVals)

    # Correct efficiency by exposures and calculate uncertainties for each dataset
    hEnr, hNat = {}, {}
    hEnrNorm, hNatNorm = {}, {}
    hEnrHi, hNatHi = {}, {}
    hEnrHiNorm, hNatHiNorm = {}, {}
    hEnrLo, hNatLo = {}, {}
    hEnrLoNorm, hNatLoNorm = {}, {}
    hEnrHi90, hNatHi90 = {}, {}
    hEnrHi90Norm, hNatHi90Norm = {}, {}
    hEnrLo90, hNatLo90 = {}, {}
    hEnrLo90Norm, hNatLo90Norm = {}, {}

    fFile = TFile(dsi.latSWDir+'/data/lat-expo-efficiency_Combined.root', 'RECREATE')
    fFile.cd()
    for ds in dsList:
        fEffEnrTemp = np.multiply(fEffEnr[::10], totalEnrEff[ds][::10])
        fEffNatTemp = np.multiply(fEffNat[::10], totalNatEff[ds][::10])

        hEnr[ds] = TH1D('hDS{}_Enr'.format(ds), 'DS{} Enr Efficiency'.format(ds), nBins, eMin, eMax)
        hNat[ds] = TH1D('hDS{}_Nat'.format(ds), 'DS{} Nat Efficiency'.format(ds), nBins, eMin, eMax)
        hEnrHi[ds] = TH1D('hDS{}_Enr_Hi'.format(ds), 'DS{} Enr Efficiency (+68% CI)'.format(ds), nBins, eMin, eMax)
        hNatHi[ds] = TH1D('hDS{}_Nat_Hi'.format(ds), 'DS{} Nat Efficiency (+68% CI)'.format(ds), nBins, eMin, eMax)
        hEnrLo[ds] = TH1D('hDS{}_Enr_Lo'.format(ds), 'DS{} Enr Efficiency (-68% CI)'.format(ds), nBins, eMin, eMax)
        hNatLo[ds] = TH1D('hDS{}_Nat_Lo'.format(ds), 'DS{} Nat Efficiency (-68% CI)'.format(ds), nBins, eMin, eMax)
        hEnrHi90[ds] = TH1D('hDS{}_Enr_Hi90'.format(ds), 'DS{} Enr Efficiency (+90% CI)'.format(ds), nBins, eMin, eMax)
        hNatHi90[ds] = TH1D('hDS{}_Nat_Hi90'.format(ds), 'DS{} Nat Efficiency (+90% CI)'.format(ds), nBins, eMin, eMax)
        hEnrLo90[ds] = TH1D('hDS{}_Enr_Lo90'.format(ds), 'DS{} Enr Efficiency (-90% CI)'.format(ds), nBins, eMin, eMax)
        hNatLo90[ds] = TH1D('hDS{}_Nat_Lo90'.format(ds), 'DS{} Nat Efficiency (-90% CI)'.format(ds), nBins, eMin, eMax)

        hEnrNorm[ds] = TH1D('hDS{}_Norm_Enr'.format(ds), 'DS{} Enr Norm Efficiency'.format(ds),  nBins, eMin, eMax)
        hNatNorm[ds] = TH1D('hDS{}_Norm_Nat'.format(ds), 'DS{} Nat Norm Efficiency'.format(ds),  nBins, eMin, eMax)
        hEnrHiNorm[ds] = TH1D('hDS{}_Norm_Enr_Hi'.format(ds), 'DS{} Enr Norm Efficiency (+68% CI)'.format(ds), nBins, eMin, eMax)
        hNatHiNorm[ds] = TH1D('hDS{}_Norm_Nat_Hi'.format(ds), 'DS{} Nat Norm Efficiency (+68% CI)'.format(ds), nBins, eMin, eMax)
        hEnrLoNorm[ds] = TH1D('hDS{}_Norm_Enr_Lo'.format(ds), 'DS{} Enr Norm Efficiency (-68% CI)'.format(ds), nBins, eMin, eMax)
        hNatLoNorm[ds] = TH1D('hDS{}_Norm_Nat_Lo'.format(ds), 'DS{} Nat Norm Efficiency (-68% CI)'.format(ds), nBins, eMin, eMax)
        hEnrHi90Norm[ds] = TH1D('hDS{}_Norm_Enr_Hi90'.format(ds), 'DS{} Enr Norm Efficiency (+90% CI)'.format(ds), nBins, eMin, eMax)
        hNatHi90Norm[ds] = TH1D('hDS{}_Norm_Nat_Hi90'.format(ds), 'DS{} Nat Norm Efficiency (+90% CI)'.format(ds), nBins, eMin, eMax)
        hEnrLo90Norm[ds] = TH1D('hDS{}_Norm_Enr_Lo90'.format(ds), 'DS{} Enr Norm Efficiency (-90% CI)'.format(ds), nBins, eMin, eMax)
        hNatLo90Norm[ds] = TH1D('hDS{}_Norm_Nat_Lo90'.format(ds), 'DS{} Nat Norm Efficiency (-90% CI)'.format(ds), nBins, eMin, eMax)

        fEffEnrHi = fEffEnrTemp + toyEnrStd*totalEnrExp[ds]
        fEffEnrLo = fEffEnrTemp - toyEnrStd*totalEnrExp[ds]
        fEffNatHi = fEffNatTemp + toyNatStd*totalNatExp[ds]
        fEffNatLo = fEffNatTemp - toyNatStd*totalNatExp[ds]

        fEffEnrHi90 = fEffEnrTemp + zVal*toyEnrStd*totalEnrExp[ds]
        fEffEnrLo90 = fEffEnrTemp - zVal*toyEnrStd*totalEnrExp[ds]
        fEffNatHi90 = fEffNatTemp + zVal*toyNatStd*totalNatExp[ds]
        fEffNatLo90 = fEffNatTemp - zVal*toyNatStd*totalNatExp[ds]

        for idx in range(len(xVals)):
            # Unnormalized (should sum up to the exposure * efficiency)
            hEnr[ds].SetBinContent(idx+1, fEffEnrTemp[idx])
            hNat[ds].SetBinContent(idx+1, fEffNatTemp[idx])

            hEnrHi[ds].SetBinContent(idx+1, fEffEnrHi[idx])
            hEnrLo[ds].SetBinContent(idx+1, fEffEnrLo[idx])
            hNatHi[ds].SetBinContent(idx+1, fEffNatHi[idx])
            hNatLo[ds].SetBinContent(idx+1, fEffNatLo[idx])

            hEnrHi90[ds].SetBinContent(idx+1, fEffEnrHi90[idx])
            hEnrLo90[ds].SetBinContent(idx+1, fEffEnrLo90[idx])
            hNatHi90[ds].SetBinContent(idx+1, fEffNatHi90[idx])
            hNatLo90[ds].SetBinContent(idx+1, fEffNatLo90[idx])

            # Normalized (efficiency only)
            hEnrNorm[ds].SetBinContent(idx+1, fEffEnrTemp[idx]/totalEnrExp[ds])
            hNatNorm[ds].SetBinContent(idx+1, fEffNatTemp[idx]/totalNatExp[ds])

            hEnrHiNorm[ds].SetBinContent(idx+1, fEffEnrHi[idx]/totalEnrExp[ds])
            hEnrLoNorm[ds].SetBinContent(idx+1, fEffEnrLo[idx]/totalEnrExp[ds])
            hNatHiNorm[ds].SetBinContent(idx+1, fEffNatHi[idx]/totalNatExp[ds])
            hNatLoNorm[ds].SetBinContent(idx+1, fEffNatLo[idx]/totalNatExp[ds])

            hEnrHi90Norm[ds].SetBinContent(idx+1, fEffEnrHi90[idx]/totalEnrExp[ds])
            hEnrLo90Norm[ds].SetBinContent(idx+1, fEffEnrLo90[idx]/totalEnrExp[ds])
            hNatHi90Norm[ds].SetBinContent(idx+1, fEffNatHi90[idx]/totalNatExp[ds])
            hNatLo90Norm[ds].SetBinContent(idx+1, fEffNatLo90[idx]/totalNatExp[ds])

        hEnr[ds].GetXaxis().SetTitle("Energy (keV)")
        hNat[ds].GetXaxis().SetTitle("Energy (keV)")
        hEnrHi[ds].GetXaxis().SetTitle("Energy (keV)")
        hNatHi[ds].GetXaxis().SetTitle("Energy (keV)")
        hEnrLo[ds].GetXaxis().SetTitle("Energy (keV)")
        hNatLo[ds].GetXaxis().SetTitle("Energy (keV)")
        hEnrHi90[ds].GetXaxis().SetTitle("Energy (keV)")
        hNatHi90[ds].GetXaxis().SetTitle("Energy (keV)")
        hEnrLo90[ds].GetXaxis().SetTitle("Energy (keV)")
        hNatLo90[ds].GetXaxis().SetTitle("Energy (keV)")

        hEnrNorm[ds].GetXaxis().SetTitle("Energy (keV)")
        hNatNorm[ds].GetXaxis().SetTitle("Energy (keV)")
        hEnrHiNorm[ds].GetXaxis().SetTitle("Energy (keV)")
        hNatHiNorm[ds].GetXaxis().SetTitle("Energy (keV)")
        hEnrLoNorm[ds].GetXaxis().SetTitle("Energy (keV)")
        hNatLoNorm[ds].GetXaxis().SetTitle("Energy (keV)")
        hEnrHi90Norm[ds].GetXaxis().SetTitle("Energy (keV)")
        hNatHi90Norm[ds].GetXaxis().SetTitle("Energy (keV)")
        hEnrLo90Norm[ds].GetXaxis().SetTitle("Energy (keV)")
        hNatLo90Norm[ds].GetXaxis().SetTitle("Energy (keV)")

        hEnr[ds].GetYaxis().SetTitle("Efficiency (kg-d)")
        hNat[ds].GetYaxis().SetTitle("Efficiency (kg-d)")
        hEnrHi[ds].GetYaxis().SetTitle("Efficiency (kg-d)")
        hNatHi[ds].GetYaxis().SetTitle("Efficiency (kg-d)")
        hEnrLo[ds].GetYaxis().SetTitle("Efficiency (kg-d)")
        hNatLo[ds].GetYaxis().SetTitle("Efficiency (kg-d)")
        hEnrHi90[ds].GetYaxis().SetTitle("Efficiency (kg-d)")
        hNatHi90[ds].GetYaxis().SetTitle("Efficiency (kg-d)")
        hEnrLo90[ds].GetYaxis().SetTitle("Efficiency (kg-d)")
        hNatLo90[ds].GetYaxis().SetTitle("Efficiency (kg-d)")

        hEnrNorm[ds].GetYaxis().SetTitle("Efficiency")
        hNatNorm[ds].GetYaxis().SetTitle("Efficiency")
        hEnrHiNorm[ds].GetYaxis().SetTitle("Efficiency")
        hNatHiNorm[ds].GetYaxis().SetTitle("Efficiency")
        hEnrLoNorm[ds].GetYaxis().SetTitle("Efficiency")
        hNatLoNorm[ds].GetYaxis().SetTitle("Efficiency")
        hEnrHi90Norm[ds].GetYaxis().SetTitle("Efficiency")
        hNatHi90Norm[ds].GetYaxis().SetTitle("Efficiency")
        hEnrLo90Norm[ds].GetYaxis().SetTitle("Efficiency")
        hNatLo90Norm[ds].GetYaxis().SetTitle("Efficiency")

        hEnr[ds].Write()
        hNat[ds].Write()
        hEnrHi[ds].Write()
        hNatHi[ds].Write()
        hEnrLo[ds].Write()
        hNatLo[ds].Write()
        hEnrHi90[ds].Write()
        hNatHi90[ds].Write()
        hEnrLo90[ds].Write()
        hNatLo90[ds].Write()

        hEnrNorm[ds].Write()
        hNatNorm[ds].Write()
        hEnrHiNorm[ds].Write()
        hNatHiNorm[ds].Write()
        hEnrLoNorm[ds].Write()
        hNatLoNorm[ds].Write()
        hEnrHi90Norm[ds].Write()
        hNatHi90Norm[ds].Write()
        hEnrLo90Norm[ds].Write()
        hNatLo90Norm[ds].Write()

    fFile.Close()


if __name__=="__main__":
    main(sys.argv[1:])
