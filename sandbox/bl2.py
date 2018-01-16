#!/usr/bin/env python3

import sys, os, imp, time, json
gStart = time.time()
sys.argv.append("-b")
ds = imp.load_source('DataSetInfo','../DataSetInfo.py')
wl = imp.load_source('waveLibs','../waveLibs.py')
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
from scipy.optimize import curve_fit
home = os.path.expanduser('~')
# print("Imports: %.2f" % (time.time()-gStart))

def main():

    # testParsers()
    # genBLFiles()
    # bkgBaselines()
    # baselinesVsTime()
    calBLFit()


def calBLFit():
    """
    key:  bl_ds[i]_idx[j]_mod[k]
    vals: {[chan]:[blLo, blHi, fitMu, fitSig, N, beta, m, fitChi2, cts]}
    """
    intMode = False
    updateDB = False

    # iterate over datasets
    # for dsNum, modNum in [(0,1),(1,1),(2,1),(3,1),(4,2),(5,1),(5,2)]:
    for dsNum, modNum in [(5,2)]:

        print("Scanning DS-%d-M%d ..." % (dsNum, modNum))

        goodList = ds.GetGoodChanListNew(dsNum)
        with open("%s/project/baselines/parseBL_ds%d_bkg.json" % (home, dsNum), 'r') as fp:
            baseDict = json.load(fp)
        baseDict = {ch : baseDict[u'%d'%ch] for ch in goodList}

        cInfo = ds.CalInfo()
        nCal = ds.getNCalIdxs(dsNum,module=modNum)
        calDicts = {}
        for calIdx in range(nCal):
            with open("%s/project/baselines/parseBL_ds%d_cal%d_m%d.json" % (home,dsNum,calIdx,modNum), 'r') as fp:
                calDict = json.load(fp)
            calDicts[calIdx] = {ch: calDict[u'%d'%ch] for ch in goodList}

        fig = plt.figure(figsize=(9,5), facecolor='w')
        p1 = plt.subplot(111)

        # loop over calIdx's
        for calIdx in range(nCal):
        # for calIdx in [1]:
            calDict = calDicts[calIdx]
            print("DS-%d, calIdx %d" % (dsNum, calIdx))

            # this is what we fill the DB with
            baseDBDict = {ch:[] for ch in goodList}

            # loop over channels (interactive to display plots)
            iCh = -1
            while True:
                iCh += 1
                if intMode==True and iCh != 0:
                    value = input()
                    if value=='q': exit(1)    # quit
                    if value=='p': iCh -= 2   # go to previous
                    if (value.isdigit()):
                        iCh = int(value)      # go to channel number
                if iCh >= len(goodList): break
                ch = goodList[iCh]

                # get calibration baseline data for this channel
                calBL = [calDict[ch][idx][0] for idx in range(len(calDict[ch]))]
                calRuns = [calDict[ch][idx][2] for idx in range(len(calDict[ch]))]
                if len(calBL)==0:
                    print("Ch %d ain't got no entries" % ch)
                    continue

                # get background baseline data matching the cal run, trim off outliers (overflow hits)
                calKey = "ds%d_m%d" % (dsNum, modNum)
                calIdx = cInfo.GetCalIdx(calKey, calRuns[0])
                runLo, runHi = cInfo.GetCalRunCoverage(calKey, calIdx)
                rList = [baseDict[ch][idx][2] for idx in range(len(baseDict[ch])) if runLo <= baseDict[ch][idx][2] <= runHi]
                bList = [baseDict[ch][idx][0] for idx in range(len(baseDict[ch])) if runLo <= baseDict[ch][idx][2] <= runHi]
                rList = [rList[idx] for idx in range(len(rList)) if abs(bList[idx]) < 8000]
                bList = [bList[idx] for idx in range(len(bList)) if abs(bList[idx]) < 8000]

                # find mode
                calInt = [int(bl) for bl in calBL]
                xMode = float(stats.mode(calInt)[0])
                if not -8000 < int(xMode) < 8000:
                    print("Ch %d, mode is bad: %d" % (ch, xMode))
                    continue

                # narrow window and estimate spread
                xLo, xHi, bpa = xMode - 25, xMode + 5, 0.2
                nBins = int((xHi-xLo)/bpa)
                h1, x1 = np.histogram(np.asarray(calBL), bins=nBins, range=(xLo, xHi))
                x1 = x1[:-1]
                hInt = [h1[0]]
                for i in range(1, len(h1)):
                    hInt.append( hInt[i-1] + h1[i] )
                hInt = np.asarray(hInt)
                idxLo = np.where(hInt >= 0.4 * hInt[-1])
                idxHi = np.where(hInt >= 0.95 * hInt[-1])
                adcLo = x1[idxLo][0]
                adcHi = x1[idxHi][0]
                sig = adcHi - adcLo

                # run fit and compute reduced chi square
                pGuess = [50., 1., 1., xMode, sig] # N,beta,m,mu,sig
                try:
                    pFit,_ = curve_fit(crystalball, x1, h1, p0=pGuess)
                except RuntimeError:
                    print("Fit failed!  Reverting to guess parameters.")
                    pFit = pGuess
                fitVals = crystalball(x1,*pFit)
                fitChi2 = np.sum(np.square(fitVals - h1)) / (len(h1) - len(pFit))
                fitMu, fitSig, N, beta, m = pFit[3], abs(pFit[4]), pFit[0], pFit[1], pFit[2]
                cts = np.sum(h1)

                # get ADC ranges
                blLo, blHi = fitMu - 5, fitMu + 20
                if len(bList) > 0:
                    if min(bList) < blLo:
                        print("Ch %d (lo %.0f) has a low-BL event: %.0f. Adjusting ..." % (ch, blLo, min(bList)))
                        blLo = min(bList)
                    if max(bList) > blHi:
                        print("Ch %d (hi %.0f) has a high-BL event: %.0f.  Adjusting ..." % (ch, blHi, max(bList)))
                        blHi = max(bList)

                # create DB entry for this channel (chop off unnecessary digits)
                results = []
                for val in [blLo, blHi, fitMu, fitSig, N, beta, m, fitChi2, cts] :
                    if isinstance(val, float):
                        results.append(float("%.2f" % val))
                    else:
                        results.append(val)
                baseDBDict[ch] = results

                print("DS-%d-M%d calIdx %d  ch %d  cts %-5d  mu %-7.2f  sig %-7.2f  chi2 %-7.2f  LO %-7.0f  HI %.0f" % (dsNum,modNum,calIdx,ch,cts,fitMu,fitSig,fitChi2,blLo,blHi))

                if not intMode: continue

                # -- make a figure --
                p1.cla()
                p1.plot(x1,h1,ls="steps",color='blue',label='data')

                # p1.axvline(xMode,color='red')
                # p1.axvline(adcLo,color='green')
                # p1.axvline(adcHi,color='green')
                p1.axvline(blLo,color='black')
                p1.axvline(blHi,color='black')
                for idx,bkgVal in enumerate(bList):
                    lbl = "bkgIdx vals" if idx == 0 else None
                    p1.axvline(bkgVal,color='cyan',alpha=0.5, label=lbl)

                guessVals = crystalball(x1,*pGuess)
                p1.plot(x1,guessVals,color='orange',label='guess',alpha=0.7)

                from scipy.interpolate import spline
                xnew = np.linspace(x1.min(),x1.max(),len(x1)*10)
                fitSmooth = spline(x1,fitVals,xnew)
                p1.plot(xnew,fitSmooth,color='red',label=r'result, $\chi^2 = %.2f$' % fitChi2)

                p1.set_title("DS-%d, calIdx %d, ch %d" % (dsNum, calIdx, ch))
                p1.set_xlabel("ADC")
                p1.set_ylabel("Counts")
                p1.legend(loc=2)

                plt.tight_layout()
                plt.pause(0.0001)

            # Update the DB for this calIdx
            if updateDB:
                dbKey = "bl_ds%d_idx%d_mod%d" % (dsNum, calIdx, modNum)
                ds.setDBRecord({"key":dbKey, "vals":baseDBDict}, False, "../calDB.json")


def crystalball(x,N,beta,m,mu,sig):
    """ Scipy's crystalball function (doesn't exist in python2.)
    https://github.com/scipy/scipy/blob/59cabc8/scipy/stats/_continuous_distns.py#L5798
    """
    from scipy._lib._util import _lazywhere
    x = (x - mu)/sig
    rhs = lambda x, beta, m: np.exp(-x**2 / 2)
    lhs = lambda x, beta, m: (m/beta)**m * np.exp(-beta**2 / 2.0) * (m/beta - beta - x)**(-m)
    return N * _lazywhere(np.atleast_1d(x > -beta), (x, beta, m), f=rhs, f2=lhs)


def baselinesVsTime():

    dsNum, modNum = 1, 1
    goodList = ds.GetGoodChanListNew(dsNum)
    with open("%s/project/baselines/parseBL_ds%d_bkg.json" % (home, dsNum), 'r') as fp:
        baseDict = json.load(fp)
    baseDict = {ch : baseDict[u'%d'%ch] for ch in goodList}

    nCal = ds.getNCalIdxs(dsNum,module=modNum)
    calDicts = {}
    for calIdx in range(nCal):
        with open("%s/project/baselines/parseBL_ds%d_cal%d_m%d.json" % (home,dsNum,calIdx,modNum), 'r') as fp:
            calDict = json.load(fp)
        calDicts[calIdx] = {ch: calDict[u'%d'%ch] for ch in goodList}

    ch = 592

    fig = plt.figure(figsize=(9,5), facecolor='w')
    p1 = plt.subplot(111)
    # plt.show(block=False)

    # loop over channels (interactive to display plots)
    iCh = -1
    while True:
        iCh += 1
        if iCh != 0:
            value = input()
            if value=='q': break
        if iCh > len(goodList): break
        ch = goodList[iCh]

        print("DS%d, ch %d" % (dsNum, ch))

        p1.cla()

        bList = [baseDict[ch][idx][0] for idx in range(len(baseDict[ch]))]
        rList = [baseDict[ch][idx][2] for idx in range(len(baseDict[ch]))]

        bMean = np.mean(bList)
        p1.set_ylim([bMean - 300, bMean + 300])

        p1.plot(rList,bList,"o",markersize=1,color='blue',label="bkg")

        for calIdx in range(nCal):
            calDict = calDicts[calIdx]
            baseListC = [calDict[ch][idx][0] for idx in range(len(calDict[ch])) if calDict[ch][idx][0] > -8000]
            runsListC = [calDict[ch][idx][2] for idx in range(len(calDict[ch])) if calDict[ch][idx][0] > -8000]

            cLabel = "Cal" if calIdx==0 else None
            p1.plot(runsListC, baseListC, "o", markersize=1, color='red', label=cLabel)


        p1.set_title("DS%d, Ch %d" % (dsNum, ch))
        p1.set_xlabel("Run")
        p1.set_ylabel("ADC value")
        p1.legend(loc='best')

        plt.tight_layout()
        plt.pause(0.00001)


def bkgBaselines():

    dsNum = 2
    goodList = ds.GetGoodChanListNew(dsNum)
    with open("../data/parseLAT2_ds%d_bkg.json" % dsNum, 'r') as fp: baseDict = json.load(fp)
    baseDict = {ch : baseDict[u'%d'%ch] for ch in goodList}

    fig = plt.figure(figsize=(9,5),facecolor='w')
    p1 = plt.subplot(111)
    cmap = plt.cm.get_cmap('hsv',len(baseDict.keys())+1)

    # 1 - baselines vs. time
    # ctr = 0
    # for ch in sorted(baseDict.keys()):
    # # for ch in [594,598,600]:
    #     baseList = [baseDict[ch][idx][0] for idx in range(len(baseDict[ch]))]
    #     timeList = [baseDict[ch][idx][1] for idx in range(len(baseDict[ch]))]
    #     runsList = [baseDict[ch][idx][2] for idx in range(len(baseDict[ch]))]
    #
    #     p1.plot(runsList,baseList,"o",markersize=1,c=cmap(ctr),label="ch%d"%ch)
    #     ctr += 1
    # p1.set_title("DS-%d Baselines" % dsNum)
    # p1.set_xlabel("Entry")
    # p1.set_ylabel("ADC value")
    # p1.legend(bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0.)
    # plt.tight_layout()
    # plt.subplots_adjust(right=0.85)
    # plt.show()
    # plt.savefig("../plots/ds%d_blHist.pdf" % dsNum)

    # 2 - baseline histograms
    bpa = 1
    adcLo, adcHi = 50, 160
    for ch in baseDict.keys():
        baseList = [baseDict[ch][idx][0] for idx in range(len(baseDict[ch]))]
        if max(baseList) > adcHi:
            adcHi = max(baseList) + 10
        elif min(baseList) < adcLo:
            adcLo = min(baseList) - 10
    nBins = int((adcHi-adcLo)/bpa)
    ctr = 0
    for ch in sorted(baseDict.keys()):
    # for ch in [594,598,600]:
        # p1.cla()
        baseList = [baseDict[ch][idx][0] for idx in range(len(baseDict[ch]))]
        p1.hist(baseList, nBins, facecolor=cmap(ctr), histtype="step", label="ch%d"%ch, range=(60,160)) # in plt.plot it's ls='steps'
        ctr += 1
    p1.legend(loc='best')
    p1.set_title("DS-%d Baselines" % dsNum)
    p1.set_xlabel("ADC")
    p1.set_ylabel("counts")
    p1.legend(bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0.)
    plt.tight_layout()
    plt.subplots_adjust(right=0.85)
        # plt.savefig("../plots/bl_ds%d_ch%d.pdf" % (dsNum,ch))
    plt.show()


def genBLFiles():
    """ Run parseLAT2 for each dataset and generate a json output file in ../data .
    Time to draw goes as the num. entries in the chains (DS0 and DS5 take ~10mins.)
    """
    # 1.
    # highECut = "channel%2==0 && trapENFCal > 100"
    # for dsNum in [0,1,2,3,4,5]:
    #     goodList = ds.GetGoodChanListNew(dsNum)
    #     parseLAT2(goodList, highECut, dsNum, None, True)
    # return

    # 2.
    # for dsNum,modNum in [(0,1),(1,1),(2,1),(3,1),(4,2),(5,1),(5,2)]:
    for dsNum,modNum in [(5,2)]:
        calCut = "channel%2==0 && trapENFCal > 50 && Entry$ < 10000" # don't need that many events for each calIdx
        goodList = ds.GetGoodChanListNew(dsNum)
        # for idx in range(ds.getNCalIdxs(dsNum, modNum)):
            # parseLAT2(goodList, calCut, dsNum, idx, True, True, modNum)

        # debug - parseBL_ds5_cal6_m2.json is missing
        parseLAT2(goodList, calCut, dsNum, 6, True, True, modNum)


def testParsers():

    dsNum = 2
    goodList = ds.GetGoodChanListNew(dsNum)

    pulserCut = "EventDC1Bits && channel%2==0 && trapENFCal > 50" # SELECTS pulsers
    highECut = "!EventDC1Bits && channel%2==0 && trapENFCal > 50"
    testCut = "channel%2==0 && trapENFCal > 50 && Entry$ < 100000" # don't need all cal events

    # functions to walk ROOT files and fill baseDict
    # baseDict = parseGAT(goodList, testCut, dsNum, 1, True)
    # baseDict1 = parseLAT1(goodList, testCut, dsNum, 1, True) # scans wfs directly
    # baseDict2 = parseLAT2(goodList, testCut, dsNum, None, True, True) # uses 'fitBL' parameter (fastest)
    # return

    # for ch in goodList:
    #     print("ch",ch)
    #     print("  avg:",baseDict1[ch])
    #     print("  fit:",baseDict2[ch])

    # load baseDict from file (very quick)
    # with open("../data/parseGAT_ds%d.json" % dsNum, 'r') as fp: baseDict = json.load(fp)
    # with open("../data/parseLAT1_ds%d.json" % dsNum, 'r') as fp: baseDict = json.load(fp)
    with open("../data/parseLAT2_ds%d.json" % dsNum, 'r') as fp: baseDict = json.load(fp)

    # convert unicode keys back to ints
    baseDict = {ch : baseDict[u'%d'%ch] for ch in goodList}

    # format: baseDict[ch] = (baseline, globalTime, run)

    # extract columns of numbers (e.g. globalTime's)

    # 1 - list comprehension (faster)
    start = time.time()
    timeList = [baseDict[608][idx][1] for idx in range(len(baseDict[608]))]
    print(timeList)
    print(time.time()-start)

    # 2 - np array, multi-axis slice
    start = time.time()
    B = np.array(baseDict[608])
    print(B[:,0].tolist())
    print(time.time()-start)


def parseGAT(goodList, theCut, dsNum, bkgIdx=None, saveMe=False):
    """ Accesses MGTWaveform objects directly and returns some arbitrary object.
    NOTE: Have to access individual TFile's because of this old stupid ROOT/MJSW bug:
    When the GetEntry command changes from file1->file2 in the loop:
        AttributeError: 'TChain' object has no attribute 'MGTWaveforms' (LAT data)
        AttributeError: 'TChain' object has no attribute 'event'        (GAT data)
    """
    from ROOT import gDirectory, TFile, TTree, TChain, MGTWaveform, MJTMSWaveform, GATDataSet

    # this is what we return
    baseDict = {ch:[] for ch in goodList}

    # create run number list
    runList = []
    bkgList = ds.bkgRunsDS[dsNum][bkgIdx]
    for idx in range(0,len(bkgList),2):
        [runList.append(run) for run in range(bkgList[idx], bkgList[idx+1]+1)]

    # get paths to data
    gds = GATDataSet()
    gatPath = gds.GetPathToRun(runList[0], GATDataSet.kGatified)
    gIdx = gatPath.find("mjd_run")
    gatPath = gatPath[:gIdx]
    bltPath = gds.GetPathToRun(runList[0], GATDataSet.kBuilt)
    bIdx = bltPath.find("OR_run")
    bltPath = bltPath[:bIdx]

    # loop over files
    for run in runList:
        start = time.time()

        # load trees
        gatFile = TFile(gatPath+"mjd_run%d.root" % run)
        gatTree = gatFile.Get("mjdTree")
        bltFile = TFile(bltPath+"OR_run%d.root" % run)
        bltTree = bltFile.Get("MGTree")
        gatTree.AddFriend(bltTree)

        # create TEntryList
        gatTree.Draw(">>elist", theCut, "entrylist")
        eList = gDirectory.Get("elist")
        gatTree.SetEntryList(eList)
        if eList.GetN()==0: continue

        # loop over entries passing cuts
        for iList in range(eList.GetN()):
            entry = gatTree.GetEntryNumber(iList);
            gatTree.LoadTree(entry)
            gatTree.GetEntry(entry)
            nChans = gatTree.channel.size()
            numPass = gatTree.Draw("channel",theCut,"GOFF",1,iList)
            chans = gatTree.GetV1()
            chanList = list(set(int(chans[n]) for n in xrange(numPass)))

            # loop over hits passing cuts (channels in good list only)
            hitList = (iH for iH in range(nChans) if gatTree.channel.at(iH) in goodList and gatTree.channel.at(iH) in chanList)
            for iH in hitList:
                if dsNum==2 or dsNum==6:
                    wf_downsampled = bltTree.event.GetWaveform(iH)
                    wf_regular = bltTree.event.GetAuxWaveform(iH)
                    wf = MJTMSWaveform(wf_downsampled,wf_regular)
                else:
                    wf = bltTree.event.GetWaveform(iH)

                # now we can pretty much save or calculate anything we want
                chan = int(gatTree.channel.at(iH))
                signal = wl.processWaveform(wf)
                baseline,_ = signal.GetBaseNoise()
                baseline = float("%.2f" % baseline) # kill precision to save on memory
                globalTime = int(latTree.globalTime.fSec)
                run = int(latTree.run)
                baseDict[chan].append(baseline,globalTime,run)

        print("Run %d, %d entries, %d passing cuts. %.2f sec." % (run, gatTree.GetEntries(), eList.GetN(), time.time()-start))

        gatFile.Close()
        bltFile.Close()

    # optionally save the output
    if saveMe:
        with open("../data/parseGAT_ds%d.json" % dsNum, 'w') as fOut:
            json.dump(baseDict, fOut)

    return baseDict


def parseLAT1(goodList, theCut, dsNum, bkgIdx=None, saveMe=False):
    """ Accesses MGTWaveform objects directly and returns some arbitrary object.
    Suffers from the same ROOT/MGSW TChain bug as parseGAT (see above).
    """
    from ROOT import gDirectory, TFile, TTree, TChain, MGTWaveform, MJTMSWaveform, GATDataSet

    # this is what we return
    baseDict = {ch:[] for ch in goodList}

    p1Start = time.time()

    # get a sequential list of file names
    latDir, latList = ds.getLATList(dsNum, bkgIdx)

    # loop over each file in the list
    for fName in latList:
        start = time.time()

        # load trees
        latFile = TFile(latDir + "/" + fName)
        latTree = latFile.Get("skimTree")

        # create TEntryList
        latTree.Draw(">>elist", theCut, "entrylist")
        eList = gDirectory.Get("elist")
        latTree.SetEntryList(eList)

        # loop over entries passing cuts
        for iList in range(eList.GetN()):
            entry = latTree.GetEntryNumber(iList);
            latTree.LoadTree(entry)
            latTree.GetEntry(entry)
            nChans = latTree.channel.size()
            numPass = latTree.Draw("channel",theCut,"GOFF",1,iList)
            chans = latTree.GetV1()
            chanList = list(set(int(chans[n]) for n in xrange(numPass)))

            # loop over hits passing cuts (channels in good list only)
            hitList = (iH for iH in range(nChans) if latTree.channel.at(iH) in goodList and latTree.channel.at(iH) in chanList)
            for iH in hitList:

                # now we can pretty much save or calculate anything we want
                wf = latTree.MGTWaveforms.at(iH)
                chan = int(latTree.channel.at(iH))
                signal = wl.processWaveform(wf)
                baseline,_ = signal.GetBaseNoise()
                baseline = float("%.2f" % baseline) # kill precision to save on memory
                globalTime = int(latTree.globalTime.fSec)
                run = int(latTree.run)
                baseDict[chan].append(baseline,globalTime,run)

        print("%s, %d entries, %d passing cuts. %.2f sec." % (fName, latTree.GetEntries(), eList.GetN(), time.time()-start))

        latFile.Close()

    # optionally save the output
    if saveMe:
        with open("../data/parseLAT1_ds%d.json" % dsNum, 'w') as fOut:
            json.dump(baseDict, fOut)

    print( "Total Elapsed:",time.time() - p1Start)
    return baseDict


def parseLAT2(goodList, theCut, dsNum, bkgIdx=None, saveMe=False, cal=False, modNum=None):
    """ Uses "fitBL" variable in LAT files to do a draw command and return an arbitrary object.
        This is able to use a TChain since it doesn't access any custom classes.
        Also it's faster.  But the TCut is global so it's maybe a little less flexible.
    """
    from ROOT import gDirectory, TFile, TTree, TChain, MGTWaveform, MJTMSWaveform, GATDataSet

    # this is what we return
    baseDict = {ch:[] for ch in goodList}
    start = time.time()

    # get a sequential list of file names
    if cal:
        calIdx = bkgIdx
        latDir = home + "/project/cal-lat"
        calList = ds.getCalFiles(dsNum, calIdx, modNum, verbose=False)
        latList = [f[f.find("latSkim"):] for f in calList ]
        if bkgIdx is None:
            print("Scanning %d cal idx's..." % ds.getNCalIdxs(dsNum,module=1))
    else:
        latDir, latList = ds.getLATList(dsNum, bkgIdx)

    # load the chain
    latChain = TChain("skimTree")
    for fName in latList: latChain.Add(latDir + "/" + fName)
    print("Loaded DS%d, %d entries.  Drawing ..." % (dsNum, latChain.GetEntries()))

    # run the draw command
    nPass = latChain.Draw("fitBL:channel:run:globalTime",theCut,"goff")
    if nPass == 0: return
    arrBL = latChain.GetV1()
    arrCH = latChain.GetV2()
    arrRN = latChain.GetV3()
    arrGT = latChain.GetV4()
    arrBL = [float("%.2f" % arrBL[idx]) for idx in range(nPass)] # save precision
    arrCH = [int(arrCH[idx]) for idx in range(nPass)]
    arrRN = [int(arrRN[idx]) for idx in range(nPass)]
    arrGT = [int(arrGT[idx]) for idx in range(nPass)]

    # fill the baseDict
    for idx in range(nPass):
        ch, bl, gt, run = arrCH[idx], arrBL[idx], arrGT[idx], arrRN[idx]
        if ch not in goodList: continue
        baseDict[ch].append((bl,gt,run))

    print("DS-%d, idx" % dsNum, bkgIdx,"%d entries passing cuts. %.2f sec." % (nPass, time.time()-start))

    # optionally save the output
    if saveMe:
        if bkgIdx is None:
            fName = "%s/project/baselines/parseBL_ds%d_bkg.json" % (home, dsNum)
        elif not cal:
            fName = "%s/project/baselines/parseBL_ds%d_%d.json" % (home, dsNum, bkgIdx)
        else:
            fName = "%s/project/baselines/parseBL_ds%d_cal%d_m%d.json" % (home, dsNum, calIdx, modNum)
        with open(fName, 'w') as fOut:
            json.dump(baseDict, fOut)
        print("Saved json file:",fName)

    return baseDict


if __name__=="__main__":
    main()