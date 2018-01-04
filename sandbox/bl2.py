#!/usr/bin/env python
import sys, imp, time, json
gStart = time.time()
sys.argv.append("-b")
ds = imp.load_source('DataSetInfo','../DataSetInfo.py')
wl = imp.load_source('waveLibs','../waveLibs.py')
from ROOT import gDirectory, TFile, TTree, TChain, MGTWaveform, MJTMSWaveform, GATDataSet
import numpy as np
import matplotlib.pyplot as plt
print "import time:",time.time()-gStart

def main():

    # testParsers()
    bkgBaselines()


def bkgBaselines():

    # generate BL files (run once)
    # highECut = "channel%2==0 && trapENFCal > 100"
    # for dsNum in [2]:
    #     goodList = ds.GetGoodChanListNew(dsNum)
    #     parseLAT2(goodList, highECut, dsNum, None, True)
    # return

    dsNum = 2
    goodList = ds.GetGoodChanListNew(dsNum)
    with open("../data/parseLAT2_ds%d.json" % dsNum, 'r') as fp: baseDict = json.load(fp)
    baseDict = {ch : baseDict[u'%d'%ch] for ch in goodList}

    # clint, you were testing the new baseDict format with DS2.
    # 1. you want to remake the pulser plots using e.g. global time
    # 2. you want to update the other to parsers to use the new format

    timeList = [baseDict[608][idx][1] for idx in range(len(baseDict[608]))]

    return

    ch = 640
    fig = plt.figure(figsize=(9,5),facecolor='w')
    p1 = plt.subplot(111)

    cmap = plt.cm.get_cmap('hsv',len(baseDict.keys())+1)

    # all baselines
    ctr = 0
    # for ch in sorted(baseDict.keys()):
    for ch in [594,598,600]:
        p1.plot(baseDict[ch],"o",markersize=1,c=cmap(ctr),label="ch%d"%ch)
        ctr += 1
    p1.set_title("DS-%d Baselines" % dsNum)
    p1.set_xlabel("Entry")
    p1.set_ylabel("ADC value")
    p1.legend(loc='best')
    plt.show(block=False)
    print "tryna save"
    plt.savefig("../plots/ds0_bljumps.pdf")

    # histograms
    # ctr = 0
    # # for ch in sorted(baseDict.keys()):
    # for ch in [594,598,600]:
    #     # p1.cla()
    #     p1.set_title("DS-%d Baselines" % dsNum)
    #     p1.set_xlabel("ADC")
    #     p1.set_ylabel("counts")
    #     adcLo, adcHi, binsPerADC = min(baseDict[ch]), max(baseDict[ch]), 0.1
    #     nBins = int((adcHi - adcLo)/binsPerADC)
    #     p1.hist(baseDict[ch], nBins, range=(40,200), facecolor=cmap(ctr), histtype="step", label="ch%d"%ch)
    #     ctr += 1
    #     p1.legend(loc='best')
    #     plt.show(block=False)
    #     plt.savefig("../plots/bl_ch%d.pdf" % ch)

    # plt.show()


def testParsers():

    dsNum = 0
    goodList = ds.GetGoodChanListNew(dsNum)

    pulserCut = "EventDC1Bits && channel%2==0 && trapENFCal > 50" # SELECTS pulsers
    highECut = "!EventDC1Bits && channel%2==0 && trapENFCal > 50"
    testCut = "channel%2==0 && trapENFCal > 50"

    # functions to walk ROOT files and fill baseDict
    # baseDict = parseGAT(goodList, testCut, dsNum, 1, True)
    # baseDict1 = parseLAT1(goodList, testCut, dsNum, 1, True) # scans wfs directly
    # baseDict2 = parseLAT2(goodList, testCut, dsNum, 1, True) # uses 'fitBL' parameter (fastest)

    # for ch in goodList:
    #     print "ch",ch
    #     print "  avg:",baseDict1[ch]
    #     print "  fit:",baseDict2[ch]

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
    print timeList
    print time.time()-start,"seconds elapsed."

    # 2 - np array, multi-axis slice
    # start = time.time()
    # B = np.array(baseDict[608])
    # print B[:,0].tolist()
    # print time.time()-start,"seconds elapsed."


def parseGAT(goodList, theCut, dsNum, bkgIdx=None, saveMe=False):
    """ Accesses MGTWaveform objects directly and returns some arbitrary object.
    NOTE: Have to access individual TFile's because of this old stupid ROOT/MJSW bug:
    When the GetEntry command changes from file1->file2 in the loop:
        AttributeError: 'TChain' object has no attribute 'MGTWaveforms' (LAT data)
        AttributeError: 'TChain' object has no attribute 'event'        (GAT data)
    """
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
                baseDict[chan].append(baseline)

        print "Run %d, %d entries, %d passing cuts. %.2f sec." % (run, gatTree.GetEntries(), eList.GetN(), time.time()-start)

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
    # this is what we return
    baseDict = {ch:[] for ch in goodList}

    p1Start = time.time()

    # get a sequential list of file names
    latDir, latList = wl.getLATList(dsNum, bkgIdx)

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
                baseDict[chan].append(baseline)

        print "%s, %d entries, %d passing cuts. %.2f sec." % (fName, latTree.GetEntries(), eList.GetN(), time.time()-start)

        latFile.Close()

    # optionally save the output
    if saveMe:
        with open("../data/parseLAT1_ds%d.json" % dsNum, 'w') as fOut:
            json.dump(baseDict, fOut)

    print "Total Elapsed:",time.time() - p1Start
    return baseDict


def parseLAT2(goodList, theCut, dsNum, bkgIdx=None, saveMe=False):
    """ Uses "fitBL" variable in LAT files to do a draw command and return an arbitrary object.
        This is able to use a TChain since it doesn't access any custom classes.
        Also it's faster.  But the TCut is global so it's maybe a little less flexible.
    """
    # this is what we return
    baseDict = {ch:[] for ch in goodList}
    start = time.time()

    # get a sequential list of file names
    latDir, latList = wl.getLATList(dsNum, bkgIdx)

    # load the chain
    latChain = TChain("skimTree")
    [latChain.Add(latDir + "/" + fName) for fName in latList]
    print "Loaded DS%d, %d entries.  Drawing ..." % (dsNum, latChain.GetEntries())

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
        baseDict[ch].append((bl,gt,run))

    print "DS-%d, bkgIdx" % dsNum, bkgIdx,"%d entries passing cuts. %.2f sec." % (nPass, time.time()-start)

    # optionally save the output
    if saveMe:
        with open("../data/parseLAT2_ds%d.json" % dsNum, 'w') as fOut:
            json.dump(baseDict, fOut)

    return baseDict


if __name__=="__main__":
    main()