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
    # for dsNum in [0,1,2,3,4,5]:
    #     goodList = ds.GetGoodChanListNew(dsNum)
    #     parseLAT2(goodList, highECut, dsNum, None, True)

    dsNum = 1
    goodList = ds.GetGoodChanListNew(dsNum)
    with open("../data/parseLAT2_ds%d.json" % dsNum, 'r') as fp: baseDict = json.load(fp)
    baseDict = {ch : baseDict[u'%d'%ch] for ch in goodList}

    # ch = 608
    # fig = plt.figure(fig)



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
    nPass = latChain.Draw("fitBL:channel",theCut,"goff")
    if nPass == 0: return
    arrBL = latChain.GetV1()
    arrCH = latChain.GetV2()
    arrBL = [float("%.2f" % arrBL[idx]) for idx in range(nPass)]
    arrCH = [float("%.2f" % arrCH[idx]) for idx in range(nPass)]

    # fill the baseDict
    for idx in range(nPass):
        bl, ch = arrBL[idx], arrCH[idx]
        baseDict[ch].append(bl)

    print "DS-%d, bkgIdx" % dsNum, bkgIdx,"%d entries passing cuts. %.2f sec." % (nPass, time.time()-start)

    # optionally save the output
    if saveMe:
        with open("../data/parseLAT2_ds%d.json" % dsNum, 'w') as fOut:
            json.dump(baseDict, fOut)

    return baseDict


if __name__=="__main__":
    main()