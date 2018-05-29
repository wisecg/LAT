#!/usr/bin/env python3
"""
===================== LAT/chan-sel.py ======================
Perform necessary run/channel selection bookkeeping for
    LAT analysis, before and after cuts are applied.
===================== C. Wiseman (USC) =====================
"""
import sys, os, time
import tinydb as db
import numpy as np

# LAT libraries
import dsi
bkg = dsi.BkgInfo()
cal = dsi.CalInfo()
det = dsi.DetInfo()
import waveLibs as wl

def main(argv):

    ds, cIdx, mod = None, None, None
    writeDB = False

    # NOTE: the options also outline a rough 'procedure' to use this code.
    for i, opt in enumerate(argv):

        # set dataset and options
        if opt=="-ds":
            ds = int(argv[i+1])
        if opt=="-cidx":
            ds = int(argv[i+1])
            cIdx = int(argv[i+2])
        if opt=="-m":
            mod = int(argv[i+1])
        if opt=="-db":
            print("Adding results to DB ...")
            writeDB = True

        # scan datasets for trapThresh and HV changes.
        if opt=="-set":
            settingsMgr(ds,cIdx,mod,writeDB)
        if opt == "-d":
            dumpSettings(ds)

        # create the settings dicts used by the DetInfo object.
        if opt=="-fillDetInfo":
            fillDetInfo()

        # since the CPD/detectorSN never changes, print a list for DetInfo.
        if opt=="-detID":
            getDetIDs()

        # make a list of sub-ranges, usable by job-panda::runAutoThresh
        if opt=="-sub":
            getSubRanges()

        # fill the DB with threshold values for each bkgIdx/sub index
        if opt=="-fillThreshDB":
            fillThreshDB()

        if opt=="-getThreshDB":
            getThreshDB()

        # check DB cuts
        if opt=="-cov":
            checkDBCutCoverage(argv[i+1],argv[i+2]) # ds, cutType


    # compareCoverage()
    # spotCheckThresh()
    # loadSubRanges()


def compareCoverage():
    # seeing some discrepancies between the old version of getSettings (just used bkg lists)
    # and the new version, which uses cal ranges.
    # w/o accessing the data itself, make sure that we're covering all the runs we THINK we're covering.

    # new method (by calIdx)
    coveredRunsNew = []
    for ds in [0,1,2,3,4,5,6]:
        for key in cal.GetKeys(ds):
            mod = -1
            if "m1" in key: mod = 1
            if "m2" in key: mod = 2
            for cIdx in range(cal.GetIdxs(key)):
                dbKeyTH = "trapThr_%s_c%d" % (key, cIdx)
                runList = []
                calLo, calHi = cal.GetCalRunCoverage(key,cIdx)
                for run in bkg.getRunList(ds):
                    if (calLo <= run <= calHi):
                        runList.append(run)
                        coveredRunsNew.append(run)
                # calList = cal.GetCalList(key,cIdx)
                # runList.append(calList[0])
                # runList = sorted(runList)

    coveredRunsOld = []
    for ds in [0,1,2,3,4,5,6]:
        runList = bkg.getRunList(ds)
        for run in runList:
            coveredRunsOld.append(run)

    # find runs that are in one list but not the other
    for run in list(set(coveredRunsOld+coveredRunsNew)):
        if run not in coveredRunsOld:
            print("run %d not in old run list!" % run)
            break
        if run not in coveredRunsNew:
            print("run %d not in new run list!" % run)
            break

    print("Run lists match!")


def spotCheckThresh():
    from ROOT import GATDataSet, TFile, MJTChannelMap, MJTChannelSettings

    run, ds, sub, cIdx = 5889, 0, 70, 23
    checkDet = 122

    gds = GATDataSet()
    runPath = gds.GetPathToRun(run,GATDataSet.kGatified)
    tf = TFile(runPath)
    chSet = tf.Get("ChannelSettings")
    chMap = tf.Get("ChannelMap")
    for ch in chSet.GetEnabledIDList():

        detName = chMap.GetString(ch, "kDetectorName")
        detCPD = chMap.GetString(ch,"kStringName")+"D"+str(chMap.GetInt(ch,"kDetectorPosition"))
        detID = ''.join(i for i in detCPD if i.isdigit())

        # access threshold and HV settings
        gretCrate = chMap.GetInt(ch,"kVME")
        gretCard = chMap.GetInt(ch,"kCardSlot")
        gretHG = chMap.GetInt(ch,"kChanHi")
        gretLG = chMap.GetInt(ch,"kChanLo")
        threshHG = chSet.GetInt("TRAP Threshold",gretCrate,gretCard,gretHG,"ORGretina4MModel")
        threshLG = chSet.GetInt("TRAP Threshold",gretCrate,gretCard,gretLG,"ORGretina4MModel")

        if detCPD=='C1P2D2':
            print(ch, detCPD, detID, threshHG)


def settingsMgr(dsIn=None, subIn=None, modIn=None, writeDB=False):
    """ ./chan-sel.py [-ds [dsNum]] [-cidx [dsNum] [cIdx]] [-m [modNum]] -set
    Manages which ds and cIdx's we call getSettings for.
    """
    global pMons, detCH # these are updated within getSettings

    # loop over datasets
    for ds in [0,1,2,3,4,5,6]:
        if dsIn is not None and ds!=dsIn:
            continue

        # these are updated in GetSettings and just printed at the end here (no need to put in DB.)
        pMons = set()
        detCH = {d:[] for d in det.allDets} # analysis channel (HG) setting

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

                # now that ds, cIdx, and module are determined, we can call getSettings
                getSettings(ds, key, mod, cIdx, writeDB)

        print("DS:%d settings found." % ds)
        np.savez("./data/ds%d_detChans.npz" % ds, detCH, pMons)


def getSettings(ds, key, mod, cIdx, writeDB=False):
    """
    Scan datasets for trapThresh and HV changes.
    Go by cIdx and write entries to calDB-v2.json.

    TRAP Thresholds:
        {"key":"trapThr_[key]_c[cIdx]", "value": {'det':[(run1,thr1),(run2,thr2)...]} }
    (det is CPD of detector, obtained from DetInfo.allDets)

    HV Thresholds:
        {"key":"hvBias_[key]_c[cIdx]", "value": {det:[(run1,thr1),(run2,thr2)...]} }
    """
    from ROOT import GATDataSet, TFile, MJTChannelMap, MJTChannelSettings
    global pMons, detCH

    # this is the data we're gonna save
    dbKeyTH = "trapThr_%s_c%d" % (key, cIdx)
    dbKeyHV = "hvBias_%s_c%d" % (key, cIdx)
    detTH = {d:[] for d in det.allDets} # trap threshold setting
    detHV = {d:[] for d in det.allDets} # HV bias setting

    # We only want to access the cal and GOOD bkg runs that fall in this run list,
    # and also assert that HV values do NOT change during calibration runs.  (why would they?)
    # These decisions are made to save processing time.
    runList = []
    calLo, calHi = cal.GetCalRunCoverage(key,cIdx)
    for run in bkg.getRunList(ds):
        if (calLo <= run <= calHi):
            runList.append(run)
    calList = cal.GetCalList(key,cIdx)
    runList.append(calList[0])
    runList = sorted(runList) # if the cal run isn't before the bkg runs, it's not the first run in this list

    print("\nDS%d M%d (%s) cIdx %d (%d runs) %s %s" % (ds,mod,key,cIdx,len(runList),dbKeyTH,dbKeyHV))

    # use GDS once just to pull out the path.
    gds = GATDataSet()
    runPath = gds.GetPathToRun(runList[0],GATDataSet.kGatified)
    filePath = '/'.join(runPath.split('/')[:-1])

    # track time per run so we can identify slowdowns on PDSF
    start = time.time()

    # begin loop over runs
    for idx, run in enumerate(runList):

        # print progress
        # f = np.fabs(100*idx/len(runList) % 10)
        # if f < 0.5:
        #     print("%d/%d (run %d) %.1f%% done." % (idx, len(runList), run, 100*idx/len(runList)))

        # make sure file exists and it's not blind before trying to load it
        fname = filePath + "/mjd_run%d.root" % run
        if not os.path.isfile(fname):
            # print("Couldn't find run",run,":",fname)
            continue
        if not os.access(fname, os.R_OK):
            # print("File is blind:",fname)
            continue
        tf = TFile(fname)

        # make sure ChannelSettings and ChannelMap exist in this file
        objs = [key.GetName() for key in tf.GetListOfKeys()]
        if "ChannelSettings" not in objs or "ChannelMap" not in objs:
            print("Settings objects not found in file:",fname)
            tf.Close()
            continue
        chSet = tf.Get("ChannelSettings")
        chMap = tf.Get("ChannelMap")

        # load pulser monitor channel list
        chPulser = chMap.GetPulserChanList()
        chPulser = [chPulser[i] for i in range(len(chPulser))]
        if ds == 1: chPulser.extend([674, 675, 677]) # 674, 675, 677 are not in the MJTChannelMap's due to a bug
        pMons.update(chPulser)

        # loop over enabled channels
        thrReset = False
        for ch in chSet.GetEnabledIDList():

            detName = chMap.GetString(ch, "kDetectorName")
            detCPD = chMap.GetString(ch,"kStringName")+"D"+str(chMap.GetInt(ch,"kDetectorPosition"))
            detID = ''.join(i for i in detCPD if i.isdigit())

            # skip pulser monitors / special channels
            if detID == '0':
                if ch not in pMons:
                    print("ch %d not in pulserMons list!  Adding it ..." % ch)
                    pMons.add(ch)
                continue

            # access threshold and HV settings
            gretCrate = chMap.GetInt(ch,"kVME")
            gretCard = chMap.GetInt(ch,"kCardSlot")
            gretHG = chMap.GetInt(ch,"kChanHi")
            gretLG = chMap.GetInt(ch,"kChanLo")
            threshHG = chSet.GetInt("TRAP Threshold",gretCrate,gretCard,gretHG,"ORGretina4MModel")
            threshLG = chSet.GetInt("TRAP Threshold",gretCrate,gretCard,gretLG,"ORGretina4MModel")
            # print(ch, detCPD, detID, gretCrate, gretCard, gretHG, gretLG)
            # if detCPD=='C1P2D2':
                # print(ch, detCPD, detID, threshHG)

            hvCrate = chMap.GetInt(ch,"kHVCrate")
            hvCard = chMap.GetInt(ch,"kHVCard")
            hvChan = chMap.GetInt(ch,"kHVChan")
            hvTarget = chMap.GetInt(ch,"kMaxVoltage")
            hvActual = chSet.GetInt("targets",hvCrate,hvCard,hvChan,"OREHS8260pModel")
            # print(ch, detCPD, detID, hvCrate, hvCard, hvChan, hvTarget, hvActual)

            # fill our data objects
            thisTH = (run, threshHG) # NOTE: HG only
            thisHV = (run, hvActual)

            if len(detTH[detID]) == 0:
                detTH[detID].append(thisTH)
            if len(detHV[detID]) == 0:
                detHV[detID].append(thisHV)

            # if we detect a difference, add a new (run,val) pair

            if len(detTH[detID]) != 0:
                prevTH = detTH[detID][-1]
                if thisTH[1] != prevTH[1]:
                    detTH[detID].append(thisTH)
                    thrReset = True

            if len(detHV[detID]) != 0:
                prevHV = detHV[detID][-1]
                if thisHV[1] != prevHV[1]:
                    detHV[detID].append(thisHV)
                    print("Found HV diff, run %d.  C%sP%sD%s, %dV -> %dV" % (run,detID[0],detID[1],detID[2],prevHV[1],thisHV[1]))

            # skip LG channels
            if ch%2!=0: continue
            thisCH = (run, ch)

            if len(detCH[detID]) == 0:
                detCH[detID].append(thisCH)

            if len(detCH[detID]) != 0:
                prevCH = detCH[detID][-1]
                if thisCH[1] != prevCH[1]:
                    detCH[detID].append(thisCH)

        if thrReset:
            print("Thresholds reset, run %d" % run)

        tf.Close()
        # return # limit to 1 run

    # note the time elapsed
    timeElapsed = time.time()-start
    print("Elapsed: %.4f sec, %.4f sec/run" % (timeElapsed, timeElapsed/len(runList)))

    # debug: print the values
    # print(dbKeyTH)
    # for val in sorted(detTH):
    #     if len(detTH[val])>0:
    #         print(val, detTH[val])

    # fill the DB
    if writeDB:
        calDB = db.TinyDB("%s/calDB-v2.json" % dsi.latSWDir)
        pars = db.Query()
        dsi.setDBRecord({"key":dbKeyTH, "vals":detTH}, forceUpdate=True, calDB=calDB, pars=pars)
        dsi.setDBRecord({"key":dbKeyHV, "vals":detHV}, forceUpdate=True, calDB=calDB, pars=pars)
        print("DB filled.")


def dumpSettings(ds):

    print("Dumping settings for DS:%d" % ds)

    # detCH and pMons
    # f = np.load("./data/ds%d_detChans.npz" % ds)
    # detCH, pMons = f['arr_0'].item(), f['arr_1']
    #
    # print("Pulser monitor channels:")
    # print(pMons)

    # print("Detector analysis channels:")
    # for key in detCH:
    #     print(key, detCH[key])

    # get HV and TF vals from DB with a regex
    calDB = db.TinyDB("%s/calDB-v2.json" % dsi.latSWDir)
    pars = db.Query()
    #
    # print("DB Threshold values:")
    # thrList = calDB.search(pars.key.matches("trapThr_ds%d" % ds))
    # for idx in range(len(thrList)):
    #     key = thrList[idx]['key']
    #     detTH = thrList[idx]['vals']
    #     print(key)
    #     # for d in detTH:
    #         # print(d, detTH[d])

    # print("DB HV values:")
    # hvList = calDB.search(pars.key.matches("hvBias_ds%d" % ds))
    # for idx in range(len(hvList)):
    #     key = hvList[idx]['key']
    #     detHV = hvList[idx]['vals']
    #     print(key)
    #     # for d in detHV:
    #         # print(d, detHV[d])

    # get a particular key
    dbKeyTH = "trapThr_ds0_m1_c23"
    dbValTH = dsi.getDBRecord(dbKeyTH,calDB=calDB,pars=pars)

    # debug: print the values
    for val in sorted(dbValTH):
        if len(dbValTH[val])>0:
            print(val, dbValTH[val])


def fillDetInfo():
    """ ./chan-sel.py -fill
    Summarize the results of getSettings in LAT/data/runSettings-v2.npz.
    Create a file accessible by the DetInfo object in dsi.py
    It contains dictionaries of TRAP threshold, HV settings,
    detector/channel mapping, and pulser monitor lists,
    that span an entire dataset (not broken down into separate calIdx's.)
    # FORMAT: {ds : {'det' : [(run1,val1),(run2,val2)...]} }
    """
    # 1. maps of analysis channel to cpd, and pulser monitor channels
    detCH, pMons = {}, {}
    for ds in [0,1,2,3,4,5,6]:
        f = np.load("%s/data/ds%d_detChans.npz" % (os.environ['LATDIR'], ds))
        detCH[ds] = f['arr_0'].item()
        pMons[ds] = f['arr_1'].item()

    # 2. maps of HV and TRAP threshold settings are stored in the DB.
    # make them global, and move them to the runSettings file.
    # FORMAT: {ds : {'det' : [(run1,val1),(run2,val2)...]} }
    detHV, detTH = {}, {}

    # load all possible values, as in settingsMgr
    detDB = db.TinyDB("%s/calDB-v2.json" % dsi.latSWDir)
    detPars = db.Query()
    cal = dsi.CalInfo()
    for ds in [0,1,2,3,4,5,6]:
    # for ds in [3]:
        print("scanning ds",ds)
        detTH[ds] = {}
        detHV[ds] = {}
        for key in cal.GetKeys(ds):
            mod = -1
            if "m1" in key: mod = 1
            if "m2" in key: mod = 2
            for cIdx in range(cal.GetIdxs(key)):

                # load the DB records
                dbKeyTH = "trapThr_%s_c%d" % (key, cIdx)
                dbValTH = dsi.getDBRecord(dbKeyTH,calDB=detDB,pars=detPars)

                dbKeyHV = "hvBias_%s_c%d" % (key, cIdx)
                dbValHV = dsi.getDBRecord(dbKeyHV,calDB=detDB,pars=detPars)

                # debug: print the record
                # for val in sorted(dbValTH):
                    # if len(dbValTH[val])>0:
                        # print(val, dbValTH[val])
                # return

                # fill the first value
                if len(detTH[ds])==0:
                    detTH[ds] = dbValTH
                    detHV[ds] = dbValHV
                    continue

                # check for new threshold values.
                for cpd in detTH[ds]:
                    nOld, nNew = len(detTH[ds][cpd]), len(dbValTH[cpd])

                    # detector just came online
                    if nOld==0 and nNew>0:
                        detTH[ds][cpd] = dbValTH[cpd]
                        continue
                    # detector still offline
                    if nOld==0 and nNew==0:
                        continue
                    # detector just went offline
                    if nOld>0 and nNew==0:
                        continue

                    # check last run/trap pair against each new one
                    prevRun, prevTH = detTH[ds][cpd][-1][0], detTH[ds][cpd][-1][1]
                    for val in dbValTH[cpd]:
                        thisRun, thisTH = val[0], val[1]
                        if thisTH != prevTH:
                            detTH[ds][cpd].append([thisRun,thisTH])
                        prevTH = thisTH

                # check for new HV values.
                for cpd in detHV[ds]:

                    nOld, nNew = len(detHV[ds][cpd]), len(dbValHV[cpd])

                    # detector just came online
                    if nOld==0 and nNew>0:
                        detHV[ds][cpd] = dbValHV[cpd]
                        continue
                    # detector still offline
                    if nOld==0 and nNew==0:
                        continue
                    # detector just went offline
                    if nOld>0 and nNew==0:
                        continue

                    # check last run/trap pair against each new one
                    prevRun, prevHV = detHV[ds][cpd][-1][0], detHV[ds][cpd][-1][1]
                    for val in dbValHV[cpd]:
                        thisRun, thisHV = val[0], val[1]
                        if thisHV != prevHV:
                            print("found HV diff.  cpd %d  prev %dV (run %d)  new %dV (run %d)" % (cpd, prevHV, prevRun, thisHV, thisRun))
                            detHV[ds][cpd].append([thisRun,thisHV])
                        prevHV = thisHV

                # return

    # # load the old file and compare
    # # GOAL: improve on this file.
    # # f = np.load("%s/data/runSettings.npz" % dsi.latSWDir)
    # # detHVOld = f['arr_0'].item()
    # # detTHOld = f['arr_1'].item()
    # # detCHOld = f['arr_2'].item()
    # # pMonsOld = f['arr_3'].item()
    #
    # ds = 3
    # print("old results, ds",ds)
    # for cpd in sorted(detTHOld[ds]):
    #     if cpd!="122":continue
    #     if len(detTHOld[ds][cpd]) > 0:
    #         print(cpd, detTHOld[ds][cpd])
    #
    # # for ds in [0,1,2,3,4,5,6]:
    # print("thresh results, ds:",ds)
    # for cpd in sorted(detTH[ds]):
    #     # if cpd!=122:continue
    #     if len(detTH[ds][cpd]) > 0:
    #         print(cpd, detTH[ds][cpd])


    np.savez("%s/data/runSettings-v2.npz" % dsi.latSWDir,detHV,detTH,detCH,pMons)


def getDetIDs():
    """ Since the detector CPD position - detector serial number
    mapping has never changed, load a channel map and print it out.
    ==> copied this output into DetInfo::detID
    """
    from ROOT import TFile, GATDataSet, MJTChannelMap
    run = bkg.getRunList(5,1)[0]
    gds = GATDataSet(run)
    chMap = gds.GetChannelMap()
    # chMap.DumpDetectorNames()

    dets = det.allDets
    for d in dets:
        detName = chMap.GetDetectorName(int(d[0]), int(d[1]), int(d[2]))

        # now match it to the IDs used in DataSetInfo.cc::Load(X)DetectorMap
        detID = '1' if detName[0]=="P" else '2'
        detID += detName[1:]
        tmp = list(detID)
        if detID[0] == '1':
            if tmp[-1] == "A": tmp[-1] = '0'
            if tmp[-1] == "B": tmp[-1] = '1'
            if tmp[-1] == "C": tmp[-1] = '2'
            detID = ''.join(tmp)

        print("'%s':%s," % (d, detID))


def getSubRanges():
    """ Detect periods where:
        1) Thresholds changed ==> identify bkgIdx & output sub-ranges s/t we can re-run auto-thresh.
        TODO: 2) HV changed - identify bkgIdx and calIdx
    Use the result to make some decisions about data ranges.
    """
    verbose = True
    writeOutput = False

    # this is what we output (for other programs)
    thRanges = [] # these will be inputs to auto-thresh (thresholds)
    hvRanges = [] # these will be inputs to PSA cut tuning (calibration)

    # loop over datasets
    # for ds in [0,1,2,3,4,"5A","5B","5C",6]:
    for ds in [4]:
        if ds in ["5A","5B","5C"]:
            dsNum = 5
        else:
            dsNum = ds
        if verbose: print("DS:",ds)

        # loop over subsets
        ranges = bkg.getRanges(ds)
        for sub in ranges:
            # print(ds,sub)
            # if sub!=37: continue
            runList = bkg.getRunList(ds,sub)

            # if verbose: print("sub:%d  runLo %d  runHi %d  nRuns %d" % (sub, runList[0], runList[-1], len(runList)))

            # continue

            # these are how we have to divide the given subset.
            thLo, thHi = runList[0], runList[-1]
            hvLo, hvHi = runList[0], runList[-1]

            # bookkeeping vars
            thReset, hvReset = False, False
            nHV, nTh = 0, 0
            prevTh = det.getTrapThreshAtRun(dsNum, runList[0]) # these interpolate between runs
            prevHV = det.getHVAtRun(dsNum, runList[0])         #
            tmpThreshRanges = []
            tmpCutTuneRanges = []

            # loop over runs
            for idx, run in enumerate(runList):
                # print(run)
                nHV += 1
                nTh += 1

                th = det.getTrapThreshAtRun(dsNum,run)
                hv = det.getHVAtRun(dsNum,run)

                # check for thresh change
                if th != prevTh:
                    nDet = 0
                    for cpd in sorted(th):
                        if prevTh[cpd] != th[cpd]:
                            nDet += 1
                            # this is REALLY verbose.  but allows you to see the actual changing values.
                            # if verbose: print("  Thr change, run %d, det %s (%s), %d -> %d" % (run,cpd,det.getCPDChan(dsNum,cpd),prevTh[cpd],th[cpd]))
                    thReset = True
                    thHi = runList[idx-1]
                    # print("  thLo %d  thHi %d  (%d runs)" % (thLo, thHi, nTh-1))
                    tmpThreshRanges.append((ds,sub,thLo,thHi,nTh-1))
                    thLo = runList[idx]
                    prevTh = th
                    nTh = 1

                # check for HV change
                if hv != prevHV:
                    m1Chg, m2Chg = False, False
                    for cpd in sorted(hv):
                        if prevHV[cpd] != hv[cpd]:
                            if verbose: print("  HV change, run %d, det %s (%s), %dV -> %dV.  nRuns: %d" % (run,cpd,det.getCPDChan(dsNum,cpd),prevHV[cpd],hv[cpd],nHV-1))
                            if cpd[0]=='1': m1Chg = True
                            if cpd[0]=='2': m2Chg = True
                    hvReset = True
                    hvHi = runList[idx-1]
                    # print("  hvLo %d  hvHi %d  (%d runs)" % (hvLo, hvHi, nHV-1))

                    # identify calIdx of this run
                    if ds not in ["5A","5B"]:
                        calKey = cal.GetKeys(dsNum)[0] # default: only 1 value
                        calIdx = cal.GetCalIdx(calKey, run)
                    else:
                        m1Cal = cal.GetCalIdx("ds%d_m1" % dsNum, run)
                        m2Cal = cal.GetCalIdx("ds%d_m2" % dsNum, run)
                        calIdx = (m1Cal, m2Cal, m1Chg, m2Chg) # 4 vals, distinguish between modules.  ugh

                    tmpCutTuneRanges.append((ds,sub,hvLo,hvHi,nHV-1,calIdx))

                    hvLo = runList[idx]
                    prevHV = hv
                    nHV = 1

            # final cleanup and output
            if hvReset:
                hvHi = runList[-1]
                # print("  hvLo %d  hvHi %d  (%d runs)" % (hvLo, hvHi, nHV))
                tmpCutTuneRanges.append((ds, sub, hvLo,hvHi,nHV,calIdx))

                if verbose:
                    print("  HV sub-ranges:")
                    for val in tmpCutTuneRanges:
                        print(" ",val)

            if thReset:
                thHi = runList[-1]
                # print("  thLo %d  thHi %d  (%d runs)" % (thLo, thHi, nTh))
                tmpThreshRanges.append((ds, sub,thLo,thHi,nTh))

                if verbose:
                    print("  thresh sub-ranges:")
                    for val in tmpThreshRanges:
                        print(" ",val)

            # extend master lists
            if len(tmpThreshRanges) > 0: thRanges.extend(tmpThreshRanges)
            if len(tmpCutTuneRanges) > 0: hvRanges.extend(tmpCutTuneRanges)

    # if verbose:
    print("----Final Ranges----")

    if len(thRanges) > 0:
        print("Auto-thresh sub-ranges:\n(ds, sub, runLo, runHi, nRuns)")
        for val in thRanges:
            print(val)

    if len(hvRanges) > 0:
        print("\nHV change sub-ranges:\n(ds, sub, runLo, runHi, nRuns, calIdx)")
        for val in hvRanges:
            print(val)

    # save to npz file
    if writeOutput:
        thRanges = np.array(thRanges,dtype=object)
        hvRanges = np.array(hvRanges,dtype=object)
        np.savez("./data/thrHV_subRanges.npz",thRanges,hvRanges)
        if verbose: print("Saved file.")


def loadSubRanges():

    f = np.load("./data/thrHV_subRanges.npz")
    thRanges = f['arr_0']
    hvRanges = f['arr_1']

    for val in hvRanges:
        print(*val)


def fillThreshDB():
    """ Fill threshold records in to LAT/calDB-v2.json.
        keys: thresh_ds[i]_bkg[j]_sub[k]
        vals: {[chan]:[50 pct mean, sigma, isGood (0 good, 1 bad)]}
    """
    from ROOT import TFile, TTree

    # Do we actually want to write new values?  Or just print stuff out?
    fillDB = True

    calDB = db.TinyDB("%s/calDB-v2.json" % dsi.latSWDir)
    pars = db.Query()
    bkg = dsi.BkgInfo()

    # loop over datasets and bkgIdx
    for ds in [0,1,2,3,4,"5A","5B","5C",6]:
        dsNum = ds if isinstance(ds, int) else 5
        goodChans = det.getGoodChanList(dsNum)

        for bkgIdx in bkg.getRanges(ds):

            # ==== loop over sub-ranges (when TF was run) ====

            rFirst, rLast = bkg.getRanges(ds)[bkgIdx][0], bkg.getRanges(ds)[bkgIdx][-1]

            subRanges = bkg.GetSubRanges(ds,bkgIdx)
            if len(subRanges) == 0: subRanges.append((rFirst, rLast))

            for subIdx, (runLo, runHi) in enumerate(subRanges):

                # Load threshold table
                if len(subRanges) > 1:
                    fname = "%s/threshDS%d_%d_%d_%d.root" % (dsi.threshDir, dsNum, bkgIdx, runLo, runHi)
                else:
                    fname = "%s/threshDS%d_%d.root" % (dsi.threshDir, dsNum, bkgIdx)
                if not os.path.isfile(fname):
                    print("Couldn't find file:",fname)
                    return
                tf = TFile(fname)
                tt = tf.Get("threshTree")
                try:
                    n = tt.GetEntries()
                except AttributeError:
                    print("skipped",fname)
                    continue
                if (n!=1):
                    print("Hmm, %d thresh table entries? %s" % (n, fname))
                    return
                tt.GetEntry(0)

                rLo, rHi = tt.runMin, tt.runMax
                key = "thresh_ds%d_bkg%d_sub%d" % (dsNum, bkgIdx, subIdx)
                vals = {}
                print("")
                print(key, runLo, runHi)
                print("chan CPD  g  threshKeV sig       ADC    err     E=0   err     nThr     nNoise   status  note")

                # loop over channels
                nGood, nTot = 0, 0
                for i in range(tt.channelList.size()):

                    # keep only HG channels, exclude pulser monitors
                    chan = tt.channelList.at(i)
                    if chan%2!=0 or chan in det.getPMon(dsNum):
                        continue

                    # load results
                    isGood = 1 if chan in goodChans else 0
                    thrKeV, thrSig = tt.threshCal.at(i), tt.sigmaCal.at(i)
                    thrADC = tt.threshADC.at(i)
                    thrADCErr = tt.threshADCErr.at(i)
                    thrStatus = tt.threshFitStatus.at(i)
                    thrEvts = tt.numTrigger.at(i)
                    sigADC = tt.sigmaADC.at(i)
                    sigADCErr = tt.sigmaADCErr.at(i)
                    sigStatus = tt.sigmaFitStatus.at(i)
                    sigEvts = tt.numNoise.at(i)
                    calScale = tt.CalOffset.at(i)
                    calOffset = tt.CalScale.at(i)

                    # Error handling
                    isBad = False
                    status = ""
                    if thrStatus==999999 or sigStatus==999999 or thrEvts < 10 or sigEvts < 10:
                        status = "Not enough events"
                        isBad = True
                    elif (0 < thrStatus < 999999) or (0 < sigStatus < 999999):
                        status = "auto-fit fail"
                        isBad = True
                    elif int(thrKeV)==99999 or int(thrSig)==99999:
                        status = "Fit fail"
                        isBad = True
                    elif int(thrKeV)==999999:
                        status = "Bad energy calibration"
                        isBad = True
                    elif thrKeV < 0.3:
                        status = "Unphysical threshold"
                        isBad = True
                    elif thrEvts < 100 or sigEvts < 100:
                        status = "Low events"
                        # this is ok
                    if not isBad: nGood += 1
                    nTot += 1

                    # pretty print the results table
                    if int(thrKeV) > 99998:
                        print("%d  %s  %d  %-7.0f s %-7.0f  %-4.3f  %-6.3f  %.2f  %-6.2f  %-8d %-8d %d %d %d %s" % (chan, det.getChanCPD(dsNum,chan), isGood, thrKeV, thrSig, thrADC, thrADCErr, sigADC, sigADCErr, thrEvts, sigEvts, thrStatus,sigStatus, int(isBad), status))
                    else:
                        print("%d  %s  %d  %-7.3f s %-7.3f  %-4.3f  %-6.3f  %.2f  %-6.2f  %-8d %-8d %d %d %d %s" % (chan, det.getChanCPD(dsNum,chan), isGood, thrKeV, thrSig, thrADC, thrADCErr, sigADC, sigADCErr, thrEvts, sigEvts, thrStatus, sigStatus, int(isBad), status))

                    # fill the dict vals
                    vals[chan] = [float("%.5f" % thrKeV), float("%.5f" % thrSig), int(isBad)]

                print("good detectors: %d/%d" % (nGood, nTot))

                # fill the DB
                if fillDB:
                    dsi.setDBRecord({"key":key, "vals":vals}, forceUpdate=True, calDB=calDB, pars=pars)

                tf.Close()
                # return


def getThreshDB():
    calDB = db.TinyDB("%s/calDB-v2.json" % dsi.latSWDir)
    pars = db.Query()
    bkg = dsi.BkgInfo()

    # loop over datasets
    for ds in [0,1,2,3,4,5,6]:
        dsNum = ds if isinstance(ds, int) else 5
        goodChans = det.getGoodChanList(dsNum)

        for bkgIdx in bkg.getRanges(ds):

            # ==== loop over sub-ranges (when TF was run) ====
            rFirst, rLast = bkg.getRanges(ds)[bkgIdx][0], bkg.getRanges(ds)[bkgIdx][-1]

            subRanges = bkg.GetSubRanges(ds,bkgIdx)
            if len(subRanges) == 0: subRanges.append((rFirst, rLast))

            for subIdx, (runLo, runHi) in enumerate(subRanges):

                key = "thresh_ds%d_bkg%d_sub%d" % (dsNum, bkgIdx, subIdx)
                print(key)


"""
    ** For BOTH of these, the DB has already been updated to reflect these changes **

    fitSlo detectors cut:  (lat2::setSloCut)
    # Detectors to cut.  This is from inspecting the 'makePlots' output.
    # Criteria: nBin must be higher than 4 (50% statistical error in the bins)
    # NOTE: these detectors could be brought back if we included the DS6 cal runs
    # to bump up the stats in M2.
    cutDets = ['211','214','221','261','274']

    fitSlo can ALSO be bad in a calIdx when fsCut and nBin are == -1.
    dbVals[ch] = [fsCut, fs200, nBin]


    riseNoise chan/calIdx cut: (lat2::badRiseChans)
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

def checkDBCutCoverage(ds, cutType):
    """ ./chan-sel.py -cov [ds] [cutType]
    Check which bkgIdx's we have good cut values for.
    This is a 5-layered loop, which is pretty badass.
    """
    # NOTE: input for DS5 must be 5A, 5B, or 5C, not 5.
    dsNum = int(ds[0]) if isinstance(ds, str) else int(ds)
    print("Getting cut run/ch vals for DS-%s (%d) ..." % (ds, dsNum))

    if cutType not in ["fr","fs","rn","thr"]:
        print("Unknown cut type:",cutType,"... exiting ...")
        return

    # set up output file.  i wish i'd made everything uppercase
    dsLabel = str(dsNum)
    if isinstance(ds, str):
        if ds=="5A": dsLabel = "5a"
        if ds=="5B": dsLabel = "5b"
        if ds=="5C": dsLabel = "5c"
    outFile = "./data/dbCut_%s_%s.txt" % (cutType,dsLabel)
    dRanges = {} # this is what we'll write to text file at the end

    calDB = db.TinyDB('%s/calDB-v2.json' % (dsi.latSWDir))
    pars = db.Query()
    dsMap = bkg.dsMap() # number of bIdx's
    bkgRanges = bkg.getRanges(ds)

    # 1. loop over modules
    mods = [1]
    if dsNum == 4: mods = [2]
    if dsNum == 5: mods = [1,2]
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
            bkgDict, calDict, bkgCov, calCov = dsi.GetDBCuts(ds,bIdx,mod,cutType,calDB,pars,False)

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

                        print("%s  bIdx %d  sbIdx %d  cIdx %d  ch %d  th %d  fs %d  rn %d  %s" % (calKey,bIdx,sbIdx,cIdx,ch,int(goodThr),int(goodSlo),int(goodRise),excludeMsg))

    # Create output suitable for ds_livetime
    with open(outFile, 'w') as f:
        for ch in sorted(dRanges):
            if len(dRanges[ch]) > 0:
                outStr = "%d" % ch
                for r in wl.niceList(dRanges[ch], "%d", "i"): outStr += " %d" % r
                outStr += "\n"
                f.write(outStr)
    print("Wrote cut file:",outFile)


if __name__=="__main__":
    main(sys.argv[1:])
