#!/usr/bin/env python3
"""
===================== LAT/chan-sel.py ======================
Perform necessary run/channel selection bookkeeping for
    LAT analysis, before and after cuts are applied.

Uses:
- Check HV settings and thresholds for all bkg/cal runs
  --> Fill objects in dsi.py

- Fill threshold results in to DB

- For LAT files: verify we don't see any bad/veto detectors

- After LAT3::ApplyChannelCuts is done:
  --> What did we cut out?
      - From outright channel selection (2 DS5 beges)
      - From no data in cut tuning
      - Files w/ no statistics after cuts (not actually a cut ...)
- Livetime:
  --> reuse chan-sel-v1.py::parseLivetimeOutput

Note: LAT/sandbox/chan-sel-v1.py has some good
routines we should re-use (incl chan. selection graphics)

===================== C. Wiseman (USC) =====================
"""
import sys, os
import tinydb as db
import numpy as np

import dsi
bkg = dsi.BkgInfo()
det = dsi.DetInfo()
cal = dsi.CalInfo()

def main(argv):

    # these have a tendency to hang.
    # if you want a live-updating log file, submit with something like
    # script -c './chan-sel.py -t -v' -f logs/chanSel.txt
    for i,opt in enumerate(argv):
        if opt == "-t": getRunSettings()
        if opt == "-v": getRunHVSettings()

    # getRunSettings()
    # getRunHVSettings()
    # dumpRunSettings()
    # getDetIDs()
    # checkSubRanges()
    # loadSubRanges()
    # fillThreshDB()
    # ds4Check()


def getRunSettings():
    """ ./chan-sel.py -t
    This takes ~2 hours to run on PDSF if /project IO is good.
    Results are used to populate dsi::DetInfo.
    """
    from ROOT import GATDataSet, TFile, MJTChannelMap, MJTChannelSettings
    print("Scanning run/threshold/channel settings ...")

    # for ds in [0,1,2,3,4,5,6]:
    for ds in [6]:

        print("Scanning DS:%d" % ds)

        pMons = set()
        detTH = {d:[] for d in det.allDets} # trap threshold setting
        detCH = {d:[] for d in det.allDets} # analysis channel (HG) setting

        # use GDS once just to pull out the path.
        gds = GATDataSet()
        runList = bkg.getRunList(ds)
        runPath = gds.GetPathToRun(runList[0],GATDataSet.kGatified)
        filePath = '/'.join(runPath.split('/')[:-1])
        print("Using path to files:",filePath)

        for idx, run in enumerate(runList):

            f = np.fabs(100*idx/len(runList) % 10)
            if f < 0.5:
                print("%d/%d, %.2f%% done." % (idx, len(runList), 100*idx/len(runList)))

            runPath = filePath + "/mjd_run%d.root" % run
            tf = TFile(runPath)
            chSet = tf.Get("ChannelSettings")
            chMap = tf.Get("ChannelMap")
            # chMap.DumpChannelMap()
            # chSet.DumpSettings() # dump to a file to see syntax for HV and TRAP

            chEnabled = chSet.GetEnabledIDList()
            chSpecial = chMap.GetSpecialChanMapString()
            chPulser = chMap.GetPulserChanList()
            chPulser = [chPulser[i] for i in range(len(chPulser))]
            if ds == 1: chPulser.extend([674, 675, 677]) # 674, 675, 677 are not in the MJTChannelMap's due to a bug
            pMons.update(chPulser)

            for ch in chEnabled:

                detName = chMap.GetString(ch, "kDetectorName")
                detCPD = chMap.GetString(ch,"kStringName")+"D"+str(chMap.GetInt(ch,"kDetectorPosition"))
                detID = ''.join(i for i in detCPD if i.isdigit())

                # skip pulser monitors / special channels
                if detID == '0':
                    if ch not in pMons:
                        print("ch %d not in pulserMons list!  Adding it ..." % ch)
                        pMons.add(ch)
                    continue

                gretCrate = chMap.GetInt(ch,"kVME")
                gretCard = chMap.GetInt(ch,"kCardSlot")
                gretHG = chMap.GetInt(ch,"kChanHi")
                gretLG = chMap.GetInt(ch,"kChanLo")
                threshHG = chSet.GetInt("TRAP Threshold",gretCrate,gretCard,gretHG,"ORGretina4MModel")
                threshLG = chSet.GetInt("TRAP Threshold",gretCrate,gretCard,gretLG,"ORGretina4MModel")
                # print(ch, detCPD, detID, gretCrate, gretCard, gretHG, gretLG)

                # fill our data objects
                thisTH = (run, threshHG) # NOTE: HG only

                if len(detTH[detID]) == 0:
                    detTH[detID].append(thisTH)

                # if we detect a difference, add a new (run,val) pair
                if len(detTH[detID]) != 0:
                    prevTH = detTH[detID][-1]
                    if thisTH[1] != prevTH[1]:
                        detTH[detID].append(thisTH)

                # skip LG channels
                if ch%2!=0: continue
                thisCH = (run, ch)

                if len(detCH[detID]) == 0:
                    detCH[detID].append(thisCH)

                if len(detCH[detID]) != 0:
                    prevCH = detCH[detID][-1]
                    if thisCH[1] != prevCH[1]:
                        detCH[detID].append(thisCH)

            tf.Close()

        # save output for dsi.py to use
        np.savez("./data/settingsBkg_ds%d.npz" % ds, detTH,detCH,pMons)


def getRunHVSettings():
    """ ./chan-sel.py -v
    Decided that we need the exact run numbers
    that the HV changes, not just the bkgIdx they correspond to.
    This is basically the same algorithm as getRunSettings
    except that it tries to access every run in the DS and only looks at HV.
    """
    from ROOT import GATDataSet, TFile, MJTChannelMap, MJTChannelSettings
    print("Scanning HV settings ...")

    # for ds in [0,1,2,3,4,5,6]:
    for ds in [6]:

        print("Scanning DS:%d" % ds)

        # fill this with [(run,val),...] pairs for each detector, each ds
        detHV = {d:[] for d in det.allDets} # high voltage setting

        # use GDS once just to pull out the path.
        gds = GATDataSet()
        runLo, runHi = bkg.dsRanges()[ds]
        runPath = gds.GetPathToRun(runLo,GATDataSet.kGatified)
        filePath = '/'.join(runPath.split('/')[:-1])

        # loop over all runs
        nRuns = runHi+1-runLo
        for idx, run in enumerate(range(runLo, runHi+1)):
            f = np.fabs(100*idx/nRuns % 10)
            if f < 0.1:
                print("%d/%d (run %d) %.2f%% done." % (idx, nRuns, run, 100*idx/nRuns))

            # make sure file exists and it's not blind
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

            for ch in chSet.GetEnabledIDList():

                detName = chMap.GetString(ch, "kDetectorName")
                detCPD = chMap.GetString(ch,"kStringName")+"D"+str(chMap.GetInt(ch,"kDetectorPosition"))
                detID = ''.join(i for i in detCPD if i.isdigit())
                if detID == '0':
                    continue

                hvCrate = chMap.GetInt(ch,"kHVCrate")
                hvCard = chMap.GetInt(ch,"kHVCard")
                hvChan = chMap.GetInt(ch,"kHVChan")
                hvTarget = chMap.GetInt(ch,"kMaxVoltage")
                hvActual = chSet.GetInt("targets",hvCrate,hvCard,hvChan,"OREHS8260pModel")
                # print(ch, detCPD, detID, hvCrate, hvCard, hvChan, hvTarget, hvActual)

                # fill our data objects
                thisHV = (run, hvActual)

                if len(detHV[detID]) == 0:
                    detHV[detID].append(thisHV)

                # if we detect a difference, add a new (run,val) pair
                if len(detHV[detID]) != 0:
                    prevHV = detHV[detID][-1]
                    if thisHV[1] != prevHV[1]:
                        detHV[detID].append(thisHV)
                        print("Found diff, run",run)

            tf.Close()

        # save output for dsi.py to use
        np.savez("./data/settingsHV_ds%d.npz" % ds,detHV)


def dumpRunSettings():

    # old file - to be deleted
    # f = np.load("./data/runSettings.npz")
    # detHV = f['arr_0'].item()
    # detTH = f['arr_1'].item()
    # detCH = f['arr_2'].item()
    # pMons = f['arr_3'].item()

    # f = np.load("./data/settingsHV.npz")
    # detHV = f['arr_0'].item()

    # print summary
    # print("HV settings:")
    # for ds in sorted(detHV):
    #     print("DS",ds)
    #     for det in sorted(detHV[ds]):
    #         print(det, detHV[ds][det])

    f = np.load("./data/settingsBkg.npz")
    detTH = f['arr_0'].item()
    detCH = f['arr_1'].item()
    pMons = f['arr_2'].item()

    print("HG Threshold settings:")
    for ds in sorted(detTH):
        print("DS",ds)
        for det in sorted(detTH[ds]):
            print(det, detTH[ds][det])

    print("HG Channel settings:")
    for ds in sorted(detCH):
        print("DS",ds)
        for det in sorted(detCH[ds]):
            print(det, detCH[ds][det])

    print("Pulser monitors:")
    print(pMons)


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


def checkSubRanges():
    """ Detect periods where:
        1) HV changed - identify bkgIdx and calIdx
        2) TH changed - identify bkgIdx
    Use the result to make some decisions about data ranges.
    """
    verbose = False

    # this is what we output (for other programs)
    thRanges = [] # these will be inputs to auto-thresh (thresholds)
    hvRanges = []    # these will be inputs to PSA cut tuning (calibration)

    # loop over datasets
    for ds in [0,1,2,3,4,"5A","5B","5C",6]:
        if ds in ["5A","5B","5C"]:
            dsNum = 5
        else:
            dsNum = ds
        if verbose: print("DS:",ds)

        # loop over subsets
        ranges = bkg.getRanges(ds)
        for sub in ranges:
            # if sub!=37: continue
            runList = bkg.getRunList(ds,sub)

            if verbose: print("sub:%d  runLo %d  runHi %d  nRuns %d" % (sub, runList[0], runList[-1], len(runList)))

            # these are how we have to divide the given subset.
            thLo, thHi = runList[0], runList[-1]
            hvLo, hvHi = runList[0], runList[-1]

            # bookkeeping vars
            thReset, hvReset = False, False
            nHV, nTh = 0, 0
            prevTh = det.getTrapThreshAtRun(dsNum, runList[0])
            prevHV = det.getHVAtRun(dsNum, runList[0])
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

    if verbose:
        print("\n----Final Ranges----")

        if len(thRanges) > 0:
            print("Auto-thresh sub-ranges:\n(ds, sub, runLo, runHi, nRuns)")
            for val in thRanges:
                print(val)

        if len(hvRanges) > 0:
            print("\nHV change sub-ranges:\n(ds, sub, runLo, runHi, nRuns, calIdx)")
            for val in hvRanges:
                print(val)

    # save to npz file
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


def ds4Check():
    """ For some reason the DS4 auto-thresh jobs didn't finish.
        Is it because the chains are super enormous?
    """
    import time
    from ROOT import GATDataSet, TFile, TTree

    ds, sub = 4, 13
    # ds, sub = 5, 47

    bkg = dsi.BkgInfo()
    for bkgIdx in bkg.getRanges(ds):
        if bkgIdx not in [sub]: continue
        runList = bkg.getRunList(ds, bkgIdx)

        print("DS:",ds,"bkgIdx:",bkgIdx)
        gds = GATDataSet()
        for run in runList:
            start = time.time()
            gPath = gds.GetPathToRun(run,GATDataSet.kGatified)
            bPath = gds.GetPathToRun(run,GATDataSet.kBuilt)
            print(run, "%.4f" % (time.time()-start))


            # f1 = TFile(gPath)
            # f2 = TFile(bPath)
            # t1 = f1.Get("mjdTree")
            # t2 = f2.Get("MGTree")
            # print("%d  gat %d  built %d" % (run, t1.GetEntries(), t2.GetEntries()))
            # f1.Close()
            # f2.Close()


def fillThreshDB():
    """ Fill threshold records in to LAT/calDB-v2.json.
        keys: thresh_ds[i]_bkg[j]_sub[k]
        vals: {[chan]:[50 pct mean, sigma, isGood (0 good, 1 bad)]}
    """
    from ROOT import TFile, TTree

    # Do we actually want to write new values?  Or just print stuff out?
    fillDB = False
    forceUpdate = False # set True to overwrite existing DB values

    calDB = db.TinyDB("%s/calDB-v2.json" % dsi.latSWDir)
    pars = db.Query()
    bkg = dsi.BkgInfo()

    # dsList = [ds] if ds is not None else [0,1,2,3,4,"5A","5B","5C",6]
    dsList = [4]

    # loop over datasets and bkgIdx
    for ds in dsList:
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
                print("chan CPD  g  threshKeV sig      ADC    err     E=0   err     nThr     nNoise   status  note")

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
                    # if int(thrKeV) > 99998:
                    #     print("%d  %s  %d  %-7.0f s %-7.0f  %-4.3f  %-6.3f  %.2f  %-6.2f  %-8d %-8d %d %d %d %s" % (chan, det.getChanCPD(dsNum,chan), isGood, thrKeV, thrSig, thrADC, thrADCErr, sigADC, sigADCErr, thrEvts, sigEvts, thrStatus,sigStatus, int(isBad), status))
                    # else:
                    #     print("%d  %s  %d  %-7.3f s %-7.3f  %-4.3f  %-6.3f  %.2f  %-6.2f  %-8d %-8d %d %d %d %s" % (chan, det.getChanCPD(dsNum,chan), isGood, thrKeV, thrSig, thrADC, thrADCErr, sigADC, sigADCErr, thrEvts, sigEvts, thrStatus, sigStatus, int(isBad), status))

                    # fill the dict vals
                    vals[chan] = [float("%.5f" % thrKeV), float("%.5f" % thrSig), int(isBad)]

                print("good detectors: %d/%d" % (nGood, nTot))

                # fill the DB
                if fillDB:
                    # for val in sorted(vals):
                        # print(val, vals[val])
                    dsi.setDBRecord({"key":key, "vals":vals}, forceUpdate, calDB=calDB, pars=pars)

                tf.Close()
                # return


if __name__=="__main__":
    main(sys.argv[1:])