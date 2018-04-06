#!/usr/bin/env python3
"""
===================== LAT/chan-sel.py ======================
Perform necessary run/channel selection bookkeeping for
    LAT analysis, before and after cuts are applied.

Uses:
- Check HV settings and thresholds for all bkg/cal runs
  --> Fill objects in dsi.py

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
import sys
import tinydb as db
import numpy as np

import dsi
bkg = dsi.BkgInfo()
det = dsi.DetInfo()
cal = dsi.CalInfo()

def main(argv):

    # get user options
    # for i,opt in enumerate(argv):
        # if opt == "-c": checkCal = True

    # getRunSettings()
    # dumpRunSettings()
    # getDetIDs()
    # checkSubRanges()
    loadSubRanges()


def getRunSettings():
    """ This takes ~2 hours to run on PDSF if projecta IO is good.
    Write a routine that identifies all active channels in a dataset,
    and outputs a plot of active channels vs. run number.
    Results are used to populate dsi::DetInfo.
    """
    from ROOT import GATDataSet, TFile, MJTChannelMap, MJTChannelSettings

    # fill these with [(run,val),...] pairs for each detector, each ds
    detHV, detTH, detCH, pMons = {}, {}, {}, {}

    for ds in [0,1,2,3,4,5,6]:
    # for ds in [2]:
        print("Scanning DS:%d" % ds)

        pMons[ds] = set()
        detHV[ds] = {d:[] for d in det.allDets} # high voltage setting
        detTH[ds] = {d:[] for d in det.allDets} # trap threshold setting
        detCH[ds] = {d:[] for d in det.allDets} # analysis channel (HG) setting

        gds = GATDataSet()
        runList = bkg.getRunList(ds)

        for idx, run in enumerate(runList):

            f = np.fabs(100*idx/len(runList) % 10)
            if f < 0.1:
                print("%d/%d, %.2f%% done." % (idx, len(runList), 100*idx/len(runList)))

            runPath = gds.GetPathToRun(run,GATDataSet.kGatified)
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
            pMons[ds].update(chPulser)

            for ch in chEnabled:

                detName = chMap.GetString(ch, "kDetectorName")
                detCPD = chMap.GetString(ch,"kStringName")+"D"+str(chMap.GetInt(ch,"kDetectorPosition"))
                detID = ''.join(i for i in detCPD if i.isdigit())

                # skip pulser monitors / special channels
                if detID == '0':
                    if ch not in pMons[ds]:
                        print("ch %d not in pulserMons list!  Adding it ..." % ch)
                        pMons[ds].add(ch)
                    continue

                gretCrate = chMap.GetInt(ch,"kVME")
                gretCard = chMap.GetInt(ch,"kCardSlot")
                gretHG = chMap.GetInt(ch,"kChanHi")
                gretLG = chMap.GetInt(ch,"kChanLo")
                threshHG = chSet.GetInt("TRAP Threshold",gretCrate,gretCard,gretHG,"ORGretina4MModel")
                threshLG = chSet.GetInt("TRAP Threshold",gretCrate,gretCard,gretLG,"ORGretina4MModel")
                # print(ch, detCPD, detID, gretCrate, gretCard, gretHG, gretLG)

                hvCrate = chMap.GetInt(ch,"kHVCrate")
                hvCard = chMap.GetInt(ch,"kHVCard")
                hvChan = chMap.GetInt(ch,"kHVChan")
                hvTarget = chMap.GetInt(ch,"kMaxVoltage")
                hvActual = chSet.GetInt("targets",hvCrate,hvCard,hvChan,"OREHS8260pModel")
                # print(ch, detCPD, detID, hvCrate, hvCard, hvChan, hvTarget, hvActual)

                # fill our data objects
                thisHV = (run, hvActual)
                thisTH = (run, threshHG) # NOTE: HG only

                if len(detHV[ds][detID]) == 0:
                    detHV[ds][detID].append(thisHV)
                if len(detTH[ds][detID]) == 0:
                    detTH[ds][detID].append(thisTH)

                # if we detect a difference, add a new (run,val) pair
                if len(detHV[ds][detID]) != 0:
                    prevHV = detHV[ds][detID][-1]
                    if thisHV[1] != prevHV[1]:
                        detHV[ds][detID].append(thisHV)

                if len(detTH[ds][detID]) != 0:
                    prevTH = detTH[ds][detID][-1]
                    if thisTH[1] != prevTH[1]:
                        detTH[ds][detID].append(thisTH)

                # skip LG channels
                if ch%2!=0: continue
                thisCH = (run, ch)

                if len(detCH[ds][detID]) == 0:
                    detCH[ds][detID].append(thisCH)

                if len(detCH[ds][detID]) != 0:
                    prevCH = detCH[ds][detID][-1]
                    if thisCH[1] != prevCH[1]:
                        detCH[ds][detID].append(thisCH)

            tf.Close()

    # save to npz file
    np.savez("./data/runSettings.npz",detHV,detTH,detCH,pMons)


def dumpRunSettings():

    f = np.load("./data/runSettings.npz")
    detHV = f['arr_0'].item()
    detTH = f['arr_1'].item()
    detCH = f['arr_2'].item()
    pMons = f['arr_3'].item()

    # print summary
    # print("HV settings:")
    # for ds in sorted(detHV):
    #     print("DS",ds)
    #     for det in sorted(detHV[ds]):
    #         print(det, detHV[ds][det])

    print("HG Threshold settings:")
    for ds in sorted(detTH):
        print("DS",ds)
        for det in sorted(detTH[ds]):
            print(det, detTH[ds][det])

    # print("HG Channel settings:")
    # for ds in sorted(detCH):
    #     print("DS",ds)
    #     for det in sorted(detCH[ds]):
    #         print(det, detCH[ds][det])
    #
    # print("Pulser monitors:")
    # print(pMons)


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


if __name__=="__main__":
    main(sys.argv[1:])