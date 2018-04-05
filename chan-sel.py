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

def main(argv):

    # get user options
    # for i,opt in enumerate(argv):
        # if opt == "-c": checkCal = True

    # getRunSettings()
    # dumpRunSettings()
    getDetIDs()

    # checkRunRanges()


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
    mapping has never changed, load a channel map and print it out. """
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

        # transcribe this into DetInfo::detID
        print("'%s':%s," % (d, detID))


def checkRunRanges():
    """ Detect periods where:
        1) HV changed - identify bkgIdx and calIdx
        2) TH changed - identify bkgIdx
    Use the result to make some decisions about
    data periods to reject.
    """
    print("hi")



if __name__=="__main__":
    main(sys.argv[1:])