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

import dsi
bkg = dsi.BkgInfo()
det = dsi.DetInfo()

def main(argv):

    # get user options
    # for i,opt in enumerate(argv):
        # if opt == "-c": checkCal = True

    getRunSettings()


def getRunSettings():
    """
    Write a routine that identifies all active channels in a dataset,
    and outputs a plot of active channels vs. run number.
    Results are used to populate dsi::DetInfo.
    """
    from ROOT import GATDataSet, TFile, MJTChannelMap, MJTChannelSettings

    # for ds in [0,1,2,3,4,5,"5A","5B","5C"]:
    ds = 1
    runList = bkg.getRunList(ds)

    # fill these with [(run,val),...] pairs for each detector
    detHV = {d:[] for d in det.allDets}
    detTH = {d:[] for d in det.allDets}

    gds = GATDataSet()
    for run in runList[:100]:
        print(run)

        runPath = gds.GetPathToRun(run,GATDataSet.kGatified)
        tf = TFile(runPath)
        chSet = tf.Get("ChannelSettings")
        chMap = tf.Get("ChannelMap")
        # chMap.DumpChannelMap()
        # chSet.DumpSettings() # dump to a file to see syntax for HV and TRAP

        # access more stuff than we really need, just to provide some examples

        chEnabled = chSet.GetEnabledIDList()
        chSpecial = chMap.GetSpecialChanMapString()
        chPulser = chMap.GetPulserChanList()
        chPulser = [chPulser[i] for i in range(len(chPulser))]
        if ds == 1: chPulser.extend([674, 675, 677]) # 674, 675, 677 are not in the MJTChannelMap's due to a bug

        for ch in chEnabled:

            detName = chMap.GetString(ch, "kDetectorName")
            detCPD = chMap.GetString(ch,"kStringName")+"D"+str(chMap.GetInt(ch,"kDetectorPosition"))
            detID = ''.join(i for i in detCPD if i.isdigit())

            # skip pulser monitors / special channels
            if detID == '0':
                if ch not in chPulser:
                    print("ch %d not in chPulser list!  Adding it ..." % ch)
                    chPulser.append(ch)
                continue

            gretCrate = chMap.GetInt(ch,"kVME")
            gretCard = chMap.GetInt(ch,"kCardSlot")
            gretHG = chMap.GetInt(ch,"kChanHi")
            gretLG = chMap.GetInt(ch,"kChanLo")
            # print(ch, detCPD, detID, gretCrate, gretCard, gretHG, gretLG)

            hvCrate = chMap.GetInt(ch,"kHVCrate")
            hvCard = chMap.GetInt(ch,"kHVCard")
            hvChan = chMap.GetInt(ch,"kHVChan")
            hvTarget = chMap.GetInt(ch,"kMaxVoltage")
            hvActual = chSet.GetInt("targets",hvCrate,hvCard,hvChan,"OREHS8260pModel")
            # print(ch, detCPD, detID, hvCrate, hvCard, hvChan, hvTarget, hvActual)

            threshHG = chSet.GetInt("TRAP Threshold",gretCrate,gretCard,gretHG,"ORGretina4MModel")
            threshLG = chSet.GetInt("TRAP Threshold",gretCrate,gretCard,gretLG,"ORGretina4MModel")
            # print(ch, detCPD, detID, threshHG, threshLG)

            # fill our data objects
            thisHV = (run, hvActual)
            thisTH = (run, threshHG)

            if len(detHV[detID]) == 0:
                detHV[detID].append(thisHV)
            if len(detTH[detID]) == 0:
                detTH[detID].append(thisTH)

            # if we detect a difference, add a new (run,val) pair
            if len(detHV[detID])!=0:
                prevHV = detHV[detID][-1]
                if thisHV[1] != prevHV[1]:
                    detHV[detID].append(thisHV)

            if len(detTH[detID])!=0:
                prevTH = detTH[detID][-1]
                if thisTH[1] != prevTH[1]:
                    detTH[detID].append(thisTH)

        tf.Close()

    # print summary
    print("HV settings:")
    for key in sorted(detHV):
        print(key, detHV[key])

    print("Threshold settings:")
    for key in sorted(detTH):
        print(key, detTH[key])

if __name__=="__main__":
    main(sys.argv[1:])