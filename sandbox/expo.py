#!/usr/bin/env python3
import numpy as np
import tinydb as db
import dsi
bkg = dsi.BkgInfo()
cal = dsi.CalInfo()
det = dsi.DetInfo()

def main():

    checkDBCuts()
    # getExposure()


def checkDBCuts():
    """ Verify that dsi.py::GetDBCuts removes the det/bIdx's it was supposed to.

    NOTE/TODO: the riseNoise cut (line896 of dsi.py) was checking for idx 3, not idx 4.
    We probably need to regenerate the fr files (and probly just delete the rn files).

    ALSO the zombie files need to be removed!

    BUT, lat2::badRiseCut DID update the riseNoise DB values correctly, and the fitSlo values are fine.
    So we don't need to edit the DB, just probably GetDBCuts.
    """
    ds = 1
    cutType = "fr"

    # spot check that we are removing detectors/bIdx's from bad fitSlo/riseNoise tuning.
    # cutDets = ['211','214','221','261','274'] lat2.py::setSloCut
    cutChans = [692] # lat2.py::badRiseChans

    dsNum = int(ds[0]) if isinstance(ds, str) else int(ds)
    calDB = db.TinyDB('%s/calDB-v2.json' % (dsi.latSWDir))
    pars = db.Query()
    dsMap = bkg.dsMap() # number of bIdx's
    bkgRanges = bkg.getRanges(ds)
    calKeys = cal.GetKeys(dsNum)

    # have to treat modules separately
    mods = [1]
    if dsNum == 4: mods = [2]
    if dsNum == 5: mods = [1,2]
    for mod in mods:

        chList = det.getGoodChanList(dsNum, mod)

        for bIdx in bkgRanges:

            # get the cut values
            bkgDict, calDict, _,_ = dsi.GetDBCuts(ds,bIdx,mod,cutType,calDB,pars,True)

            for ch in chList:
                cpd = det.getChanCPD(dsNum,ch)

                thData = True if ch in bkgDict.keys() else False
                fsData = True if ch in calDict.keys() and "fitSlo" in calDict[ch] else False
                rnData = True if ch in calDict.keys() and "riseNoise" in calDict[ch] else False

                goodCh = True if thData and fsData and rnData else False

                # if cpd in cutDets:
                # if ch in cutChans:
                # if not rnData:
                print("DS%d  bIdx %d  ch %d  cpd %s  th %d  fs %d  rn %d  good? %d" % (dsNum,bIdx,ch,cpd,int(thData),int(fsData),int(rnData), int(goodCh)))



def getExposure():
    """ Similar to LAT2:ApplyCuts . """
    from ROOT import TFile, TTree

    ds = 0
    cutType = "fr"


    dsNum = int(ds[0]) if isinstance(ds, str) else int(ds)
    calDB = db.TinyDB('%s/calDB-v2.json' % (dsi.latSWDir))
    pars = db.Query()
    dsMap = bkg.dsMap() # number of bIdx's
    bkgRanges = bkg.getRanges(ds)
    calKeys = cal.GetKeys(dsNum)

    # load ds_livetime output
    tl = TFile("../data/ds_%s_livetime.root" % str(ds))
    lt = tl.Get("dsTree")

    # have to treat modules separately
    mods = [1]
    if dsNum == 4: mods = [2]
    if dsNum == 5: mods = [1,2]
    for mod in mods:

        chList = det.getGoodChanList(dsNum, mod)

        for bIdx in bkgRanges:

            # get the livetime and exposure of each channel in this bkgIdx
            # live = {ch:0 for ch in chList} # days
            # expo = {ch:0 for ch in chList} # kg-days
            #
            # n = lt.Draw("run:channel:livetime","run>=%d && run<=%d" % (bkgRanges[bIdx][0], bkgRanges[bIdx][1]), 'goff')
            # ltRun, ltChan, ltLive = lt.GetV1(), lt.GetV2(), lt.GetV3()
            # for i in range(n):
            #     ch = ltChan[i]
            #     if ch not in chList: continue # this skips M2 channels when we are looking at M1 and vice versa
            #     cpd = det.getChanCPD(dsNum,ch)
            #     detID = det.getDetIDChan(dsNum,ch)
            #     aMass = det.allActiveMasses[detID]
            #     live[ch] += ltLive[i]/86400
            #     expo[ch] += ltLive[i]*aMass/86400/1000

            # for ch in chList:
            #     cpd = det.getChanCPD(dsNum,ch)
            #     print("%s  %d  %s  %d  lt %.2f e %.2f" % (str(ds), bIdx, cpd, ch, live[ch], expo[ch]))
            # continue

            # get the cut values
            bkgDict, calDict, _,_ = dsi.GetDBCuts(ds,bIdx,mod,cutType,calDB,pars,False)
            for ch in chList:
                cpd = det.getChanCPD(dsNum,ch)

                # To make a "fr" cut file, LAT2 requires ALL THREE [thData, fsData, rnData] to be True.
                # If any bkg or cal sub-index was missing, the whole bIdx (for this channel) gets thrown out.
                # This was easier than being really nitpicky about the run boundaries for a minimal exposure increase.

                thData = True if ch in bkgDict.keys() else False
                fsData = True if ch in calDict.keys() and "fitSlo" in calDict[ch] else False
                rnData = True if ch in calDict.keys() and "riseNoise" in calDict[ch] else False

                goodCh = True if thData and fsData and rnData else False

                # if cpd in cutDets:
                # if ch in cutChans:
                print("DS%d  bIdx %d  ch %d  cpd %s  th %d  fs %d  rn %d  good? %d" % (dsNum,bIdx,ch,cpd,int(thData),int(fsData),int(rnData), int(goodCh)))





if __name__=="__main__":
    main()