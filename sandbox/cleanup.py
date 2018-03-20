#!/usr/bin/env python3
import sys, os, glob
import numpy as np
import dsi

def main():

    # cleanUpLogs()
    # print("heyy")
    checkFiles()

def cleanUpLogs():

    skimFails = [
        "skim-ds5-116.txt","skim-ds5-121.txt","skim-ds5-115.txt","skim-ds5-120.txt","skim-ds5-113.txt","skim-ds5-114.txt",
        "skim-ds5-117.txt","skim-ds5-118.txt","skim-ds5-119.txt",
        ]

    # delete bad files?

    # make a new job sub list

    # ./skim_mjd_data 1 11 -l -g -d -t 0.7 /global/projecta/projectdirs/majorana/users/wisecg/bkg/skim >& ./logs/skim-ds1-11.txt &


    f = open("../jobLists/test.ls","w")

    for logFile in sorted(skimFails):
        tmp = logFile.split("-")
        dsNum = int(tmp[1][-1])
        subNum = int(tmp[2].split(".")[0])
        dub = "-d" if dsNum < 6 else ""

        # yep, no & at the end.
        # cmd = """./skim_mjd_data %d %d -l -g %s -t 0.7 /global/projecta/projectdirs/majorana/users/wisecg/bkg/skim >& ./logs/skim-ds%d-%d.txt\n""" % (dsNum, subNum, dub, dsNum, subNum)
        # f.write(cmd)

        if os.path.isfile("../logs/"+logFile):
            os.remove("../logs/"+logFile)
            print("deleted",logFile)

    f.close()


def checkFiles():
    from ROOT import TFile, TTree

    # non-vectors
    br1 = [
        "skimgatrev", "gatrev", "run", "iEvent", "startTime_s", "startTime0_s", "stopTime_s", "startClockTime_s",
        "clockTime_s", "localTime_s", "dtPulserGlobal", "sumEH", "sumEL", "sumEHClean", "sumELClean", "mH", "mL",
        "mHClean", "mLClean", "isLNFill1", "isLNFill2", "EventDC1Bits", "muType", "muTUnc", "muVeto",
        "mAct_M1Total_kg", "mAct_M1enr_kg", "mAct_M1nat_kg", "mAct_M2Total_kg", "mAct_M2enr_kg", "mAct_M2nat_kg",
        "mVeto_M1Total_kg", "mVeto_M2Total_kg"
        ]
    # vectors
    br2 = [
        "rawRun", "index", "iHit", "channel", "P", "D", "C", "gain", "mageID", "detID", "detName", "isEnr", "isNat",
        "isGood", "c0Channels", "globalTime", "tOffset", "triggerTrapt0", "dtPulserCard", "trapENMSample", "blrwfFMR50",
        "trapENFCal", "trapENMCal", "trapENFCalC", "trapECal", "onBoardE", "trapENF", "trapENM", "avse", "nlcblrwfSlope",
        "dcr99", "dcr95", "dcr90", "dcrctc90", "trapETailMin", "kvorrT", "dcr85", "dcr98", "dcr995", "dcr999", "RawWFblSlope",
        "RawWFblChi2", "wfDCBits", "d2wfnoiseTagNorm", "nX", "nFlippedBits", "d2wf5MHzTo30MHzPower",
        "d2wf30MHzTo35MHzPower", "d2wf0MHzTo50MHzPower", "threshKeV", "threshSigma", "dtmu_s", "mAct_g", "MGTWaveforms"
        ]

    dsNum = 5
    subLo, subHi = 113, 121

    for i in range(subLo,subHi+1):
        f = TFile("%s/waveSkimDS%d_%d.root" % (dsi.waveDir,dsNum,i))
        t = f.Get("skimTree")
        print(f,t.GetEntries())
        # for br in t.GetListOfBranches():
            # print(br)
            # print(br.GetClassName())

        # check ~1000 entries
        n = t.GetEntries()
        entList = np.arange(1,500,1) + np.arange(n-500,n-1,1)
        for iEnt in entList:
            t.GetEntry(iEnt)
            for b in br1:
                val = getattr(t,b)
                # print(b,val)
            for b in br2:
                try:
                    nCh = getattr(t,b).size()
                except AttributeError:
                    # print("error:",b)
                    pass
                for j in range(nCh):
                    val = getattr(t,b)[j]
                    # print(b,nCh,val)

        f.Close()
        # return


if __name__=="__main__":
    main()