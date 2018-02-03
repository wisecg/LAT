#!/usr/bin/env python3
import sys, os, imp, glob
sys.argv.append("-b")
ds = imp.load_source('DataSetInfo',os.environ['LATDIR']+'/DataSetInfo.py')
wl = imp.load_source('waveLibs',os.environ['LATDIR']+'/waveLibs.py')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np
import tinydb as db


def main():
    riseTime()
    # avgSlo()
    # getEff()
    # getEffMultiChan()
    # compareRunTypes()


calInfo = ds.CalInfo()
extPulserInfo = {
        # Test 1 - attenuation, --donotbuild was used
        "7a": [[5942, 5945], [10,14,18,22], 200, 674], # att, rt, chan
        7:  [[5947, 5960], [10,14,18,22,26,30,34,38,42,46,50,50,54,58,62], 190, 674],
        8:  [[5964, 5977], [10,14,18,22,26,30,34,38,42,46,50,50,54,58,62], 190, 624],
        9:  [[5979, 5992], [10,14,18,22,26,30,34,38,42,46,50,50,54,58,62], 190, 688],
        10: [[6191, 6204], [10,14,18,22,26,30,34,38,42,46,50,50,54,58,62], 190, 662],
        11: [[6206, 6219], [10,14,18,22,26,30,34,38,42,46,50,50,54,58,62], 190, 608],
        # Test 2 - rise time
        12: [[6934, 6944], [140,145,150,155,160,165,170,175,180,185,190], 18, 674], # rt, att, chan
        13: [[6964, 6970], [4354,1257,1296,654,1278,1278,1278],0,[674,624,688,662,608,608,608]], # adc,att,chan
        14: [[6971, 6976], [140,145,150,155,160,165,170,175,180,185,190], 18, 614],
        15: [[6977, 6982], [140,145,150,155,160,165,170,175,180,185,190], 18, 624],
        16: [[7002, 7007], [140,145,150,155,160,165,170,175,180,185,190], 18, 688],
        17: [[7008, 7013], [140,145,150,155,160,165,170,175,180,185,190], 18, 662],
        # Test 3 - attenuation
        18: [[7219, 7233], [14,18,22,26,30,999,30,34,38,42,46,50,54,58,62], 155, 674], # att, rt, chan
        19: [[7234, 7246], [14,18,22,26,30,34,38,42,46,50,54,58,62], 164, 624],
        20: [[7247, 7259], [14,18,22,26,30,34,38,42,46,50,54,58,62], 146, 688],
        21: [[7260, 7272], [14,18,22,26,30,34,38,42,46,50,54,58,62], 138, 662],
        22: [[13168, 13181], [14,18,22,26,30,34,38,42,46,50,54,58,62], 999, 690]
    }
def skipRuns(runList):
    return [run for run in runList if run not in [6936,6937,6940,6942,6944,6974,6977,7224] and run < 7266]


def riseTime():
    """ fitSlo vs. rise time """
    from ROOT import TChain

    pIdx = 14
    extChan = extPulserInfo[pIdx][-1]
    syncChan = wl.getChan(0,10,0) # 672
    runList = skipRuns(calInfo.GetSpecialRuns("extPulser",pIdx))

    for run in runList:
        fileList = ds.getLATRunList([run],"%s/lat" % (ds.specialDir))
        latChain = TChain("skimTree")
        for f in fileList:
            print("Adding",f)
            latChain.Add("%s/lat/%s" % (ds.specialDir,f))

        tNames = ["Entry$","mH","channel","trapENFCal","fitSlo","den90","den10"]
        theCut = "(channel==%d || channel==%d) && mH==2" % (syncChan, extChan) # enforce correct sync
        # theCut = "Entry$ < 100 && trapENFCal > 10 && gain==0"
        tVals = wl.GetVX(latChain,tNames,theCut)
        nPass = len(tVals["Entry$"])
        if nPass == 0: continue

        # for i in range(nPass):
        #     ent = tVals["Entry$"][i]
        #     mH = tVals["mH"][i]
        #     chan = tVals["channel"][i]
        #     ene = tVals["trapENFCal"][i]
        #     slo = tVals["fitSlo"][i]
        #     d90 = tVals["den90"][i]
        #     d10 = tVals["den10"][i]
        #     print("%d  %d  %d  e%-10.2f  s%-8.2f  %-8.2f" % (ent,chan,mH,ene,slo,d90-d10))

        enfArr = [tVals["trapENFCal"][i] for i in range(nPass) if tVals["channel"][i]==extChan]
        sloArr = [tVals["fitSlo"][i] for i in range(nPass) if tVals["channel"][i]==extChan]
        rtArr =  [tVals["den90"][i] - tVals["den10"][i] for i in range(nPass) if tVals["channel"][i]==extChan]
        nSync = len(enfArr)

        for idx in range(len(enfArr)):
            print("%.2f  %.2f  %.2f" % (enfArr[idx],sloArr[idx],rtArr[idx]))

        # enfArr, sloArr = np.asarray(enfArr), np.asarray(sloArr)



def avgSlo():
    """ for each RT, get fitSlo's mu, sigma.  plot it for each channel """


def getEff():
    """ Efficiency vs. energy, one channel """


def getEffMultiChan():
    """ Efficiency vs. energy, comparing two different rise times.
        Requires that Test 1 be matched w/ the sync channel properly.
        Maybe try manually grouping the events within 4us of each other?
    """

def compareRunTypes():
    """ Plot fitSlo for extPulser, bkg, and cal runs for one channel. """


def wfFitter():
    """ fitSlo vs energy for different wf fit methods.
        Probably should save this for ext3.py
    """


if __name__=="__main__":
    main()
