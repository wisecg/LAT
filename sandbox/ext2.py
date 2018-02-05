#!/usr/bin/env python3
import sys, os, imp, glob
sys.argv.append("-b")
ds = imp.load_source('DataSetInfo',os.environ['LATDIR']+'/DataSetInfo.py')
wl = imp.load_source('waveLibs',os.environ['LATDIR']+'/waveLibs.py')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np
import tinydb as db
calInfo = ds.CalInfo()

def main():
    riseTime()
    # avgSlo()
    # getEff()
    # getEffMultiChan()
    # compareRunTypes()




def riseTime():
    """ fitSlo vs. rise time """
    from ROOT import TChain

    extPulserInfo = calInfo.GetSpecialList()["extPulserInfo"]
    syncChan = wl.getChan(0,10,0) # 672

    pIdx = 13
    runList = calInfo.GetSpecialRuns("extPulser",pIdx)
    noFiles = [6936,6937,6940,6942,6944,6965,6968,6969,6974,6977,7224,7267,7268,13168]

    attList = extPulserInfo[pIdx][0]
    extChan = extPulserInfo[pIdx][-1]

    fig = plt.figure(figsize=(10,5),facecolor='w')
    cmap = plt.cm.get_cmap('hsv',len(runList)+1)
    xLo, xHi = 65, 80

    enfArr, sloArr, rtArr = [], [], []
    for i, run in enumerate(runList):
        if run in noFiles:
            continue

        fileList = ds.getLATRunList([run],"%s/lat" % (ds.specialDir))
        if len(fileList)==0:
            print("pIdx %d, run %d.  No files found." % (pIdx, run))
            continue

        latChain = TChain("skimTree")
        for f in fileList:
            latChain.Add("%s/lat/%s" % (ds.specialDir,f))

        tNames = ["Entry$","mH","channel","trapENFCal","fitSlo","den90","den10"]
        theCut = "(channel==%d || channel==%d) && mH==2" % (syncChan, extChan) # enforce correct sync
        tVals = wl.GetVX(latChain,tNames,theCut)
        nPass = len(tVals["Entry$"])
        if nPass == 0: continue

        enfTmp = [tVals["trapENFCal"][i] for i in range(nPass) if tVals["channel"][i]==extChan]
        sloTmp = [tVals["fitSlo"][i] for i in range(nPass) if tVals["channel"][i]==extChan]
        rtTmp  = [tVals["den90"][i] - tVals["den10"][i] for i in range(nPass) if tVals["channel"][i]==extChan]

        attRun = extPulserInfo[pIdx][0][i]
        attVal = attList[i]
        print("Run %d, atten %d" % (run,attVal))

        enfTmp, sloTmp, rtTmp = np.asarray(enfTmp), np.asarray(sloTmp), np.asarray(rtTmp)

        idx = np.where((enfTmp > xLo) & (enfTmp < xHi))

        plt.plot(sloTmp[idx], rtTmp[idx], ".", markersize=1, c=cmap(i), label="rt %d" % attVal)

        # enfArr.extend(enfTmp)
        # sloArr.extend(sloTmp)
        # rtArr.extend(rtTmp)

    # enfArr, sloArr, rtArr = np.asarray(enfArr), np.asarray(sloArr), np.asarray(rtArr)
    # print("enf",min(enfArr),max(enfArr),"slo",min(sloArr),max(sloArr),"rt",min(rtArr),max(rtArr))

    # xArr, yArr = enfArr, sloArr
    # xLo, xHi, bpX, yLo, yHi, bpY = 65, 80, 0.2, 50, 100, 0.1

    # xArr, yArr = enfArr, rtArr
    # xLo, xHi, bpX, yLo, yHi, bpY = 65, 80, 0.2, 320, 380, 0.1

    # xArr, yArr = sloArr, rtArr
    # xLo, xHi, bpX, yLo, yHi, bpY = 50, 100, 0.2, 335, 375, 0.2

    # nBY, nBX = int((yHi-yLo)/bpY), int((xHi-xLo)/bpY)
    # plt.hist2d(xArr, yArr, bins=[nBX,nBY], range=[[xLo,xHi],[yLo,yHi]], norm=LogNorm())
    # plt.colorbar()

    plt.xlim(50, 100)
    plt.title("pIdx %d, channel %d" % (pIdx, extChan))
    plt.xlabel("fitSlo",horizontalalignment='right',x=1.0)
    plt.ylabel("riseTime (ns)",horizontalalignment='right',y=1.0)
    plt.legend(loc="best")
    plt.savefig("../plots/rtStudy_idx%d.pdf" % (pIdx))
    # plt.show()


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
