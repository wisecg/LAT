#!/usr/bin/env python3
import sys, os, imp, glob
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

sys.argv.append("-b")
ds = imp.load_source('DataSetInfo',os.environ['LATDIR']+'/DataSetInfo.py')
wl = imp.load_source('waveLibs',os.environ['LATDIR']+'/waveLibs.py')

def main():

    getSumSpec()


def getSumE(tree, theCut):
    """ get sum energy for all events passing cuts """
    sumArr = []
    tNames = ["Entry$","mH","channel","trapENFCal","sumEH","gain","isGood"]
    tVals = wl.GetVX(tree, tNames, theCut, False)
    nPass = len(tVals["Entry$"])
    prevEnt, sumE = -1, 0.
    for idx in range(nPass):
        ent = tVals["Entry$"][idx]
        mH = tVals["mH"][idx]
        chan = tVals["channel"][idx]
        tmp = tVals["sumEH"][idx]
        enf = tVals["trapENFCal"][idx]
        if enf > 99999: enf = 0
        if ent!=prevEnt and idx!=0:
            # print("%.2f" % sumE)
            sumArr.append(sumE)
            sumE = 0
            sumE += enf
        prevEnt = ent
        # print("%d  %d  %d  %.2f  %2f" % (ent,chan,mH,enf,sumE))
    # print("%.2f" % sumE)
    sumArr.append(sumE)
    sumArr = np.asarray(sumArr)
    return sumArr


def getSumSpec():
    from ROOT import TFile, TTree

    fileList = []
    calInfo = ds.CalInfo()
    runList = calInfo.GetSpecialRuns("longCal",5) # 54 total
    runList = runList[:25] # truncate
    for run in runList: fileList.extend(ds.getLATRunList([run],"%s/lat" % (ds.specialDir)))
    nFiles = len(fileList)

    xLo, xHi, xpb = 0., 4000., 2.
    nb = int((xHi-xLo)/xpb)
    sum1, sum2, sum3, sum4, sum5 = np.zeros(nb), np.zeros(nb), np.zeros(nb), np.zeros(nb), np.zeros(nb)

    for iFile, f in enumerate(fileList):
        print("%d/%d %s" % (iFile,nFiles,f))
        tf = TFile("%s/lat/%s" % (ds.specialDir,f))
        latTree = tf.Get("skimTree")

        sumArr = getSumE(latTree, "mH==4 && gain==0 && isGood")
        y, x = np.histogram(sumArr, bins=nbx, range=(xLo,xHi))
        sum4 = np.add(sum4, y)

        tf.Close()


if __name__=="__main__":
    main()
