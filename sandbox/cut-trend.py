#!/usr/bin/env python
"""
===================== cut-trend.py =====================
Draws trend of cuts from database

Usage:
    ./cut-trend.py

"""
import sys, os, glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import waveLibs as wl
import DataSetInfo as ds
sns.set_style('whitegrid')
sns.set_context('talk')

def main(argv):
    cInfo = ds.CalInfo()
    dsNum, modNum = 1, 1
    inDir, outDir = '.', '.'
    chNum, dfList = -1, []
    parList, chList = [], []
    perc = 90
    if len(argv) == 0:
        return
    for i, opt in enumerate(argv):
        # -- Cut tuning options --
        if opt == "-all":
            parList.append('bcMax')
            parList.append('noiseWeight')
            parList.append('bcTime')
            parList.append('pol2')
            parList.append('pol3')
            parList.append('fitSlo')
            parList.append('riseNoise')
            print ("Drawing all cuts")
        if opt == "-bcMax":
            parList.append('bcMax')
            print ("Drawing bcMax")
        if opt == "-noiseWeight":
            parList.append('noiseWeight')
            print ("Drawing noiseWeight")
        if opt == "-bcTime":
            parList.append('bcTime')
            print ("Drawing bcTime")
        if opt == "-tailSlope":
            parList.append('pol2')
            parList.append('pol3')
            print ("Drawing tailSlope")
        if opt == "-fitSlo":
            parList.append('fitSlo')
            print ("Drawing fitSlo")
        if opt == "-riseNoise":
            parList.append('riseNoise')
            print ("Drawing riseNoise")
        if opt == "-wfstd":
            parList.append('wfstd')
            print ("Drawing wfstd")

        # -- Input/output options --
        if opt == "-s":
            dsNum, modNum = int(argv[i+1]), int(argv[i+2])
            print ("Drawing DS-%d Module-%d"%(dsNum, modNum))
        if opt == "-d":
            inDir, outDir = argv[i+1], argv[i+2]
            print ("Custom paths: Input %s" % (inDir))
        if opt == "-ch":
            chNum = int(argv[i+1])
            print ("Drawing specific channel %d" % (chNum))
        if opt == "-p":
            perc = str(argv[i+1])
            print ("Drawing for percentage %s" % (perc))

    # -- Load channel list --
    if chNum == -1:
        chList = ds.GetGoodChanList(dsNum)
        if dsNum==5 and modNum == 1: # remove 692 and 1232
            chList = [584, 592, 598, 608, 610, 614, 624, 626, 628, 632, 640, 648, 658, 660, 662, 672, 678, 680, 688, 690, 694]
        if dsNum==5 and modNum == 2:
            chList = [1106, 1110, 1120, 1124, 1128, 1170, 1172, 1174, 1176, 1204, 1208, 1298, 1302, 1330, 1332]
    else:
        chList = [chNum]

    # Everything done by DB now!
    fig = plt.figure(figsize=(12,7))
    ax = fig.add_subplot(111)
    ch_idx = np.linspace(0, 1, len(chList))
    parValList = {}
    for ch in chList: parValList[ch] = []
    for par in parList:
        ax.cla()
        for subNum in cInfo.master["ds%d_m%d"%(dsNum,modNum)].keys():
            # fsD = wl.getDBCalRecord("%s_ds%d_idx%d_m%d_Peak"%(par,dsNum,subNum,modNum))
            fsD = wl.getDBCalRecord("%s_ds%d_idx%d_m%d"%(par,dsNum,subNum,modNum))
            for idx2, ch in enumerate(chList):
                # if parValList[idx2]:
                parValList[ch].append(fsD[ch][2]) # 90% value
                # else:
                    # parValList[idx2] = [fsD[ch][2]]

        print parValList
        for idx,ch in zip(ch_idx, chList):
            # Grab the array for the channel + cut combo
            # parVals0 = dfCut.loc[par, str(ch)].values
            # Fill any potential zeros
            # parVals = FillZeros(parVals0)
            # ax.plot(np.linspace(0, len(parVals), len(parVals)), parVals, marker='o', label='Ch %s'%(ch), color = plt.cm.tab20c(idx))
            ax.plot(np.linspace(0, len(parValList[ch]), len(parValList[ch])), parValList[ch], marker='o', label='Ch %s'%(ch), color = plt.cm.tab20c(idx))
        # Make stuff pretty
        ax.set_title('DS%d %s Trend (%s)'%(dsNum, par, perc))
        ax.set_xlabel('SubDS')
        ax.set_ylabel(par)
        ax.xaxis.set_ticks(np.arange(0, len(parValList[chList[0]])+1, 2) )
        plt.tight_layout()
        # Shrink current axis and put legend on the side
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width*0.88, box.height])
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), ncol=1)
        fig.savefig('%s/%s_ds%d_m%d_p%s.png'%(outDir, par, dsNum, modNum, perc))

def FillZeros(parVals):
    for idx,val in enumerate(parVals):
        # Fill with previous value
        if val == 0 and idx == 0:
            parVals[idx] = parVals[idx+1]
        elif val == 0:
            parVals[idx] = parVals[idx-1]
        else: continue
    return parVals

if __name__ == '__main__':
    main(sys.argv[1:])
