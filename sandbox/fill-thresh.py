#!/usr/bin/env python
"""
    This script fills in DB with threshold values
    The key is: thresh_ds#_bkgidx#
    In the dictionary, the values are ch:[threshold, sigma]

    Notes:
        - Sigma can be negative (it's a Gaussian!), use abs(Sigma) if a positive value is necessary
        - For DS1: ch 672 is off for bkgidx 32, 33, 43, 44, and 45
        - For DS5: bkgidx 12 is missing all of M2, M2 is probably off for that subDS

"""

import ROOT
import DataSetInfo as ds
import waveLibs as wl

if __name__ == "__main__":
    dsNum = 1
    threshTree = ROOT.TChain("threshTree")
    # Note: for DS0, the original threshold generation had 76 subDS instead of 75, so include 1 more for DS0!
    for idx in range(ds.dsMap[dsNum]+1):
        threshTree.Add("/projecta/projectdirs/majorana/users/bxyzhu/thresholds/threshDS%d_%d.root"%(dsNum, idx))

    chList = ds.GetGoodChanList(dsNum)
    if dsNum==5:
        chList = [584, 592, 598, 608, 610, 614, 624, 626, 628, 632, 640, 648, 658, 660, 662, 672, 678, 680, 688, 690, 694, 1106, 1110, 1120, 1124, 1128, 1170, 1172, 1174, 1176, 1204, 1208, 1298, 1302, 1330, 1332]

    # Loop over threshold files
    keyList = []
    dictList = []
    prevDict = {}
    iList = -1
    while True:
        iList += 1
        if iList >= threshTree.GetEntries(): break
        threshDict = {}
        threshTree.GetEntry(iList)
        # Check if min or max run ranges have an idx associated with it
        # These thresholds were calculated a long time ago so minor run selection was done
        bIDX1 = ds.GetBkgIdx(dsNum, threshTree.runMin)
        bIDX2 = ds.GetBkgIdx(dsNum, threshTree.runMax)
        # It's not a problem unless both values are -1 -- then some super serious run selection was done
        if bIDX1 == -1 and bIDX2 == -1:
            print "Uh oh, we have a problem with idx", iList
        key = "thresh_ds%d_bkgidx%d"%(dsNum, max(bIDX1, bIDX2))
        nChans = threshTree.channelList.size()
        for iH in range(nChans):
            chan = threshTree.channelList.at(iH)
            # Keep only high gain
            if chan%2 != 0 or chan not in chList: continue
            threshDict[chan] = [threshTree.threshCal.at(iH), threshTree.sigmaCal.at(iH)]

        # Print message if a channel is missing in the keys -- this means the channel was off
        # The channel will NOT be included in the dictionary (we can change later)
        for ch in chList:
            if ch not in threshDict.keys():
                print "Channel %d not in %s"%(ch, key)

        # Add shit to lists
        keyList.append(key)
        dictList.append(threshDict)


    # This way is super ghetto but it's careful and we don't make any mistakes
    # Solve issues with missing parameters by looping through lists
    for idx, tDict in enumerate(dictList):
        # Check if a value is 99999.0 or 999999.0 (low statistics or failed fit)
        for ch, values in tDict.iteritems():
            if idx == 0 and (99999.0 in values or 999999.0 in values):
                print "Filling idx%d ch%d with"%(idx, ch), dictList[idx+1][ch]
                dictList[idx][ch] = dictList[idx+1][ch]
            elif idx < len(dictList)-1 and idx!= 0 and (99999.0 in values or 999999.0 in values):
                if ch in dictList[idx-1].keys() and (99999.0 not in dictList[idx-1][ch] or 999999.0 not in dictList[idx-1][ch]):
                    print "Filling idx%d ch%d with"%(idx, ch), dictList[idx-1][ch]
                    dictList[idx][ch] = dictList[idx-1][ch]
                elif ch in dictList[idx+1].keys() and (99999.0 not in dictList[idx+1][ch] or 999999.0 not in dictList[idx+1][ch]):
                    print "Filling idx%d ch%d with"%(idx, ch), dictList[idx+1][ch]
                    dictList[idx][ch] = dictList[idx+1][ch]
            elif idx == len(dictList)-1 and (99999.0 in values or 999999.0 in values):
                print "Filling idx%d ch%d with"%(idx, ch), dictList[idx-1][ch]
                dictList[idx][ch] = dictList[idx-1][ch]


    # Loop again to check for negative threshold values
    print "Second Loop for Negative"
    for idx, tDict in enumerate(dictList):
        for ch, values in tDict.iteritems():
            if idx == 0 and values[0] < 0:
                print "Filling Negative idx%d ch%d with"%(idx, ch), dictList[idx+1][ch]
                dictList[idx][ch] = dictList[idx+1][ch]
            elif idx < len(dictList)-1 and idx!=0 and values[0] < 0:
                if ch in dictList[idx-1].keys() and dictList[idx-1][ch][0] > 0:
                    print "Filling Negative idx%d ch%d with"%(idx, ch), dictList[idx-1][ch]
                    dictList[idx][ch] = dictList[idx-1][ch]
                elif ch in dictList[idx+1].keys() and dictList[idx+1][ch][0] > 0:
                    print "Filling Negative idx%d ch%d with"%(idx, ch), dictList[idx+1][ch]
                    dictList[idx][ch] = dictList[idx+1][ch]
            elif idx == len(dictList)-1 and values[0] < 0:
                print "Filling Negative idx%d ch%d with"%(idx, ch), dictList[idx-1][ch]
                dictList[idx][ch] = dictList[idx-1][ch]

    # For DS5 -- channel 1302 had bad values for a series of runs
    # Manually filling in here with a good value
    # dictList[13][1302] = dictList[16][1302]
    # dictList[14][1302] = dictList[16][1302]
    # dictList[15][1302] = dictList[16][1302]

    # Final check
    print "Third Loop for final checks"
    for idx, tDict in enumerate(dictList):
        for ch, values in tDict.iteritems():
            if 99999.0 in values or 999999.0 in values:
                print idx, ch
                print dictList[idx][ch]
            if values[0] < 0:
                print "NEGATIVE!", idx, ch, values

        # Finally fill DB here if all values look good
        # wl.setDBCalRecord({"key":keyList[idx],"vals":dictList[idx]}, forceUpdate=False)
        # wl.setDBCalRecord({"key":keyList[idx],"vals":dictList[idx]}, forceUpdate=True)
