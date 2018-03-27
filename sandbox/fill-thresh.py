#!/usr/bin/env python
"""
    This script fills in DB with threshold values
    The key is: thresh_ds#_bkgidx#
    The main code outputs a list of dictionaries, one per subDS
    In the dictionary, the values are {ch:[threshold, sigma]}

    Notes:
        - Sigma can be negative (it's a Gaussian!), use abs(Sigma) if a positive value is necessary
        - For DS1: ch 672 is off for bkgidx 32, 33, 43, 44, and 45
        - For DS5: bkgidx 12 is missing all of M2, M2 is probably off for that subDS

"""

import os, imp, ROOT
ds = imp.load_source('DataSetInfo',os.environ['LATDIR']+'/sandbox/DataSetInfo.py')
wl = imp.load_source('waveLibs',os.environ['LATDIR']+'/waveLibs.py')

def main():
    dsNum = 5
    inDir = '/mnt/mjdDisk1/Majorana/data/sandbox/thresholds'
    threshDictList = generateThreshDict(dsNum, inDir)
    print(len(threshDictList))
    print(sorted(threshDictList[0].keys()))


def findGood(threshDictList, dsNum, subDS, ch, searchStyle = 0):
    """
        This search function returns a valid threshold idx
        Search style: 0 = alternate between forward one idx and backward one idx
    """
    searchList, goodIdx = [], -1
    searchMax = min(5, ds.dsMap[dsNum]-subDS)
    searchMin = min(5, subDS)
    if searchStyle == 0:
        # Build alternating idx list
        searchList = [None]*(searchMax + searchMin)
        searchList[::2] = [subDS+i for i in range(searchMax)]
        searchList[1::2] = [subDS-i for i in range(searchMin)]
    elif searchStyle == -1:
        # Backward only
        searchList = [subDS-i for i in range(searchMin)]
    elif searchStyle == 1:
        # Forward only
        searchList = [subDS+i for i in range(searchMax)]
    else:
        print('Search style not valid!')
        return

    for idx in searchList:
        try:
            values = threshDictList[idx][ch]
            if 99999.0 in values or 999999.0 in values:
                continue
            elif values[0] < 0:
                continue
            else:
                goodIdx = subDS+idx
        # If ch doesn't exist in the dictionary keys, continue to next
        except:
            continue
    return goodIdx


def generateThreshDict(dsNum = 1, inDir = '', fillDb = False):

    threshTree = ROOT.TChain("threshTree")
    # Note: for DS0, the original threshold generation had 76 subDS instead of 75, so include 1 more for DS0!
    for idx in range(ds.dsMap[dsNum]+1):
        threshTree.Add("{}/threshDS{}_{}.root".format(inDir, dsNum, idx))

    chList = ds.GetGoodChanList(dsNum)
    if dsNum == 5: # remove 692 and 1232 (both beges, so who cares)
        if 692 in chList: chList.remove(692)
        if 1232 in chList: chList.remove(1232)

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
            print ("Uh oh, we have a problem with the background idx:", iList)
        key = "thresh_ds{:d}_bkgidx{:d}".format(dsNum, max(bIDX1, bIDX2))
        nChans = threshTree.channelList.size()
        for iH in range(nChans):
            chan = threshTree.channelList.at(iH)
            # Keep only high gain
            if chan%2 != 0 or chan not in chList: continue
            threshDict[chan] = [threshTree.threshCal.at(iH), abs(threshTree.sigmaCal.at(iH))]

        # Print message if a channel is missing in the keys -- this means the channel was off
        # The channel will NOT be included in the dictionary (we can change later)
        for ch in chList:
            if ch not in threshDict.keys():
                print ("Channel {} not in {}".format(ch, key))

        # Add shit to lists
        keyList.append(key)
        dictList.append(threshDict)


    # This way is slow but it's careful and we don't make any mistakes
    # Solve issues with missing parameters by looping through lists
    for idx, tDict in enumerate(dictList):
        # Check if a value is 99999.0 or 999999.0 (typically bad calibration)
        for ch, values in tDict.items():
            # If the index is 0, can't fill with the previous index so must fill with the next
            if idx == 0 and (99999.0 in values or 999999.0 in values or values[0] < 0):
                goodIdx = findGood(dictList, dsNum, idx, ch, 1)
                print "Filling idx{:d} ch{:d} with".format(idx, ch), dictList[goodIdx][ch]
                dictList[idx][ch] = dictList[goodIdx][ch]
                # print "Filling idx{:d} ch{:d} with".format(idx, ch), dictList[idx+1][ch]
                # dictList[idx][ch] = dictList[idx+1][ch]
            elif idx < len(dictList)-1 and idx!= 0 and (99999.0 in values or 999999.0 in values or values[0] < 0):
                goodIdx = findGood(dictList, dsNum, idx, ch, 0)
                print "Filling idx{:d} ch{:d} with".format(idx, ch), dictList[goodIdx][ch]
                dictList[idx][ch] = dictList[goodIdx][ch]

                # if ch in dictList[idx-1].keys() and (99999.0 not in dictList[idx-1][ch] or 999999.0 not in dictList[idx-1][ch]):
                #     print "Filling idx{:d} ch{:d} with".format(idx, ch), dictList[idx-1][ch]
                #     dictList[idx][ch] = dictList[idx-1][ch]
                # elif ch in dictList[idx+1].keys() and (99999.0 not in dictList[idx+1][ch] or 999999.0 not in dictList[idx+1][ch]):
                #     print "Filling idx{:d} ch{:d} with".format(idx, ch), dictList[idx+1][ch]
                #     dictList[idx][ch] = dictList[idx+1][ch]
            # If index is the last one, must fill with previous
            elif idx == len(dictList)-1 and (99999.0 in values or 999999.0 in values or values[0] < 0):
                goodIdx = findGood(dictList, dsNum, idx, ch, -1)
                print "Filling idx{:d} ch{:d} with".format(idx, ch), dictList[goodIdx][ch]
                dictList[idx][ch] = dictList[goodIdx][ch]

                # print "Filling idx{:d} ch{:d} with".format(idx, ch), dictList[idx-1][ch]
                # dictList[idx][ch] = dictList[idx-1][ch]


    # Loop again to check for negative threshold values
    # print "Second Loop for Negative values"
    # for idx, tDict in enumerate(dictList):
    #     for ch, values in tDict.items():
    #         if idx == 0 and values[0] < 0:
    #             print "Filling Negative idx{:d} ch{:d} with".format(idx, ch), dictList[idx+1][ch]
    #             dictList[idx][ch] = dictList[idx+1][ch]
    #         elif idx < len(dictList)-1 and idx!=0 and values[0] < 0:
    #             if ch in dictList[idx-1].keys() and dictList[idx-1][ch][0] > 0:
    #                 print "Filling Negative idx{:d} ch{:d} with".format(idx, ch), dictList[idx-1][ch]
    #                 dictList[idx][ch] = dictList[idx-1][ch]
    #             elif ch in dictList[idx+1].keys() and dictList[idx+1][ch][0] > 0:
    #                 print "Filling Negative idx{:d} ch{:d} with".format(idx, ch), dictList[idx+1][ch]
    #                 dictList[idx][ch] = dictList[idx+1][ch]
    #         elif idx == len(dictList)-1 and values[0] < 0:
    #             print "Filling Negative idx{:d} ch{:d} with".format(idx, ch), dictList[idx-1][ch]
    #             dictList[idx][ch] = dictList[idx-1][ch]

    # For DS5 -- channel 1302 had bad values for a series of subDS early on
    # Manually filling in here with a good value
    # if dsNum == 5:
    #     print("DS5, manually filling in some values because of bad runs")
    #     dictList[13][1302] = dictList[16][1302]
    #     dictList[14][1302] = dictList[16][1302]
    #     dictList[15][1302] = dictList[16][1302]

    # Final check
    print "Third Loop for final checks"
    for idx, tDict in enumerate(dictList):
        for ch, values in tDict.items():
            if 99999.0 in values or 999999.0 in values:
                print "Warning, bad fit found:", idx, ch
                print dictList[idx][ch]
            if values[0] < 0:
                print "Warning, negative value found:", idx, ch, values

        # Finally fill DB here if all values look good
        if fillDb:
            wl.setDBRecord({"key":keyList[idx],"vals":dictList[idx]}, forceUpdate=False)
            wl.setDBRecord({"key":keyList[idx],"vals":dictList[idx]}, forceUpdate=True)
    return dictList


if __name__ == "__main__":
    main()
