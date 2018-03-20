#!/usr/bin/env python3
""" 'dsi.py': DataSetInfo for LAT.
    C. Wiseman, 18 March 2018
"""
import os, json, glob
import numpy as np

latSWDir    = os.environ['LATDIR']
dataDir     = "/global/projecta/projectdirs/majorana/users/wisecg"
bkgDir      = dataDir+"/bkg"     # subfolders: skim waves split lat
calDir      = dataDir+"/cal"     # subfolders: skim waves split lat
specialDir  = dataDir+"/special" # subfolders: skim waves split lat
skimDir     = bkgDir+"/skim"
waveDir     = bkgDir+"/waves"
splitDir    = bkgDir+"/split"
latDir      = bkgDir+"/lat"
calSkimDir  = calDir+"/skim"
calWaveDir  = calDir+"/waves"
calSplitDir = calDir+"/split"
calLatDir   = calDir+"/lat"
pandaDir    = dataDir+"/pandas"

class BkgInfo:
    def __init__(self):
        with open("%s/data/runsBkg.json" % latSWDir) as f:
            self.master = scrubDict(json.load(f))

    def dsRanges(self):
        """ Manual edit.  Must cover BG list and cal runs. """
        return {
            0:[2571,7614],
            1:[9407,14502],
            2:[14699,15892],
            3:[16797,18589],
            4:[60000791,60002394],
            5:[18623,23958],
            6:[25672,100000]
        }

    def dsMap(self):
        """{ds:lastRange}"""
        lastRange = {}
        for key in sorted(list(self.master)):
            last = int(list(self.master[key].keys())[-1])
            lastRange[int(key)] = last
        return lastRange

    def GetDSNum(self,run):
        ranges = self.dsRanges()
        for ids in range(len(ranges)):
            if ranges[ids][0] <= run <= ranges[ids][1]:
                dsNum = ids
        return dsNum

    def GetBkgIdx(self, dsNum, runNum):
        """ Finds the bkgIdx of a given run.  Must be IN the dataset! """
        bkgRuns = self.master[dsNum]

        for bkgIdx in bkgRuns:
            runCov = bkgRuns[bkgIdx]
            runList = []
            for idx in range(0,len(runCov),2):
                runList.extend(list(range(runCov[idx],runCov[idx+1]+1)))
            if runNum in runList:
                return bkgIdx
        return -1


class CalInfo:
    def __init__(self):
        with open("%s/data/runsCal.json" % latSWDir) as f:
            self.master = scrubNotes(json.load(f))
        with open("%s/data/runsSpecial.json" % latSWDir) as f:
            self.special = scrubNotes(json.load(f))

        # Track all the 'hi' run coverage numbers for fast run range lookups
        self.covIdx = {}
        for key in self.master:
            tmp = []
            for idx in self.master[key]:
                tmp.append(self.master[key][idx][2])
            self.covIdx[key] = np.asarray(tmp)

    def GetMasterList(self):
        return self.master

    def GetCovArr(self,key):
        return self.covIdx[key]

    def GetIdxs(self,key):
        return len(self.covIdx[key])

    def GetKeys(self,dsNum=None):
        keyList = sorted(self.master.keys())
        if dsNum==None:
            return keyList
        else:
            thisDSList = []
            for key in keyList:
                if "ds%d" % dsNum in key: thisDSList.append(key)
            return thisDSList

    def GetCalIdx(self,key,run):
        """ Look up the calibration index corresponding to a particular run. """
        if key not in self.covIdx:
            print("Key %s not found in master list!" % key)
            return None

        idx = np.searchsorted(self.covIdx[key], run)
        if idx not in self.master[key]:
            print("Run %d out of range of key %s.  calIdx was %d" % (run, key, idx))
            return None
        lst = self.master[key][idx]
        lo, hi = lst[1], lst[2]
        if lo <= run <= hi:
            return idx
        else:
            print("Run %d not found with key %s, lo=%d hi=%d" % (run,key,lo,hi))
            return None

    def GetNCalIdxs(self,dsNum,module):
        """ Get the number of calIdx's in a given dataset. """
        calKeys = self.GetKeys(dsNum)
        for key in calKeys:
            if "m%d" % module in key:
                return self.GetIdxs(key)
        return 0

    def GetCalList(self,key,idx,runLimit=None):
        """ Generate a list of runs for a given calibration index. """
        if key not in self.master:
            print("Key %s not found in master list!" % key)
            return None

        runList = []
        if idx not in self.master[key]:
            return None
        lst = self.master[key][idx][0]
        for i in range(0,len(lst),2):
            lo, hi = lst[i], lst[i+1]
            runList += range(lo, hi+1)
        if runLimit is not None:
            del runList[runLimit:]
        return runList

    def GetCalRunCoverage(self,key,idx):
        """ Return the (runLo, runHi) coverage of a particular calIdx"""
        if key not in self.master:
            print("Key %s not found in master list!" % key)
            return None
        return self.master[key][idx][1], self.master[key][idx][2]

    def GetCalFiles(self, dsNum, calIdx=None, modNum=None, verbose=False, cDir=calDir):
        """ Get a list of all files for a particular dsNum+calIdx.
            This will match the cut record entries in the DB.
        """
        calKeys = self.GetKeys(dsNum)

        fList = []
        for key in calKeys:
            if modNum is not None and str(modNum) not in key:
                continue
            if verbose: print(key)
            nIdx = self.GetIdxs(key)  # number of cal subsets

            # get the runs in each calIdx
            runList = []
            if calIdx!=None:
                runList = self.GetCalList(key, calIdx, 10)
                if verbose: print(runList)
            else:
                for idx in range(nIdx):
                    tmp = self.GetCalList(key, idx, 10)
                    if verbose: print(tmp)
                    runList += tmp

            # make a list of the actual file paths
            for run in runList:
                fPath = "%s/latSkimDS%d_run%d*.root" % (calDir, dsNum, run)
                fList += glob.glob(fPath)

        # for f in fList: print(f)
        return fList

    def GetSpecialKeys(self):
        return self.special.keys()

    def GetSpecialNIdxs(self,key):
        return len(self.special[key])

    def GetSpecialRuns(self,key,idx=None):

        noFiles = [6936,6937,6940,6942,6944,6965,6968,6969,6974,6977,7224,7267,7268,7269,7270,7271,7272,13168]

        if idx is not None:
            runLo, runHi = self.special[key][idx][0], self.special[key][idx][1]
            runList = [run for run in range(runLo, runHi+1) if run not in noFiles]
            return runList

        runList = []
        for idx in self.special[key].keys():
            runLo, runHi = self.special[key][idx][0], self.special[key][idx][1]
            runList.extend([run for run in range(runLo, runHi+1) if run not in noFiles])
        return runList

    def GetSpecialList(self):
        return self.special


class SimInfo:
    """ Adapted from ~mjdsim/analysisScriptsV2/analysisUtilities.py """
    dets = {}
    dets["M1"] = ['1010101', '1010102', '1010103', '1010104',
        '1010201', '1010202', '1010203', '1010204',
        '1010301', '1010302', '1010303', '1010304',
        '1010401', '1010402', '1010403', '1010404', '1010405',
        '1010501', '1010502', '1010503', '1010504',
        '1010601', '1010602', '1010603', '1010604',
        '1010701', '1010702', '1010703', '1010704']
    dets["M2"] = ['1020101', '1020102', '1020103', '1020104',
        '1020201', '1020202', '1020203', '1020204', '1020205',
        '1020301', '1020302', '1020303',
        '1020401', '1020402', '1020403', '1020404', '1020405',
        '1020501', '1020502', '1020503', '1020504',
        '1020601', '1020602', '1020603', '1020604',
        '1020701', '1020702', '1020703', '1020704']
    detectors = dets["M1"] + dets["M2"]

    activeDets = {}
    activeDets["M1"] = {
              # C1P1         C1P2         C1P3         C1P4            C1P5         C1P6         C1P7
        'All':[	1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1],
        'DS0':[	1, 1, 1, 1,  0, 1, 1, 0,  0, 0, 0, 1,  1, 1, 1, 1, 1,  1, 1, 1, 1,  0, 1, 1, 0,  1, 1, 1, 0],
        'DS1':[	0, 1, 1, 1,  1, 1, 1, 0,  0, 1, 1, 1,  0, 0, 0, 0, 0,  0, 0, 1, 0,  1, 0, 1, 1,  1, 1, 1, 1],
        'DS2':[	0, 1, 1, 1,  1, 1, 1, 0,  0, 1, 1, 1,  0, 0, 0, 0, 0,  0, 0, 1, 0,  1, 0, 1, 1,  1, 1, 0, 1],
        'DS3':[	0, 1, 1, 1,  1, 1, 1, 0,  0, 1, 1, 1,  1, 0, 0, 1, 1,  0, 1, 1, 0,  1, 0, 1, 1,  1, 1, 1, 1],
        'DS4':[	0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0],
        'DS5':[	0, 1, 1, 1,  1, 1, 1, 0,  0, 1, 1, 1,  1, 1, 1, 1, 1,  0, 1, 1, 0,  1, 0, 1, 1,  1, 1, 1, 1]
        }
    activeDets["M2"] = {
        'All':[	1, 1, 1, 1,  1, 1, 1, 1, 1,  1, 1, 1,  1, 1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1],
        'DS0':[	0, 0, 0, 0,  0, 0, 0, 0, 0,  0, 0, 0,  0, 0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0],
        'DS1':[	0, 0, 0, 0,  0, 0, 0, 0, 0,  0, 0, 0,  0, 0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0],
        'DS2':[	0, 0, 0, 0,  0, 0, 0, 0, 0,  0, 0, 0,  0, 0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0],
        'DS3':[	0, 0, 0, 0,  0, 0, 0, 0, 0,  0, 0, 0,  0, 0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0],
        'DS4':[	0, 0, 0, 1,  1, 1, 1, 0, 0,  1, 1, 0,  1, 0, 0, 1, 0,  1, 0, 1, 0,  1, 1, 0, 0,  0, 0, 1, 1],
        'DS5':[	1, 0, 0, 1,  1, 1, 1, 0, 0,  1, 1, 0,  1, 1, 0, 1, 0,  1, 0, 1, 1,  0, 1, 0, 0,  0, 0, 1, 1]
        }

    dtCutoffs = {}  #  P1           P2            P3            P4              P5           P6           P7
    dtCutoffs["M1"] = [7, 7, 6, 6,  8, 8, 7, 6,   6, 6, 8, 8,   6, 6, 6, 6, 7,  6, 8, 7, 8,  6, 6, 6, 6,  7, 7, 6, 9]
    dtCutoffs["M2"] = [6, 7, 5, 6,  6, 6, 6, 6, 6,   7, 8, 8,   6, 6, 6, 5, 6,  6, 7, 6, 6,  7, 9, 6, 6,  6, 7, 7, 5]

    def __init__(self, config):
        self.config = config

    def GetDetectorList(self, module=None):
        return (self.detectors if module is None else self.dets[module])

    def GetActiveDets(self, config, module):
        detList = []
        for iD, det in enumerate(self.GetDetectorList(module)):
            if self.activeDets[module][config][iD] == 1:
                detList.append(det)
        return detList

    def GetDTCutoff(self, module, detector):
        iD = self.dets[module].index(detector)
        return self.dtCutoffs[module][iD]


def scrubDict(myDict):
    """ Give dict imported by json files integer keys, and scrub out notes"""
    for key in list(myDict):
        if "note" in key:
            del myDict[key]
            continue
        for key2 in list(myDict[key]):
            if "note" in key2:
                del myDict[key][key2]
    makeIntKeys = {int(key):{int(key2):myDict[key][key2] for key2 in myDict[key]} for key in myDict}
    return makeIntKeys


def GetExposureDict(dsNum, modNum, dPath="%s/data" % latSWDir, verbose=False):

    chList = GetGoodChanList(dsNum)
    if dsNum==5 and modNum==1: chList = [ch for ch in chList if ch < 1000 and ch!=692]
    if dsNum==5 and modNum==2: chList = [ch for ch in chList if ch > 1000 and ch!=1232]

    expDict = {ch:[] for ch in chList}
    tmpDict, bkgIdx, prevBkgIdx = {}, -1, -1

    with open("%s/expos_ds%d.txt" % (dPath, dsNum), "r") as f:
        table = f.readlines()

    for idx, line in enumerate(table):
        tmp = (line.rstrip()).split(" ")
        if len(tmp)==0: continue

        if bkgIdx != prevBkgIdx:
            for ch in chList:
                if ch in tmpDict.keys():
                    expDict[ch].append(tmpDict[ch])
                else:
                    expDict[ch].append(0.)
            tmpDict = {}
            prevBkgIdx = bkgIdx

        if tmp[0] == "bkgIdx":
            bkgIdx = tmp[1]

        if len(tmp) > 1 and tmp[1] == ":" and tmp[0].isdigit() and int(tmp[0]) in chList:
            ch, exp = int(tmp[0]), float(tmp[2])
            tmpDict[ch] = exp

        if line == "All-channel summary: \n":
            summaryIdx = idx

    # get last bkgIdx
    for ch in chList:
        if ch in tmpDict.keys():
            expDict[ch].append(tmpDict[ch])
        else:
            expDict[ch].append(0.)

    # knock off the first element (it's 0).  Now expDict is done
    for ch in expDict:
        if expDict[ch][0] > 0:
            print("ERROR, WTF")
            exit(1)
        expDict[ch].pop(0)

    # now get the all-channel summary for HG channels
    summaryDict = {ch:[] for ch in chList}
    for line in table[summaryIdx+2:]:
        tmp = (line.rstrip()).split()
        ch, detID, aMass, runTime, expo = int(tmp[0]), int(tmp[1]), float(tmp[2]), float(tmp[3]), float(tmp[4])
        summaryDict[ch] = [expo, aMass]

    # now a final cross check
    if verbose:
        print("DS%d, M%d" % (dsNum, modNum))
    for ch in chList:

        if sum(expDict[ch]) > 0 and len(summaryDict[ch]) == 0:
            print("That ain't shoulda happened")
            exit(1)
        elif len(summaryDict[ch]) == 0:
            continue;

        mySum, ltResult, aMass = sum(expDict[ch]), summaryDict[ch][0], summaryDict[ch][1]
        diff = ((ltResult-mySum)/aMass) * 86400

        if verbose:
            print("%d   %.4f   %-8.4f    %-8.4f    %-8.4f" % (ch, aMass, mySum, ltResult, diff))

    return expDict


def test():
    print("testing...")

    # cal = CalInfo()
    # runsCal = cal.GetSpecialList()
    # print(runsCal)

    ds = BkgInfo()
    # print(ds.dsMap())
    # print(ds.dsRanges())
    # print(ds.GetDSNum(18588))
    print(ds.GetBkgIdx(6,27065))


if __name__=="__main__":
    test()