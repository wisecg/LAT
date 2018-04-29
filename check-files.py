#!/usr/bin/env python
import sys, os
sys.argv.append("-b")
import tinydb as db
import numpy as np

import waveLibs as wl
import dsi
bkg = dsi.BkgInfo()
cal = dsi.CalInfo()
det = dsi.DetInfo()

from ROOT import TFile, TTree, MGTWaveform

# switches
fLimit = None    # set to None to run over everything
skipDS6Cal = True
verbose = True
testMode = False

def main(argv):
    """ NOTE: don't use globs when getting files.
    Manually make sure everything is here.
    Can submit these commands as separate batch jobs:
        ./check-files.py -all
        ./check-files.py -c -all
    """
    global checkCal
    checkCal = False
    if checkCal: print("Skip DS6 cal?",skipDS6Cal)
    if testMode: print("Test mode active.")

    for i,opt in enumerate(argv):
        if opt == "-c": checkCal = True
        if opt == "-s": checkSkim()
        if opt == "-w": checkWave()
        if opt == "-p": checkSplit()
        if opt == "-l": checkLAT()
        if opt == "-all":
            checkSkim()
            checkWave()
            checkSplit()
            checkLAT()
        if opt == "-m": makeJobList()
        if opt == "-t": testDraw()


def checkSkim():
    """ ./check-files.py [-c] -s """
    print("Checking skims.  Cal?", checkCal)

    # make file list
    fileList = []

    # bkg
    if not checkCal:
        dsMap = bkg.dsMap()
        for ds in dsMap:
            for sub in range(dsMap[ds]+1):
                fname = "%s/skimDS%d_%d_low.root" % (dsi.skimDir, ds, sub)
                if os.path.isfile(fname):
                    fileList.append(fname)
                else:
                    print("File not found:",fname)
                    continue
    # cal
    else:
        for key in cal.GetKeys():
            ds = int(key[2])
            if skipDS6Cal and ds==6: continue
            for cIdx in range(cal.GetIdxs(key)): # 0-indexed
                cRuns = cal.GetCalList(key, cIdx)
                for run in cRuns:
                    fname = "%s/skimDS%d_run%d_low.root" % (dsi.calSkimDir, ds, run)
                    if os.path.isfile(fname):
                        fileList.append(fname)
                    else:
                        print("File not found:",fname)
                        continue
    if testMode:
        fileList = []
        ds = 1
        dsMap = bkg.dsMap()
        for sub in range(dsMap[ds]+1):
            fname = "%s/skimDS%d_%d_low.root" % (dsi.skimDir, ds, sub)
            fileList.append(fname)


    # check first and last few entries of every file
    for idx, fname in enumerate(fileList[:fLimit]):
        f = TFile(fname)
        t = f.Get("skimTree")
        n = t.GetEntries()

        if n==0:
            print("No entries in file",fname)
            continue

        brSingle, brVector = [], []
        for br in t.GetListOfBranches():
            if "vector" in br.GetClassName():
                brVector.append(br.GetName())
            else:
                brSingle.append(br.GetName())

        if n < 20:
            eList = np.arange(0,n,1)
        else:
            eList = np.arange(0,20,1)
            eList = np.append(eList, np.arange(n-20,n,1))

        for i in eList:
            t.GetEntry(i)

            # make sure individual entries are accessible (no segfaults)
            for br in brSingle:
                val = getattr(t,br)
                # print(val)
            for br in brVector:
                try:
                    nCh = getattr(t,br).size()
                except AttributeError:
                    print("Error:",br)
                for j in range(nCh):
                    val = getattr(t,br)[j]
                    # print(br,nCh,val)

        if verbose:
            print("%d/%d  %s  nEnt %d" % (idx, len(fileList[:fLimit]), fname.split("/")[-1], t.GetEntries()))

        f.Close()


def checkWave():
    """ ./check-files.py [-c] -w """
    print("Checking waves.  Cal?", checkCal)

    fileList = []

    # bkg
    if not checkCal:
        dsMap = bkg.dsMap()
        for ds in dsMap:
            for sub in range(dsMap[ds]+1):
                fname = "%s/waveSkimDS%d_%d.root" % (dsi.waveDir, ds, sub)
                if os.path.isfile(fname):
                    fileList.append(fname)
                else:
                    print("File not found:",fname)
                    continue
    # cal
    else:
        for key in cal.GetKeys():
            ds = int(key[2])
            if skipDS6Cal and ds==6: continue
            for cIdx in range(cal.GetIdxs(key)): # 0-indexed
                cRuns = cal.GetCalList(key, cIdx)
                for run in cRuns:
                    fname = "%s/waveSkimDS%d_run%d.root" % (dsi.calWaveDir, ds, run)
                    if os.path.isfile(fname):
                        fileList.append(fname)
                    else:
                        print("File not found:",fname)
                        continue
    if testMode:
        fileList = []
        ds = 1
        dsMap = bkg.dsMap()
        for sub in range(dsMap[ds]+1):
            fname = "%s/waveSkimDS%d_%d.root" % (dsi.waveDir, ds, sub)
            fileList.append(fname)


    # check first and last few entries of every file
    for idx, fname in enumerate(fileList[:fLimit]):
        f = TFile(fname)
        t = f.Get("skimTree")
        n = t.GetEntries()

        if n==0:
            print("No entries in file",fname)
            continue

        brSingle, brVector = [], []
        for br in t.GetListOfBranches():
            if "vector" in br.GetClassName():
                brVector.append(br.GetName())
            else:
                brSingle.append(br.GetName())

        if n < 20:
            eList = np.arange(0,n,1)
        else:
            eList = np.arange(0,20,1)
            eList = np.append(eList, np.arange(n-20,n,1))

        for i in eList:
            t.GetEntry(i)

            # make sure individual entries are accessible (no segfaults)
            for br in brSingle:
                val = getattr(t,br)
                # print(val)
            for br in brVector:
                try:
                    nCh = getattr(t,br).size()
                except AttributeError:
                    print("Error:",br)
                for j in range(nCh):
                    val = getattr(t,br)[j]
                    # print(br,nCh,val)

            # make sure we can process waveforms
            for j in range(nCh):
                wf = t.MGTWaveforms.at(j)
                ch = t.channel.at(j)

                # be absolutely sure you're matching the right waveform to this hit
                if wf.GetID() != ch:
                    print("ERROR -- Vector matching failed. entry %d, file: %s" % (i, fname))

                # run the LAT routine to convert into numpy arrays
                truncLo, truncHi = 0, 2
                if ds==6 or ds==2: truncLo = 4
                signal = wl.processWaveform(wf, truncLo, truncHi)

        if verbose:
            print("%d/%d  %s  nEnt %d" % (idx, len(fileList[:fLimit]), fname.split("/")[-1], t.GetEntries()))

        f.Close()


def checkSplit():
    """ ./check-files.py [-c] -p """
    # make sure every file has a TCut
    # do a sum file size check - make sure sum(size(splits)) = size(waves)
    print("Checking splits.  Cal?", checkCal)

    # bkg
    if not checkCal:
        dsMap = bkg.dsMap()
        for ds in dsMap:
            for sub in range(dsMap[ds]+1)[:fLimit]:

                # dict: {'DSX_X_X':filePath} to list[file paths]
                tmpList = dsi.getSplitList("%s/splitSkimDS%d_%d*" % (dsi.splitDir, ds, sub), sub)
                sList = [f for idx, f in sorted(tmpList.items())]

                # check file existence
                wFile = "%s/waveSkimDS%d_%d.root" % (dsi.waveDir, ds, sub)
                if len(sList)==0 or not os.path.isfile(wFile):
                    print("Files not found! len(split): %d, wave file:" % len(sList), wFile)
                    continue

                # check file sizes - result: wave files are slightly larger.  weird.
                # I guess ROOT has some kinda overhead for larger files?
                wSize = os.stat(wFile).st_size/1e6
                sSize = sum([os.stat(f).st_size for f in sList])/1e6

                # critical check - number of entries match, and TCut is the same
                f1 = TFile(wFile)
                t1 = f1.Get("skimTree")
                wCut = f1.Get("theCut").GetTitle()
                nEnt = t1.GetEntries()
                f1.Close()

                nEntS = 0
                for f in sList:
                    f2 = TFile(f)
                    t2 = f2.Get("skimTree")
                    sCut = f2.Get("theCut").GetTitle()
                    nEntS += t2.GetEntries()
                    f2.Close()

                if nEnt != nEntS:
                    print("Error: nEntries don't match!  ds %d  sub %d  size(wave) %.2f  size(split) %.2f  nEnt(wave) %d  nEnt(split) %d" % (ds,sub,wSize,sSize,nEnt,nEntS))

                if wCut != sCut:
                    print("Error: TCut doesn't match!\nwaveFile:%s\nsplitFile:%s" % (wCut,sCut))

                if verbose:
                    print("%d/%d  DS-%d-%d nEnt %d" % (sub, dsMap[ds]+1, ds,sub,nEnt))

    # cal
    else:
        dsMap = cal.GetKeys()
        for key in dsMap:
            ds = int(key[2])
            if skipDS6Cal and ds==6: continue
            for sub in range(cal.GetIdxs(key))[:fLimit]:
                cRuns = cal.GetCalList(key, sub)
                for run in cRuns:
                    tmpList = dsi.getSplitList("%s/splitSkimDS%d_run%d*" % (dsi.calSplitDir, ds, run), run)
                    sList = [f for idx, f in sorted(tmpList.items())]

                    # check file existence
                    wFile = "%s/waveSkimDS%d_run%d.root" % (dsi.calWaveDir, ds, run)
                    if len(sList)==0 or not os.path.isfile(wFile):
                        print("Files not found! len(split): %d, wave file:" % len(sList), wFile)
                        continue

                    # check file sizes - result: wave files are slightly larger.  weird.
                    # i guess ROOT has some kinda overhead for larger files?
                    wSize = os.stat(wFile).st_size/1e6
                    sSize = sum([os.stat(f).st_size for f in sList])/1e6

                    # critical check - number of entries match, and TCut is the same
                    f1 = TFile(wFile)
                    t1 = f1.Get("skimTree")
                    wCut = f1.Get("theCut").GetTitle()
                    nEnt = t1.GetEntries()
                    f1.Close()

                    nEntS = 0
                    for f in sList:
                        f2 = TFile(f)
                        t2 = f2.Get("skimTree")
                        sCut = f2.Get("theCut").GetTitle()
                        nEntS += t2.GetEntries()
                        f2.Close()

                    if nEnt != nEntS:
                        print("Error: nEntries don't match!  ds %d  sub %d  size(wave) %.2f  size(split) %.2f  nEnt(wave) %d  nEnt(split) %d" % (ds,sub,wSize,sSize,nEnt,nEntS))

                    if wCut != sCut:
                        print("Error: TCut doesn't match!\nwaveFile:%s\nsplitFile:%s" % (wCut,sCut))

                    if verbose:
                        print("%d/%d  DS-%d-c%d  run %d  nEnt %d" % (sub, cal.GetIdxs(key), ds, sub, run, nEnt))


def checkLAT():
    """ ./check-files.py [-c] -l """
    # make sure nLat matches nSplit files
    # repeat the checkWave checks
    print("Checking lats.  Cal?", checkCal)

    fileList = []

    # bkg
    if not checkCal:
        dsMap = bkg.dsMap()
        for ds in dsMap:
            for sub in range(dsMap[ds]+1):

                sList = dsi.getSplitList("%s/splitSkimDS%d_%d*" % (dsi.splitDir, ds, sub), sub)
                latList = dsi.getSplitList("%s/latSkimDS%d_%d*" % (dsi.latDir, ds, sub), sub)
                if len(sList) != len(latList):
                    print("Error: ds %d sub %d.  Found %d split files but %d lat files." % (ds,sub,len(sList),len(latList)))

                tmpList = [f for idx, f in sorted(latList.items())]
                fileList.extend(tmpList)
    # cal
    else:
        dsMap = cal.GetKeys()
        for key in dsMap:
            ds = int(key[2])
            if skipDS6Cal and ds==6: continue
            for sub in range(cal.GetIdxs(key)):
                cRuns = cal.GetCalList(key, sub)
                for run in cRuns:

                    if not (ds==5 and run>=25506): continue

                    sList = dsi.getSplitList("%s/splitSkimDS%d_run%d*" % (dsi.calSplitDir, ds, run), run)
                    latList = dsi.getSplitList("%s/latSkimDS%d_run%d*" % (dsi.calLatDir, ds, run), run)
                    if len(sList) != len(latList):
                        print("Error: ds %d  sub %d  run %d.  Found %d split files but %d lat files." % (ds,sub,run,len(sList),len(latList)))

                    tmpList = [f for idx, f in sorted(latList.items())]
                    fileList.extend(tmpList)
    if testMode:
        fileList = []
        ds = 1
        # dsMap = bkg.dsMap()
        # for sub in range(dsMap[ds]+1):
        #     latList = dsi.getSplitList("%s/latSkimDS%d_%d*" % (dsi.latDir, ds, sub), sub)
        #     tmpList = [f for idx, f in sorted(latList.items())]
        #     fileList.extend(tmpList)
        dsMap = cal.GetKeys()
        key = dsMap[1]
        for sub in range(cal.GetIdxs(key)):
            cRuns = cal.GetCalList(key, sub)
            for run in cRuns:
                latList = dsi.getSplitList("%s/latSkimDS%d_run%d*" % (dsi.calLatDir, ds, run), run)
                fileList.extend([f for idx, f in sorted(latList.items())])


    # loop over files, repeating the same checks in checkWave
    for idx, fname in enumerate(fileList[:fLimit]):
        f = TFile(fname)
        t = f.Get("skimTree")
        n = t.GetEntries()

        if n==0:
            print("No entries in file",fname)
            continue

        brSingle, brVector, brNames = [], [], []
        for br in t.GetListOfBranches():
            if "vector" in br.GetClassName():
                brVector.append(br.GetName())
            else:
                brSingle.append(br.GetName())
            brNames.append(br.GetName())

        # make sure a typical LAT branch exists
        if "fitSlo" not in brNames:
            print("fitSlo branch not found in:",fname)

        if n < 20:
            eList = np.arange(0,n,1)
        else:
            eList = np.arange(0,20,1)
            eList = np.append(eList, np.arange(n-20,n,1))
        for i in eList:
            t.GetEntry(i)

            # make sure individual entries are accessible (no segfaults)
            for br in brSingle:
                val = getattr(t,br)
                # print(val)
            for br in brVector:
                try:
                    nCh = getattr(t,br).size()
                except AttributeError:
                    print("Error:",br)
                for j in range(nCh):
                    val = getattr(t,br)[j]
                    # print(br,nCh,val)

            # make sure we can process waveforms
            for j in range(nCh):
                wf = t.MGTWaveforms.at(j)
                ch = t.channel.at(j)

                # be absolutely sure you're matching the right waveform to this hit
                if wf.GetID() != ch:
                    print("ERROR -- Vector matching failed. entry %d, file: %s" % (i, fname))

                # run the LAT routine to convert into numpy arrays
                truncLo, truncHi = 0, 2
                if ds==6 or ds==2: truncLo = 4
                signal = wl.processWaveform(wf, truncLo, truncHi)

        if verbose:
            print("%d/%d  %s  nEnt %d" % (idx, len(fileList[:fLimit]), fname.split("/")[-1], t.GetEntries()))

        f.Close()


def unpackFileName(splitList):
    """ Takes output of dsi.getSplitList and returns a list [ds,sub1,sub2]. """
    tmpS = []
    for idx, f in sorted(splitList.items()):
        tmp = f.split("DS")[-1].split("_")
        tmp[-1] = tmp[-1].split(".root")[0]
        if "run" in tmp[1]: tmp[1] = tmp[1].split("run")[-1]
        tmp = [int(i) for i in tmp]
        if len(tmp)==2: tmp.append(0)
        tmpS.append(tuple(tmp))
    return tmpS


def makeJobList():
    """ ./check-files.py -m """

    # bkg
    if not checkCal:

        # open the job list to compare
        with open("./jobs/bkgLAT.ls") as f:
            latJobs = f.readlines()
        latJobs = [j.rstrip() for j in latJobs]

        # write failed jobs to a cleanup list
        fOut = open("./jobs/bkgLAT_cleanup.ls","w")

        dsMap = bkg.dsMap()
        for ds in dsMap:
            for sub in range(dsMap[ds]+1):
                splList = dsi.getSplitList("%s/splitSkimDS%d_%d*" % (dsi.splitDir, ds, sub), sub)
                latList = dsi.getSplitList("%s/latSkimDS%d_%d*" % (dsi.latDir, ds, sub), sub)

                if len(splList) != len(latList):
                    print("Error: ds %d sub %d.  Found %d split files but %d lat files." % (ds,sub,len(splList),len(latList)))
                    spl = unpackFileName(splList)
                    lat = unpackFileName(latList)
                    diffs = list(set(spl).symmetric_difference(set(lat)))
                    for diff in diffs:
                        fname = "latSkimDS%d_%d_%d.root" % (diff[0],diff[1],diff[2])
                        job = [s for s in latJobs if fname in s][0]
                        fOut.write(job+"\n")
        fOut.close()

    # cal
    else:

        # open the job lists to compare
        with open("./jobs/calLAT.ls") as f:
            latJobs = f.readlines()
        latJobs = [j.rstrip() for j in latJobs]
        with open("./jobs/calLAT_ds5c.ls") as f:
            latJobs2 = f.readlines()
        latJobs.extend([j.rstrip() for j in latJobs2])

        # write failed jobs to a cleanup list
        fOut = open("./jobs/calLAT_cleanup.ls","w")

        dsMap = cal.GetKeys()
        for key in dsMap:
            ds = int(key[2])
            if skipDS6Cal and ds==6: continue
            for sub in range(cal.GetIdxs(key)):
                cRuns = cal.GetCalList(key, sub)
                for run in cRuns:

                    splList = dsi.getSplitList("%s/splitSkimDS%d_run%d*" % (dsi.calSplitDir, ds, run), run)
                    latList = dsi.getSplitList("%s/latSkimDS%d_run%d*" % (dsi.calLatDir, ds, run), run)

                    if len(splList) != len(latList):
                        print("Error: ds %d run %d.  Found %d split files but %d lat files." % (ds,run,len(splList),len(latList)))
                        spl = unpackFileName(splList)
                        lat = unpackFileName(latList)
                        diffs = list(set(spl).symmetric_difference(set(lat)))
                        for diff in diffs:
                            fname = "latSkimDS%d_run%d_%d.root" % (diff[0],diff[1],diff[2])
                            job = [s for s in latJobs if fname in s][0]
                            fOut.write(job+"\n")
        fOut.close()


def testDraw():
    """ ./check-files.py -t
    Make sure we can actually do a Draw on the final product! """
    from ROOT import TChain

    # don't have this pkg in the shifter image
    # import psutil
    # process = psutil.Process(os.getpid())

    # METHOD 1: load massive chainz dawg.
    # this segfaults on Draw() if the chain's boot too big (not enough mem available)
    # the default batch options can't handle this for most of the datasets. (login node can ...)
    # ds = 1
    # dsMap = bkg.dsMap()
    # ch = TChain("skimTree")
    # for sub in range(dsMap[ds]+1):
    #     latList = dsi.getSplitList("%s/latSkimDS%d_%d*" % (dsi.latDir, ds, sub), sub)
    #     tmpList = [f for idx, f in sorted(latList.items())]
    #     for f in tmpList: ch.Add(f)
    # nEnt = ch.GetEntries()
    # print("nEnt",nEnt)
    # n = ch.Draw("trapENFCal:riseNoise:fitSlo","trapENFCal < 250","goff")
    # t1, t2, t3 = ch.GetV1(), ch.GetV2(), ch.GetV3()
    # t1 = np.asarray([t1[i] for i in range(n)])
    # t2 = np.asarray([t2[i] for i in range(n)])
    # print("nEnt %-8d nDraw %-8d  mem: %d MB" % (nEnt, n, process.memory_info().rss/1e6))

    # METHOD 2: file-by-file, to save memory
    fileList = []
    dsMap = bkg.dsMap()
    # for ds in dsMap:
    for ds in [1]:
        for sub in range(dsMap[ds]+1):
            latList = dsi.getSplitList("%s/latSkimDS%d_%d*" % (dsi.latDir, ds, sub), sub)
            tmpList = [f for idx, f in sorted(latList.items())]
            fileList.extend(tmpList)
    for f in fileList:
        tf = TFile(f)
        tt = tf.Get("skimTree")
        nEnt = tt.GetEntries()
        n = tt.Draw("trapENFCal:riseNoise:fitSlo","trapENFCal < 250","goff")
        t1, t2, t3 = tt.GetV1(), tt.GetV2(), tt.GetV3()
        t1 = np.asarray([t1[i] for i in range(n)])
        t2 = np.asarray([t2[i] for i in range(n)])
        # print("%s: nEnt %-8d nDraw %-8d  mem: %d MB" % (f.split("/")[-1], nEnt, n, process.memory_info().rss/1e6))
        print("%s: nEnt %-8d nDraw %-8d" % (f.split("/")[-1], nEnt, n))
        tf.Close()


if __name__=="__main__":
    main(sys.argv[1:])