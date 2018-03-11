#!/usr/bin/env python3
"""
========================= job-panda.py ========================
A cute and adorable way to do various processing tasks.

The functions are arranged mostly sequentially, i.e. this file
documents the "procedure" necessary to produce LAT data,
starting from running skim_mjd_data.

====================== C. Wiseman, B. Zhu =====================
"""
import sys, shlex, glob, os, re, time
import subprocess as sp
import DataSetInfo as ds

nc = {"pdsf":[30,33], "cori":[60,66], "edison":[48,52]}
jobQueue = ds.latSWDir+"/job.queue"

# =============================================================
def main(argv):

    # defaults
    global jobStr, useJobQueue
    # jobStr = "sbatch slurm-job.sh" # SLURM mode
    jobStr = "sbatch pdsf.slr" # SLURM+Shifter mode
    dsNum, subNum, runNum, argString, calList, useJobQueue = None, None, None, None, [], False

    # loop over user args
    for i,opt in enumerate(argv):

        # set queue mode
        if opt == "-q":
            useJobQueue = True
            jobStr = ""

        # set ds, sub, cal
        if opt == "-ds":  dsNum = int(argv[i+1])
        if opt == "-sub": dsNum, subNum = int(argv[i+1]), int(argv[i+2])
        if opt == "-run": dsNum, runNum = int(argv[i+1]), int(argv[i+2])
        if opt == "-cal": calList = getCalRunList(dsNum,subNum,runNum)

        # main skim routines
        if opt == "-skim":      runSkimmer(dsNum, subNum, runNum, calList=calList)
        if opt == "-wave":      runWaveSkim(dsNum, subNum, runNum, calList=calList)
        if opt == "-qsubSplit": qsubSplit(dsNum, subNum, runNum, calList=calList)
        if opt == "-writeCut":  writeCut(dsNum, subNum, runNum, calList=calList)
        if opt == "-lat":       runLAT(dsNum, subNum, runNum, calList=calList)
        if opt == "-pandify":   pandifySkim(dsNum, subNum, runNum, calList=calList)

        # mega modes
        if opt == "-mskim":  [runSkimmer(i) for i in range(0,5+1)]
        if opt == "-mwave":  [runWaveSkim(i) for i in range(0,5+1)]
        if opt == "-msplit": [qsubSplit(i) for i in range(0,5+1)]
        if opt == "-mlat":   [runLAT(i) for i in range(0,5+1)]
        if opt == "-mcut":   [writeCut(i) for i in range(0,5+1)]

        # misc
        if opt == "-split":     splitTree(dsNum, subNum, runNum)
        if opt == "-skimLAT":   skimLAT(argv[i+1],argv[i+2],argv[i+3])
        if opt == "-tuneCuts":  tuneCuts(argv[i+1],dsNum)
        if opt == "-lat3":      applyCuts(int(argv[i+1]),argv[i+2])
        if opt == "-cron":      cronJobs()
        if opt == "-shifter":   shifterTest()
        if opt == "-test":      quickTest()
        if opt == "-b":         runBatch()

        # process systematic study runs
        if opt == "-sskim":  specialSkim()
        if opt == "-swave":  specialWave()
        if opt == "-ssplit": specialSplit()
        if opt == "-splitf": splitFile(argv[i+1],argv[i+2])
        if opt == "-swrite": specialWrite()
        if opt == "-sdel":   specialDelete()
        if opt == "-slat":   specialLAT()
        if opt == "-scheck": specialCheck()
        if opt == "-sbuild": specialBuild()


# =============================================================

def sh(cmd):
    """ Either call the shell or put the command in our job queue.
        Uses the global bool 'useJobQueue' and the global path 'jobQueue'.
    """
    if not useJobQueue:
        sp.call(shlex.split(cmd))
        return

    with open(jobQueue, 'a+') as f:
        if cmd not in open(jobQueue).read():
            print("+job.queue: %s" % (cmd))
            f.write(cmd + "\n")


def getSBatch(opt):
    """
    http://www.nersc.gov/users/computational-systems/cori/running-jobs/batch-jobs/
    http://www.nersc.gov/users/computational-systems/edison/running-jobs/batch-jobs/
    http://www.nersc.gov/users/computational-systems/pdsf/using-slurm-pdsf-batch-sbatch/
    """
    batchOpts = {
        "chos": [
            "--workdir=%s" % (ds.latSWDir),
            "--output=%s/logs/chos-%%j.txt" % (ds.latSWDir),
            "-p shared-chos",
            "-t 24:00:00",
            "--mem=3400",
        ],
        "pdsf-single": [
            "--workdir=%s" % (ds.latSWDir),
            "--output=%s/logs/pdsf-%%j.txt" % (ds.latSWDir),
            "--image=wisecg/mjsw:v2",
            "-p shared",
            "-t 24:00:00",
        ],
        "pdsf-pump": [
            "--workdir=%s" % (ds.latSWDir),
            "--output=%s/logs/pdsf-%%j.txt" % (ds.latSWDir),
            "--image=wisecg/mjsw:v2",
            "-p shared",
            "-n%d" % nc["pdsf"][0],
            "-t 24:00:00",
        ],
        "cori": [
            "--workdir=%s" % (ds.latSWDir),
            "--output=%s/logs/cori-%%j.txt" % (ds.latSWDir),
            "--image=wisecg/mjsw:v2",
            "-C haswell",
            "-N 1",
            "--qos=debug",
            # "--qos=regular",
            # "-t 24:00:00",
            "-t 00:10:00"
        ],
        "edison": [
            "--workdir=%s" % (ds.latSWDir),
            "--output=%s/logs/edison-%%j.txt" % (ds.latSWDir),
            "--image=wisecg/mjsw:v2",
            "-t 10:00:00",
            # "-t 00:10:00",
            # "--qos=debug"
            "--qos=regular"
        ]
    }
    return "sbatch "+' '.join(batchOpts[opt])


def runBatch():
    """ ./job-panda.py -b
    http://www.nersc.gov/users/computational-systems/cori/running-jobs/queues-and-policies/
    http://www.nersc.gov/users/computational-systems/edison/running-jobs/queues-and-policies/
    http://www.nersc.gov/users/accounts/user-accounts/how-usage-is-charged/
    """
    # EX. 1 - single job -- chos
    # cmd = "./skim_mjd_data -f 22513 -l -t 0.7 %s/skim" % (ds.specialDir)
    # sh("%s slurm-job.sh %s" % (getSBatch("chos"), cmd))

    # EX. 2 - run 'make clean && make' in a new env (pdsf-single, cori, edison)
    # sh("make clean")
    # sh("%s slurm.slr make" % (getSBatch("edison"),cmd))

    # EX. 3 - single job -- pdsf-single / edison / cori
    # cmd = "./skim_mjd_data -f 22513 -l -t 0.7 %s/skim" % (ds.specialDir)
    # sh("%s slurm.slr %s" % (getSBatch("pdsf-single"),cmd))

    # EX. 4 - run a job pump -- pdsf-single / edison / cori
    # sbatch = "sbatch "+' '.join(batchOpts["cori"])
    # sh("%s slurm.slr './job-pump.sh jobLists/test.list skim_mjd_data %d %d'" % (getSBatch("cori"), nC["cori"][0], nC["cori"][1]))

    # EX. 5 - skim long cal on edison
    # nCores, peakLoad = nC["edison"][0], nC["edison"][1]
    # sh("%s slurm.slr './job-pump.sh jobLists/skimLongCal.list skim_mjd_data %d %d'" % (getSBatch("edison"), nCores, peakLoad))
    # sh("%s slurm.slr './job-pump.sh jobLists/waveLongCal.list wave-skim %d %d'" % (getSBatch("edison"), nCores, peakLoad))
    # sh("%s slurm.slr './job-pump.sh jobLists/splitLongCal.list python3 %d %d'" % (getSBatch("edison"), nCores, peakLoad))
    # sh("%s slurm.slr './job-pump.sh jobLists/test.list python3 %d %d'" % (getSBatch("edison"), nCores, peakLoad))
    # sh("%s slurm.slr './job-pump.sh jobLists/latLongCal.list python3 %d %d'" % (getSBatch("edison"), nCores, peakLoad))
    # sh("%s slurm.slr './job-pump.sh jobLists/latLongCal_cleanup.list python3 %d %d'" % (getSBatch("edison"), nCores, peakLoad))

    nCores, peakLoad = nc["pdsf"][0], nc["pdsf"][1]

    # EX. 6 - bkg file sequence
    # sh("%s slurm.slr './job-pump.sh jobLists/skimBkg.list skim_mjd_data %d %d'" % (getSBatch("pdsf-pump"), nCores, peakLoad))
    # sh("%s slurm.slr './job-pump.sh jobLists/waveBkg.list wave-skim %d %d'" % (getSBatch("pdsf-pump"), nCores, peakLoad))

    # EX. 7 - cal file sequence
    # sh("%s slurm.slr './job-pump.sh jobLists/skimCal.list skim_mjd_data %d %d'" % (getSBatch("pdsf-pump"), nCores, peakLoad))
    # sh("%s slurm.slr './job-pump.sh jobLists/waveCal.list wave-skim %d %d'" % (getSBatch("pdsf-pump"), nCores, peakLoad))
    # sh("%s slurm.slr './job-pump.sh jobLists/test.list wave-skim %d %d'" % (getSBatch("pdsf-pump"), nCores, peakLoad))

    # EX. 8 - ext pulser
    sh("%s slurm.slr './job-pump.sh jobLists/extPulserLAT.list python3 %d %d'" % (getSBatch("edison"), nCores, peakLoad))


def getCalRunList(dsNum=None,subNum=None,runNum=None):
    """ ./job-panda.py -cal (-ds [dsNum] -sub [dsNum] [calIdx] -run [runNum])
        Create a calibration run list, using the CalInfo object in DataSetInfo.py .
        Note that the -sub option is re-defined here to mean a calibration range idx.
        Note that running with -cal alone will create a list for all datasets (mega mode).
    """
    runLimit = 15 # yeah I'm hardcoding this, sue me.
    calList = []
    calInfo = ds.CalInfo()
    calKeys = calInfo.GetKeys(dsNum)

    # single-run mode
    if runNum!=None:
        calList.append(runNum)
        print(calList)
        return calList

    # multi-run mode:
    for key in calKeys:
        print("key:",key)

        # -cal (mega mode)
        if dsNum==None:
            for idx in range(calInfo.GetIdxs(key)):
                lst = calInfo.GetCalList(key,idx,runLimit)
                print(lst)
                calList += lst
        # -ds
        elif subNum==None:
            for idx in range(calInfo.GetIdxs(key)):
                lst = calInfo.GetCalList(key,idx,runLimit)
                print(lst)
                calList += lst
        # -sub
        else:
            lst = calInfo.GetCalList(key,subNum,runLimit)
            if lst==None: continue
            print(lst)
            calList += lst

    # remove any duplicates, but there probably aren't any
    calList = sorted(list(set(calList)))

    return calList


def runSkimmer(dsNum, subNum=None, runNum=None, calList=[]):
    """ ./job-panda.py [-q] -skim (-ds dsNum) (-sub dsNum subNum) (-run dsNum subNum) [-cal]
        Submit skim_mjd_data jobs.
    """
    # bkg
    if not calList:
        # -ds
        if subNum==None and runNum==None:
            for i in range(ds.dsMap[dsNum]+1):
                job = "./skim_mjd_data %d %d -l -g -t 0.7 %s" % (dsNum, i, ds.skimDir)
                if useJobQueue: sh("%s >& ./logs/skim-ds%d-%d.txt" % (job, dsNum, i))
                else: sh("%s '%s'" % (jobStr, job))
        # -sub
        elif runNum==None:
            job = "./skim_mjd_data %d %d -l -g -t 0.7 %s" % (dsNum, subNum, ds.skimDir)
            if useJobQueue: sh("%s >& ./logs/skim-ds%d-%d.txt" % (job, dsNum, subNum))
            else: sh("%s '%s'" % (jobStr, job))
        # -run
        elif subNum==None:
            job = "./skim_mjd_data -f %d -l -g -t 0.7 %s" % (runNum, ds.skimDir)
            if useJobQueue: sh("%s >& ./logs/skim-ds%d-run%d.txt" % (job, ds.GetDSNum(run), run))
            else: sh("%s '%s'" % (jobStr, job))
    # cal
    else:
        for run in calList:
            job = "./skim_mjd_data -f %d -l -g -t 0.7 %s" % (run, ds.calSkimDir)
            if useJobQueue: sh("%s >& ./logs/skim-ds%d-run%d.txt" % (job, ds.GetDSNum(run),run))
            else: sh("%s '%s'" % (jobStr, job))


def runWaveSkim(dsNum, subNum=None, runNum=None, calList=[]):
    """ ./job-panda.py [-q] -wave (-ds dsNum) (-sub dsNum subNum) (-run dsNum subNum) [-cal]
        Submit wave-skim jobs.
    """
    # bkg
    if not calList:
        # -ds
        if subNum==None and runNum==None:
            for i in range(ds.dsMap[dsNum]+1):
                job = "./wave-skim -n -r %d %d -p %s %s" % (dsNum, i, ds.skimDir, ds.waveDir)
                if useJobQueue: sh("%s >& ./logs/wave-ds%d-%d.txt" % (job, dsNum, i))
                else: sh("%s '%s'" % (jobStr, job))
        # -sub
        elif runNum==None:
            job = "./wave-skim -n -r %d %d -p %s %s" % (dsNum, subNum, ds.skimDir, ds.waveDir)
            if useJobQueue: sh("%s >& ./logs/wave-ds%d-%d.txt" % (job, dsNum, subNum))
            else: sh("%s '%s'" % (jobStr, job))
        # -run
        elif subNum==None:
            job = "./wave-skim -n -f %d %d -p %s %s" % (dsNum, dsNum, runNum, ds.skimDir, ds.waveDir)
            if useJobQueue: sh("%s >& ./logs/wave-ds%d-run%d.txt" % (job, dsNum, runNum))
            else: sh("%s '%s'" % (jobStr, job))
    # cal
    else:
        for run in calList:
            job = "./wave-skim -n -c -f %d %d -p %s %s" % (ds.GetDSNum(run), run, ds.calSkimDir, ds.calWaveDir)
            if useJobQueue: sh("%s >& ./logs/wave-ds%d-run%d.txt" % (job, ds.GetDSNum(run), run))
            else: sh("%s '%s'" % (jobStr, job))


def getFileList(filePathRegexString, subNum, uniqueKey=False, dsNum=None):
    """ Creates a dict of files w/ the format {'DSX_X_X':filePath.}
        Used to combine and split apart files during the LAT processing.
        Used by: splitTree, writeCut, runLAT.
    """
    files = {}
    for fl in glob.glob(filePathRegexString):
        int(re.search(r'\d+',fl).group())
        ints = map(int, re.findall(r'\d+',fl))
        if (ints[1]==subNum):
            if (len(ints)==2):
                ints.append(0)
            if not uniqueKey:
                files[ints[2]] = fl # zero index
            else:
                files["DS%d_%d_%d" % (dsNum,subNum,ints[2])] = fl
    return files


def splitTree(dsNum, subNum=None, runNum=None):
    """ ./job-panda.py -split (-sub dsNum subNum) (-run dsNum runNum)

        Split a SINGLE waveSkim file into small (~50MB) files to speed up LAT parallel processing.
        Can call 'qsubSplit' instead to submit each run in the list as a job, splitting the files in parallel.
        NOTE: The cut written into the first file is NOT copied into the additional files
              (I couldnt get it to work within this function -- kept getting "file not closed" errors.)
              To clean up, do that with the 'writeCut' function below, potentially AFTER a big parallel job.
    """
    from ROOT import TFile, TTree, gDirectory, TEntryList, TNamed, TObject, gROOT

    # Set input and output paths.  Clear out any files from a previous
    # try before you attempt a copy (avoid the double underscore)
    inPath, outPath = "", ""
    if runNum==None:
        # bg mode
        inPath = "%s/waveSkimDS%d_%d.root" % (ds.waveDir,dsNum,subNum)
        outPath = "%s/split/splitSkimDS%d_%d.root" % (ds.waveDir,dsNum,subNum)
        fileList = getFileList("%s/split/splitSkimDS%d_%d*.root" % (ds.waveDir,dsNum,subNum),subNum)
        for key in fileList: os.remove(fileList[key])
    elif subNum==None:
        # cal mode
        inPath = "%s/waveSkimDS%d_run%d.root" % (ds.calWaveDir,dsNum,runNum)
        outPath = "%s/split/splitSkimDS%d_run%d.root" % (ds.calWaveDir,dsNum,runNum)
        fileList = getFileList("%s/split/splitSkimDS%d_run%d*.root" % (ds.calWaveDir,dsNum,runNum),runNum)
        for key in fileList: os.remove(fileList[key])

    inFile = TFile(inPath)
    bigTree = inFile.Get("skimTree")
    theCut = inFile.Get("theCut").GetTitle()
    bigTree.Draw(">>elist",theCut,"entrylist")
    elist = gDirectory.Get("elist")
    bigTree.SetEntryList(elist)
    nList = elist.GetN()

    outFile = TFile(outPath,"RECREATE")
    lilTree = TTree()
    lilTree.SetMaxTreeSize(50000000) # 50 MB
    thisCut = TNamed("theCut",theCut)
    thisCut.Write("",TObject.kOverwrite)
    lilTree = bigTree.CopyTree("") # this does NOT write the cut into the extra files
    lilTree.Write("",TObject.kOverwrite)


def qsubSplit(dsNum, subNum=None, runNum=None, calList=[]):
    """ ./job-panda.py -qsubSplit (-ds dsNum) (-sub dsNum subNum) (-run dsNum subNum) [-cal]
        Submit jobs that call splitTree for each run, splitting files into small (~100MB) chunks.
        NOTE: The data cleaning cut is NOT written into the output files and the
              function 'writeCut' must be called after these jobs are done.
    """
    from shutil import copyfile

    # bg
    if not calList:
        # -ds
        if subNum==None and runNum==None:
            for i in range(ds.dsMap[dsNum]+1):
                inPath = "%s/waveSkimDS%d_%d.root" % (ds.waveDir,dsNum,i)
                if not os.path.isfile(inPath):
                    print("File",inPath,"not found. Continuing ...")
                    continue
                if (os.path.getsize(inPath)/1e6 < 45):
                    copyfile(inPath, "%s/splitSkimDS%d_%d.root" % (ds.splitDir, dsNum, i))
                else:
                    sh("""%s './job-panda.py -split -sub %d %d'""" % (jobStr, dsNum, i))
        # -sub
        elif runNum==None:
            inPath = "%s/waveSkimDS%d_%d.root" % (ds.waveDir,dsNum,subNum)
            if not os.path.isfile(inPath):
                print("File",inPath,"not found.")
                return
            if (os.path.getsize(inPath)/1e6 < 45):
                copyfile(inPath, "%s/splitSkimDS%d_%d.root" % (ds.splitDir, dsNum, subNum))
            else:
                sh("""%s './job-panda.py -split -sub %d %d'""" % (jobStr, dsNum, subNum))
        # -run
        elif subNum==None:
            inPath = "%s/waveSkimDS%d_run%d.root" % (ds.waveDir,dsNum,runNum)
            if not os.path.isfile(inPath):
                print("File",inPath,"not found.")
                return
            if (os.path.getsize(inPath)/1e6 < 45):
                copyfile(inPath, "%s/splitSkimDS%d_%d.root" % (ds.splitDir, dsNum, runNum))
            else:
                sh("""%s './job-panda.py -split -run %d %d'""" % (jobStr, dsNum, runNum))
    # cal
    else:
        for run in calList:
            for key in ds.dsRanges:
                if ds.dsRanges[key][0] <= run <= ds.dsRanges[key][1]:
                    dsNum=key
            inPath = "%s/waveSkimDS%d_run%d.root" % (ds.calWaveDir,dsNum,run)
            if not os.path.isfile(inPath):
                print("File",inPath,"not found. Continuing ...")
                continue
            if (os.path.getsize(inPath)/1e6 < 45):
                copyfile(inPath, "%s/splitSkimDS%d_run%d.root" % (ds.calSplitDir, dsNum, run))
            else:
                sh("""%s './job-panda.py -split -run %d %d'""" % (jobStr, dsNum, run))


def writeCut(dsNum, subNum=None, runNum=None, calList=[]):
    """ ./job-panda.py -writeCut (-ds dsNum) (-sub dsNum subNum) (-run dsNum subNum) [-cal]
        Assumes the cut used in the FIRST file (even in the whole DS) should be applied
        to ALL files.  This should be a relatively safe assumption.
    """
    from ROOT import TFile, TNamed, TObject
    mainList = {}

    # bg
    if not calList:
        # -ds
        if subNum==None and runNum==None:
            for i in range(ds.dsMap[dsNum]+1):
                inPath = "%s/splitSkimDS%d_%d*" % (ds.splitDir,dsNum,i)
                fileList = getFileList(inPath,i,True,dsNum)
                mainList.update(fileList)
        # -sub
        elif runNum==None:
            inPath = "%s/splitSkimDS%d_%d*" % (ds.splitDir,dsNum,subNum)
            fileList = getFileList(inPath,subNum,True,dsNum)
            mainList.update(fileList)
        # -run
        elif subNum==None:
            inPath = "%s/splitSkimDS%d_run%d*" % (ds.splitDir,dsNum,runNum)
            fileList = getFileList(inPath,runNum,True,dsNum)
            mainList.update(fileList)
    # cal
    else:
        for run in calList:
            for key in ds.dsRanges:
                if ds.dsRanges[key][0] <= run <= ds.dsRanges[key][1]:
                    dsNum=key
            inPath = "%s/splitSkimDS%d_run%d*" % (ds.calSplitDir,dsNum,run)
            fileList = getFileList(inPath,run,True,dsNum)
            mainList.update(fileList)

    # Pull the cut off the FIRST file and add it to the sub-files
    if len(mainList) <= 1:
        print("No files found!  Exiting...")
        exit(1)
    theCut = ""
    foundFirst = False
    for key, inFile in sorted(mainList.items()):
        if not foundFirst:
            firstFile = TFile(mainList[key])
            theCut = firstFile.Get("theCut").GetTitle()
            print("Applying this cut:\n",theCut)
            foundFirst = True
        print(key, inFile)
        subRangeFile = TFile(inFile,"UPDATE")
        thisCut = TNamed("theCut",theCut)
        thisCut.Write("",TObject.kOverwrite)


def runLAT(dsNum, subNum=None, runNum=None, calList=[]):
    """ ./job-panda.py -lat (-ds dsNum) (-sub dsNum subNum) (-run dsNum subNum) [-cal]
        Runs LAT on splitSkim output.  Does not combine output files back together.
    """
    # bg
    if not calList:
        # -ds
        if subNum==None and runNum==None:
            for i in range(ds.dsMap[dsNum]+1):
                files = getFileList("%s/splitSkimDS%d_%d*" % (ds.splitDir,dsNum,i),i)
                for idx, inFile in sorted(files.items()):
                    outFile = "%s/latSkimDS%d_%d_%d.root" % (ds.latDir,dsNum,i,idx)
                    sh("""%s './lat.py -b -r %d %d -p %s %s'""" % (jobStr,dsNum,i,inFile,outFile))
        # -sub
        elif runNum==None:
            files = getFileList("%s/splitSkimDS%d_%d*" % (ds.splitDir,dsNum,subNum),subNum)
            for idx, inFile in sorted(files.items()):
                outFile = "%s/latSkimDS%d_%d_%d.root" % (ds.latDir,dsNum,subNum,idx)
                sh("""%s './lat.py -b -r %d %d -p %s %s'""" % (jobStr,dsNum,subNum,inFile,outFile))
        # -run
        elif subNum==None:
            files = getFileList("%s/splitSkimDS%d_run%d*" % (ds.splitDir,dsNum,runNum),runNum)
            for idx, inFile in sorted(files.items()):
                outFile = "%s/latSkimDS%d_run%d_%d.root" % (ds.latDir,dsNum,runNum,idx)
                sh("""%s './lat.py -b -f %d %d -p %s %s'""" % (jobStr,dsNum,runNum,inFile,outFile))
    # cal
    else:
        for run in calList:
            for key in ds.dsRanges:
                if ds.dsRanges[key][0] <= run <= ds.dsRanges[key][1]:
                    dsNum=key
            files = getFileList("%s/splitSkimDS%d_run%d*" % (ds.calSplitDir,dsNum,run),run)
            for idx, inFile in sorted(files.items()):
                outFile = "%s/latSkimDS%d_run%d_%d.root" % (ds.calLatDir,dsNum,run,idx)
                sh("""%s './lat.py -b -f %d %d -p %s %s'""" % (jobStr,dsNum,run,inFile,outFile))


def mergeLAT():
    """ It seems like a good idea, right?
        Merging all the LAT files back together after splitting?
    """
    print("hey")


def pandifySkim(dsNum, subNum=None, runNum=None, calList=[]):
    """ ./job-panda.py -pandify (-ds dsNum) (-sub dsNum subNum) (-run dsNum subNum) [-cal]
        Run ROOTtoPandas jobs.
    """
    # bg
    if not calList:
        # -ds
        if subNum==None and runNum==None:
            for i in range(ds.dsMap[dsNum]+1):
                sh("""%s 'python3 ./sandbox/ROOTtoPandas.py -ws %d %d -p -d %s %s'""" % (jobStr, dsNum, i, ds.waveDir, ds.pandaDir))
        # -sub
        elif runNum==None:
            sh("""%s 'python3 ./sandbox/ROOTtoPandas.py -ws %d %d -p -d %s %s'""" % (jobStr, dsNum, subNum, ds.waveDir, ds.pandaDir))
        # -run
        elif subNum==None:
            sh("""%s 'python3 ./sandbox/ROOTtoPandas.py -f %d %d -p -d %s %s'""" % (jobStr, dsNum, runNum, ds.waveDir, ds.pandaDir))
    # cal
    else:
        for i in calList:
            sh("""%s 'python3 ./sandbox/ROOTtoPandas.py -f %d %d -p -d %s %s'""" % (jobStr, dsNum, i, ds.calWaveDir, ds.pandaDir))


def tuneCuts(argString, dsNum=None):
    """ ./job-panda.py -tuneCuts '[argString]' -- run over all ds's
        ./job-panda.py -ds [dsNum] -tuneCuts '[argString]' -- just one DS

    Submit a bunch of lat3.py jobs to the queues.
    NOTE:
        1) If processing individual dataset, the -ds option MUST come before -tuneCuts.
        2) Make sure to put argString in quotes.
        3) argString may be multiple options separated by spaces

    Options for argString:
        -all, -bcMax, -noiseWeight, -bcTime, -tailSlope, -fitSlo, -riseNoise
    """
    calInfo = ds.CalInfo()
    if dsNum==None:
        for i in ds.dsMap.keys():
            if i == 6: continue
            for mod in [1,2]:
                try:
                    for j in range(calInfo.GetIdxs("ds%d_m%d"%(i, mod))):
                        print("%s './lat3.py -db -tune %s -s %d %d %d %s" % (jobStr, ds.calLatDir, i, j, mod, argString))
                        sh("""%s './lat3.py -db -tune %s -s %d %d %d %s '""" % (jobStr, ds.calLatDir, i, j, mod, argString))
                except: continue
    # -ds
    else:
        for mod in [1,2]:
            try:
                for j in range(calInfo.GetIdxs("ds%d_m%d"%(dsNum, mod))):
                    print("%s './lat3.py -db -tune %s -s %d %d %d %s" % (jobStr, ds.calLatDir, dsNum, j, mod, argString))
                    sh("""%s './lat3.py -db -tune %s -s %d %d %d %s '""" % (jobStr, ds.calLatDir, dsNum, j, mod, argString))
            except: continue


def applyCuts(dsNum, cutType):
    """ ./job-panda.py -lat3 [dsNum] [cutType]"""

    if dsNum==-1:
        for ds in range(6):
            sh("""%s './lat3.py -cut %d %s'""" % (jobStr, ds, cutType))
    else:
        sh("""%s './lat3.py -cut %d %s'""" % (jobStr, dsNum, cutType))


def cronJobs():
    """ ./job-panda.py -cron
    Uses the global string 'jobQueue'.
    Crontab should contain the following lines (crontab -e):
    SHELL=/bin/bash
    MAILTO="" # can put in some address here if you LOVE emails
    #*/10 * * * * source ~/env/EnvBatch.sh; ~/lat/job-panda.py -cron >> ~/lat/cron.log 2>&1
    """
    os.chdir(home+"/lat/")
    print("Cron:",time.strftime('%X %x %Z'),"cwd:",os.getcwd())

    nMaxRun, nMaxPend = 15, 200

    with open(jobQueue) as f:
        jobList = [line.rstrip('\n') for line in f]
    nList = len(jobList)

    status = os.popen('slusers | grep wisecg').read()
    status = status.split()
    nRun = int(status[0]) if len(status) > 0 else 0  # Rjob Rcpu Rcpu*h PDjob PDcpu user:account:partition
    nPend = int(status[3]) if len(status) > 0 else 0

    nSubmit = (nMaxRun-nRun) if nRun < nMaxRun else 0
    nSubmit = nList if nList < nSubmit else nSubmit
    nSubmit = 0 if nPend >= nMaxPend else nSubmit

    print("   nRun %d  (max %d)  nPend %d (max %d)  nList %d  nSubmit %d" % (nRun,nMaxRun,nPend,nMaxPend,nList,nSubmit))

    with open(jobQueue, 'w') as f:
        for idx, job in enumerate(jobList):
            if idx < nSubmit:
                print("Submitted:",job)
                sh(job)
            else:
                # print("Waiting:", job)
                f.write(job + "\n")


def shifterTest():
    """ ./job-panda.py -shifter """
    print("Shifter:",time.strftime('%X %x %Z'),"cwd:",os.getcwd())
    sh("""sbatch shifter.slr 'python sandbox/bl2.py'""")


def specialSkim():
    """ ./job-panda.py [-q (use job queue)] -sskim """
    cal = ds.CalInfo()
    # runList = cal.GetSpecialRuns("extPulser")
    # runList = cal.GetSpecialRuns("delayedTrigger")
    runList = cal.GetSpecialRuns("longCal",5)
    for run in runList:
        if useJobQueue:
            sh("""./skim_mjd_data -f %d -l -t 0.7 %s/skim >& ./logs/specialSkim-DS%d-%d.txt""" % (run,ds.specialDir,ds.GetDSNum(run),run))
        else:
            sh("""%s './skim_mjd_data -f %d -l -t 0.7 %s/skim >& ./logs/specialSkim-DS%d-%d.txt""" % (jobStr,run,ds.specialDir,ds.GetDSNum(run),run))


def specialWave():
    """ ./job-panda.py [-q (use queue)] -swave """
    cal = ds.CalInfo()
    # runList = cal.GetSpecialRuns("extPulser")
    runList = cal.GetSpecialRuns("longCal",5)
    for run in runList:
        if useJobQueue:
            sh("""./wave-skim -l -n -f %d %d -p %s/skim %s/waves >& ./logs/wave-ds%d-%d.txt""" % (ds.GetDSNum(run),run,ds.specialDir,ds.specialDir,ds.GetDSNum(run),run))
        else:
            sh("""%s './wave-skim -x -n -f %d %d -p %s/skim %s/waves'""" % (jobStr, ds.GetDSNum(run), run, ds.specialDir, ds.specialDir) )


def specialSplit():
    """ ./job-panda.py [-q] -ssplit
    External pulser runs have no data cleaning cut.
    Has a memory leak (can't close both TFiles, damn you, ROOT); submit each run as a batch job.
    """
    cal = ds.CalInfo()
    # runList = cal.GetSpecialRuns("extPulser")
    runList = cal.GetSpecialRuns("longCal",5)
    for run in runList:
        # print(run)
        # if run <= 4592: continue

        inPath = "%s/waves/waveSkimDS%d_run%d.root" % (ds.specialDir, ds.GetDSNum(run), run)
        outPath = "%s/split/splitSkimDS%d_run%d.root" % (ds.specialDir, ds.GetDSNum(run), run)

        outFiles = glob.glob("%s/split/splitSkimDS%d_run%d*.root" % (ds.specialDir, ds.GetDSNum(run), run))
        for filename in outFiles:
            try:
                os.remove(filename)
            except OSError:
                pass

        if useJobQueue:
            sh("""./job-panda.py -splitf %s %s""" % (inPath,outPath))
        else:
            sh("""%s './job-panda.py -splitf %s %s'""" % (jobStr,inPath, outPath))


def splitFile(inPath,outPath):
    """ ./job-panda.py -splitf [inPath] [outPath]
    Used by specialSplit. """
    from ROOT import TFile, TTree, TObject
    inFile = TFile(inPath)
    bigTree = inFile.Get("skimTree")
    outFile = TFile(outPath, "RECREATE")
    lilTree = TTree()
    lilTree.SetMaxTreeSize(30000000) # 30MB
    lilTree = bigTree.CopyTree("")
    lilTree.Write("",TObject.kOverwrite)


def specialWrite():
    """ ./job-panda.py -swrite
    Write TCuts from waveSkim files into splitSkim files.
    """
    from ROOT import TFile, TNamed, TObject
    cal = ds.CalInfo()
    runList = cal.GetSpecialRuns("longCal",5)

    for run in runList:
        dsNum = ds.GetDSNum(run)
        wavePath = "%s/waves/waveSkimDS%d_run%d.root" % (ds.specialDir, dsNum, run)
        waveFile = TFile(wavePath)
        theCut = waveFile.Get("theCut").GetTitle()
        print(wavePath)

        splitFiles = glob.glob("%s/split/splitSkimDS%d_run%d*.root" % (ds.specialDir, dsNum, run))
        for idx in range(len(splitFiles)):
            if idx==0:
                splitPath = "%s/split/splitSkimDS%d_run%d.root" % (ds.specialDir, dsNum, run)
            else:
                splitPath = "%s/split/splitSkimDS%d_run%d_%d.root" % (ds.specialDir, dsNum, run, idx)
            if not os.path.isfile(splitPath) :
                print("File doesn't exist:",splitPath)
                return

            splitFile = TFile(splitPath,"UPDATE")
            thisCut = TNamed("theCut",theCut)
            thisCut.Write("",TObject.kOverwrite)
        print(splitFiles[-1])


def specialDelete():
    """./job-panda.py -sdel"""

    # remove all files for specific run numbers
    removeList = [5940, 5941, 5946, 5961, 5962, 5963, 5978, 6205]
    for run in removeList:
        outFiles = glob.glob("%s/split/splitSkimDS%d_run%d*.root" % (ds.specialDir, ds.GetDSNum(run), run))
        outFiles.extend(["%s/skim/skimDS%d_run%d_low.root" % (ds.specialDir, ds.GetDSNum(run), run)])
        outFiles.extend(["%s/waves/waveSkimDS%d_run%d.root" % (ds.specialDir, ds.GetDSNum(run), run)])
        outFiles.extend(glob.glob("%s/lat/latSkimDS%d_run%d_*.root" % (ds.specialDir, ds.GetDSNum(run), run)))
        for filename in outFiles:
            try:
                os.remove(filename)
                print(filename)
            except OSError:
                pass

    # remove all files from ext pulser range
    # cal = ds.CalInfo()
    # for idx in [6]:
    #     runList = cal.GetSpecialRuns("extPulser",idx)
    #     for run in runList:
    #         outFiles = glob.glob("%s/split/splitSkimDS%d_run%d*.root" % (ds.specialDir, ds.GetDSNum(run), run))
    #         outFiles.extend(["%s/skim/skimDS%d_run%d_low.root" % (ds.specialDir, ds.GetDSNum(run), run)])
    #         outFiles.extend(["%s/waves/waveSkimDS%d_run%d.root" % (ds.specialDir, ds.GetDSNum(run), run)])
    #         outFiles.extend(glob.glob("%s/lat/latSkimDS%d_run%d_*.root" % (ds.specialDir, ds.GetDSNum(run), run)))
    #         for filename in outFiles:
    #             print(filename)
    #             try:
    #                 os.remove(filename)
    #             except OSError:
    #                 pass

    # remove lat files without the _X.root
    # import datetime
    # cal = ds.CalInfo()
    # runList = cal.GetSpecialRuns("extPulser")
    # for run in runList:
    #     outFile = "%s/lat/latSkimDS%d_run%d.root" % (ds.specialDir, ds.GetDSNum(run), run)
    #     try:
    #         modDate = os.path.getmtime(outFile)
    #         modDate = datetime.datetime.fromtimestamp(int(modDate)).strftime('%Y-%m-%d %H:%M:%S')
    #         print(outFile, modDate)
    #         os.remove(outFile)
    #     except OSError:
    #         pass


def specialLAT():
    """ ./job-panda.py [-q (use job queue)] -slat"""
    cal = ds.CalInfo()
    runList = cal.GetSpecialRuns("extPulser")
    # runList = cal.GetSpecialRuns("longCal",5)

    # deal with unsplit files
    # run = runList[0]
    # dsNum = ds.GetDSNum(run)
    # inFile = "%s/waves/waveSkimDS%d_run%d.root" % (ds.specialDir,dsNum,run)
    # outFile = "%s/lat/latSkimDS%d_run%d.root" % (ds.specialDir,dsNum,run)
    # sh("""./lat.py -x -b -f %d %d -p %s %s""" % (dsNum,run,inFile,outFile))

    # deal with split files
    for run in runList:

        dsNum = ds.GetDSNum(run)
        inFiles = glob.glob("%s/split/splitSkimDS%d_run%d*.root" % (ds.specialDir, dsNum, run))
        for idx in range(len(inFiles)):
            if idx==0:
                inFile = "%s/split/splitSkimDS%d_run%d.root" % (ds.specialDir, dsNum, run)
            else:
                inFile = "%s/split/splitSkimDS%d_run%d_%d.root" % (ds.specialDir, dsNum, run, idx)
            if not os.path.isfile(inFile) :
                print("File doesn't exist:",inFile)
                return
            outFile = "%s/lat/latSkimDS%d_run%d_%d.root" % (ds.specialDir, dsNum, run, idx)

            if useJobQueue:
                # this is what you would want for a normal cron queue
                # sh("""./lat.py -x -b -f %d %d -p %s %s""" % (dsNum, run, inFile, outFile))

                # this is what i need for a 1-node job pump
                # sh("""./lat.py -x -b -f %d %d -p %s %s >& ./logs/extPulser-%d-%d.txt""" % (dsNum, run, inFile, outFile, run, idx))
                sh("""./lat.py -b -f %d %d -p %s %s >& ./logs/longCal-%d-%d.txt""" % (dsNum, run, inFile, outFile, run, idx))
            else:
                sh("""%s './lat.py -x -b -f %d %d -p %s %s' """ % (jobStr, dsNum, run, inFile, outFile))


def specialCheck():
    """./job-panda.py -scheck
    A next step could be to 'hadd' split files back together, but we'll wait for now.
    """
    from ROOT import TFile, TTree
    cal = ds.CalInfo()
    runList = cal.GetSpecialRuns("extPulser")
    for run in runList:
        fileList = glob.glob("%s/lat/latSkimDS%d_run%d_*.root" % (ds.specialDir, ds.GetDSNum(run), run))
        for f in fileList:
            tf = TFile(f)
            tr = tf.Get("skimTree")
            print(f)
            print(tr.GetEntries())
            tr.GetEntry(0)
            tf.Close()


def specialBuild():
    """ ./job-panda.py -sbuild
    Set 1 of the external pulser runs were built with --donotbuild.
    Let's manually build them.

    Example:
    - Make sure MkCookie has been run recently
    majorcaroot --nomultisampling --setspecialchanmap /global/homes/m/mjd/production/P3JDYspecchanmap.txt --donotbuild 2015-9-3-P3JDY_Run5942
    process_mjd_cal OR_*.root
    pulsertag 5942
    auto-thresh 5942
    process_mjd_gat OR_*.root
    """
    cal = ds.CalInfo()
    rawDir = "/global/project/projectdirs/majorana/data/mjd/surfmjd/data/raw/P3JDY/Data"
    buildDir = ds.dataDir + "/mjddatadir"
    os.chdir(buildDir)

    runList = []
    for pIdx in [7,8,9,10,11,12]:
        runList.extend(cal.GetSpecialRuns("extPulser",pIdx))

    for run in runList:

        os.chdir(buildDir + "/raw")
        rawFile = glob.glob("%s/*%d*" % (rawDir, run))[0]
        rawName = rawFile.rsplit('/',1)[1]
        rawName = rawName.rsplit('.')[0]
        # sp.call("""zcat %s > %s""" % (rawFile, rawName), shell=True)

        os.chdir(buildDir + "/built")
        sp.call("""majorcaroot --nomultisampling --setspecialchanmap /global/homes/m/mjdproduction/P3JDYspecchanmap.txt ../raw/%s""" % (rawName), shell=True)

        # ... i kinda lost steam here.  Set 1 runs are not worth this much hassle.

        return


def quickTest():
    """./job-panda.py -test """
    # from ROOT import MGTWaveform, GATDataSet, TChain
    # import datetime, time
    # print("Sleeping 10 sec ...")
    # time.sleep(10)
    # now = datetime.datetime.now()
    # print("Done. Date: ",str(now))

    cal = ds.CalInfo()
    runList = cal.GetSpecialRuns("extPulser")
    print(runList)


if __name__ == "__main__":
    main(sys.argv[1:])

