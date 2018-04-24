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
import dsi
jobQueue = dsi.latSWDir+"/job.q"

# =============================================================
def main(argv):

    # defaults
    global jobStr, useJobQueue
    # jobStr = "sbatch slurm-job.sh" # SLURM mode
    jobStr = "sbatch pdsf.slr" # SLURM + Shifter mode
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
        #if opt == "-cal": calList = getCalRunList(dsNum,subNum,runNum)
        if opt == "-cal": calList = getCalRunList(dsNum)
            print 'calibration list ', calList

        # main skim routines
        if opt == "-skim":      runSkimmer(dsNum, subNum, runNum, calList=calList)
        if opt == "-wave":      runWaveSkim(dsNum, subNum, runNum, calList=calList)
        if opt == "-thresh":    runAutoThresh(dsNum)
        if opt == "-batchSplit": batchSplit(dsNum, subNum, runNum, calList=calList)
        if opt == "-writeCut":  writeCut(dsNum, subNum, runNum, calList=calList)
        if opt == "-lat":       runLAT(dsNum, subNum, runNum, calList=calList)
        if opt == "-pandify":   pandifySkim(dsNum, subNum, runNum, calList=calList)

        # mega modes
        if opt == "-mskim":  [runSkimmer(i) for i in range(0,6+1)]
        if opt == "-mwave":  [runWaveSkim(i) for i in range(0,6+1)]
        if opt == "-msplit": [batchSplit(i) for i in range(0,6+1)]
        if opt == "-mcut":   [writeCut(i) for i in range(0,6+1)]
        if opt == "-mlat":   [runLAT(i) for i in range(0,6+1)]

        # special run routines
        if opt == "-sskim":  specialSkim()
        if opt == "-swave":  specialWave()
        if opt == "-ssplit": specialSplit()
        if opt == "-splitf": splitFile(argv[i+1],argv[i+2])
        if opt == "-swrite": specialWrite()
        if opt == "-sdel":   specialDelete()
        if opt == "-slat":   specialLAT()
        if opt == "-scheck": specialCheck()
        if opt == "-sbuild": specialBuild()

        # misc
        if opt == "-split":     splitTree(dsNum, subNum, runNum)
        if opt == "-skimLAT":   skimLAT(argv[i+1],argv[i+2],argv[i+3])
        if opt == "-tuneCuts":  tuneCuts(argv[i+1],dsNum)
        if opt == "-lat3":      applyCuts(int(argv[i+1]),argv[i+2])
        if opt == "-cron":      cronJobs()
        if opt == "-shifter":   shifterTest()
        if opt == "-test":      quickTest()
        if opt == "-b":         runBatch()
        if opt == "-chunk":     chunkJobList()

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
            print("+%s: %s" % (jobQueue.split("/")[-1],cmd))
            f.write(cmd + "\n")


def getSBatch(opt, getCores=True, nArr=1):
    """
    http://www.nersc.gov/users/computational-systems/cori/running-jobs/batch-jobs/
    http://www.nersc.gov/users/computational-systems/edison/running-jobs/batch-jobs/
    http://www.nersc.gov/users/computational-systems/pdsf/using-slurm-pdsf-batch-sbatch/
    """
    batchOpts = {
        "chos": [
            "--workdir=%s" % (dsi.latSWDir),
            "--output=%s/logs/chos-%%j.txt" % (dsi.latSWDir),
            "-p shared-chos",
            "-t 24:00:00",
            "--mem=3400",
        ],
        "pdsf-single": [
            "--workdir=%s" % (dsi.latSWDir),
            "--output=%s/logs/pdsf-%%j.txt" % (dsi.latSWDir),
            "--image=wisecg/mjsw:v2",
            "-p shared",
            "-t 24:00:00",
        ],
        "pdsf-pump": [
            "--workdir=%s" % (dsi.latSWDir),
            "--output=%s/logs/pdsf-%%j.txt" % (dsi.latSWDir),
            "--image=wisecg/mjsw:v2",
            "-p shared",
            "-n30", # match nCores for pdsf below
            "-t 24:00:00",
        ],
        "pdsf-test": [
            "--workdir=%s" % (dsi.latSWDir),
            "--output=%s/logs/pdsf-%%j.txt" % (dsi.latSWDir),
            "--image=wisecg/mjsw:v2",
            "-p shared",
            # "-n30", # match nCores for pdsf below
            "-t 00:20:00",
        ],
        "cori": [
            "--workdir=%s" % (dsi.latSWDir),
            "--output=%s/logs/cori-%%j.txt" % (dsi.latSWDir),
            "--image=wisecg/mjsw:v2",
            "-C haswell",
            "-N 1",
            "--qos=debug",
            "-t 00:10:00"
            # "--qos=regular",
            # "-t 24:00:00",
        ],
        "edison": [
            "--workdir=%s" % (dsi.latSWDir),
            "--output=%s/logs/edison-%%j.txt" % (dsi.latSWDir),
            "--image=wisecg/mjsw:v2",
            "-t 24:00:00",
            # "-t 00:10:00",
            # "--qos=debug"
            "--qos=shared"
        ],
        "edison-arr": [
            "--workdir=%s" % (dsi.latSWDir),
            "--output=%s/logs/edison-%%j.txt" % (dsi.latSWDir),
            "--image=wisecg/mjsw:v2",
            # "-t 00:10:00",
            # "--qos=debug"
            " --array=0-%d" % nArr,
            "-t 24:00:00",
            "--qos=regular",
            "-N 1"
        ],
        "pdsf-arr": [
            "--workdir=%s" % (dsi.latSWDir),
            "--output=%s/logs/pdsf-%%j.txt" % (dsi.latSWDir),
            "--image=wisecg/mjsw:v2",
            # "-p shared",
            "-p short",
            "-t 00:20:00",
            "--array=0-%d" % nArr
        ]
    }
    sbStr = "sbatch "+' '.join(batchOpts[opt])

    if getCores:
        if "pdsf" in opt: return sbStr, 30, 33 # sbatch cmd, nCores, peakLoad
        if "cori" in opt: return sbStr, 60, 66
        if "edison" in opt: return sbStr, 48, 52
    else:
        return sbStr


def runBatch():
    """ ./job-panda.py -b
    http://www.nersc.gov/users/computational-systems/cori/running-jobs/queues-and-policies/
    http://www.nersc.gov/users/computational-systems/edison/running-jobs/queues-and-policies/
    http://www.nersc.gov/users/accounts/user-accounts/how-usage-is-charged/

    Did you make sure your job lists have a newline at the end ??
    """
    # EX. 1: Run a single job under chos
    # cmd = "./skim_mjd_data -f 22513 -l -t 0.7 %s/skim" % (dsi.specialDir)
    # sh("%s slurm-job.sh %s" % (getSBatch("chos",False), cmd))

    # EX. 2 Run 'make clean && make' in a new env (pdsf-single, cori, edison)
    # sh("make clean")
    # sh("%s slurm.slr make" % (getSBatch("edison",False),cmd))

    # EX. 3: Single job
    # cmd = "./skim_mjd_data -f 22513 -l -t 0.7 %s/skim" % (dsi.specialDir)
    # cmd = "./skim_mjd_data 1 11 -l -g -d -t 0.7"
    # sh("%s slurm.slr %s" % (getSBatch("pdsf-single",False),cmd))

    # EX. 4: Run a job pump
    # sh("%s slurm.slr './job-pump.sh jobs/test.ls skim_mjd_data %d %d'" % getSBatch("cori"))

    # EX. 5: Run a job array (note: slurm task id matches the integer given by the --array option)
    # sh("%s slurm.slr 'eval ./job-pump.sh jobs/test/test_${SLURM_ARRAY_TASK_ID}.ls python3 %d %d'" % getSBatch("edison-arr",nArr=2))
    # sh("%s slurm.slr 'eval ./job-pump.sh jobs/test/test_${SLURM_ARRAY_TASK_ID}.ls python3 %d %d'" % getSBatch("pdsf-arr",nArr=2))

    # EX. 5: Long calibration
    # sh("%s slurm.slr './job-pump.sh jobs/longSkim.ls skim_mjd_data %d %d'" % getSBatch("edison"))
    # sh("%s slurm.slr './job-pump.sh jobs/longWave.ls wave-skim %d %d'" % getSBatch("edison"))
    # sh("%s slurm.slr './job-pump.sh jobs/longSplit.ls python3 %d %d'" % getSBatch("edison"))
    # sh("%s slurm.slr './job-pump.sh jobs/longLAT.ls python3 %d %d'" % getSBatch("edison"))
    # sh("%s slurm.slr './job-pump.sh jobs/longLAT_2.ls python3 %d %d'" % getSBatch("edison"))

    # EX. 6: External pulser
    # sh("%s slurm.slr './job-pump.sh jobs/extPulserLAT.ls python3 %d %d'" % getSBatch("edison"))

    # EX. 7: Force trigger
    # sh("%s slurm.slr './job-pump.sh jobs/trigSkim.ls skim_mjd_data %d %d'" % getSBatch("pdsf-pump"))
    # sh("%s slurm.slr './job-pump.sh jobs/trigWave.ls wave-skim %d %d'" % getSBatch("pdsf-pump"))
    # sh("%s slurm.slr './job-pump.sh jobs/trigSplit.ls python3 %d %d'" % getSBatch("pdsf-pump"))
    # sh("%s slurm.slr './job-pump.sh jobs/trigLAT.ls python3 %d %d'" % getSBatch("edison"))

    # EX. 8: PROCESS BKG DATA
    # sh("%s slurm.slr './job-pump.sh jobs/bkgSkim.ls skim_mjd_data %d %d'" % getSBatch("pdsf-pump"))
    # sh("%s slurm.slr './job-pump.sh jobs/bkgWave.ls wave-skim %d %d'" % getSBatch("pdsf-pump"))
    # sh("%s slurm.slr './job-pump.sh jobs/bkgWave_ds5c6.ls wave-skim %d %d'" % getSBatch("pdsf-pump"))
    # sh("%s slurm.slr './job-pump.sh jobs/bkgSplit.ls python3 %d %d'" % getSBatch("pdsf-pump"))
    # sh("%s slurm.slr './job-pump.sh jobs/bkgSplit_ds5c6.ls python3 %d %d'" % getSBatch("pdsf-pump"))
    # sh("%s slurm.slr './job-pump.sh jobs/test.ls python3 %d %d'" % getSBatch("pdsf-pump"))
    # sh("%s slurm.slr './job-pump.sh jobs/bkgLAT_ds4.ls python3 %d %d'" % getSBatch("pdsf-pump"))
    # sh("%s slurm.slr 'eval ./job-pump.sh jobs/bkgLAT/bkgLAT_${SLURM_ARRAY_TASK_ID}.ls python3 %d %d'" % getSBatch("edison-arr", nArr=49))
    # sh("%s slurm.slr 'eval ./job-pump.sh jobs/bkgLAT_ds5c6/bkgLAT_${SLURM_ARRAY_TASK_ID}.ls python3 %d %d'" % getSBatch("edison-arr", nArr=14))
    # sh("%s slurm.slr './job-pump.sh jobs/bkgLAT_cleanup.ls python3 %d %d'" % getSBatch("pdsf-pump"))
    # sh("%s slurm.slr 'eval ./job-pump.sh jobs/bkgThresh/bkgThresh_${SLURM_ARRAY_TASK_ID}.ls auto-thresh %d %d'" % getSBatch("edison-arr", nArr=10))
    # sh("%s slurm.slr './job-pump.sh jobs/bkgThresh_cleanup.ls auto-thresh %d %d'" % getSBatch("edison"))

    # EX. 9: PROCESS CAL DATA
    # sh("%s slurm.slr './job-pump.sh jobs/calSkim.ls skim_mjd_data %d %d'" % getSBatch("pdsf-pump"))
    # sh("%s slurm.slr './job-pump.sh jobs/calSkim_ds6.ls skim_mjd_data %d %d'" % getSBatch("pdsf-pump"))
    # sh("%s slurm.slr './job-pump.sh jobs/calWave.ls wave-skim %d %d'" % getSBatch("pdsf-pump"))
    # sh("%s slurm.slr './job-pump.sh jobs/calWave_ds5c.ls wave-skim %d %d'" % getSBatch("pdsf-pump"))
    # sh("%s slurm.slr './job-pump.sh jobs/calWave_ds6.ls wave-skim %d %d'" % getSBatch("pdsf-pump"))
    # sh("%s slurm.slr './job-pump.sh jobs/calSplit.ls python3 %d %d'" % getSBatch("pdsf-pump"))
    # sh("%s slurm.slr './job-pump.sh jobs/calSplit_2.ls python3 %d %d'" % getSBatch("pdsf-pump"))
    # sh("%s slurm.slr './job-pump.sh jobs/calSplit_ds5c.ls python3 %d %d'" % getSBatch("pdsf-pump"))
    # sh("%s slurm.slr 'eval ./job-pump.sh jobs/calLAT/calLAT_${SLURM_ARRAY_TASK_ID}.ls python3 %d %d'" % getSBatch("edison-arr",nArr=99))
    # sh("%s slurm.slr 'eval ./job-pump.sh jobs/calLAT_ds5c/calLAT_${SLURM_ARRAY_TASK_ID}.ls python3 %d %d'" % getSBatch("edison-arr",nArr=11))
    # sh("%s slurm.slr './job-pump.sh jobs/calLAT_cleanup.ls python3 %d %d'" % getSBatch("pdsf-pump"))

    # EX. 10: file integrity checks
    # sh("%s slurm.slr %s" % (getSBatch("pdsf-single",False),"./check-files.py -all"))
    # sh("%s slurm.slr %s" % (getSBatch("pdsf-single",False),"./check-files.py -c -all"))

    # EX. 11: raw threshold, channel, and HV settings
    # sh("%s slurm.slr %s" % (getSBatch("edison",False),"./chan-sel.py -t -v"))
    # sh("%s slurm.slr './job-pump.sh jobs/test.ls python3 %d %d'" % getSBatch("pdsf-test"))


def getCalRunList(dsNum=None,subNum=None,runNum=None):
    """ ./job-panda.py -cal (-ds [dsNum] -sub [dsNum] [calIdx] -run [runNum])
        Create a calibration run list, using the CalInfo object in DataSetInfo.py .
        Note that the -sub option is re-defined here to mean a calibration range idx.
        Note that running with -cal alone will create a list for all datasets (mega mode).
    """
    runLimit = None # need all the stats we can get ...
    calList = []
    calInfo = dsi.CalInfo()
    calKeys = calInfo.GetKeys(dsNum)

    # single-run mode
    if runNum!=None:
        calList.append(runNum)
        print(calList)
        return calList

    # multi-run mode:
    for key in calKeys:
        print("key:",key)
        if key=="ds6": continue # comment this out and do -cal -ds 6 to get the ds6 list

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

    print("Total Runs:",len(calList))

    return calList


def runSkimmer(dsNum, subNum=None, runNum=None, calList=[]):
    """ ./job-panda.py [-q] -skim (-ds dsNum) (-sub dsNum subNum) (-run dsNum subNum) [-cal]
        Submit skim_mjd_data jobs.
    """
    bkg = dsi.BkgInfo()

    # bkg
    if not calList:
        dub = "-d" if dsNum < 6 else ""
        dsMap = bkg.dsMap()
        # -ds
        if subNum==None and runNum==None:
            for i in range(dsMap[dsNum]+1):
                if dsNum==5 and i >= 113: dub = ""
                job = "./skim_mjd_data %d %d -l -g %s -t 0.7 %s" % (dsNum, i, dub, dsi.skimDir)
                if useJobQueue: sh("%s >& ./logs/skim-ds%d-%d.txt" % (job, dsNum, i))
                else: sh("%s '%s'" % (jobStr, job))
        # -sub
        elif runNum==None:
            if dsNum==5 and i >= 113: dub = ""
            job = "./skim_mjd_data %d %d -l -g %s -t 0.7 %s" % (dsNum, subNum, dub, dsi.skimDir)
            if useJobQueue: sh("%s >& ./logs/skim-ds%d-%d.txt" % (job, dsNum, subNum))
            else: sh("%s '%s'" % (jobStr, job))
        # -run
        elif subNum==None:
            dub = "-d" if run < 23959 or run > 6000000 else ""
            job = "./skim_mjd_data -f %d -l -g %s -t 0.7 %s" % (runNum, dub, dsi.skimDir)
            if useJobQueue: sh("%s >& ./logs/skim-ds%d-run%d.txt" % (job, bkg.GetDSNum(run), run))
            else: sh("%s '%s'" % (jobStr, job))
    # cal
    else:
        for run in calList:
            dub = "-d" if run < 23959 or run > 6000000 else ""
            job = "./skim_mjd_data -f %d -l -g %s -t 0.7 %s" % (run, dub, dsi.calSkimDir)
            if useJobQueue: sh("%s >& ./logs/skim-ds%d-run%d.txt" % (job, bkg.GetDSNum(run),run))
            else: sh("%s '%s'" % (jobStr, job))


def runWaveSkim(dsNum, subNum=None, runNum=None, calList=[]):
    """ ./job-panda.py [-q] -wave (-ds dsNum) (-sub dsNum subNum) (-run dsNum subNum) [-cal]
        Submit wave-skim jobs.
    """
    bkg = dsi.BkgInfo()

    # bkg
    if not calList:
        dub = "-d" if dsNum < 6 else ""
        dsMap = bkg.dsMap()
        # -ds
        if subNum==None and runNum==None:
            for i in range(dsMap[dsNum]+1):
                if dsNum==5 and i >= 113: dub = ""
                job = "./wave-skim -n %s -r %d %d -p %s %s" % (dub, dsNum, i, dsi.skimDir, dsi.waveDir)
                if useJobQueue: sh("%s >& ./logs/wave-ds%d-%d.txt" % (job, dsNum, i))
                else: sh("%s '%s'" % (jobStr, job))
        # -sub
        elif runNum==None:
            if dsNum==5 and i >= 113: dub = ""
            job = "./wave-skim -n %s -r %d %d -p %s %s" % (dub, dsNum, subNum, dsi.skimDir, dsi.waveDir)
            if useJobQueue: sh("%s >& ./logs/wave-ds%d-%d.txt" % (job, dsNum, subNum))
            else: sh("%s '%s'" % (jobStr, job))
        # -run
        elif subNum==None:
            dub = "-d" if run < 23959 or run > 6000000 else ""
            job = "./wave-skim -n %s -f %d %d -p %s %s" % (dub, dsNum, runNum, dsi.skimDir, dsi.waveDir)
            if useJobQueue: sh("%s >& ./logs/wave-ds%d-run%d.txt" % (job, dsNum, runNum))
            else: sh("%s '%s'" % (jobStr, job))
    # cal
    else:
        for run in calList:
            dub = "-d" if run < 23959 or run > 6000000 else ""
            job = "./wave-skim -n %s -c -f %d %d -p %s %s" % (dub, bkg.GetDSNum(run), run, dsi.calSkimDir, dsi.calWaveDir)
            if useJobQueue: sh("%s >& ./logs/wave-ds%d-run%d.txt" % (job, bkg.GetDSNum(run), run))
            else: sh("%s '%s'" % (jobStr, job))


def runAutoThresh(ds=None):
    """./job-panda.py [-q] [-ds dsNum] -thresh
    Generates auto-thresh jobs.  Default is to do it for all datasets.
    """
    bkg = dsi.BkgInfo()
    dsList = [ds] if ds is not None else [0,1,2,3,4,"5A","5B","5C",6]
    for ds in dsList:
        dsNum = ds if isinstance(ds, int) else 5
        dub = "-d" if dsNum < 6 else ""
        for sub in bkg.getRanges(ds):
            if dsNum==5 and sub >= 113: dub = ""
            subr = bkg.GetSubRanges(ds,sub)
            if len(subr) > 0:
                for runLo, runHi in subr:
                    job = "./auto-thresh %d %d -s %d %d %s -o %s" % (dsNum, sub, runLo, runHi, dub, dsi.threshDir)
                    if useJobQueue: sh("%s >& ./logs/thresh-ds%d-%d-%d-%d.txt" % (job, dsNum, sub, runLo, runHi))
                    else: sh("%s '%s'" % (jobStr, job))
            else:
                job = "./auto-thresh %d %d %s -o %s" % (dsNum, sub, dub, dsi.threshDir)
                if useJobQueue: sh("%s >& ./logs/thresh-ds%d-%d.txt" % (job, dsNum, sub))
                else: sh("%s '%s'" % (jobStr, job))


def splitTree(dsNum, subNum=None, runNum=None):
    """ ./job-panda.py -split (-sub dsNum subNum) (-run dsNum runNum)

        Split a SINGLE waveSkim file into small (~50MB) files to speed up LAT parallel processing.
        Can call 'batchSplit' instead to submit each run in the list as a job, splitting the files in parallel.
        NOTE: The cut written into the first file is NOT copied into the additional files
              (I couldn't get it to work within this function -- kept getting "file not closed" errors.)
              To clean up, do that with the 'writeCut' function below, potentially AFTER a big parallel job.
    """
    from ROOT import TFile, TTree, gDirectory, TEntryList, TNamed, TObject, gROOT

    print("Splitting tree.  dsNum:",dsNum,"subNum:",subNum,"runNum:",runNum)

    # Set input and output paths.  Clear out any files from a previous
    # try before you attempt a copy (avoid the double underscore)
    inPath, outPath = "", ""
    if runNum==None:
        # bg mode
        inPath = "%s/waveSkimDS%d_%d.root" % (dsi.waveDir,dsNum,subNum)
        outPath = "%s/splitSkimDS%d_%d.root" % (dsi.splitDir,dsNum,subNum)
        fileList = sorted(glob.glob("%s/splitSkimDS%d_%d*.root" % (dsi.splitDir,dsNum, subNum)))
        for f in fileList: os.remove(f)
    elif subNum==None:
        # cal mode
        inPath = "%s/waveSkimDS%d_run%d.root" % (dsi.calWaveDir,dsNum,runNum)
        outPath = "%s/splitSkimDS%d_run%d.root" % (dsi.calSplitDir,dsNum,runNum)
        fileList = sorted(glob.glob("%s/splitSkimDS%d_run%d*.root" % (dsi.calSplitDir,dsNum,runNum)))
        for f in fileList: os.remove(f)

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


def batchSplit(dsNum, subNum=None, runNum=None, calList=[]):
    """ ./job-panda.py [-q] [-cal] -batchSplit (-ds dsNum) (-sub dsNum subNum) (-run dsNum subNum)
        Submit jobs that call splitTree for each run, splitting files into small chunks.
        NOTE: The data cleaning cut is NOT written into the output files and the
              function 'writeCut' must be called after these jobs are done.
    """
    bkg = dsi.BkgInfo()

    # bg
    if not calList:

        dsMap = bkg.dsMap()
        # -ds
        if subNum==None and runNum==None:
            for i in range(dsMap[dsNum]+1):
                inPath = "%s/waveSkimDS%d_%d.root" % (dsi.waveDir,dsNum,i)
                if not os.path.isfile(inPath):
                    print("File",inPath,"not found. Continuing ...")
                    continue
                else:
                    job = "./job-panda.py -sub %d %d -split" % (dsNum, i)
                    if useJobQueue: sh("%s >& ./logs/split-ds%d-%d.txt" % (job, dsNum, i))
                    else: sh("""%s '%s'""" % (jobStr, job))
        # -sub
        elif runNum==None:
            inPath = "%s/waveSkimDS%d_%d.root" % (dsi.waveDir,dsNum,subNum)
            if not os.path.isfile(inPath):
                print("File",inPath,"not found.")
                return
            else:
                job = "./job-panda.py -sub %d %d -split" % (dsNum, i)
                if useJobQueue: sh("%s >& ./logs/split-ds%d-%d.txt" % (job, dsNum, i))
                else: sh("""%s '%s'""" % (jobStr, job))
        # -run
        elif subNum==None:
            inPath = "%s/waveSkimDS%d_run%d.root" % (dsi.waveDir,dsNum,runNum)
            if not os.path.isfile(inPath):
                print("File",inPath,"not found.")
                return
            else:
                job = "./job-panda.py -run %d %d -split" % (dsNum, runNum)
                if useJobQueue: sh("%s >& ./logs/split-ds%d-run%d.txt" % (job, dsNum, runNum))
                else: sh("""%s '%s'""" % (jobStr, job))


    # cal
    else:
        for run in calList:
            dsRanges = bkg.dsRanges()
            for key in dsRanges:
                if dsRanges[key][0] <= run <= dsRanges[key][1]:
                    dsNum=key
            inPath = "%s/waveSkimDS%d_run%d.root" % (dsi.calWaveDir,dsNum,run)
            if not os.path.isfile(inPath):
                print("File",inPath,"not found. Continuing ...")
                continue
            else:
                job = "./job-panda.py -run %d %d -split" % (dsNum, run)
                if useJobQueue: sh("%s >& ./logs/split-ds%d-run%d.txt" % (job, dsNum, run))
                else: sh("""%s '%s'""" % (jobStr, job))


def writeCut(dsNum, subNum=None, runNum=None, calList=[]):
    """ ./job-panda.py -writeCut (-ds dsNum) (-sub dsNum subNum) (-run dsNum subNum) [-cal]
        Assumes the cut used in the FIRST file (even in the whole DS) should be applied
        to ALL files.  This should be a relatively safe assumption.
    """
    from ROOT import TFile, TNamed, TObject
    fileList = []
    bkg = dsi.BkgInfo()

    # bg
    if not calList:
        dsMap = bkg.dsMap()
        # -ds
        if subNum==None and runNum==None:
            for i in range(dsMap[dsNum]+1):
                inPath = "%s/splitSkimDS%d_%d*" % (dsi.splitDir,dsNum,i)
                fileList.extend(sorted(glob.glob(inPath)))
        # -sub
        elif runNum==None:
            inPath = "%s/splitSkimDS%d_%d*" % (dsi.splitDir,dsNum,subNum)
            fileList.extend(sorted(glob.glob(inPath)))
        # -run
        elif subNum==None:
            inPath = "%s/splitSkimDS%d_run%d*" % (dsi.splitDir,dsNum,runNum)
            fileList.extend(sorted(glob.glob(inPath)))
    # cal
    else:
        dsRanges = bkg.dsRanges()
        for run in calList:
            for key in dsRanges:
                if dsRanges[key][0] <= run <= dsRanges[key][1]:
                    dsNum=key
            inPath = "%s/splitSkimDS%d_run%d*" % (dsi.calSplitDir,dsNum,run)
            fileList.extend(sorted(glob.glob(inPath)))

    # Pull the cut off the FIRST file and add it to the sub-files
    if len(fileList) <= 1:
        print("No files found!  Exiting...")
        exit(1)

    firstFile = TFile(fileList[0])
    theCut = firstFile.Get("theCut").GetTitle()
    print("Applying this cut:\n",theCut)
    for f in fileList:
        print(f)
        subRangeFile = TFile(f,"UPDATE")
        thisCut = TNamed("theCut",theCut)
        thisCut.Write("",TObject.kOverwrite)


def runLAT(dsNum, subNum=None, runNum=None, calList=[]):
    """ ./job-panda.py [-q] -lat (-ds dsNum) (-sub dsNum subNum) (-run dsNum subNum) [-cal]
        Runs LAT on splitSkim output.  Does not combine output files back together.
    """
    bkg = dsi.BkgInfo()

    # bg
    if not calList:
        dsMap = bkg.dsMap()
        # -ds
        if subNum==None and runNum==None:
            for subNum in range(dsMap[dsNum]+1):
                files = dsi.getSplitList("%s/splitSkimDS%d_%d*" % (dsi.splitDir,dsNum,subNum),subNum)
                for idx, inFile in sorted(files.items()):
                    outFile = "%s/latSkimDS%d_%d_%d.root" % (dsi.latDir,dsNum,subNum,idx)
                    job = "./lat.py -b -r %d %d -p %s %s" % (dsNum,subNum,inFile,outFile)

                    # jspl = job.split() # make SUPER sure stuff is matched
                    # print(jspl[3],jspl[4],jspl[6].split("/")[-1],jspl[7].split("/")[-1])

                    if useJobQueue: sh("%s >& ./logs/lat-ds%d-%d-%d.txt" % (job, dsNum, subNum, idx))
                    else: sh("""%s '%s'""" % (jobStr, job))
        # -sub
        elif runNum==None:
            files = dsi.getSplitList("%s/splitSkimDS%d_%d*" % (dsi.splitDir,dsNum,subNum),subNum)
            for idx, inFile in sorted(files.items()):
                outFile = "%s/latSkimDS%d_%d_%d.root" % (dsi.latDir,dsNum,subNum,idx)
                job = "./lat.py -b -r %d %d -p %s %s" % (dsNum,subNum,inFile,outFile)
                if useJobQueue: sh("%s >& ./logs/lat-ds%d-run%d-%d.txt" % (job, dsNum, subNum, idx))
                else: sh("""%s '%s'""" % (jobStr, job))
        # -run
        elif subNum==None:
            files = dsi.getSplitList("%s/splitSkimDS%d_run%d*" % (dsi.splitDir,dsNum,runNum),runNum)
            for idx, inFile in sorted(files.items()):
                outFile = "%s/latSkimDS%d_run%d_%d.root" % (dsi.latDir,dsNum,runNum,idx)
                job = "./lat.py -b -r %d %d -p %s %s" % (dsNum,runNum,inFile,outFile)
                if useJobQueue: sh("%s >& ./logs/lat-ds%d-run%d-%d.txt" % (job, dsNum, runNum, idx))
                else: sh("""%s '%s'""" % (jobStr, job))
    # cal
    else:
        dsRanges = bkg.dsRanges()
        for run in calList:
            for key in dsRanges:
                if dsRanges[key][0] <= run <= dsRanges[key][1]:
                    dsNum=key
            files = dsi.getSplitList("%s/splitSkimDS%d_run%d*" % (dsi.calSplitDir,dsNum,run),run)
            for idx, inFile in sorted(files.items()):
                outFile = "%s/latSkimDS%d_run%d_%d.root" % (dsi.calLatDir,dsNum,run,idx)
                job = "./lat.py -b -f %d %d -p %s %s" % (dsNum,run,inFile,outFile)
                if useJobQueue: sh("%s >& ./logs/lat-ds%d-run%d-%d.txt" % (job, dsNum, run, idx))
                else: sh("""%s '%s'""" % (jobStr, job))


def mergeLAT():
    """ It seems like a good idea, right?
        Merging all the LAT files back together after splitting?
        Ehhh, maybe they're too large to handle.
    """
    print("hey")


def pandifySkim(dsNum, subNum=None, runNum=None, calList=[]):
    """ ./job-panda.py -pandify (-ds dsNum) (-sub dsNum subNum) (-run dsNum subNum) [-cal]
        Run ROOTtoPandas jobs.
    """
    # bg
    if not calList:
        bkg = dsi.BkgInfo()
        dsMap = bkg.dsMap()
        # -ds
        if subNum==None and runNum==None:
            for i in range(dsMap[dsNum]+1):
                sh("""%s 'python3 ./sandbox/ROOTtoPandas.py -ws %d %d -p -d %s %s'""" % (jobStr, dsNum, i, dsi.waveDir, dsi.pandaDir))
        # -sub
        elif runNum==None:
            sh("""%s 'python3 ./sandbox/ROOTtoPandas.py -ws %d %d -p -d %s %s'""" % (jobStr, dsNum, subNum, dsi.waveDir, dsi.pandaDir))
        # -run
        elif subNum==None:
            sh("""%s 'python3 ./sandbox/ROOTtoPandas.py -f %d %d -p -d %s %s'""" % (jobStr, dsNum, runNum, dsi.waveDir, dsi.pandaDir))
    # cal
    else:
        for i in calList:
            sh("""%s 'python3 ./sandbox/ROOTtoPandas.py -f %d %d -p -d %s %s'""" % (jobStr, dsNum, i, dsi.calWaveDir, dsi.pandaDir))


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
    calInfo = dsi.CalInfo()
    if dsNum==None:
        for i in dsi.dsMap.keys():
            if i == 6: continue
            for mod in [1,2]:
                try:
                    for j in range(calInfo.GetIdxs("ds%d_m%d"%(i, mod))):
                        print("%s './lat3.py -db -tune %s -s %d %d %d %s" % (jobStr, dsi.calLatDir, i, j, mod, argString))
                        sh("""%s './lat3.py -db -tune %s -s %d %d %d %s '""" % (jobStr, dsi.calLatDir, i, j, mod, argString))
                except: continue
    # -ds
    else:
        for mod in [1,2]:
            try:
                for j in range(calInfo.GetIdxs("ds%d_m%d"%(dsNum, mod))):
                    print("%s './lat3.py -db -tune %s -s %d %d %d %s" % (jobStr, dsi.calLatDir, dsNum, j, mod, argString))
                    sh("""%s './lat3.py -db -tune %s -s %d %d %d %s '""" % (jobStr, dsi.calLatDir, dsNum, j, mod, argString))
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
    cal = dsi.CalInfo()
    # runList = cal.GetSpecialRuns("extPulser")
    # runList = cal.GetSpecialRuns("delayedTrigger")
    # runList = cal.GetSpecialRuns("longCal",5)
    runList = cal.GetSpecialRuns("forcedAcq",8)
    for run in runList:
        if useJobQueue:
            # sh("""./skim_mjd_data -f %d -l -t 0.7 %s/skim >& ./logs/specialSkim-DS%d-%d.txt""" % (run,dsi.specialDir,dsi.GetDSNum(run),run))
            sh("""./skim_mjd_data -f %d -x -l %s/skim >& ./logs/specialSkim-DS%d-%d.txt""" % (run,dsi.specialDir,dsi.GetDSNum(run),run))
        else:
            sh("""%s './skim_mjd_data -f %d -x -l %s/skim'""" % (jobStr,run,dsi.specialDir))


def specialWave():
    """ ./job-panda.py [-q (use queue)] -swave """
    cal = dsi.CalInfo()
    # runList = cal.GetSpecialRuns("extPulser")
    # runList = cal.GetSpecialRuns("longCal",5)
    runList = cal.GetSpecialRuns("forcedAcq",8)
    for run in runList:
        if useJobQueue:
            # sh("""./wave-skim -l -n -f %d %d -p %s/skim %s/waves >& ./logs/wave-ds%d-%d.txt""" % (dsi.GetDSNum(run),run,dsi.specialDir,dsi.specialDir,dsi.GetDSNum(run),run))
            sh("""./wave-skim -x -n -f %d %d -p %s/skim %s/waves >& ./logs/wave-ds%d-%d.txt""" % (dsi.GetDSNum(run),run,dsi.specialDir,dsi.specialDir,dsi.GetDSNum(run),run))
        else:
            sh("""%s './wave-skim -x -n -f %d %d -p %s/skim %s/waves'""" % (jobStr, dsi.GetDSNum(run), run, dsi.specialDir, dsi.specialDir) )


def specialSplit():
    """ ./job-panda.py [-q] -ssplit
    External pulser runs have no data cleaning cut.
    Has a memory leak (can't close both TFiles, damn you, ROOT); submit each run as a batch job.
    """
    cal = dsi.CalInfo()
    # runList = cal.GetSpecialRuns("extPulser")
    # runList = cal.GetSpecialRuns("longCal",5)
    runList = cal.GetSpecialRuns("forcedAcq",8)
    for run in runList:

        inPath = "%s/waves/waveSkimDS%d_run%d.root" % (dsi.specialDir, dsi.GetDSNum(run), run)
        outPath = "%s/split/splitSkimDS%d_run%d.root" % (dsi.specialDir, dsi.GetDSNum(run), run)

        outFiles = glob.glob("%s/split/splitSkimDS%d_run%d*.root" % (dsi.specialDir, dsi.GetDSNum(run), run))
        for filename in outFiles:
            try:
                os.remove(filename)
            except OSError:
                pass

        if useJobQueue:
            sh("""./job-panda.py -splitf %s %s""" % (inPath,outPath))
        else:
            sh("""%s './job-panda.py -splitf %s %s'""" % (jobStr,inPath,outPath))


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
    cal = dsi.CalInfo()
    runList = cal.GetSpecialRuns("longCal",5)

    for run in runList:
        dsNum = dsi.GetDSNum(run)
        wavePath = "%s/waves/waveSkimDS%d_run%d.root" % (dsi.specialDir, dsNum, run)
        waveFile = TFile(wavePath)
        theCut = waveFile.Get("theCut").GetTitle()
        print(wavePath)

        splitFiles = glob.glob("%s/split/splitSkimDS%d_run%d*.root" % (dsi.specialDir, dsNum, run))
        for idx in range(len(splitFiles)):
            if idx==0:
                splitPath = "%s/split/splitSkimDS%d_run%d.root" % (dsi.specialDir, dsNum, run)
            else:
                splitPath = "%s/split/splitSkimDS%d_run%d_%d.root" % (dsi.specialDir, dsNum, run, idx)
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
        outFiles = glob.glob("%s/split/splitSkimDS%d_run%d*.root" % (dsi.specialDir, dsi.GetDSNum(run), run))
        outFiles.extend(["%s/skim/skimDS%d_run%d_low.root" % (dsi.specialDir, dsi.GetDSNum(run), run)])
        outFiles.extend(["%s/waves/waveSkimDS%d_run%d.root" % (dsi.specialDir, dsi.GetDSNum(run), run)])
        outFiles.extend(glob.glob("%s/lat/latSkimDS%d_run%d_*.root" % (dsi.specialDir, dsi.GetDSNum(run), run)))
        for filename in outFiles:
            try:
                os.remove(filename)
                print(filename)
            except OSError:
                pass

    # remove all files from ext pulser range
    # cal = dsi.CalInfo()
    # for idx in [6]:
    #     runList = cal.GetSpecialRuns("extPulser",idx)
    #     for run in runList:
    #         outFiles = glob.glob("%s/split/splitSkimDS%d_run%d*.root" % (dsi.specialDir, dsi.GetDSNum(run), run))
    #         outFiles.extend(["%s/skim/skimDS%d_run%d_low.root" % (dsi.specialDir, dsi.GetDSNum(run), run)])
    #         outFiles.extend(["%s/waves/waveSkimDS%d_run%d.root" % (dsi.specialDir, dsi.GetDSNum(run), run)])
    #         outFiles.extend(glob.glob("%s/lat/latSkimDS%d_run%d_*.root" % (dsi.specialDir, dsi.GetDSNum(run), run)))
    #         for filename in outFiles:
    #             print(filename)
    #             try:
    #                 os.remove(filename)
    #             except OSError:
    #                 pass

    # remove lat files without the _X.root
    # import datetime
    # cal = dsi.CalInfo()
    # runList = cal.GetSpecialRuns("extPulser")
    # for run in runList:
    #     outFile = "%s/lat/latSkimDS%d_run%d.root" % (dsi.specialDir, dsi.GetDSNum(run), run)
    #     try:
    #         modDate = os.path.getmtime(outFile)
    #         modDate = datetime.datetime.fromtimestamp(int(modDate)).strftime('%Y-%m-%d %H:%M:%S')
    #         print(outFile, modDate)
    #         os.remove(outFile)
    #     except OSError:
    #         pass


def specialLAT():
    """ ./job-panda.py [-q (use job queue)] -slat"""
    cal = dsi.CalInfo()
    # runList = cal.GetSpecialRuns("extPulser")
    # runList = cal.GetSpecialRuns("longCal",5)
    runList = cal.GetSpecialRuns("forcedAcq",8)

    # deal with unsplit files
    # run = runList[0]
    # dsNum = dsi.GetDSNum(run)
    # inFile = "%s/waves/waveSkimDS%d_run%d.root" % (dsi.specialDir,dsNum,run)
    # outFile = "%s/lat/latSkimDS%d_run%d.root" % (dsi.specialDir,dsNum,run)
    # sh("""./lat.py -x -b -f %d %d -p %s %s""" % (dsNum,run,inFile,outFile))

    # deal with split files
    for run in runList:

        dsNum = dsi.GetDSNum(run)
        inFiles = glob.glob("%s/split/splitSkimDS%d_run%d*.root" % (dsi.specialDir, dsNum, run))
        for idx in range(len(inFiles)):
            if idx==0:
                inFile = "%s/split/splitSkimDS%d_run%d.root" % (dsi.specialDir, dsNum, run)
            else:
                inFile = "%s/split/splitSkimDS%d_run%d_%d.root" % (dsi.specialDir, dsNum, run, idx)
            if not os.path.isfile(inFile) :
                print("File doesn't exist:",inFile)
                return
            outFile = "%s/lat/latSkimDS%d_run%d_%d.root" % (dsi.specialDir, dsNum, run, idx)

            if useJobQueue:
                # this is what you would want for a normal cron queue
                # sh("""./lat.py -x -b -f %d %d -p %s %s""" % (dsNum, run, inFile, outFile))

                # this is what i need for a 1-node job pump
                # sh("""./lat.py -x -b -f %d %d -p %s %s >& ./logs/extPulser-%d-%d.txt""" % (dsNum, run, inFile, outFile, run, idx))
                # sh("""./lat.py -b -f %d %d -p %s %s >& ./logs/longCal-%d-%d.txt""" % (dsNum, run, inFile, outFile, run, idx))
                sh("""./lat.py -x -b -f %d %d -p %s %s >& ./logs/forceAcq-%d-%d.txt""" % (dsNum, run, inFile, outFile, run, idx))
            else:
                sh("""%s './lat.py -x -b -f %d %d -p %s %s' """ % (jobStr, dsNum, run, inFile, outFile))


def specialCheck():
    """./job-panda.py -scheck
    A next step could be to 'hadd' split files back together, but we'll wait for now.
    """
    from ROOT import TFile, TTree
    cal = dsi.CalInfo()
    runList = cal.GetSpecialRuns("extPulser")
    for run in runList:
        fileList = glob.glob("%s/lat/latSkimDS%d_run%d_*.root" % (dsi.specialDir, dsi.GetDSNum(run), run))
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
    cal = dsi.CalInfo()
    rawDir = "/global/project/projectdirs/majorana/data/mjd/surfmjd/data/raw/P3JDY/Data"
    buildDir = dsi.dataDir + "/mjddatadir"
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


def chunkJobList():
    """ ./job-panda.py -chunk
    Split huge job lists into chunks.
    NOTE:  Edison has 48 cores/node, PDSF has 30, Cori has 60.
    It's probably good to have at least that many jobs in a chunk.
    """
    # with open("%s/jobs/bkgLAT.ls" % dsi.latSWDir) as f:
    # with open("%s/jobs/bkgLAT_3.ls" % dsi.latSWDir) as f:
    # with open("%s/jobs/calLAT.ls" % dsi.latSWDir) as f:
    # with open("%s/jobs/calLAT_ds5c.ls" % dsi.latSWDir) as f:
    with open("%s/jobs/bkgThresh.ls" % dsi.latSWDir) as f:
        jobList = [line.rstrip('\n') for line in f]

    # nChunks = 50 # full lat BG process
    # nChunks = 15 # ds5c & ds6
    # nChunks = 100 # DS0-5c cal process
    # nChunks = 11 # DS5c cal
    nChunks = 10 # auto-thresh calculation
    jobsInChunk = int(len(jobList)/nChunks)
    print("nChunks %d  jobs in chunk %d" % (nChunks, jobsInChunk))

    for iCh in range(nChunks):
        cLo, cHi = iCh*jobsInChunk, (iCh+1) * jobsInChunk - 1
        if iCh == nChunks-1: cHi = len(jobList)-1
        print("cLo, cHi:",cLo, cHi)
        # jobFile = "%s/jobs/bkgLAT/bkgLAT_%d.ls" % (dsi.latSWDir, iCh)
        # jobFile = "%s/jobs/bkgLAT_ds5c6/bkgLAT_%d.ls" % (dsi.latSWDir, iCh)
        # jobFile = "%s/jobs/calLAT/calLAT_%d.ls" % (dsi.latSWDir, iCh)
        # jobFile = "%s/jobs/calLAT_ds5c/calLAT_%d.ls" % (dsi.latSWDir, iCh)
        jobFile = "%s/jobs/bkgThresh/bkgThresh_%d.ls" % (dsi.latSWDir, iCh)
        with open(jobFile,"w") as f:
            print(jobFile)
            for idx, job in enumerate(jobList[cLo:cHi+1]):
                f.write(job+"\n")
                print(idx, job)


def quickTest():
    """./job-panda.py -test """
    from ROOT import MGTWaveform, GATDataSet, TChain
    import datetime, time
    print("Sleeping 10 sec ...")
    time.sleep(10)
    now = datetime.datetime.now()
    print("Done. Date: ",str(now))

    # from ROOT import TFile, TTree
    # cal = dsi.CalInfo()
    # runList = cal.GetSpecialRuns("forcedAcq",8)
    # for run in runList:
    #     dsNum = dsi.GetDSNum(run)
    #     inFiles = glob.glob("%s/split/splitSkimDS%d_run%d*.root" % (dsi.specialDir, dsNum, run))
    #     # for idx in range(len(inFiles)):
    #     for f in inFiles:
    #         print(f)
    #         # tf = TFile()



if __name__ == "__main__":
    main(sys.argv[1:])

