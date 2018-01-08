#!/usr/bin/env python
"""
========================= job-panda.py ========================
A cute and adorable way to do various low-energy processing tasks.

The functions are arranged mostly sequentially, i.e. this file
documents the "procedure" necessary to produce LAT data,
starting from running skim_mjd_data.

Dependencies:
- skim_mjd_data.cc & DataSetInfo.hh
- wave-skim.cc
- lat.py, waveModel.py, waveLibs.py

=================== C. Wiseman, 2 June 2017 ===================
"""
import sys, shlex, glob, os, re, time
import subprocess as sp
import DataSetInfo as ds

home    = os.path.expanduser('~')
skimDir    = home+"/project/bg-skim"
waveDir    = home+"/project/bg-waves" # split files are in waveDir/split
latDir     = home+"/project/bg-lat"
calSkimDir = home+"/project/cal-skim"
calWaveDir = home+"/project/cal-waves"
calLatDir  = home+"/project/cal-lat"
pandaDir   = home+"/project/panda-skim"

# qsubStr = "qsub -l h_vmem=2G qsub-job.sh" # SGE mode
qsubStr = "sbatch slurm-job.sh" # SLURM mode
# qsubStr = "sbatch shifter.slr" # SLURM+Shifter mode

cronFile = home + "/lat/cron.queue"


# =============================================================
def main(argv):

    if len(argv)==0:
        print "Read the main function, broh!"
        return

    global useJobQueue
    useJobQueue = False

    # get some margs
    dsNum, subNum, runNum, argString = None, None, None, None

    f = dict.fromkeys(shlex.split('a b c d e f g h i j k l m n p q r s t'),False) # make a bunch of bools
    for i,opt in enumerate(argv):

        if opt == "-q": useJobQueue = True

        if opt == "-ds": dsNum = int(argv[i+1])
        if opt == "-sub": dsNum, subNum = int(argv[i+1]), int(argv[i+2])
        if opt == "-run": dsNum, runNum = int(argv[i+1]), int(argv[i+2])

        if opt == "-skim":      f['a'] = True
        if opt == "-wave":      f['b'] = True
        if opt == "-qsubSplit": f['d'] = True
        if opt == "-writeCut":  f['e'] = True
        if opt == "-split":     f['f'] = True
        if opt == "-lat":       f['g'] = True
        if opt == "-merge":     f['h'] = True

        if opt == "-megaLAT":   f['i'] = True
        if opt == "-megaSkim":  f['j'] = True
        if opt == "-megaWave":  f['k'] = True
        if opt == "-megaCut":   f['l'] = True
        if opt == "-megaSplit": f['m'] = True
        if opt == "-pandify":   f['p'] = True

        if opt == "-cal":       f['c'] = True
        if opt == "-lat2":      f['n'] = True
        if opt == "-force":     f['q'] = True
        if opt == "-up2":       f['r'] = True
        if opt == "-c":         f['s'] = True

        # one-offs
        if opt == "-purge": purgeLogs()
        if opt == "-makeScript": makeScript()
        if opt == "-makeSlurm": makeSlurm()
        if opt == "-makeShifter": makeShifter()
        if opt == "-checkLogs": checkLogErrors()
        if opt == "-checkLogs2": checkLogErrors2()
        if opt == "-checkFiles": checkFiles()
        if opt == "-skimLAT": skimLAT(argv[i+1],argv[i+2],argv[i+3])
        if opt == "-getEff": getEff()
        if opt == "-threshCut" : threshCut()
        if opt == "-push" : pushOneRun()
        if opt == "-tuneCuts": tuneCuts(argv[i+1],dsNum)
        if opt == "-applyChannelCut" : applyChannelCut(int(argv[i+1]),int(argv[i+2]))
        if opt == "-applyCuts": applyCuts(int(argv[i+1]))
        if opt == "-cleanUpCuts": cleanUpCuts(int(argv[i+1]))
        if opt == "-lat3": lat3ApplyCuts(int(argv[i+1]),argv[i+2])
        if opt == "-cron": cronJobs()
        if opt == "-shifter": shifterTest()


    # -- calibration stuff --
    calList = []
    if f['c']: calList = getCalRunList(dsNum,subNum,runNum)
    if f['n']: runLAT2Cal(dsNum, subNum, f['q'])
    if f['r']: updateLAT2(dsNum, f['s'])

    # -- go running --
    if f['a']: runSkimmer(dsNum, subNum, runNum, calList=calList)
    if f['b']: runWaveSkim(dsNum, subNum, runNum, calList=calList)
    if f['d']: qsubSplit(dsNum, subNum, runNum, calList=calList)
    if f['e']: writeCut(dsNum, subNum, runNum, calList=calList)
    if f['f']: splitTree(dsNum, subNum, runNum)
    if f['g']: runLAT(dsNum, subNum, runNum, calList=calList)
    if f['p']: pandifySkim(dsNum, subNum, runNum, calList=calList)

    # -- mega modes --
    if f['i']:
        for i in range(0,5+1): runLAT(i)
    if f['j']:
        for i in range(0,5+1): runSkimmer(i)
    if f['k']:
        for i in range(0,5+1): runWaveSkim(i)
    if f['m']:
        for i in range(0,5+1): qsubSplit(i)
    if f['l']:
        for i in range(0,5+1): writeCut(i)

    # print "Done! Job Panda loves you."
# =============================================================

def sh(cmd):
    """ Either call the shell or put the command in our cron queue.
        Uses the global bool 'useJobQueue' and the global path 'cronFile'.
    """
    if not useJobQueue:
        sp.call(shlex.split(cmd))
        return
    with open(cronFile,"a+") as f:
        print "Adding to cronfile (%s): %s" % (cronFile, cmd)
        f.write(cmd + "\n")


def makeSlurm():
    """ ./job-panda.py -makeSlurm
        Makes a SLURM submission script.
        Only need to do this once. Don't forget to change the user-specific things at the top!
    """
    from textwrap import dedent
    print "Generating SLURM submission script, slurm-job.sh"

    outFile = open('slurm-job.sh','w+')
    slurm_file_text = """
    #!/bin/bash -l
    #SBATCH -t 24:00:00  --ntasks=1
    #SBATCH --mem 3400
    #SBATCH --account=majorana
    #SBATCH --workdir=/global/homes/w/wisecg/lat
    #SBATCH --output=/global/homes/w/wisecg/lat/logs/slurm-%j.txt

    # THIS FILE WAS AUTOGENERATED BY JOB-PANDA.

    echo "Job Start:"
    date
    echo "Node(s):  "$SLURM_JOB_NODELIST
    echo "Job ID:  "$SLURM_JOB_ID

    # This runs whatever commands job-panda.py passes to it.
    echo "${@}"
    ${@}

    echo "Job Complete:"
    date
    """
    # remove tabs and first empty line
    script = dedent(slurm_file_text).split("\n",1)[1]
    outFile.write(script)
    outFile.close()


def makeScript():
    """ ./job-panda.py -makeScript
        DEPRECATED -- use makeSlurm!
        Makes a qsub submission script that job-panda can use.
        Only need to do this once. Don't forget to change the user-specific things at the top!
    """
    from textwrap import dedent
    print "Generating qsub submission script, 'qsub-job.sh'"

    outFile = open('qsub-job.sh','w+')
    qsub_file_text = """
    #!/bin/bash
    #$ -cwd
    #$ -j y
    #$ -o /global/homes/w/wisecg/lat/logs/
    #$ -P majorana
    source /global/homes/w/wisecg/env/EnvBatch.sh # can also comment this out and run with qsub -V
    cd /global/homes/w/wisecg/lat

    # THIS FILE WAS AUTOGENERATED BY JOB-PANDA.

    echo "Job Start:"
    date
    echo "Node:  "$HOSTNAME
    echo "Job ID:  "$JOB_ID

    # This runs whatever commands job-panda.py passes to it.
    echo "${@}"
    ${@}

    echo "Job Complete:"
    date
    """
    # remove tabs and first empty line
    script = dedent(qsub_file_text).split("\n",1)[1]
    outFile.write(script)
    outFile.close()


def makeShifter():
    """ ./job-panda.py -makeShifter
        Makes SLURM+Shifter submission scripts (it needs two)
        Only need to run this function once.
    """
    from textwrap import dedent
    print "Generating SLURM+shifter scripts ..."

    outFile = open('shifter.slr','w+')
    shifter_file_text = """
    #!/bin/bash
    #SBATCH --workdir=/global/homes/w/wisecg/lat
    #SBATCH --output=/global/homes/w/wisecg/lat/logs/shifter-%j.txt
    #SBATCH -p shared --image=custom:pdsf-chos-sl64:v4
    shifter --volume=/global/project:/project /bin/bash shifter-job.sh ${@}
    """
    script = dedent(shifter_file_text).split("\n",1)[1]
    outFile.write(script)
    outFile.close()

    outFile = open('shifter-job.sh','w+')
    job_file_text = """
    #!/bin/bash
    echo "Job Start:"
    date
    echo "Node(s):  "$SLURM_JOB_NODELIST
    echo "Job ID:  "$SLURM_JOB_ID
    echo inShifter:`env|grep  SHIFTER_RUNTIME`
    echo $CHOS
    echo `pwd`
    echo $SHELL
    echo $MJSWDIR
    echo "homedir is:"$HOMEDIR

    # This runs whatever commands we pass to it.
    echo "${@}"
    ${@}

    echo "Job Complete:"
    date
    """
    script = dedent(job_file_text).split("\n",1)[1]
    outFile.write(script)
    outFile.close()


def purgeLogs():
    print "Purging logs ..."
    for fl in glob.glob("./logs/*"): os.remove(fl)


def getCalRunList(dsNum=None,subNum=None,runNum=None):
    """ ./job-panda.py -cal (-ds [dsNum] -sub [dsNum] [calIdx] -run [runNum])
        Create a calibration run list, using the CalInfo object in DataSetInfo.py .
        Note that the -sub option is re-defined here to mean a calibration range idx.
        Note that running with -cal alone will create a list for all datasets (mega mode).
    """
    runLimit = 10 # yeah I'm hardcoding this, sue me.
    calList = []
    calInfo = ds.CalInfo()
    calKeys = calInfo.GetKeys(dsNum)

    # single-run mode
    if runNum!=None:
        calList.append(runNum)
        print calList
        return calList

    # multi-run mode:
    for key in calKeys:
        print "key:",key

        # -cal (mega mode)
        if dsNum==None:
            for idx in range(calInfo.GetIdxs(key)):
                lst = calInfo.GetCalList(key,idx,runLimit)
                print lst
                calList += lst
        # -ds
        elif subNum==None:
            for idx in range(calInfo.GetIdxs(key)):
                lst = calInfo.GetCalList(key,idx,runLimit)
                print lst
                calList += lst
        # -sub
        else:
            lst = calInfo.GetCalList(key,subNum,runLimit)
            if lst==None: continue
            print lst
            calList += lst

    # remove any duplicates, but there probably aren't any
    calList = sorted(list(set(calList)))

    return calList


def runSkimmer(dsNum, subNum=None, runNum=None, calList=[]):
    """ ./job-panda.py -skim (-ds dsNum) (-sub dsNum subNum) (-run dsNum subNum) [-cal]
        Submit skim_mjd_data jobs.
    """
    # bg
    if not calList:
        # -ds
        if subNum==None and runNum==None:
            for i in range(ds.dsMap[dsNum]+1):
                sh("""%s './skim_mjd_data %d %d -n -l -t 0.7 %s'""" % (qsubStr, dsNum, i, skimDir))
        # -sub
        elif runNum==None:
            sh("""%s './skim_mjd_data %d %d -n -l -t 0.7 %s'""" % (qsubStr, dsNum, subNum, skimDir))
        # -run
        elif subNum==None:
            sh("""%s './skim_mjd_data -f %d -n -l -t 0.7 %s'""" % (qsubStr, runNum, skimDir))
    # cal
    else:
        for run in calList:
            sh("""%s './skim_mjd_data -f %d -n -l -t 0.7 %s'""" % (qsubStr, run, calSkimDir))


def runWaveSkim(dsNum, subNum=None, runNum=None, calList=[]):
    """ ./job-panda.py -wave (-ds dsNum) (-sub dsNum subNum) (-run dsNum subNum) [-cal]
        Submit wave-skim jobs.
    """
    # bg
    if not calList:
        # -ds
        if subNum==None and runNum==None:
            for i in range(ds.dsMap[dsNum]+1):
                sh("""%s './wave-skim -n -r %d %d -p %s %s'""" % (qsubStr, dsNum, i, skimDir, waveDir) )
        # -sub
        elif runNum==None:
            sh("""%s './wave-skim -n -r %d %d -p %s %s'""" % (qsubStr, dsNum, subNum, skimDir, waveDir) )
        # -run
        elif subNum==None:
            sh("""%s './wave-skim -n -f %d %d -p %s %s""" % (qsubStr, dsNum, runNum, skimDir, waveDir) )
    # cal
    else:
        for run in calList:
            for key in ds.dsRanges:
                if ds.dsRanges[key][0] <= run <= ds.dsRanges[key][1]:
                    dsNum=key
            sh("""%s './wave-skim -n -c -f %d %d -p %s %s'""" % (qsubStr, dsNum, run, calSkimDir, calWaveDir) )


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
        inPath = "%s/waveSkimDS%d_%d.root" % (waveDir,dsNum,subNum)
        outPath = "%s/split/splitSkimDS%d_%d.root" % (waveDir,dsNum,subNum)
        fileList = getFileList("%s/split/splitSkimDS%d_%d*.root" % (waveDir,dsNum,subNum),subNum)
        for key in fileList: os.remove(fileList[key])
    elif subNum==None:
        # cal mode
        inPath = "%s/waveSkimDS%d_run%d.root" % (calWaveDir,dsNum,runNum)
        outPath = "%s/split/splitSkimDS%d_run%d.root" % (calWaveDir,dsNum,runNum)
        fileList = getFileList("%s/split/splitSkimDS%d_run%d*.root" % (calWaveDir,dsNum,runNum),runNum)
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
                inPath = "%s/waveSkimDS%d_%d.root" % (waveDir,dsNum,i)
                if not os.path.isfile(inPath):
                    print "File",inPath,"not found. Continuing ..."
                    continue
                if (os.path.getsize(inPath)/1e6 < 45):
                    copyfile(inPath, "%s/split/splitSkimDS%d_%d.root" % (waveDir, dsNum, i))
                else:
                    sh("""%s './job-panda.py -split -sub %d %d'""" % (qsubStr, dsNum, i))
        # -sub
        elif runNum==None:
            inPath = "%s/waveSkimDS%d_%d.root" % (waveDir,dsNum,subNum)
            if not os.path.isfile(inPath):
                print "File",inPath,"not found."
                return
            if (os.path.getsize(inPath)/1e6 < 45):
                copyfile(inPath, "%s/split/splitSkimDS%d_%d.root" % (waveDir, dsNum, subNum))
            else:
                sh("""%s './job-panda.py -split -sub %d %d'""" % (qsubStr, dsNum, subNum))
        # -run
        elif subNum==None:
            inPath = "%s/waveSkimDS%d_run%d.root" % (waveDir,dsNum,runNum)
            if not os.path.isfile(inPath):
                print "File",inPath,"not found."
                return
            if (os.path.getsize(inPath)/1e6 < 45):
                copyfile(inPath, "%s/split/splitSkimDS%d_%d.root" % (waveDir, dsNum, runNum))
            else:
                sh("""%s './job-panda.py -split -run %d %d'""" % (qsubStr, dsNum, runNum))
    # cal
    else:
        for run in calList:
            for key in ds.dsRanges:
                if ds.dsRanges[key][0] <= run <= ds.dsRanges[key][1]:
                    dsNum=key
            inPath = "%s/waveSkimDS%d_run%d.root" % (calWaveDir,dsNum,run)
            if not os.path.isfile(inPath):
                print "File",inPath,"not found. Continuing ..."
                continue
            if (os.path.getsize(inPath)/1e6 < 45):
                copyfile(inPath, "%s/split/splitSkimDS%d_run%d.root" % (calWaveDir, dsNum, run))
            else:
                sh("""%s './job-panda.py -split -run %d %d'""" % (qsubStr, dsNum, run))


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
                inPath = "%s/split/splitSkimDS%d_%d*" % (waveDir,dsNum,i)
                fileList = getFileList(inPath,i,True,dsNum)
                mainList.update(fileList)
        # -sub
        elif runNum==None:
            inPath = "%s/split/splitSkimDS%d_%d*" % (waveDir,dsNum,subNum)
            fileList = getFileList(inPath,subNum,True,dsNum)
            mainList.update(fileList)
        # -run
        elif subNum==None:
            inPath = "%s/split/splitSkimDS%d_run%d*" % (waveDir,dsNum,runNum)
            fileList = getFileList(inPath,runNum,True,dsNum)
            mainList.update(fileList)
    # cal
    else:
        for run in calList:
            for key in ds.dsRanges:
                if ds.dsRanges[key][0] <= run <= ds.dsRanges[key][1]:
                    dsNum=key
            inPath = "%s/split/splitSkimDS%d_run%d*" % (calWaveDir,dsNum,run)
            fileList = getFileList(inPath,run,True,dsNum)
            mainList.update(fileList)

    # Pull the cut off the FIRST file and add it to the sub-files
    if len(mainList) <= 1:
        print "No files found!  Exiting..."
        exit(1)
    theCut = ""
    foundFirst = False
    for key, inFile in sorted(mainList.iteritems()):
        if not foundFirst:
            firstFile = TFile(mainList[key])
            theCut = firstFile.Get("theCut").GetTitle()
            print "Applying this cut:\n",theCut
            foundFirst = True
        print key, inFile
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
                files = getFileList("%s/split/splitSkimDS%d_%d*" % (waveDir,dsNum,i),i)
                for idx, inFile in sorted(files.iteritems()):
                    outFile = "%s/latSkimDS%d_%d_%d.root" % (latDir,dsNum,i,idx)
                    sh("""%s './lat.py -b -r %d %d -p %s %s'""" % (qsubStr,dsNum,i,inFile,outFile))
        # -sub
        elif runNum==None:
            files = getFileList("%s/split/splitSkimDS%d_%d*" % (waveDir,dsNum,subNum),subNum)
            for idx, inFile in sorted(files.iteritems()):
                outFile = "%s/latSkimDS%d_%d_%d.root" % (latDir,dsNum,subNum,idx)
                sh("""%s './lat.py -b -r %d %d -p %s %s'""" % (qsubStr,dsNum,subNum,inFile,outFile))
        # -run
        elif subNum==None:
            files = getFileList("%s/split/splitSkimDS%d_run%d*" % (waveDir,dsNum,runNum),runNum)
            for idx, inFile in sorted(files.iteritems()):
                outFile = "%s/latSkimDS%d_run%d_%d.root" % (latDir,dsNum,runNum,idx)
                sh("""%s './lat.py -b -f %d %d -p %s %s'""" % (qsubStr,dsNum,runNum,inFile,outFile))
    # cal
    else:
        for run in calList:
            for key in ds.dsRanges:
                if ds.dsRanges[key][0] <= run <= ds.dsRanges[key][1]:
                    dsNum=key
            files = getFileList("%s/split/splitSkimDS%d_run%d*" % (calWaveDir,dsNum,run),run)
            for idx, inFile in sorted(files.iteritems()):
                outFile = "%s/latSkimDS%d_run%d_%d.root" % (calLatDir,dsNum,run,idx)
                sh("""%s './lat.py -b -f %d %d -p %s %s'""" % (qsubStr,dsNum,run,inFile,outFile))


def mergeLAT():
    """ It seems like a good idea, right?
        Merging all the LAT files back together after splitting?
    """
    print "hey"


def checkLogErrors():
    """ ./job-panda.py -checkLogs
        This isn't really complete. but you get the idea.  Error checking via bash inside python is kind of a PITA.
        Maybe it would be better to just have python look at the files directly.
    """

    # Shell commands
    c1 = "ls -F ./logs/ | grep -v / | wc -l" # count total files
    c2 = "grep -rIl \"Done! Job Panda\" ./logs/ | wc -l" # count completed files
    c3 = "grep -rL \"Done! Job Panda\" ./logs/"  # negative file matching, can also count # fails.  gives a file list
    c4 = "grep -rIl \"Segmentation\" ./logs/" # segfaults
    c5 = "grep -rIl \"bad_alloc\" ./logs/" # memory errors

    # using sp to deal with a pipe is kind of annoying
    p1 = sp.Popen('ls -F ./logs/'.split(), stdout=sp.PIPE)
    p2 = sp.Popen('grep -v /'.split(), stdin=p1.stdout, stdout=sp.PIPE)
    p3 = sp.Popen('wc -l'.split(), stdin=p2.stdout,stdout=sp.PIPE)
    output = p3.communicate()[0]
    num = int(output.strip('\n'))
    print num

    # make a dummy bash script that runs all the shell commands.  who knows if this is smart or not
    outFile = open('logCheck.sh','w+')
    dummyScript = "#!/bin/bash \n %s \n %s \n %s \n %s \n %s \n" % (c1,c2,c3,c4,c5)
    outFile.write(dummyScript)
    outFile.close()
    sh('chmod a+x logCheck.sh')
    sh('./logCheck.sh')
    os.remove('logCheck.sh')


def checkLogErrors2():
    """ Usage: ./job-panda -checkLogs2
        Globs together log files and then searches for "Error", returning the failed ROOT files.
    """
    print "Checking log errors ..."

    ErrList = []
    for fl in glob.glob("./logs/*"):
        fErr = open(fl,'r').read()
        if 'Error' in open(fl, 'r').read():
            print ErrList.append(fl)

    for errFile in ErrList:
        fErr = open(errFile,'r')
        for lineErr in fErr:
            if '/lat.py -b' in lineErr:
                print 'Error from: ', lineErr
            if 'Error' in lineErr:
                print lineErr


def checkFiles():
    """ ./job-panda.py -checkFiles """
    from ROOT import TFile, TTree
    import os.path, imp
    import DataSetInfo as ds

    # Check BG skim and waveskim files
    # dsMap = {0:75,1:51,2:7,3:24,4:18,5:112}
    dsMap = {2:7}
    for dsNum in dsMap:
        for sub in range(dsMap[dsNum]+1):

            # check skims
            fileName = "/global/homes/w/wisecg/project/bg-skim/skimDS%d_%d_low.root" % (dsNum,sub)
            if not os.path.isfile(fileName):
                print "file not found! name:", fileName
                continue
            f1 = TFile(fileName)
            t1 = f1.Get("skimTree")
            n1 = t1.GetEntriesFast()
            print "DS %d  sub %d  skim entries %d" % (dsNum, sub, n1)
            if n1==0:
                print "no skim entries found! file:", fileName
                continue

            # check waveskims
            fileName = "/global/homes/w/wisecg/project/bg-waves/waveSkimDS%d_%d.root" % (dsNum,sub)
            if not os.path.isfile(fileName):
                print "file not found! name:", fileName
                continue
            f2 = TFile(fileName)
            t2 = f2.Get("skimTree")
            n2 = t2.GetEntriesFast()
            print "DS %d  sub %d  wave entries %d" % (dsNum, sub, n2)
            if n2==0:
                print "no waveskim entries found! file:", fileName
                continue

    # Check CAL skim and waveskim files
    calList = getCalRunList(dsNum=2) # None checks all ds's
    for run in calList:
        dsNum=-1
        for key in ds.dsRanges:
            if ds.dsRanges[key][0] <= run <= ds.dsRanges[key][1]:
                dsNum=key

        # check skims
        fileName = "/global/homes/w/wisecg/project/cal-skim/skimDS%d_run%d_low.root" % (dsNum,run)
        if not os.path.isfile(fileName):
            print "file not found! name:", fileName
            continue
        f1 = TFile(fileName)
        t1 = f1.Get("skimTree")
        n1 = t1.GetEntriesFast()
        print "DS %d  run %d  skim entries %d" % (dsNum, run, n1)
        if n1==0:
            print "no skim entries found! file:", fileName
            continue

        # check waveskims
        fileName = "/global/homes/w/wisecg/project/cal-waves/waveSkimDS%d_run%d.root" % (dsNum,run)
        if not os.path.isfile(fileName):
            print "file not found! name:", fileName
            continue
        f2 = TFile(fileName)
        t2 = f2.Get("skimTree")
        n2 = t2.GetEntriesFast()
        print "DS %d  run %d  wave entries %d" % (dsNum, run, n2)
        if n2==0:
            print "no waveskim entries found! file:", fileName
            continue


def updateLAT2(dsNum, cal=False):
    """ ./job-panda.py -up2 -ds [dsNum] (-c) """

    if not cal:
        print "Submitting DS%d BKG" % (dsNum)
        sh("""%s './lat2.py -upd -d %d' """ % (qsubStr,dsNum))
    else:
        print "Submitting DS%d CAL" % (dsNum)
        sh("""%s './lat2.py -upd -d %d -c' """ % (qsubStr,dsNum))


def runLAT2Cal(dsNum, calIdx=None, forceUpdate=False):
    """ ./job-panda.py -cal (-ds [dsNum]) or (-sub [dsNum] [calIdx])  -force (optional)
        Run LAT2 in cal mode.
    """
    import waveLibs as wl
    if forceUpdate: print "Force updating DB entries."

    # get calIdx's for this dataset from the DB
    calTable = wl.getDBCalTable(dsNum)
    calIdxs = calTable.keys()

    # -ds
    if calIdx==None:
        for idx in calIdxs:
            print "======================================"
            rec = wl.getDBRecord( "ds%d_idx%d" % (dsNum, idx) )

            if rec==None or forceUpdate:
                if forceUpdate:
                    # sh("""./lat2.py -cal -b -p -s %d %d -force""" % (dsNum, idx) )
                    sh("""qsub -l h_vmem=2G qsub-job.sh './lat2.py -cal -b -p -s %d %d -force'""" % (dsNum, idx))
                else:
                    # sh("""./lat2.py -cal -b -p -s %d %d""" % (dsNum, idx) )
                    sh("""qsub -l h_vmem=2G qsub-job.sh './lat2.py -cal -b -p -s %d %d'""" % (dsNum, idx))
    # -sub
    else:
        if calIdx not in calIdxs:
            print "calIdx %d doesn't exist for DS-%d.  Exiting ..."
            return

        rec = wl.getDBRecord( "ds%d_idx%d" % (dsNum, calIdx) )
        if rec==None or forceUpdate:
            if forceUpdate:
                # sh("""./lat2.py -cal -b -p -s %d %d -force""" % (dsNum, calIdx) )
                sh("""qsub -l h_vmem=2G qsub-job.sh './lat2.py -cal -b -p -s %d %d -force'""" % (dsNum, calIdx))
            else:
                # sh("""./lat2.py -cal -b -p -s %d %d""" % (dsNum, calIdx) )
                sh("""qsub -l h_vmem=2G qsub-job.sh './lat2.py -cal -b -p -s %d %d'""" % (dsNum, calIdx))


def skimLAT(inPath,outPath,thisCut):
    """ ./job-panda.py -skimLat [inDir] [outFile] [custom TCut]
        Chains together a given set of LAT files and makes an output file w/ a custom TCut.
        This is handy for quick processing of a particular interesting subset,
        but NOT intended to do any final data production.
        EXAMPLE:
            ./job-panda.py -skimLAT '/Users/wisecg/project/lat/latSkimDS1*.root' "ds1NoisyRuns.root" " && (run==12736||run==12767||run==13005)"
        Note that you can't use a tilde at the beginning.
    """
    from ROOT import TFile, TChain, TTree, TEntryList, TNamed, TObject, gDirectory

    fileList = glob.glob(inPath)
    cutFile = TFile(fileList[0])
    fileCut = cutFile.Get("theCut").GetTitle()
    theCut = fileCut + thisCut
    print "Skimming %d files from %s, with this cut:\n%s" % (len(fileList),inPath,theCut)

    bigChain = TChain("skimTree")
    bigChain.Add(inPath)
    bigChain.Draw(">>elist",theCut,"entrylist")
    elist = gDirectory.Get("elist")
    bigChain.SetEntryList(elist)

    outFile = TFile(outPath,"RECREATE")
    lilTree = TTree()
    # lilTree.SetMaxTreeSize(50000000) # 50 MB
    lilTree = bigChain.CopyTree("") # this does NOT write the cut into the extra files
    lilTree.Write("",TObject.kOverwrite)
    thisCut = TNamed("theCut",theCut)
    thisCut.Write("",TObject.kOverwrite)
    print "Wrote",lilTree.GetEntries(),"entries to the cut tree."

    outFile.Close()


def pushOneRun():
    """ ./job-panda.py -push
    IDEA:  can you make one node do everything necessary to LAT-process 1 run?
    """
    # dsNum = 5
    # rlo, rhi = 21970, 21998
    # for run in range(rlo,rhi+1):
    #     jobStr = ""
    #     # jobStr += "'./skim_mjd_data -n -l -t 0.8 -f %d %s'" % (run,skimDir)
    #     # sh("""qsub -l h_vmem=2G qsub-job.sh './skim_mjd_data -n -l -t 0.8 -f %d %s' """ % (run,homePath+"/project/cal-skim/"))
    #     inPath = homePath + "/project/cal-skim/"
    #     outPath = homePath + "/project/cal-waveskim/"
    #     jobStr += "'./wave-skim -f %d %d -p %s %s'" % (dsNum,run,inPath,outPath)
    #     sh(""" qsub -l h_vmem=2G qsub-job.sh '%s'""" % (jobStr))


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
                        print "%s './lat3.py -db -tune %s -s %d %d %d %s" % (qsubStr, calLatDir, i, j, mod, argString)
                        sh("""%s './lat3.py -db -tune %s -s %d %d %d %s '""" % (qsubStr, calLatDir, i, j, mod, argString))
                except: continue
    # -ds
    else:
        for mod in [1,2]:
            try:
                for j in range(calInfo.GetIdxs("ds%d_m%d"%(dsNum, mod))):
                    print "%s './lat3.py -db -tune %s -s %d %d %d %s" % (qsubStr, calLatDir, dsNum, j, mod, argString)
                    sh("""%s './lat3.py -db -tune %s -s %d %d %d %s '""" % (qsubStr, calLatDir, dsNum, j, mod, argString))
            except: continue


def lat3ApplyCuts(dsNum, cutType):
    """ ./job-panda.py -lat3 [dsNum] [cutType]"""

    if dsNum==-1:
        for ds in range(6):
            sh("""%s './lat3.py -cut %d %s'""" % (qsubStr, ds, cutType))
    else:
        sh("""%s './lat3.py -cut %d %s'""" % (qsubStr, dsNum, cutType))


def getEff():
    """ ./job-panda.py -getEff

    METHOD:
    open up the latskim file for each channel.
    loop over the good run ranges.
    for each good range, make an energy histogram.
    then calculate the efficiency curve based on the sigma value
    and convolve it with the histogram points.
    """
    import numpy as np
    import waveLibs as wl
    import scipy.special as spec
    import matplotlib.pyplot as plt
    from ROOT import TFile, TTree, TH1D, TF1, TCanvas, gROOT
    import ROOT, random

    gROOT.ProcessLine(".x ~/env/MJDClintPlotStyle.C")
    # gROOT.ProcessLine("gErrorIgnoreLevel = 3001;") # suppress ROOT messages

    bins, xlo, xhi = 50,0,15  # set it just high enough that the first bin center isn't negative

    hSumCorr = TH1D("hSumCorr","hSumCorr",bins,xlo,xhi)
    hSumUncr = TH1D("hSumUncr","hSumUncr",bins,xlo,xhi)

    dsNum = 1
    # ch = 578
    for ch in ds.GetGoodChanList(dsNum):

        inFile = TFile(homePath+"/project/latskim/latSkimDS%d_ch%d.root" % (dsNum,ch))
        tree = inFile.Get("skimTree")
        fileCut = inFile.Get("theCut").GetTitle()

        _,_,goodRunErfs = ds.GetThreshDicts(dsNum)

        hUnc = wl.H1D(tree,bins,xlo,xhi,"trapENFCal",fileCut)
        hSumUncr.Add(hUnc)

        for erfs in goodRunErfs[ch]:

            runCut = " && run >= %d && run <= %d" % (erfs[0],erfs[1])
            theCut = fileCut + runCut

            h1 = wl.H1D(tree,bins,xlo,xhi,"trapENFCal",theCut)
            h1x,h1y = wl.npTH1D(h1)

            # calculate efficiency curve
            thisErf = TF1("thisErf","0.5*(1+TMath::Erf((x-[0])/(TMath::Sqrt(2)*[1]) ))")
            thisErf.SetParameter(0,erfs[2]) # mu
            thisErf.SetParameter(1,erfs[3]) # sigma
            h1.Divide(thisErf)

            # thisErf = 0.5 * (1 + spec.erf( (h1x - mu) / (np.sqrt(2) * sig) ))
            # h1yScaled = h1y / thisErf
            # nameStr = str(random.uniform(1.,2.))
            # h2 = TH1D(nameStr,nameStr,bins,xlo,xhi)
            # for i in range(bins):
            #     h2.SetBinContent(i,h1yScaled[i])

            hSumCorr.Add(h1)


    # eff-corrected spectrum.
    c = TCanvas("c","c",800,600)
    c.SetLogy(1)

    hSumCorr.SetLineColor(ROOT.kBlue)
    hSumCorr.Draw("hist")

    hSumUncr.SetLineColor(ROOT.kRed)
    hSumUncr.Draw("hist same")

    # l1 = TLegend

    c.Print("./plots/effWeight/eff_DS%d.pdf" % dsNum)


def threshCut():
    """ ./job-panda.py -threshCut
        Applies

    """
    import numpy as np
    import pandas as pd
    import os

    threshCut = 0.9

    for dsNum in range(0,6):

        if dsNum==2: continue

        print "dataset",dsNum
        df = pd.read_hdf("./data/ThreshDS%d_Processed.h5" % dsNum, 'threshTree')

        goodRuns, badRuns, goodRunErfs = {}, {}, {}
        for chan in df:
            col = int(chan)
            goodRuns[col], badRuns[col], goodRunErfs[col] = [], [], []

            for idx, vals in enumerate(df.loc[:,chan]):

                if np.isnan(vals).any(): continue # skip NaN run ranges where data wasn't collected for the channel

                thr, sigma, hi, lo = vals[0], vals[1], int(vals[2]), int(vals[3])

                if vals[0] <= threshCut:
                    goodRuns[col].append([hi,lo])
                    goodRunErfs[col].append([hi,lo,thr,sigma])
                else:
                    badRuns[col].append([hi,lo])


        # Now make a list of channels to cut in each bad run.

        runList = []
        for ch in badRuns:
            for pair in badRuns[ch]:
                for run in range(pair[0],pair[1]+1):
                    runList.append(run)
        runSet = set(runList)
        runList = sorted(list(runSet))

        runChanPairs = {}
        for run in runList:

            if run not in runChanPairs:
                runChanPairs[run] = []

            for ch in badRuns:
                for pair in badRuns[ch]:
                    if run >= pair[0] and run <= pair[1]:
                        runChanPairs[run].append(ch)
                        # print run,ch

        text_file = open("threshCut_v1.txt", "a")
        for run in sorted(runChanPairs):
            sarr = [str(a) for a in runChanPairs[run]]
            text_file.write(str(run) + ' ' + (' '.join(sarr)) + "\n")
        text_file.close()


def pandifySkim(dsNum, subNum=None, runNum=None, calList=[]):
    """ ./job-panda.py -pandify (-ds dsNum) (-sub dsNum subNum) (-run dsNum subNum) [-cal]
        Run ROOTtoPandas jobs.
    """
    # bg
    if not calList:
        # -ds
        if subNum==None and runNum==None:
            for i in range(ds.dsMap[dsNum]+1):
                sh("""%s './ROOTtoPandas.py -ws %d %d -p -d %s %s'""" % (qsubStr, dsNum, i, waveDir, pandaDir))
        # -sub
        elif runNum==None:
            sh("""%s './ROOTtoPandas.py -ws %d %d -p -d %s %s'""" % (qsubStr, dsNum, subNum, waveDir, pandaDir))
        # -run
        elif subNum==None:
            sh("""%s './ROOTtoPandas.py -f %d %d -p -d %s %s'""" % (qsubStr, dsNum, runNum, waveDir, pandaDir))
    # cal
    else:
        for i in calList:
            sh("""%s './ROOTtoPandas.py -f %d %d -p -d %s %s'""" % (qsubStr, dsNum, i, calWaveDir, pandaDir))


def cronJobs():
    """ ./job-panda.py -cron
    Uses the global string 'cronFile'.
    Crontab should contain the following lines (crontab -e):
    SHELL=/bin/bash
    MAILTO="" # can put in some address here if you LOVE emails
    #*/10 * * * * source ~/env/EnvBatch.sh; ~/lat/job-panda.py -cron >> ~/lat/cron.log 2>&1
    """
    os.chdir(home+"/lat/")
    print "Cron:",time.strftime('%X %x %Z'),"cwd:",os.getcwd()

    nMaxRun, nMaxPend = 15, 200

    with open(cronFile) as f:
        jobList = [line.rstrip('\n') for line in f]
    nList = len(jobList)

    status = os.popen('slusers | grep wisecg').read()
    status = status.split()
    nRun = int(status[0]) if len(status) > 0 else 0  # Rjob Rcpu Rcpu*h PDjob PDcpu user:account:partition
    nPend = int(status[3]) if len(status) > 0 else 0

    nSubmit = (nMaxRun-nRun) if nRun < nMaxRun else 0
    nSubmit = nList if nList < nSubmit else nSubmit
    nSubmit = 0 if nPend >= nMaxPend else nSubmit

    print "   nRun %d  (max %d)  nPend %d (max %d)  nList %d  nSubmit %d" % (nRun,nMaxRun,nPend,nMaxPend,nList,nSubmit)

    with open(cronFile, 'w') as f:
        for idx, job in enumerate(jobList):
            if idx < nSubmit:
                print "Submitted:",job
                sh(job)
            else:
                # print "Waiting:", job
                f.write(job + "\n")


def shifterTest():
    """ ./job-panda.py -shifter """

    # os.chdir(home+"/lat/cron")
    print "Shifter:",time.strftime('%X %x %Z'),"cwd:",os.getcwd()
    sh("""sbatch shifter.slr './lat3.py -cut 2 fs'""")


if __name__ == "__main__":
    main(sys.argv[1:])
