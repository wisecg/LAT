#!/usr/common/usg/software/python/2.7.9/bin/python
#!/usr/local/bin/python
"""
========================= job-panda.py ========================
A cute and adorable way to do various low-energy tasks.

Dependencies:
- skim_mjd_data.cc & DataSetInfo.hh
  (very slightly modified from this revision:
  https://github.com/mppmu/GAT/commit/3bad1ebe7b9d0d8984d3861764f8d34e9955f284)
- wave-skim.cc
- lat.py, waveModel.py, waveLibs.py

=================== C. Wiseman, 2 June 2017 ===================
"""
Usage = """
./job-panda.py [options]
[opt1]: ... eh. TBD.
"""

import sys, shlex, glob, os, re
import subprocess as sp
import DataSetInfo as ds
from textwrap import dedent
from shutil import copyfile

homePath = os.path.expanduser('~')
skimDir = "/global/homes/w/wisecg/project/skim"
# waveDir = "/global/homes/w/wisecg/project/waveskim"
# latDir = "/global/homes/w/wisecg/project/lat"
waveDir = "."
latDir = "."


# =============================================================
def main(argv):

    if len(argv)==0:
        print Usage
        return

    # hey, let's get some margs
    dsNum, subNum, runNum, argString = None, None, None, None
    f = dict.fromkeys(shlex.split('a b c d e f g h i j k l m n'),False) # make a bunch of bools
    for i,opt in enumerate(argv):

        if opt == "-clean": cleanUp()
        if opt == "-purge": purgeLogs()
        if opt == "-makeScript": makeScript()
        if opt == "-checkLogs": checkLogErrors()
        if opt == "-skimLAT": skimLAT(argv[i+1],argv[i+2],argv[i+3])
        if opt == "-getEff": getEff()
        if opt == "-cleanlat" : cleanlat()
        if opt == "-threshCut" : threshCut()
        if opt == "-push" : pushOneRun()

        if opt == "-skim":      f['a'] = True
        if opt == "-wave":      f['b'] = True
        if opt == "-cal":       f['c'] = True
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

        if opt == "-ds": dsNum = int(argv[i+1])
        if opt == "-sub": dsNum, subNum = int(argv[i+1]), int(argv[i+2])
        if opt == "-run": dsNum, runNum = int(argv[i+1]), int(argv[i+2])

        if opt == "-tuneCuts": tuneCuts(argv[i+1],dsNum)
        if opt == "-applyChannelCut" : applyChannelCut(int(argv[i+1]),int(argv[i+2]))
        if opt == "-applyCuts": applyCuts(int(argv[i+1]))
        if opt == "-cleanUpCuts": cleanUpCuts(int(argv[i+1]))

    # -- go running --
    if f['a']: runSkimmer(dsNum, subNum, runNum, cal=f['c'])
    if f['b']: runWaveSkim(dsNum, subNum, runNum, cal=f['c'])
    if f['d']: qsubSplit(dsNum, subNum, runNum)
    if f['e']: writeCut(dsNum, subNum, runNum)
    if f['f']: splitTree(dsNum, subNum, runNum)
    if f['g']: runLAT(dsNum, subNum, runNum)
    # mega modes
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

    # if c: combineLatOutput()
    # TODO: if you're going to submit a c++ job, do a make first?

    # print "Done! Job Panda loves you."
# =============================================================


def makeScript():
    """ Makes a qsub submission script that job-panda can use.
    Only need to do this once.
    Don't forget to change the user-specific things at the top!"""

    print "Generating qsub submission script, 'qsub-job.sh'"

    outFile = open('qsub-job.sh','w+')
    qsub_file_text = """
    #!/bin/bash
    #$ -cwd
    #$ -j y
    #$ -o /global/homes/w/wisecg/skim-clean/logs/
    #$ -P majorana
    source /global/homes/w/wisecg/env/EnvBatch.sh # can also comment this out and run with qsub -V
    cd /global/homes/w/wisecg/skim-clean

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


def sh(cmd):
    sp.call(shlex.split(cmd))

def runSkimmer(dsNum, subNum=None, runNum=None, cal=False):
    """ ./job-panda.py -skim (-ds dsNum) (-sub dsNum subNum) (-run dsNum subNum)
        Run skim_mjd_data.
    """
    # Calibration
    if cal:
        if subNum==None and runNum==None: # -ds
            for i in range(0,len(ds.calRanges[dsNum])):
                numRuns = ds.calRanges[dsNum][i][1] - ds.calRanges[dsNum][i][0]
                if numRuns > 4:
                    for rx in range(ds.calRanges[dsNum][i][0]+numRuns/2-2, ds.calRanges[dsNum][i][0]+numRuns/2+2):
                        print "ds %i  run %i" % (dsNum,rx)
                        sh("""qsub -l h_vmem=2G qsub-job.sh './skim_mjd_data -f %d -n -l -t 0.7 %s'""" % (rx, skimDir))
                else:
                    for rx in range(ds.calRanges[dsNum][i][2], ds.calRanges[dsNum][i][1]+1):
                        print "ds %i  run %i" % (dsNum,rx)
                        sh("""qsub -l h_vmem=2G qsub-job.sh './skim_mjd_data -f %d -n -l -t 0.7 %s'""" % (rx, skimDir))
        elif runNum==None: # -sub
            numRuns = ds.calRanges[dsNum][subNum][1] - ds.calRanges[dsNum][subNum][0]
            if numRuns > 4:
                for rx in range(ds.calRanges[dsNum][subNum][0]+numRuns/2-2, ds.calRanges[dsNum][subNum][0]+numRuns/2+2):
                    print "ds %i  run %i" % (dsNum,rx)
                    sh("""qsub -l h_vmem=2G qsub-job.sh './skim_mjd_data -f %d -n -l -t 0.7 %s'""" % (rx, skimDir))
            else:
                for rx in range(ds.calRanges[dsNum][subNum][2], ds.calRanges[dsNum][i][1]+1):
                    print "ds %i  run %i" % (dsNum,rx)
                    sh("""qsub -l h_vmem=2G qsub-job.sh './skim_mjd_data -f %d -n -l -t 0.7 %s'""" % (rx, skimDir))

    # Background
    else:
        if subNum==None and runNum==None: # -ds
            for i in range(ds.dsMap[dsNum]+1):
                print "ds %i  sub %i" % (dsNum,i)
                sh("""qsub -l h_vmem=2G qsub-job.sh './skim_mjd_data %d %d -n -l -t 0.7 %s'""" % (dsNum, i, skimDir))
        elif runNum==None: # -sub
            print "ds %i  sub %i" % (dsNum,subNum)
            sh("""qsub -l h_vmem=2G qsub-job.sh './skim_mjd_data %d %d -n -l -t 0.7 %s'""" % (dsNum, subNum, skimDir))
        elif subNum==None: # -run
            print "ds %i  run %i" % (dsNum,runNum)
            sh("""qsub -l h_vmem=2G qsub-job.sh './skim_mjd_data -f %d -n -l -t 0.7 %s'""" % (runNum, skimDir))


def runWaveSkim(dsNum, subNum=None, runNum=None, cal=False):
    """ ./job-panda.py -wave (-ds dsNum) (-sub dsNum subNum) (-run dsNum subNum [-cal])
        Run wave-skim.  Optionally set a calibration cut.
    """
    # Calibration
    if cal:
        if subNum==None and runNum==None: # -ds
            for i in range(0,len(ds.calRanges[dsNum])):
                numRuns = ds.calRanges[dsNum][i][1] - ds.calRanges[dsNum][i][0]
                if numRuns > 4:
                    for rx in range(ds.calRanges[dsNum][i][0]+numRuns/2-2, ds.calRanges[dsNum][i][0]+numRuns/2+2):
                        print "ds %i  run %i" % (dsNum,rx)
                        sh("""qsub -l h_vmem=2G qsub-job.sh './wave-skim -n -c -f %d %d -p %s %s'""" % (dsNum, rx, skimDir, waveDir) )
                else:
                    for rx in range(ds.calRanges[dsNum][i][2], ds.calRanges[dsNum][i][1]+1):
                        print "ds %i  run %i" % (dsNum,rx)
                        sh("""qsub -l h_vmem=2G qsub-job.sh './wave-skim -n -c -f %d %d -p %s %s'""" % (dsNum, rx, skimDir, waveDir) )

        elif runNum==None: # -sub
            numRuns = ds.calRanges[dsNum][subNum][1] - ds.calRanges[dsNum][subNum][0]
            if numRuns > 4:
                for rx in range(ds.calRanges[dsNum][subNum][0]+numRuns/2-2, ds.calRanges[dsNum][subNum][0]+numRuns/2+2):
                    print "ds %i  run %i" % (dsNum,rx)
                    sh("""qsub -l h_vmem=2G qsub-job.sh './wave-skim -n -c -f %d %d -p %s %s'""" % (dsNum, rx, skimDir, waveDir) )
            else:
                for rx in range(ds.calRanges[dsNum][subNum][2], ds.calRanges[dsNum][subNum][1]+1):
                    print "ds %i  run %i" % (dsNum,rx)
                    sh("""qsub -l h_vmem=2G qsub-job.sh './wave-skim -n -c -f %d %d -p %s %s'""" % (dsNum, rx, skimDir, waveDir) )        

        elif subNum==None: # -run
            print "ds %i  run %i" % (dsNum,runNum)
            sh("""qsub -l h_vmem=2G qsub-job.sh './wave-skim -n -c -f %d %d -p %s %s""" % (dsNum, runNum, skimDir, waveDir) )
    # Background
    else:
        if subNum==None and runNum==None: # -ds
            for i in range(ds.dsMap[dsNum]+1):
                print "ds %i  sub %i" % (dsNum,i)
                sh("""qsub -l h_vmem=2G qsub-job.sh './wave-skim -n -r %d %d -p %s %s'""" % (dsNum, i, skimDir, waveDir) )
        elif runNum==None: # -sub
            print "ds %i  sub %i" % (dsNum,subNum)
            sh("""qsub -l h_vmem=2G qsub-job.sh './wave-skim -n -r %d %d -p %s %s'""" % (dsNum, subNum, skimDir, waveDir) )
        elif subNum==None: # -run
            print "ds %i  run %i" % (dsNum,runNum)
            sh("""qsub -l h_vmem=2G qsub-job.sh './wave-skim -n -f %d %d -p %s %s""" % (dsNum, runNum, skimDir, waveDir) )


def qsubSplit(dsNum, subNum=None, runNum=None):
    """ ./job-panda.py -qsubSplit (-ds dsNum) (-sub dsNum subNum) (-run dsNum subNum)

    Submit jobs to the cluster that call splitTree for each run.
    """
    if subNum==None and runNum==None: # -ds
        for i in range(ds.dsMap[dsNum]+1):

            inPath = "%s/waveSkimDS%d_%d.root" % (waveDir,dsNum,i)
            fileSize = os.path.getsize(inPath)/1e6 # mb
            if (fileSize < 45):
                copyfile(inPath, "%s/splitSkimDS%d_%d.root" % (waveDir, dsNum, i))
            else:
                print "ds %i  sub %i" % (dsNum,i)
                sh("""qsub -l h_vmem=2G qsub-job.sh './job-panda.py -split -sub %d %d'""" % (dsNum, i))

    elif runNum==None: # -sub

        inPath = "%s/waveSkimDS%d_%d.root" % (waveDir,dsNum,subNum)
        fileSize = os.path.getsize(inPath)/1e6 # mb
        if (fileSize < 45):
            copyfile(inPath, "%s/splitSkimDS%d_%d.root" % (waveDir, dsNum, subNum))
        else:
            print "ds %i  sub %i" % (dsNum,i)
            sh("""qsub -l h_vmem=2G qsub-job.sh './job-panda.py -split -sub %d %d'""" % (dsNum, subNum))

    elif subNum==None: # -run

        inPath = "%s/waveSkimDS%d_%d.root" % (waveDir,dsNum,runNum)
        fileSize = os.path.getsize(inPath)/1e6 # mb
        if (fileSize < 45):
            copyfile(inPath, "%s/splitSkimDS%d_%d.root" % (waveDir, dsNum, runNum))
        else:
            print "ds %i  sub %i" % (dsNum,i)
            sh("""qsub -l h_vmem=2G qsub-job.sh './job-panda.py -split -run %d %d'""" % (dsNum, runNum))


def splitTree(dsNum, subNum=None, runNum=None):
    """ ./job-panda.py -split (-sub dsNum subNum) (-run dsNum subNum)

    NOTE: For a waveSkim input file of 120k entries, and 150MB file size,
    testing of LAT indicates it will do about ~6000 entries/hr.
    So a 50MB waveSkim input file will take 6.66 hrs to process.
    If the input file is 2GB (worst case), then it will be split into
    312 output files.  This will make the logs hard to grep for errors,
    but not impossible.

    NOTE: This is relatively slow (has to loop over entries).  If you call it with qsubSplit,
    each run will be submitted to the grid, and job-panda will split the files for you in parallel.
    So this function is only intended to operate on ONE file at a time.

    NOTE: I learned later (farther down in job-panda) that TTree::CopyTree can also be used to pass a cut.
    I don't know if it would play nice with the file splitting that's going on here though.
    """
    from ROOT import TFile, TTree, gDirectory, TEntryList, TNamed, TObject, gROOT

    # Set input and output paths.  Clear out any files from a previous
    # try before you attempt a copy (avoid the double underscore)
    inPath, outPath = "", ""
    if runNum==None:
        inPath = "%s/waveSkimDS%d_%d.root" % (waveDir,dsNum,subNum)
        outPath = "%s/splitSkimDS%d_%d.root" % (waveDir,dsNum,subNum)

        fileList = getFileList("%s/splitSkimDS%d_%d*.root" % (waveDir,dsNum,subNum),subNum)
        for key in fileList: os.remove(fileList[key])

    elif subNum==None:
        inPath = "%s/waveSkimDS%d_run%d.root" % (waveDir,dsNum,runNum)
        outPath = "%s/splitSkimDS%d_run%d.root" % (waveDir,dsNum,runNum)

        fileList = getFileList("%s/splitSkimDS%d_run%d*.root" % (waveDir,dsNum,subNum),subNum)
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

    # can't get this to work within this function -- get "file not closed" errors.
    # Have to do it separately in the "writeCut" function below.
    # gROOT.Reset()
    # files2 = getFileList("%s/splitSkimDS%d_%d*" % (waveDir,dsNum,subNum),subNum)
    # for idx, inFile in sorted(files2.iteritems()):
    #     print "inFile:",inFile
    #     outFile2 = TFile(inFile,"UPDATE")
    #     thisCut = TNamed("theCut",theCut)
    #     thisCut.Write("",TObject.kOverwrite)


def writeCut(dsNum, subNum=None, runNum=None):
    """ ./job-panda.py -writeCut (-ds dsNum) (-sub dsNum subNum) (-run dsNum subNum)

    Assumes the cut used in the FIRST file (even in the whole DS) should be applied
    to ALL files.  This should be a relatively safe assumption for now (except maybe
    if we decide to split DS5 ...)
    """
    from ROOT import TFile, TNamed, TObject

    mainList = {}

    if subNum==None and runNum==None: # -ds
        for i in range(ds.dsMap[dsNum]+1):
            inPath = "%s/splitSkimDS%d_%d*" % (waveDir,dsNum,i)
            fileList = getFileList(inPath,i,True,dsNum)
            mainList.update(fileList)

    elif runNum==None: # -sub
        inPath = "%s/splitSkimDS%d_%d*" % (waveDir,dsNum,subNum)
        fileList = getFileList(inPath,subNum,True,dsNum)
        mainList.update(fileList)

    elif subNum==None: # -run
        inPath = "%s/splitSkimDS%d_run%d*" % (waveDir,dsNum,runNum)
        fileList = getFileList(inPath,runNum,True,dsNum)
        mainList.update(fileList)

    # Pull the cut off the FIRST file and add it to the sub-files
    if len(mainList) <= 1: return

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


def runLAT(dsNum, subNum=None, runNum=None):
    """ ./job-panda.py -lat (-ds dsNum) (-sub dsNum subNum) (-run dsNum subNum) """

    if subNum==None and runNum==None: # -ds
        for i in range(ds.dsMap[dsNum]+1):
            files = getFileList("%s/splitSkimDS%d_%d*" % (waveDir,dsNum,i),i)
            for idx, inFile in sorted(files.iteritems()):
                outFile = "%s/latSkimDS%d_%d_%d.root" % (latDir,dsNum,i,idx)
                # print """qsub -l h_vmem=2G qsub-job.sh './lat.py -b -r %d %d -p %s %s'""" % (dsNum,i,inFile,outFile)
                sh("""qsub -l h_vmem=2G qsub-job.sh './lat.py -b -r %d %d -p %s %s'""" % (dsNum,i,inFile,outFile))

    elif runNum==None: # -sub
        files = getFileList("%s/splitSkimDS%d_%d*" % (waveDir,dsNum,subNum),subNum)
        for idx, inFile in sorted(files.iteritems()):
            outFile = "%s/latSkimDS%d_%d_%d.root" % (latDir,dsNum,subNum,idx)
            sh("""qsub -l h_vmem=2G qsub-job.sh './lat.py -b -r %d %d -p %s %s'""" % (dsNum,subNum,inFile,outFile))

    elif subNum==None: # -run
        print "hi"
        files = getFileList("%s/splitSkimDS%d_run%d*" % (waveDir,dsNum,runNum),runNum)
        for idx, inFile in sorted(files.iteritems()):
            outFile = "%s/latSkimDS%d_run%d_%d.root" % (latDir,dsNum,runNum,idx)
            sh("""qsub -l h_vmem=2G qsub-job.sh './lat.py -b -r %d %d -p %s %s'""" % (dsNum,runNum,inFile,outFile))


def getFileList(filePathRegexString, subNum, uniqueKey=False, dsNum=None):
    files = {}
    for fl in glob.glob(filePathRegexString):
        int(re.search(r'\d+',fl).group())
        ints = map(int, re.findall(r'\d+',fl))
        if (ints[1]==subNum):
            if (len(ints)==2): ints.append(0)

            if not uniqueKey:
                files[ints[2]] = fl # zero index
            else:
                files["DS%d_%d_%d" % (dsNum,subNum,ints[2])] = fl



    return files


def cleanUp():
    print "Cleaning up ..."
    sh("make clean"); sh("make -s")


def purgeLogs():
    print "Purging logs ..."
    for fl in glob.glob("./logs/*"): os.remove(fl)


def checkLogErrors():
    print "Checking log errors ..."

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

    # make a dummy bash script
    outFile = open('logCheck.sh','w+')
    dummyScript = "#!/bin/bash \n %s \n %s \n %s \n %s \n %s \n" % (c1,c2,c3,c4,c5)
    outFile.write(dummyScript)
    outFile.close()
    sh('chmod a+x logCheck.sh')
    sh('./logCheck.sh')
    os.remove('logCheck.sh')

    # This isn't really complete. but you get the idea.  error checking via bash inside python is kind of a PITA.
    # Maybe it would be better to just have python look at the files directly.


def pushOneRun():
    """ ./job-panda.py -push

    # Idea: Run a single run through the entire skim-clean routine.
    # Right now: must be two stages, because we have to split the wave-skim tree into multiple files to make LAT run fast.
    """
    dsNum = 5
    rlo, rhi = 21970, 21998

    for run in range(rlo,rhi+1):

        jobStr = ""

        # jobStr += "'./skim_mjd_data -n -l -t 0.8 -f %d %s'" % (run,skimDir)
        # sh("""qsub -l h_vmem=2G qsub-job.sh './skim_mjd_data -n -l -t 0.8 -f %d %s' """ % (run,homePath+"/project/cal-skim/"))

        inPath = homePath + "/project/cal-skim/"
        outPath = homePath + "/project/cal-waveskim/"
        jobStr += "'./wave-skim -f %d %d -p %s %s'" % (dsNum,run,inPath,outPath)
        sh(""" qsub -l h_vmem=2G qsub-job.sh '%s'""" % (jobStr))


def printCalFileSizes():
    for dsNum in ds.calMap:
        totalSize = 0
        runList = []

        for i in xrange(0,len(ds.calMap[dsNum]),2):
            lower = ds.calMap[dsNum][i]
            upper = ds.calMap[dsNum][i+1]
            for j in xrange(lower,upper+1):
                runList.append(j)

        for run in runList:
            inPath = "/global/homes/w/wisecg/project/v1-cal-waveskim/waveSkimDS%d_run%d.root" % (dsNum, run)
            inFile = TFile()
            try:
                infile = TFile(inPath)
                inFile = TFile(inPath)
                gatTree = inFile.Get("skimTree")
                nEnt = gatTree.GetEntriesFast()
                fileSize = os.path.getsize(inPath)/1e6 # mb
                totalSize += fileSize
                print "DS%d, run %d: entries %d  size %.0f" % (dsNum,run,nEnt,fileSize)
            except:
                AttributeError
                print "file ain't exist"

        print "Total Size: %.0f\n" % totalSize


def printFileSizes():

    for dsNum, subDS in ds.dsMap.iteritems():
        totalSize = 0

        for subNum in range(subDS+1):
            waveSkimPath = "/global/homes/w/wisecg/project/waveskim/waveSkimDS%d_%d.root" % (dsNum, subNum)
            skimPath = "/global/homes/w/wisecg/project/skim/skimDS%d_%d_low.root" % (dsNum, subNum)

            wFile = TFile(waveSkimPath)
            theCut = wFile.Get("theCut").GetTitle()
            wTree = wFile.Get("skimTree")
            wEnt = wTree.GetEntriesFast()

            sFile = TFile(skimPath)
            sTree = sFile.Get("skimTree")
            # sTree.Draw(">>elist",theCut,"entrylist")
            # elist = gDirectory.Get("elist")
            # sTree.SetEntryList(elist)
            # nList = elist.GetN()
            sEnt = sTree.GetEntriesFast()

            sSize = os.path.getsize(skimPath)/1e6 # mb
            wSize = os.path.getsize(waveSkimPath)/1e6

            if wEnt > 100000:
                print "DS%d-%d:  sEnt %-8d  wEnt %-8d  s/w %-4.1f  sSize %-4d  wSize %-4d" % (dsNum,subNum,sEnt,wEnt,sEnt/wEnt,sSize,wSize)

            totalSize += fileSize

            if nEnt > 100000:
                print "DS%d-%d: entries %d  size %.0f" % (dsNum,subNum,nEnt,fileSize)
                # print "%d %.1f" % (nEnt,fileSize)
                # print "%.1f," % fileSize

        print "Total Size: %.0f\n" % totalSize


def skimLAT(inPath,outPath,thisCut):
    """ EXAMPLE:
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


def tuneCuts(argString,dsNum=None):
    """ ./job-panda.py -tuneCuts '[argString]' -- run over all ds's
        ./job-panda.py -ds [dsNum] -tuneCuts '[argString]' -- just one DS

    Submit a bunch of tuneCuts.py jobs to the queues.
    Make sure to put argString in quotes.
    Options for argString:
        -f (fastMode), -bcMax, -noiseWeight, -bcTime, -tailSlope, -fitSlo
    """
    if dsNum == None:
        for i in range(6):
            sh("""qsub -l h_vmem=2G qsub-job.sh './tuneCuts.py %d %s'""" % (i, argString))
    else:
        sh("""qsub -l h_vmem=2G qsub-job.sh './tuneCuts.py %d %s'""" %(dsNum,argString))


def applyChannelCut(dsNum,ch):
    """ ./job-panda.py -applyChannelCut [dsNum] [ch]

        Creates ROOT files for each channel in a dataset.
    The "mega cut" is saved into the file, though it **appears**
    that it's not necessary to re-apply the cut in a draw command.

    This is painfully slow, so the next function is a wrapper script to submit these jobs.
    """
    from ROOT import TFile, TChain, TTree, TEntryList, TNamed, TObject, gDirectory

    # -- Load this DS --
    inPath = homePath + "/project/lat/latSkimDS%d*.root" % dsNum

    fileList = glob.glob(inPath)
    cutFile = TFile(fileList[0])
    fileCut = cutFile.Get("theCut").GetTitle()

    bigChain = TChain("skimTree")
    bigChain.Add(inPath)

    # -- Get channel info --
    goodRuns,badRuns,goodRunSigmas = ds.GetThreshDicts(dsNum)

    # -- Make the mega TCut for this channel --
    channelCut = " && channel==%d" % ch
    bandTimeCut = " && bandTime-tOffset-1100 < 11000"
    bcMaxCut = " && bcMax < %.2e" % ds.bcMax[dsNum][ch][2] # 90% value
    noiseWt = "(waveS4-waveS1)/bcMax/trapENFCal"
    noiseWtCut = " && %s < %.2e && %s > %.2e" % (noiseWt,ds.noiseWt[dsNum][ch][0],noiseWt,ds.noiseWt[dsNum][ch][1])
    tailSlopeCut = " && pol2 > %.2e && pol2 < %.2e && pol3 > %.2e && pol3 < %.2e" % (ds.pol2[dsNum][ch][1],ds.pol2[dsNum][ch][0],ds.pol3[dsNum][ch][1],ds.pol3[dsNum][ch][0])
    bcTimeCut = " && (bandTime-tOffset-1100)/(matchTime-tOffset) > %.2e" % ds.bcTime[dsNum][ch]
    fitSloCut = " && fitSlo <= %.2f" % (ds.fitSloCutDep[dsNum][ch][2]) # 95% value
    megaCut = fileCut + bandTimeCut + channelCut + bcMaxCut + noiseWtCut + tailSlopeCut + bcTimeCut + fitSloCut

    # -- Apply the threshold run cut as a series of Draw's --
    outPaths = []
    for idx,runRange in enumerate(goodRuns[ch]):

        runCut = " && (run >= %d && run <= %d)" % (runRange[0],runRange[1])
        channelCut = megaCut + runCut

        outPath = homePath + "/project/latskim/latSkimDS%d_ch%d_%d.root" % (dsNum,ch,idx)
        outFile = TFile(outPath,"RECREATE")

        lilTree = TTree()
        # lilTree.SetMaxTreeSize(50000000) # 50 MB - do I want this?
        lilTree = bigChain.CopyTree(channelCut)
        lilTree.Write("",TObject.kOverwrite)
        nWrote = lilTree.GetEntries()
        print ch, idx, runRange, "nWrote:",nWrote,"file:",outPath
        outPaths.append(outPath)
        outFile.Close()

    # The next function, cleanUpCuts, merges the output files and writes the cut.
    # The PDSF nodes can't seem to handle the hadd, but the login nodes do it fine.


def applyCuts(dsNum):
    """ ./job-panda.py -applyCuts [dsNum]
    Submits a bunch of applyChannelCut jobs to the queue.
    """
    for ch in ds.GetGoodChanList(dsNum):
        sh("""qsub -l h_vmem=2G qsub-job.sh './job-panda.py -applyChannelCut %d %d'""" % (dsNum,ch))


def cleanUpCuts(dsNum):
    """ ./job-panda.py -cleanUpCuts [dsNum]
    applyCuts made a big mess on the queues.  Stupid nodes can't handle it.
    """
    import re
    from ROOT import TFile, TChain, TTree, TEntryList, TNamed, TObject, gDirectory

    print "Scanning DS-%d" % dsNum
    chList = ds.GetGoodChanList(dsNum)
    goodRuns,badRuns,goodRunSigmas = ds.GetThreshDicts(dsNum)

    for ch in chList:
        inPath = homePath + "/project/latskim/latSkimDS%d_ch%d*" % (dsNum,ch)
        fileList = glob.glob(inPath)
        nFiles = len(fileList)

        if nFiles == 0:
            print "Merge failed -- No files found for ch.",ch,".  path:",inPath
            continue
        if nFiles == 1:
            thisFile = TFile(fileList[0])
            thisTree = thisFile.Get("skimTree")
            if not isinstance(thisTree, TTree):
                print "Merge failed -- channel",ch,"has a corrupted output TTree."
            else:
                print "Channel",ch,":Found",thisTree.GetEntries(),"entries."
            continue

        # How many files were there supposed to be?
        # Does copyTree still make a file with 0 entries? (I think it does.)
        for idx,runRange in enumerate(goodRuns[ch]):

            filePath = homePath + "/project/latskim/latSkimDS%d_ch%d_%d.root" % (dsNum,ch,idx)
            if not os.path.isfile(filePath):
                print "DS-%d, Ch %d is incomplete, needs re-run." % (dsNum,ch)
                break

        # Check for the hadd output file
        finalPath = homePath + "/project/latskim/latSkimDS%d_ch%d.root" % (dsNum,ch)
        if os.path.isfile(finalPath):
            print "found a broken final file:",finalPath
            os.remove(finalPath) # if it already exists, hadd will mess up.
            print "...removed."

        fileList = glob.glob(inPath)
        nFiles = len(fileList)
        nRanges = len(goodRuns[ch])
        print "Found %d files for ch %d.  (Expected %d files)" % (nFiles,ch,nRanges)
        if (nFiles != nRanges): print "WARNING: DIDN'T GET ALL THE FILES\n\n"

        # Now sort the files into the proper order
        sortList, sortFiles = [], ""
        for f in fileList:
            import re
            nums = re.findall(r'\d+',f)
            if len(nums) < 3:
                print "skipping file",f
                continue
            sortList.append([int(nums[2]),f])
        sortList = sorted(sortList)
        for f in sortList: sortFiles += f[1] + " "

        print "Do you want to delete input files and go on with merge? (n: exit, s: delete files and skip to next channel)"
        value = raw_input()
        if value=='n': return
        if value!='s':
            sh("hadd %s %s" % (finalPath,sortFiles))
            sh("rm %s" % sortFiles)
        if value=='s':
            sh("rm %s" % sortFiles)
            continue

        # Write the mega cut (sans the run cut) into the file
        origFiles = glob.glob(homePath + "/project/lat/latSkimDS%d*.root" % (dsNum))
        cutFile = TFile(origFiles[0])
        fileCut = cutFile.Get("theCut").GetTitle()
        channelCut = " && channel==%d" % ch
        bandTimeCut = " && bandTime-tOffset-1100 < 11000"
        bcMaxCut = " && bcMax < %.2e" % ds.bcMax[dsNum][ch][2] # 90% value
        noiseWt = "(waveS4-waveS1)/bcMax/trapENFCal"
        noiseWtCut = " && %s < %.2e && %s > %.2e" % (noiseWt,ds.noiseWt[dsNum][ch][0],noiseWt,ds.noiseWt[dsNum][ch][1])
        tailSlopeCut = " && pol2 > %.2e && pol2 < %.2e && pol3 > %.2e && pol3 < %.2e" % (ds.pol2[dsNum][ch][1],ds.pol2[dsNum][ch][0],ds.pol3[dsNum][ch][1],ds.pol3[dsNum][ch][0])
        bcTimeCut = " && (bandTime-tOffset-1100)/(matchTime-tOffset) > %.2e" % ds.bcTime[dsNum][ch]
        fitSloCut = " && fitSlo <= %.2f" % (ds.fitSloCutDep[dsNum][ch][2]) # 95% value
        megaCut = fileCut + bandTimeCut + channelCut + bcMaxCut + noiseWtCut + tailSlopeCut + bcTimeCut + fitSloCut

        finalFile = TFile(finalPath,"UPDATE")
        thisCut = TNamed("theCut",megaCut)
        thisCut.Write("",TObject.kOverwrite)
        finalFile.Close()


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
    """ ./job-panda.py -threshCut """

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


if __name__ == "__main__":
    main(sys.argv[1:])












