#!/usr/bin/env python
"""
===================== LAT2.py =====================

Calibrate energy parameters in LAT data.

Modes:
    -cal  : Scans a calibration range and updates parameters in a TinyDB file.
    -upd  : Updates an input file with data from the calibration database file.
    -test : Database stuff.

v1: 07 Aug 2017

========= C. Wiseman (USC), B. Zhu (LANL) =========
"""
import sys, time, os, glob, ROOT, signal
import numpy as np
import tinydb as db
import DataSetInfo as ds
import waveLibs as wl
from ROOT import gROOT, TFile, TChain
from scipy.optimize import curve_fit
from scipy.stats import linregress

import matplotlib
matplotlib.use('Agg') # noninteractive backend
import matplotlib.pyplot as plt

homePath = os.path.expanduser('~')
bgDir = homePath + "/project/bg-lat"
calDir = homePath + "/project/cal-lat"


def main(argv):

    print "========================================"
    print "LAT2 started:",time.strftime('%X %x %Z')
    startT = time.clock()
    gROOT.ProcessLine("gErrorIgnoreLevel = 3001;") # suppress ROOT error messages

    dsNum, subNum, runNum = None,None,None
    fCal, fUpd, fBat, fFor, fPlt, fRaw, fCalPath = 0,0,0,0,0,0,0
    fPaths = [".","."]
    if len(argv)==0: return
    for i,opt in enumerate(argv):
        if opt == "-cal":
            fCal = True
            print "Calibration mode."
        if opt == "-upd":
            fUpd = True
            print "File update mode."
        if opt == "-p":
            fPlt = True
            print "Writing plots mode."
        if opt == "-debug":
            debugFiles()
            stopT = time.clock()
            print "Stopped:",time.strftime('%X %x %Z'),"\nProcess time (min):",(stopT - startT)/60
            exit()
        if opt == "-fitRN":
            fitDBRiseNoise()
            # stopT = time.clock()
            # print "Stopped:",time.strftime('%X %x %Z'),"\nProcess time (min):",(stopT - startT)/60
            exit()
        if opt == "-force":
            fFor = True
            print "Force DB update mode."
        if opt == "-raw":
            fRaw = True
            print "Printing raw peakfinder plots."
        if opt == "-d":
            dsNum = int(argv[i+1])
            print "Scanning DS-%d" % (dsNum)
        if opt == "-c":
            fCalPath = True
            print "pointing to cal path"
        if opt == "-f":
            print "Scanning DS-%d, run %d" % (dsNum, runNum)
        if opt == "-s":
            dsNum, subNum = int(argv[i+1]), int(argv[i+2])
            print "Scanning DS-%d sub-range %d" % (dsNum, subNum)
        if opt == "-p":
            fPaths = [argv[i+1], argv[i+2]]
            print "Manual I/O paths set."
        if opt == "-test":
            testDB()
        if opt == "-b":
            fBat = True
            import matplotlib
            if os.environ.get('DISPLAY','') == '':
                print "No display found. Using non-interactive Agg backend"
                matplotlib.use('Agg')
            print "Batch mode selected."
    global plt
    import matplotlib.pyplot as plt

    if fCal:
        rec = calibrateRuns(dsNum,subNum,fPaths,fBat,fPlt,fRaw)
        ds.setDBRecord(rec,fFor)

    if fUpd:
        updateFile(dsNum, fCalPath)

    stopT = time.clock()
    print "Stopped:",time.strftime('%X %x %Z'),"\nProcess time (min):",(stopT - startT)/60


def calibrateRuns(dsNum,calIdx,fPaths,batMode,saveFig,plotRaw):
    """ ./lat2.py -cal [options]
    Source: http://nucleardata.nuclear.lu.se/toi/radSearch.asp
    Pb212 : 238.6322, Ra224 : 240.9866
    """
    latTree = TChain("skimTree")

    # chain together all available cal runs for this calIdx
    # calTable = ds.getDBCalTable(dsNum) # TODO: this function is obsolete, it should use the CalInfo object
    calTable = None
    calRuns = [ calTable[calIdx][0], calTable[calIdx][1] ]
    inPath = "/projecta/projectdirs/majorana/users/bxyzhu/cal-lat/latSkimDS%d*" % dsNum
    fList = glob.glob(inPath)
    fPass = []
    for f in fList:
        f0 = f.split('/')
        f1 = f0[len(f0)-1]
        f2 = f1[f1.find('run')+3:f.find('.root')-5].split('_')
        runNum = int(f2[0])
        if calRuns[0] <= runNum <= calRuns[1]:
            latTree.Add(f)
            fPass.append(f)
    if latTree.GetNtrees()==0:
        print "No cal runs found for DS-%d, calIdx %d, run range [%d, %d].  Exiting..." % (dsNum, calIdx, calRuns[0], calRuns[1])
        exit(1)
    print "Found %d files for calIdx %d, run range [%d, %d]:" % (latTree.GetNtrees(), calIdx, calRuns[0], calRuns[1])
    for f in fPass: print f

    f1 = TFile(fPass[0])
    theCut = f1.Get("theCut").GetTitle()

    # this is the record we'll return
    key = "ds%d_idx%d" % (dsNum,calIdx)
    calRecord = {}

    fig = plt.figure(figsize=(10,7), facecolor='w')
    # p0 = plt.subplot(111)
    p0 = plt.subplot(221)
    p1 = plt.subplot(222)
    p2 = plt.subplot(223)
    p3 = plt.subplot(224)
    if not batMode: plt.show(block=False)

    # loop over channels
    counter = 0
    cpd = ds.CPD[dsNum]
    for ch in sorted(ds.DetID[dsNum]):

        # skip dumb stuff and get det pos
        if ch%2==1: continue
        if ch in ds.PMon[dsNum]: continue
        if ds.CPD[dsNum][ch] == 1 or ds.CPD[dsNum][ch] == 2: continue
        pos = "C%sP%sD%s" % (str(cpd[ch])[0],str(cpd[ch])[1],str(cpd[ch])[2])

        calRecord[ch] = [-1,-1,-1,-1]

        n = latTree.Draw("trapENF:fitAmp:latAF:latAFC", theCut + " && fitAmp > 50 && fitAmp < 1000 && channel==%d && avse >-1" % ch, "GOFF")
        if n==0:
            print "%s  Chan %-4d - totCts 0" % (pos, ch)
            continue

        if not batMode:
            if counter!=0:
                value = raw_input()
                if value=='q': exit(1)
            counter += 1

        # do fitting and pull parameters.
        fit2 = True
        pk1 = peakFit(latTree.GetV1(), n, ch, dsNum, "trapENF", plotRaw)
        pk2 = peakFit(latTree.GetV2(), n, ch, dsNum, "fitAmp", plotRaw)
        pk3 = peakFit(latTree.GetV3(), n, ch, dsNum, "latAF", plotRaw)
        pk4 = peakFit(latTree.GetV4(), n, ch, dsNum, "latAFC", plotRaw)

        r1 = pk1.GetResults()
        r2 = pk2.GetResults()
        r3 = pk3.GetResults()
        r4 = pk4.GetResults()

        calRecord[ch] = [r1[0],r2[0],r3[0],r4[0]]

        if pk1.fail() or pk2.fail() or pk3.fail() or pk4.fail():
            print "%s  Chan %-4d - totCts %-5d  Fails: trapENF %d  fitAmp %d  latAF %d  latAFC %d" % (pos,ch,n,pk1.fail(),pk2.fail(),pk3.fail(),pk4.fail())
        else:
            print "%s  Chan %-4d - totCts %-5d  modes (%d %d %d %d)  calRecord:[%.3f,%.3f,%.3f,%.3f]" % (pos,ch,n,r1[7],r2[7],r3[7],r4[7],r1[0],r2[0],r3[0],r4[0])

        # plotting
        if batMode and not saveFig: continue
        p0.cla()
        pk1.GetPlot(p0)
        p1.cla()
        pk2.GetPlot(p1)
        p2.cla()
        pk3.GetPlot(p2)
        p3.cla()
        pk4.GetPlot(p3)
        plt.tight_layout(rect=[0,0,1,0.97])
        plt.suptitle("DS-%d  %s (%d)  calIdx %d" % (dsNum,pos,ch,calIdx))
        if not batMode: plt.pause(0.00001)
        if saveFig:
            plt.savefig("./plots/ds%d_idx%d_ch%d.pdf" % (dsNum, calIdx, ch))

    return {"key":key,"vals":calRecord}


class peakFit:
    """ For the fitting routine. """
    def __init__(self, vec, n, ch, dsNum, name, plotRaw=False):
        self.Name = name
        res = self.find238Peak(vec, n, ch, dsNum, plotRaw)
        self.fitFail, self.ch, self.nTot, self.nPkCts, self.chi2, self.rawX, self.rawH, self.rawMu, self.rawThr, self.pks, self.fineX, self.fineH, self.fitFunc, self.fitMu, self.fitFWHM, self.calConst, self.fitMode = res

    def fail(self): return self.fitFail
    def GetResults(self): return [self.calConst, self.fitMu, self.fitFWHM, self.nTot, self.nPkCts, self.chi2, self.rawMu, self.fitMode]
    def GetRawArrays(self): return [self.rawX, self.rawH]
    def GetFitArrays(self): return [self.fineX, self.fineH, self.fitFunc]
    def GetPlot(self, ax=None):
        """ If the fit fails, show the raw histogram. """
        if ax is None:
            ax = plt.gca()
        const, mu, fwhm, nTot, nPk, chi2, rawMu, fitMode = self.GetResults()
        if self.fitFail:
            xvals, hist = self.GetRawArrays()
            ax.margins(x=0.01,y=0.15)
            ax.plot(xvals, hist, color='red', ls='steps-post', label='%s:%dcts\npkThr:%.1f' % (self.Name,nTot,self.rawThr))
            for pk in self.pks:
                ax.axvline(pk,label="pk @ %.1f" % pk)
            ax.legend(loc='best',fontsize=10)
            return ax
        else:
            xvals, hist, fit = self.GetFitArrays()
            ax.margins(x=0.01,y=0.1)
            ax.plot(xvals, hist, color='blue', ls='steps-post', label='%s (%d nPk)' % (self.Name,nPk))
            ax.axvline(rawMu, color='orange', label='muGuess: %.3f' % rawMu)
            ax.axvline(mu, color='magenta', label='muFit: %.3f' % mu)
            ax.plot(xvals, fit, color='red', label='bestfit\nchi2: %.3f\nconst: %.3f\ncalFWHM %.2f keV' % (chi2,const,const*fwhm))
            ax.legend(loc='best',fontsize=10)
            return ax

    def find238Peak(self, vec, nTot, ch, dsNum, plotRaw=False):
        """ Assumes that the 238 peak is the highest energy peak in the spectrum.
            This means you have to feed it a cal file truncated at 250 kev.
            Return list:
            [fitFail,ch,nTot,nPkCts,chi2,rawX,rawH,rawMu,rawThr,pks,fineX,fineH,fitFunc,fitMu,fitFWHM,calConst,fitMode]
        """
        # get det position
        cpd = ds.CPD[dsNum]
        pos = "C%sP%sD%s" % (str(cpd[ch])[0],str(cpd[ch])[1],str(cpd[ch])[2])

        # read in raw data from a tree.GetVX() object
        lst = []
        for i in range(nTot): lst.append(vec[i])
        arr = np.asarray(lst)

        # coarse histogram & peak detection w/ dynamic threshold
        nBins = 70
        rawH, rawX = np.histogram(arr,bins=nBins)
        rawX = rawX[:-1]
        idx = np.where(rawH>10) # eliminate outliers
        avgCts = np.sum(rawH[idx][10:20])/10.
        rawThr = avgCts * 0.5
        if rawThr <= 0:
            print "Bad peakdet rawThreshold: rawThr %.1f  avgCts %d" % (rawThr,avgCts)
            rawThr = 50
        pks, cts = wl.GetPeaks(rawH, rawX, rawThr)
        rawMu = 0
        if len(pks > 0): rawMu = pks[-1]
        if plotRaw or len(pks) < 1:
            return [True,ch,nTot,0,-1,rawX,rawH,rawMu,rawThr,pks,np.zeros(1),np.zeros(1),np.zeros(1),0,0,-1,0]

        # fine histogram.  Assume the last peak in pks is the 238 peak.
        # adjust the first guess into the bin center
        pkGuess = pks[-1] + (rawX[-1]-rawX[0])/float(nBins)/2.
        lo, hi = pks[-1] - 0.02*pks[-1], pks[-1] + 0.03*pks[-1]
        fineH, fineX = np.histogram(arr, bins=100, range=(lo,hi))
        fineX = fineX[:-1]

        # fit step.  Try fitting w/ a couple different models before giving up.
        rawMu = fineX[np.argmax(fineH[1:-2])]
        rawCt = fineH[np.argmax(fineH[1:-2])]

        pars, guess = [],[]

        fitMode = 0
        try:
            # a1, c0, mu, sig, c1, tau, c2, a2, b
            pars, guess = [0.]*10, [rawCt, 1., rawMu, 0.8, 10., 100., 0.1, 0.1*rawCt, 5.]
            pars,_ = curve_fit(wl.peakModel238240, fineX, fineH, p0=guess)
            fitMode = 1
        except RuntimeError: pass
        if fitMode==0:
            try:
                # a1,c0,mu,sig,c1,tau,c2,b
                pars, guess = [0.]*8, [rawCt, 1., rawMu, 0.8, 10., 100., 0.1, 5.]
                pars,_ = curve_fit(wl.peakModel238_2, fineX, fineH, p0=guess)
                fitMode = 2
            except RuntimeError: pass
        if fitMode==0:
            try:
                # a, mu, sig, tau
                pars, guess = [0.]*4, [rawCt, rawMu, 0.8, 2.]
                pars,_ = curve_fit(wl.peakModel238, fineX, fineH, p0=guess)
                fitMode = 3
            except ValueError:
                print "%s  Chan. %-4d - ValueError. ydata or xdata contain nan's." % (pos,ch)
                return [True,ch,nTot,0,-1,rawX,rawH,rawMu,rawThr,pks,np.zeros(1),np.zeros(1),np.zeros(1),0,0,-1,0]
            except RuntimeError:
                print "%s  Chan. %-4d - RuntimeError.  Leastsq minimization failed." % (pos,ch)
                return [True,ch,nTot,0,-1,rawX,rawH,rawMu,rawThr,pks,np.zeros(1),np.zeros(1),np.zeros(1),0,0,-1,0]

        # calculate stuff
        fitFunc = np.zeros(1)
        fitMu, sig, fitFWHM = 0.,0.,0.

        if fitMode==1:
            fitMu, sig, fitFWHM = pars[2], pars[3], 2.3548 * pars[3]
            fitFunc = wl.peakModel238240(fineX, *pars)
        elif fitMode==2:
            fitMu, sig, fitFWHM = pars[2], pars[3], 2.3548 * pars[3]
            fitFunc = wl.peakModel238_2(fineX, *pars)
        elif fitMode==3:
            fitMu, sig, fitFWHM = pars[1], pars[2], pars[2] * 2.3548
            fitFunc = wl.peakModel238(fineX, *pars)

        if fitMu==0: print "fitmu is zero.  mode is ",fitMode
        calConst = 238.6322 / fitMu

        idx = np.where((fineX >= fitMu-3.*sig) & (fineX <= fitMu+3.*sig))
        nPkCts = np.sum(fineH[idx])
        denom = len(fitFunc[idx]) - 4.
        if denom == 0: denom = 1
        chi2 = np.sum( np.square(fineH[idx] - fitFunc[idx]) / fitFunc[idx] ) / denom

        # consistency checks
        # TODO: a more verbose error message?
        fitFail = False
        err = self.Name + " Error: "
        if fitMu < 400: err += "fitMu=%.1f  " % fitMu
        if nPkCts < 100: err += "nPkCts=%d  " % nPkCts
        if chi2 < 0: err += "chi2=%.1f  " % chi2
        if fitMu < 400 or nPkCts < 100 or chi2 < 0:
            print err
            fitFail = True

        return [fitFail,ch,nTot,nPkCts,chi2,rawX,rawH,rawMu,rawThr,pks,fineX,fineH,fitFunc,fitMu,fitFWHM,calConst,fitMode]


def testDB():
    """ ./lat2.py -test
    Do database stuff. """
    # ds.getDBKeys()
    # ds.getDBRecord("ds1_idx0")
    # ds.delDBRecord("ds1_idx0")

    cal = ds.CalInfo()

    # get a cal index for a run
    # key, run = "ds1_m1",10770
    # print cal.GetCalIdx(key,run)

    # generate a list of cal runs for a given index
    # key, idx = "ds3_m1",4
    # print cal.GetCalList(key,idx,runLimit=10)

    # generate cal runs for a given dataset
    # key = "ds2_m1"
    # for idx in range(cal.GetIdxs(key)):
        # print cal.GetCalList(key,idx,runLimit=10)

    # generate all possible cal runs
    # for key in cal.GetKeys():
        # print key
        # for idx in range(cal.GetIdxs(key)):
            # print cal.GetCalList(key,idx,runLimit=10)


def updateFile(dsNum, cal):
    """ ./lat2.py -upd [options]

    Cruel twist of fate: Not gonna use the calibration stuff above rn.
    Am gonna use this to update each dataset's LAT files with Ralph's
    sigma parameter.

    11/16/2017 - This has big problems when used on the PDSF queue.  It randomly corrupts files and
    causes segfaults when the tree is read. For some reason it works fine when run on login nodes.
    Short term fix for wfStd is to do it with the login nodes.
    If we ever implement LAT2 energy calibration, we will have to deal with this issue.
    """
    import ROOT
    from ROOT import std

    filePath = bgDir + "/latSkimDS%d*.root" % dsNum
    if cal==True:
        filePath = calDir + "/latSkimDS%d*.root" % dsNum

    print "Globbing",filePath
    files = glob.glob(filePath)
    # print files

    # return

    # comment this block out for normal running
    # skip ahead to where these jobs died last
    # smallList = []
    # foundLeftOff = False
    # for fileName in sorted(files):
    #     if dsNum==0:
    #         if foundLeftOff:
    #             smallList.append(fileName)
    #         if fileName == calDir + "/latSkimDS0_run3439_0.root":
    #             foundLeftOff = True
    #     elif dsNum==1:
    #         if foundLeftOff:
    #             smallList.append(fileName)
    #         if fileName == calDir + "/latSkimDS1_run12729_2.root":
    #             foundLeftOff = True
    # print len(files), len(smallList)
    # files = smallList

    # run over just one file
    # files = ["/global/homes/w/wisecg/project/cal-lat/latSkimDS0_run4841_0.root"]

    for fileName in sorted(files):

        start = time.clock()
        print "Now scanning",fileName

        f = TFile(fileName,"UPDATE")
        tree = f.Get("skimTree")
        nEnt = tree.GetEntries()

        # wipe the wfstd branch if it already exists before creating a new one
        b0 = tree.GetListOfBranches().FindObject("wfstd")
        if isinstance(b0, ROOT.TBranchElement):
            tree.GetListOfBranches().Remove(b0)
            tree.Write()

        wfstd = std.vector("double")()
        b1 = tree.Branch("wfstd",wfstd)

        # loop over events
        for iList in range(nEnt):
            tree.GetEntry(iList)
            nChans = tree.channel.size()
            nWFs = tree.MGTWaveforms.size()
            if (nChans != nWFs):
                print "Wrong num entries.  Bailing!"
                exit(1)
            wfstd.assign(nChans,-88888)

            # loop over hits
            for iH in range(nWFs):
                run = tree.run
                chan = tree.channel.at(iH)
                energy = tree.trapENFCal.at(iH)
                wf = tree.MGTWaveforms.at(iH)
                signal = wl.processWaveform(wf,0,0)
                waveRaw = signal.GetWaveRaw()
                waveTS = signal.GetTS()
                # print "%d / %d  Run %d  nCh %d  chan %d  trapENF %.1f" % (iList,nList,run,nChans,chan,energy)

                wfstd[iH] = np.std(waveRaw[5:-5])
                maxAdc = max(waveRaw[5:-5])
                minAdc = min(waveRaw[5:-5])
                nBins = int(maxAdc-minAdc)

            # End loop over hits, fill branches
            b1.Fill()

        # End loop over events
        tree.Write("",ROOT.TObject.kOverwrite)
        print " - Tree entries: %d  Branch entries: %d  Time (min): %.2f" % (tree.GetEntries(),b1.GetEntries(),(time.clock()-start)/60.)

        f.Close()


def debugFiles():
    """ ./lat2.py -debug

    When we tried running this on the SLURM queue in Oct. 2017, it caused a bunch of
    files to segfault.  We've repaired the damage now.
    This works by adjusting the "glob" string to scan a group of files, and attempts
    to access the wfStd and trapENFCalC branches.  If it doesn't segfault, then
    the file is OK for reading.
    """
    # try to catch ROOT segfaults, but this just hangs forever
    # def sig_handler(signum, frame):
    #     print "segfault"
    #     return None
    # signal.signal(signal.SIGSEGV, sig_handler)

    startHere = False
    # fList = glob.glob("/global/homes/w/wisecg/project/cal-lat/latSkimDS*")
    fList = glob.glob("/global/homes/w/wisecg/project/bg-lat/latSkimDS*")
    fList = sorted(fList)

    # fList = ["/global/homes/w/wisecg/project/bg-lat/latSkimDS3_3_4.root"] # check one file
    for idx, fName in enumerate(fList):
        # if fName == "/global/homes/w/wisecg/project/cal-lat/latSkimDS1_run12729_2.root":
            # startHere = True
        # if not startHere: continue

        # be careful with this!
        # print "going to delete", fName
        # os.remove(fName)

        # standard segfault trap
        print "Current: ",fName
        # if idx < len(fList)-1: print "   (next):",fList[idx+1]
        f = TFile(fName)
        t = f.Get("skimTree")
        b1 = t.GetBranch("trapENFCalC")
        b2 = t.GetBranch("wfstd")
        print t.GetEntries(), b1.GetEntries(), b2.GetEntries()
        if b1.GetEntries() != b2.GetEntries(): print "WTF"
        t.GetEntry(0) # this catches it
        print t.trapENFCalC.at(0), t.wfstd.at(0)
        f.Close()

    # f = TFile("/global/homes/w/wisecg/project/bg-lat/latSkimDS3_3_4.root","UPDATE")
    # t = f.Get("skimTree")

    # wipe the wfstd branch
    # b0 = t.GetListOfBranches().FindObject("wfstd")
    # b0 = t.GetBranch("wfstd")
    # print type(b0)
    # if isinstance(b0, ROOT.TBranchElement):

    # t.GetListOfBranches().Remove(b0)
    # t.Write()

    # ROOT.gInterpreter.ProcessLine("gDebug=2") # this prints out a crazy amount of output

    # b1 = t.GetBranch("trapENFCalC")
    # b2 = t.GetBranch("wfstd")
    # print t.GetEntries(), b1.GetEntries() #, b2.GetEntries()
    # if b1.GetEntries() != b2.GetEntries(): print "WTF"
    # t.GetEntry(0)
    # print t.trapENFCalC.at(0) #, t.wfstd.at(0)
    # f.Close()


def fitDBRiseNoise():
    """ ./lat2.py -fitRN

    For each dataset, for every channel, for every calIdx:
    pull the 95% riseNoise value for 50-90, 90-130, 130-170, 170-210.
    Then do a simple linear fit and save the channel value.
    """

    dsNum = 1

    # parse database
    # { recordName : { chan : 95% value } }  - recordName is idxN-e50-90, etc

    dbDict = {}
    calDB = db.TinyDB('calDB.json')
    for item in calDB:
        d = dict(item)
        key = d["key"]
        vals = d["vals"]
        tmp = key.split("_")
        tmp = [str(t) for t in tmp]

        # ex - ['riseNoise', 'ds1', 'idx12', 'm1', '130', '170']
        if tmp[0]=="riseNoise" and tmp[1]=="ds%d" % dsNum and tmp[4]!="Peak" and tmp[4]!="Continuum":
            calIdx = int(tmp[2][3:])
            eLo, eHi = int(tmp[4]), int(tmp[5])
            recName = "idx%s-%d-%d" % (calIdx, eLo, eHi)
            recDict = {}
            for ch in vals:
                chan = int(ch)
                val95 = vals[ch][3] # ex - [1%, 5%, 90%, 95%, 99%]
                recDict[chan] = val95
            dbDict[recName] = recDict


    dsNum, module = 1, 1
    if dsNum == 4: module = 2

    nIdx = ds.getNCalIdxs(dsNum, module)

    # for plots
    # { chan : { calIdx : [vals] } }
    # plotDict = {}

    # for ch in ds.GetGoodChanList(dsNum):
    ch = 578

    calIdx = 0

    chEne = [70,110,150,190]

    chVals = []
    chVals.append(dbDict["idx%d-%d-%d" % (calIdx,50,90)][ch])
    chVals.append(dbDict["idx%d-%d-%d" % (calIdx,90,130)][ch])
    chVals.append(dbDict["idx%d-%d-%d" % (calIdx,130,170)][ch])
    chVals.append(dbDict["idx%d-%d-%d" % (calIdx,170,210)][ch])

    chEne, chVals = np.array(chEne), np.array(chVals)

    slope, yint, r_val, p_val, std_err = linregress(chEne, chVals)

    # plt.figure()
    # plt.plot(chEne, chVals, 'o', label='riseNoise vals')
    # plt.plot(chEne, yint + slope*chEne, 'r', label='fit')
    # plt.legend(loc="best")
    # plt.xlabel("Energy (keV)")
    # plt.ylabel("riseNoise 95% val")
    # plt.show()
    # plt.savefig("./plots/linFit.png")

    print "ch",ch,"vals:",chVals,slope,yint,r_val,p_val,std_err

    # build a db record
    # key: [Name]_ds[i]_idx[j]_module[k]_[descriptor]
    # vals: {[chan] : [slope, yint]} --> need a value for every channel.
    # fk.  the linear fit is going to overcut. i don't even want to use it.

    rnKey = "riseNoise_ds%d_idx%d_m%d_idx%d-ch%d"

    return


if __name__ == "__main__":
    main(sys.argv[1:])