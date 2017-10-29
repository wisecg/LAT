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
import sys, time, os, glob, imp
import numpy as np
import tinydb as db
sys.argv.append("-b") # kill all interactive crap
ds = imp.load_source('DataSetInfo','../DataSetInfo.py')
wl = imp.load_source('waveLibs','../waveLibs.py')
from ROOT import gROOT, TFile, TChain
from scipy.optimize import curve_fit

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
        wl.setDBCalRecord(rec,fFor)

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
    calTable = wl.getDBCalTable(dsNum)
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
    # wl.setDBCalTable()
    # wl.getDBKeys()
    # wl.getDBCalRecord("ds1_idx0")
    # wl.getDBCalRecord("ds1_calIdx")
    # wl.delDBRecord("ds1_idx0")
    # wl.getDBCalTable(5)
    # wl.getDBRunCoverage(1,9999)

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
    """
    import ROOT
    from ROOT import std

    filePath = bgDir + "/latSkimDS%d*.root" % dsNum
    if cal==True:
        filePath = calDir + "/latSkimDS%d*.root" % dsNum

    files = glob.glob(filePath)
    for fileName in files:

        start = time.clock()
        print "Now scanning",fileName

        f = TFile(fileName,"UPDATE")
        tree = f.Get("skimTree")
        nEnt = tree.GetEntries()

        # wipe the wfstd branch if it already exists before creating a new one
        b0 = tree.GetListOfBranches().FindObject("wfstd")
        if isinstance(b0, ROOT.TBranchElement): tree.GetListOfBranches().Remove(b0)

        wfstd = std.vector("double")()
        b1 = tree.Branch("wfstd",wfstd)

        print "entries:", nEnt

        # loop over events
        for iList in range(nEnt):
            print 1
            tree.GetEntry(iList)
            print 2
            nChans = tree.channel.size()
            print 3
            nWFs = tree.MGTWaveforms.size()
            if (nChans != nWFs):
                print "Wrong num entries.  Bailing!"
                exit(1)
            wfstd.assign(nChans,-88888)

            continue

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



if __name__ == "__main__":
    main(sys.argv[1:])