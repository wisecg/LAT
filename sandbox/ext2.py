#!/usr/bin/env python3
import sys, os, imp, glob
sys.argv.append("-b")
ds = imp.load_source('DataSetInfo',os.environ['LATDIR']+'/DataSetInfo.py')
wl = imp.load_source('waveLibs',os.environ['LATDIR']+'/waveLibs.py')
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np
import tinydb as db
from scipy.optimize import curve_fit
calInfo = ds.CalInfo()

def main():
    # riseTime()
    # getEff()
    # compareRunTypes()
    manualBuild()
    # getEffMultiChan()


def riseTime():
    """ fitSlo vs. rise time study """
    from ROOT import TChain

    extPulserInfo = calInfo.GetSpecialList()["extPulserInfo"]
    syncChan = wl.getChan(0,10,0) # 672

    rtfsSlope, rtfsInt = [], []

    pIdxs = [13,16,17,18]
    for pIdx in pIdxs:

        runList = calInfo.GetSpecialRuns("extPulser",pIdx)
        attList = extPulserInfo[pIdx][0]
        extChan = extPulserInfo[pIdx][-1]

        print("pIdx",pIdx,"Chan",extChan)

        fig = plt.figure(figsize=(8,8),facecolor='w')
        cmap = plt.cm.get_cmap('hsv', len(runList)+1)
        p1 = plt.subplot(211)
        p2 = plt.subplot(212)

        sloVals, rtVals, attVals = [], [], []
        for i, run in enumerate(runList):

            fileList = ds.getLATRunList([run],"%s/lat" % (ds.specialDir))
            latChain = TChain("skimTree")
            for f in fileList:
                latChain.Add("%s/lat/%s" % (ds.specialDir,f))

            tNames = ["Entry$","mH","channel","trapENFCal","fitSlo","den90","den10"]
            theCut = "(channel==%d || channel==%d) && mH==2" % (syncChan, extChan) # enforce correct sync
            theCut += " && abs(fitSlo) < 2000 && abs(den90-den10) < 1000 && fitSlo > 60" # make nice plots
            tVals = wl.GetVX(latChain,tNames,theCut)
            nPass = len(tVals["Entry$"])

            enfArr = np.asarray([tVals["trapENFCal"][i] for i in range(nPass) if tVals["channel"][i]==extChan])
            sloArr = np.asarray([tVals["fitSlo"][i] for i in range(nPass) if tVals["channel"][i]==extChan])
            rtArr  = np.asarray([tVals["den90"][i] - tVals["den10"][i] for i in range(nPass) if tVals["channel"][i]==extChan])

            if len(enfArr)==0:
                print("No hits in channel %d found.  Continuing ..." % extChan)
                continue

            eAvg, attVal = np.mean(enfArr), attList[i]
            attVals.append(attVal)
            sloVals.append((np.mean(sloArr),np.std(sloArr)))
            rtVals.append((np.mean(rtArr), np.std(rtArr)))

            p1.plot(sloArr, rtArr, ".", markersize=2., c=cmap(i), label="rt %d" % (attVal), alpha=0.3)
            print("Run %d, rt %d, eAvg %.2f" % (run,attVal,eAvg))

        xdata, xe = np.asarray([val[0] for val in sloVals]), np.asarray([val[1] for val in sloVals])
        ydata, ye = np.asarray([val[0] for val in rtVals]), np.asarray([val[1] for val in rtVals])

        p1.errorbar(xdata, ydata, yerr=ye, xerr=xe, fmt="--o", ecolor='black')

        def linear(x,m,b):
            return m*x + b
        popt,_ = curve_fit(linear, xdata, ydata)
        p1.plot(xdata, linear(xdata, *popt), 'r-', label="fit\nm=%.2f\nb=%.2f" % tuple(popt))

        print("RT = %.2f * FS + %.2f" % tuple(popt))

        rtfsSlope.append(popt[0])
        rtfsInt.append(popt[1])

        p1.set_title("pIdx %d  channel %d  avgE %.2f" % (pIdx, extChan, eAvg))
        p1.set_xlabel("fitSlo",horizontalalignment='right',x=1.0)
        p1.set_ylabel("riseTime (ns)",horizontalalignment='right',y=1.0)
        p1.legend(numpoints=1,loc="best")

        xdata = np.asarray(attVals)
        popt,_ = curve_fit(linear,xdata,ydata)
        p2.plot(xdata, ydata, "o", c='b', markersize=5, label='data')
        p2.plot(xdata, linear(xdata, *popt), 'r-', label="fit\nm=%.2f\nb=%.2f" % tuple(popt))
        p2.set_xlabel("Waveform Gen. setting (ns)", horizontalalignment='right',x=1.0)
        p2.set_ylabel("t90-t10 riseTime (ns)", horizontalalignment='right',y=1.0)
        p2.set_xlim(min(xdata)*0.98, max(xdata)*1.02)
        p2.set_ylim(min(ydata)*0.98, max(ydata)*1.02)
        p2.legend(loc='best')

        plt.tight_layout()
        plt.savefig("../plots/rtStudy_idx%d.pdf" % (pIdx))

        print("RT = %.2f * FGen + %.2f" % tuple(popt))

    rtfsSlope, rtfsInt = np.asarray(rtfsSlope), np.asarray(rtfsInt)
    print("For pIdxs:",pIdxs)
    print("RT/FS slope: %.3f  stdev: %.3f" % (np.mean(rtfsSlope), np.std(rtfsSlope)))
    print("RT/FS x-int: %.3f  stdev: %.3f" % (np.mean(rtfsInt), np.std(rtfsInt)))


def getEff():
    """ Efficiency vs. energy, each detector in Test 3"""
    from ROOT import TChain, GATDataSet

    extPulserInfo = calInfo.GetSpecialList()["extPulserInfo"]
    syncChan = wl.getChan(0,10,0) # 672

    dsNum, modNum, calIdx = 0, 1, 33
    calDB = db.TinyDB('../calDB.json')
    pars = db.Query()
    fsD = ds.getDBRecord("fitSlo_ds%d_idx%d_m%d_Peak" % (dsNum, calIdx, modNum), False, calDB, pars)

    bkgIdx = 75 # runs 6887-6963
    thD = ds.getDBRecord("thresh_ds%d_bkgidx%d" % (dsNum,bkgIdx), False, calDB, pars)

    for pIdx in [19,20,21]:
    # for pIdx in [19]:

        runList = calInfo.GetSpecialRuns("extPulser",pIdx)
        attList = extPulserInfo[pIdx][0]
        extChan = extPulserInfo[pIdx][-1]
        fsCut = fsD[extChan][2] # 90% value (used in LAT3)

        effVals, threshVals, trigVals = [], [], []
        eneVals, sloVals, rtVals = [], [], []
        for i, run in enumerate(runList):
            if run in [7225, 7233]:
                continue

            # elogs: "20 Hz, 150 second runs"
            gds = GATDataSet(run)
            runTime = gds.GetRunTime() # sec
            pulseRate = 20 # Hz

            fileList = ds.getLATRunList([run],"%s/lat" % (ds.specialDir))
            latChain = TChain("skimTree")
            for f in fileList:
                latChain.Add("%s/lat/%s" % (ds.specialDir,f))

            tNames = ["Entry$","mH","channel","trapENFCal","fitSlo","den90","den10","threshKeV","threshSigma"]
            theCut = "(channel==%d || channel==%d) && mH==2" % (syncChan, extChan) # enforce correct sync
            tVals = wl.GetVX(latChain,tNames,theCut)
            nPass = len(tVals["Entry$"])

            enfArr = [tVals["trapENFCal"][i] for i in range(nPass) if tVals["channel"][i]==extChan]
            sloArr = [tVals["fitSlo"][i] for i in range(nPass) if tVals["channel"][i]==extChan]
            rtArr  = [tVals["den90"][i] - tVals["den10"][i] for i in range(nPass) if tVals["channel"][i]==extChan]

            if len(enfArr)==0:
                print("Run %d, No hits in channel %d found.  Continuing ..." % (run,extChan))
                continue

            eneVals.extend(enfArr)
            sloVals.extend(sloArr)
            rtVals.extend(rtArr)

            thr = [tVals["threshKeV"][i] for i in range(nPass) if tVals["channel"][i]==extChan][0]
            sig = [tVals["threshSigma"][i] for i in range(nPass) if tVals["channel"][i]==extChan][0]
            if thr < 99999 and sig < 99999:
                threshVals.append((thr,sig))

            muE, stdE = np.mean(np.asarray(enfArr)), np.std(np.asarray(enfArr))
            muF, stdF = np.mean(np.asarray(sloArr)), np.std(np.asarray(sloArr))
            nTot = len(sloArr)
            nAcc = len([fs for fs in sloArr if fs < fsCut])
            eff = (nAcc/nTot)
            effVals.append((muE,eff))

            nHits = len(enfArr)
            expHits = runTime * pulseRate
            trigEff = nHits / expHits
            trigVals.append((muE,trigEff))

            print("pIdx %d  run %d  chan %d  nHits %d  (exp %d) muE %.2f  muFS %.2f  eff %.2f  trigEff %.2f" % (pIdx,run,extChan,nHits,expHits,muE,muF,eff,trigEff))

        eneVals, sloVals, rtVals = np.asarray(eneVals), np.asarray(sloVals), np.asarray(rtVals)

        fig = plt.figure(figsize=(8,8),facecolor='w')
        p1 = plt.subplot(211)
        p2 = plt.subplot(212)

        xLo, xHi, bpX = 0, 50, 0.2
        yLo, yHi, bpY = 0, 300, 1.
        nbX, nbY = int((xHi-xLo)/bpX), int((yHi-yLo)/bpY)

        _,_,_,im = p1.hist2d(eneVals, sloVals, bins=[nbX, nbY], range=[[xLo,xHi],[yLo,yHi]], norm=LogNorm(), cmap='jet')
        fig.colorbar(im, ax=p1)
        p1.axhline(fsCut, color='black', linewidth=3)

        p1.set_title("pIdx %d  channel %d  fsCut %.2f" % (pIdx, extChan, fsCut))
        p1.set_xlabel("trapENFCal (keV)", horizontalalignment='right',x=1.0)
        p1.set_ylabel("fitSlo", horizontalalignment='right',y=1.0)

        xvals = np.asarray([val[0] for val in sorted(effVals) if 1. < val[0] < 10.])
        yvals = np.asarray([val[1] for val in sorted(effVals) if 1. < val[0] < 10.])
        p2.plot(xvals, yvals, "o", c='b', markersize=5, label='data')

        popt1kev,_ = curve_fit(threshFunc, xvals, yvals)
        xnew = np.arange(0, max(xvals), 0.1)
        p2.plot(xnew, threshFunc(xnew, *popt1kev), 'r-', label="fit mu=%.2f\nfit sig=%.2f" % tuple(popt1kev))

        p2.set_title("fitSlo efficiency vs. trapENFCal")
        p2.set_xlabel("trapENFCal (keV)", horizontalalignment='right',x=1.0)
        p2.set_ylabel("fitSlo Efficiency (%)", horizontalalignment='right',y=1.0)
        p2.legend(loc='best')

        plt.tight_layout()
        plt.savefig("../plots/efficiency_idx%d.pdf" % pIdx)


        # plot 3 - show how the trigger efficiency is hurting us
        # compare run by run avg, db vals, and measured trig. efficiency
        fig2 = plt.figure(figsize=(9,7),facecolor='w')

        xvals = np.asarray([val[0] for val in sorted(effVals) if val[0] < 10])
        yvals = np.asarray([val[1] for val in sorted(effVals) if val[0] < 10])
        plt.plot(xvals, yvals, "o", c='b', markersize=10, label='data')
        plt.plot(xnew, threshFunc(xnew, *popt1kev), 'b-', label="efficiency > 1 keV\nmu=%.2f sig=%.2f" % tuple(popt1kev))

        # threshMu = np.mean(np.asarray([val[0] for val in threshVals]))
        # threshSig = np.mean(np.asarray([val[1] for val in threshVals]))
        threshMu = thD[extChan][0]
        threshSig = thD[extChan][1]
        ytrigDB = threshFunc(xnew,threshMu,threshSig)
        plt.plot(xnew, ytrigDB, '-', color='gray',
            label="DB trigger efficiency\nmu %.2f sig %.2f" % (threshMu, threshSig))

        trigVals = np.asarray([val[1] for val in sorted(trigVals) if val[0] < 10])
        print("trigVals:",trigVals)
        plt.plot(xvals, trigVals, marker='o', linestyle='-', color='black', label="Meas. trigger efficiency")

        ycorr = []
        for idx in range(len(xvals)):
            corr = threshFunc(xvals[idx], threshMu, threshSig)
            print("kev %.2f  eff %.2f  corr %.2f  corrected eff %.2f" % (xvals[idx], yvals[idx], corr, yvals[idx]/corr))
            ycorr.append(yvals[idx]/corr)
        ycorr = np.asarray(ycorr)
        plt.plot(xvals, ycorr, "o", c='r', markersize=7, label='trigger eff. corrected')

        poptPt7kev,_ = curve_fit(threshFunc, xvals, ycorr)
        plt.plot(xnew, threshFunc(xnew, *poptPt7kev), 'r-', label="efficiency fit to corrected\nmu %.2f sig %.2f" % tuple(poptPt7kev))

        plt.title("fitSlo efficiency vs. trapENFCal")
        plt.xlabel("trapENFCal (keV)", horizontalalignment='right',x=1.0)
        plt.ylabel("fitSlo Efficiency (%)", horizontalalignment='right',y=1.0)
        plt.legend(loc='best')
        plt.savefig("../plots/efficiency_idx%d_corr.pdf" % pIdx)


def threshFunc(x,mu,sig):
    from scipy.special import erf
    # return erf((x-mu)/sig)
    return 0.5*(1 + erf( (x - mu)/(sig*np.sqrt(2) ) ))
    # from scipy.special import expit
    # return expit((x-mu)/sig)


def compareRunTypes():
    """ Plot fitSlo for extPulser, bkg, and cal runs for one channel. """
    from ROOT import TChain, GATDataSet, TFile

    extPulserInfo = calInfo.GetSpecialList()["extPulserInfo"]
    syncChan = wl.getChan(0,10,0) # 672

    # for pIdx in [19,20,21]:
    for pIdx in [19]:

        runList = calInfo.GetSpecialRuns("extPulser",pIdx)
        attList = extPulserInfo[pIdx][0]
        extChan = extPulserInfo[pIdx][-1]

        eLo, eHi = 1., 100
        fsLo, fsHi = 0, 300

        # ext pulser
        extEne, extSlo = [], []
        for i, run in enumerate(runList):
            fileList = ds.getLATRunList([run],"%s/lat" % (ds.specialDir))
            extChain = TChain("skimTree")
            for f in fileList:
                extChain.Add("%s/lat/%s" % (ds.specialDir,f))
            tNames = ["Entry$","mH","channel","trapENFCal","fitSlo"]
            theCut = "(channel==%d || channel==%d) && mH==2" % (syncChan, extChan) # enforce correct sync
            theCut += " && trapENFCal > %.2f && trapENFCal < %.2f && fitSlo > %.2f && fitSlo < %.2f" % (eLo, eHi, fsLo, fsHi)
            tVals = wl.GetVX(extChain,tNames,theCut)
            nPass = len(tVals["Entry$"])
            enfArr = [tVals["trapENFCal"][i] for i in range(nPass) if tVals["channel"][i]==extChan]
            sloArr = [tVals["fitSlo"][i] for i in range(nPass) if tVals["channel"][i]==extChan]
            if len(enfArr)==0:
                print("Run %d, No hits in channel %d found.  Continuing ..." % (run,extChan))
                continue
            extEne.extend(enfArr)
            extSlo.extend(sloArr)
        extEne, extSlo = np.asarray(extEne), np.asarray(extSlo)

        # cal
        dsNum, calIdx = 0, 33
        fileList = ds.getCalFiles(dsNum, calIdx, 1, False)
        calChain = TChain("skimTree")
        for f in fileList: calChain.Add(f)
        cutFile = TFile(f)
        theCut = cutFile.Get("theCut").GetTitle()
        theCut += " && channel==%d && trapENFCal > %.2f && trapENFCal < %.2f && fitSlo > %.2f && fitSlo < %.2f" % (extChan,eLo,eHi,fsLo,fsHi)
        tNames = ["channel","trapENFCal","fitSlo"]
        tVals = wl.GetVX(calChain,tNames,theCut)
        nPass = len(tVals["channel"])
        calEne = np.asarray([tVals["trapENFCal"][i] for i in range(nPass)])
        calSlo = np.asarray([tVals["fitSlo"][i] for i in range(nPass)])

        # TODO: don't do a bkg chain, try looping over files
        # bkg - ds1, 0.7 keV threshold (because ds0 files are enormous)
        dsNum = 1
        # thisChan = extChan # keep DS0 setting - will take FOREVER to draw
        thisChan = 648 # P2D2 in DS1
        bkgChain = TChain("skimTree")
        for bkgIdx in range(ds.dsMap[dsNum]):
        # for bkgIdx in range(10):
            fileList = glob.glob("%s/bg-lat/latSkimDS%d_%d_*.root" % (ds.dataDir,dsNum,bkgIdx))
            for f in fileList:
                bkgChain.Add(f)
        cutFile = TFile(f)
        theCut = cutFile.Get("theCut").GetTitle()
        theCut += " && channel==%d && trapENFCal > %.2f && trapENFCal < %.2f && fitSlo > %.2f && fitSlo < %.2f" % (thisChan,eLo,eHi,fsLo,fsHi)
        tNames = ["channel","trapENFCal","fitSlo"]
        tVals = wl.GetVX(bkgChain,tNames,theCut)
        nPass = len(tVals["channel"])
        bkgEne = np.asarray([tVals["trapENFCal"][i] for i in range(nPass)])
        bkgSlo = np.asarray([tVals["fitSlo"][i] for i in range(nPass)])

        fig = plt.figure(figsize=(9,6),facecolor='w')

        plt.plot(calEne, calSlo, ',', c='blue', markersize=0.5, label="cal", alpha=0.5)
        plt.plot(extEne, extSlo, '.', c='black', markersize=0.5, label="extPulser")
        plt.plot(bkgEne, bkgSlo, '.', c='red', markersize=2, label="ds-1 bkg")

        plt.xlabel("trapENFCal (keV)", horizontalalignment='right',x=1.0)
        plt.ylabel("fitSlo", horizontalalignment='right',y=1.0)
        plt.legend(loc='best')
        plt.savefig("../plots/comparison_idx%d.pdf" % pIdx)


def manualBuild():
    from ROOT import TFile, TTree

    f = TFile("/global/homes/w/wisecg/project/mjddatadir/built/mjd_run5942.root")
    t = f.Get("mjdTree")

    for idx in range(t.GetEntries()):
        t.GetEntry(idx)

        sct = t.timeinfo.startClockTime/1e9
        ct = t.timeinfo.clockTime/1e9

        nHit = t.channel.size()
        for iHit in range(nHit):
            toff = t.tOffset.at(iHit)/1e9
            chan = t.channel.at(iHit)

            print("%d  %d  %d  %.9f  %.2e" % (idx, iHit, chan, ct-sct, toff))


        # if idx > 100:
        #     print("    ")
        #     break




def getEffMultiChan():
    """ Efficiency vs. energy, comparing two different rise times.
        Requires that Test 1 be matched w/ the sync channel properly.
        Try manually building the events within 4us of each other.
        NOTE: The timestamps of card 11 (channels 688-697) have an incorrect starting value.
            So try to deal with that shit separately.
            Can probably use the sync channel.
    """
    from ROOT import TChain

    extPulserInfo = calInfo.GetSpecialList()["extPulserInfo"]
    syncChan = wl.getChan(0,10,0) # 672

    # for pIdx in [7,8,9,10,11,12]:
    # for pIdx in [10]:
    for pIdx in [11]:
        print("PIDX", pIdx)

        runList = calInfo.GetSpecialRuns("extPulser",pIdx)
        attList = extPulserInfo[pIdx][0]
        extChan = extPulserInfo[pIdx][-1]

        eneVals, sloVals = [], []
        for i, run in enumerate(runList):
            print(run, extChan)

            latChain = TChain("skimTree")
            fileList = ds.getLATRunList([run],"%s/lat" % (ds.specialDir))
            for f in fileList: latChain.Add("%s/lat/%s" % (ds.specialDir,f))

            iEvt = 0
            for idx in range(latChain.GetEntries()):
                latChain.GetEvent(idx)
                if latChain.channel.size() > 1:
                    print("Error, run %d idx %d, channel size is %d." % (run,idx,latChain.channel.size()))
                    return
                chan = latChain.channel.at(0)
                ene = latChain.trapENFCal.at(0)

                sct = latChain.startClockTime_s
                ct = latChain.clockTime_s

                # if chan==syncChan or chan==extChan:
                print(idx,chan,ct-sct)

                if idx > 100:
                    print("    ")
                    break


def wfFitter():
    """ fitSlo vs energy for different wf fit methods.
        Probably should save this for ext3.py
    """


if __name__=="__main__":
    main()
