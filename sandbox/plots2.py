#!/usr/bin/env python
import sys, imp, glob, os, time
sys.argv.append("-b") # kill interactivity before loading ROOT (PDSF)
ds = imp.load_source('DataSetInfo','../DataSetInfo.py')
wl = imp.load_source('waveLibs','../waveLibs.py')
from ROOT import gStyle, gROOT
from ROOT import TFile, TChain, TTree, TNamed, TCanvas, TLegend
import ROOT

def main(argv):

    gStyle.SetOptStat(0)
    gROOT.ProcessLine("gErrorIgnoreLevel = 3001;")
    gStyle.SetPalette(ROOT.kBlueRedYellow)

    for i,opt in enumerate(argv):
        if opt == "-gen":  genRawHists(int(argv[i+1]))
        if opt == "-read": readRawHists(int(argv[i+1]))
        if opt == "-raw":
            for dsNum in [0,1,2,3,4,5]:
                genRawHists(dsNum)
                readRawHists(dsNum)
        if opt == "-ord": optimizeCutOrder(int(argv[i+1]))
        if opt == "-cut": makeCutSpectra(int(argv[i+1]))
        if opt == "-all":
            [makeCutSpectra(dsNum) for dsNum in [0,1,2,3,4,5]]
        if opt == "-thr": applyThresholdCut(int(argv[i+1]))


def genRawHists(dsNum):
    """ ./plots2.py -gen [dsNum]
        Create "raw" LAT histograms for each dataset.
        Use a super fine binning and huge energy range for versatility.
        Don't distinguish between enriched and natural detectors.
        1. sum histogram 'hSum'
        2. cpd vs energy 'hCPDE'
        3. cpd vs run    'hCPDrun'
        Could also do channel-by-channel, bkgIdx by bkgIdx,
            but hopefully the 2D CPD plots will be enough for diagnostics (can examine their bin contents).
            Plus we can loop over the bins in the raw histo instead of re-drawing everything.
    """
    print "Plots2 started:",time.strftime('%X %x %Z')
    startT = time.clock()

    kpb = 0.01
    eLo, eHi = 0., 15000.
    nBins = int((eHi-eLo)/kpb)

    # build the file list sequentially
    fList = []
    dsPath = "/global/homes/w/wisecg/project/bg-lat"
    for bkgIdx in range(ds.dsMap[dsNum]+1):
        tmp = glob.glob("%s/latSkimDS%d_%d_*.root" % (dsPath, dsNum, bkgIdx))
        nFiles = len(tmp)
        for subIdx in range(nFiles):
            fName = "%s/latSkimDS%d_%d_%d.root" % (dsPath, dsNum, bkgIdx, subIdx)
            if os.path.isfile(fName) == True:
                fList.append(fName)
            else:
                print "File doesn't exist, exiting ...",fName
                return

    # fList = fList[:2] # limit fList (debug only)

    fCut = TFile(fList[0])
    dcCut = fCut.Get("theCut").GetTitle()
    dcCut += "&& gain==0"

    ch = TChain("skimTree")
    [ch.Add(f) for f in fList]
    print "Generating DS-%d (%d entries)." % (dsNum,ch.GetEntries())

    cpdLo, cpdHi = 111, 175
    if dsNum==4: cpdLo, cpdHi = 211, 275
    if dsNum==5: cpdLo, cpdHi = 111, 275
    nCPD = cpdHi - cpdLo

    ch.GetEntry(0)
    runLo = ch.run
    ch.GetEntry(ch.GetEntries()-1)
    runHi = ch.run
    nRuns = runHi - runLo + 1
    print "First run:",runLo,"Last run:",runHi

    fOut = TFile("../data/latRaw_DS%d.root" % dsNum,"RECREATE")

    print "Generating sum spectrum.  %.2f mins elapsed. %s" % ((time.clock()-startT)/60., time.strftime('%X %x %Z'))
    h1 = wl.H1D(ch,nBins,eLo,eHi,"trapENFCal",dcCut,"trapENFCal","Counts","raw sum LAT spectrum","hSum")
    h1.Write()

    print "Generating CPD vs. E.     %.2f mins elapsed. %s" % ((time.clock()-startT)/60., time.strftime('%X %x %Z'))
    h2 = wl.H2D(ch,nBins,eLo,eHi,nCPD,cpdLo,cpdHi,"C*100+P*10+D:trapENFCal",dcCut,"trapENFCal","CPD","CPD vs. E","hCPDE")
    h2.Write()

    print "Generating CPD vs. run.   %.2f mins elapsed. %s" % ((time.clock()-startT)/60., time.strftime('%X %x %Z'))
    h3 = wl.H2D(ch,nRuns,runLo,runHi,nCPD,cpdLo,cpdHi,"C*100+P*10+D:run",dcCut,"run","CPD","CPD vs. run","hCPDrun")
    h3.Write()

    fOut.Close()
    print "Time Elapsed:",(time.clock()-startT)/60.,"minutes."

    print "Done at",time.strftime('%X %x %Z')


def readRawHists(dsNum):
    """ ./plots2.py -read [dsNum] """

    fIn = TFile("../data/latRaw_DS%d.root" % dsNum)

    print "Reading raw hists ..."

    c = TCanvas("c","Bob Ross's Canvas",800,600)
    c.SetLogy(1)

    h1 = fIn.Get("hSum")

    # convert to a more reasonable binning and energy range
    kpb = 0.1
    rebinFactor = int(kpb/0.01)
    h1 = h1.Rebin(rebinFactor) # this creates a new histogram
    h1.GetXaxis().SetRangeUser(0.,12.)

    h1.SetTitle("")
    h1.Draw("hist")
    c.Print("../plots/latRawDS%d.pdf" % dsNum)

    c.SetLogy(0)
    # c.SetLogz(1)
    h2 = fIn.Get("hCPDE")
    h2 = h2.RebinX(rebinFactor)
    h2.GetXaxis().SetRangeUser(0.,10.)
    h2.SetTitle("")
    h2.Draw("COLZ")
    c.Print("../plots/latRawDS%d_cpde.pdf" % dsNum)

    h3 = fIn.Get("hCPDrun")
    h3.SetTitle("")
    h3.Draw("COLZ")
    c.Print("../plots/latRawDS%d_cpdRun.pdf" % dsNum)


def getIntegralCounts(hist,int2,int5,int10,int50,int250):
    """ Used by optimizeCutOrder. """
    ax = hist.GetXaxis()
    int2.append(hist.Integral(ax.FindBin(0.1), ax.FindBin(2.),"width"))
    int5.append(hist.Integral(ax.FindBin(2.), ax.FindBin(5.),"width"))
    int10.append(hist.Integral(ax.FindBin(5.), ax.FindBin(10.),"width"))
    int50.append(hist.Integral(ax.FindBin(10.), ax.FindBin(50.),"width"))
    int250.append(hist.Integral(ax.FindBin(50), ax.FindBin(100.),"width"))


def optimizeCutOrder(dsNum):
    """ ./plots2.py -ord [dsNum]

    Optimize the order of reduction:
        fs, rn, wf, fs+rn, fs+wf, rn+wf, fs+rn+wf

    now I'm also thinking about applying the threshold cut to every file ...
    ... that will tell you how good the cuts are when we're above/below the detector threshold ... do I want that?
    anyway, try it first without.  (simpler.)

    Results:
         [raw, fs, rn, wf, fs+rn, fs+wf, fs+rn, rn+wf, fs+rn+wf]
    DS1 outputs:
    Int2: [330357.60, 606.30, 314160.80, 1359.80, 52.80, 21.70, 1358.30, 21.50]
    Int5: [4362.90, 52.90, 1018.90, 54.90, 52.60, 30.50, 54.40, 30.20]
    Int10: [1991.90, 75.20, 98.00, 63.40, 74.80, 52.10, 63.40, 52.10]
    Int50: [94.20, 75.00, 90.70, 63.40, 74.90, 51.50, 63.40, 51.50]
    Int250: [51.60, 42.50, 49.80, 39.20, 42.40, 35.80, 39.10, 35.70]

    Reduction order (least to most, <2 kev): raw, rn, wf, rn+wf, fs, fs+rn, fs+wf, fs+rn+wf
    Order to plot cut reductions:
        raw, fs, fs+wf, fs+rn+wf.
    Clearly fs is our best noise discriminator.
    """
    int2, int5, int10, int50, int250 = [], [], [], [], []

    kpb = 0.1
    eLo, eHi = 0., 250.
    nBins = int((eHi-eLo)/kpb)

    fIn = TFile("../data/latRaw_DS%d.root" % dsNum)
    hRaw = fIn.Get("hSum")
    rebinFactor = int(kpb/0.01)
    hRaw = hRaw.Rebin(rebinFactor)
    getIntegralCounts(hRaw, int2, int5, int10, int50, int250)

    fileFS = sorted(glob.glob("/global/homes/w/wisecg/project/cuts/fs/fitSlo-DS%d-*.root" % dsNum))
    chainFS = TChain("skimTree")
    [chainFS.Add(f) for f in fileFS]
    hFS = wl.H1D(chainFS,nBins,eLo,eHi,"trapENFCal","","trapENFCal","","hFS")
    getIntegralCounts(hFS, int2, int5, int10, int50, int250)

    fileRN = sorted(glob.glob("/global/homes/w/wisecg/project/cuts/rn/riseNoise-DS%d-*.root" % dsNum))
    chainRN = TChain("skimTree")
    [chainRN.Add(f) for f in fileRN]
    hRN = wl.H1D(chainRN,nBins,eLo,eHi,"trapENFCal","","trapENFCal","","hRN")
    getIntegralCounts(hRN, int2, int5, int10, int50, int250)

    fileWF = sorted(glob.glob("/global/homes/w/wisecg/project/cuts/wf/wfstd-DS%d-*.root" % dsNum))
    chainWF = TChain("skimTree")
    [chainWF.Add(f) for f in fileWF]
    hWF = wl.H1D(chainWF,nBins,eLo,eHi,"trapENFCal","","trapENFCal","","hWF")
    getIntegralCounts(hWF, int2, int5, int10, int50, int250)

    fileFSRN = sorted(glob.glob("/global/homes/w/wisecg/project/cuts/fs_rn/fs_rn-DS%d-*.root" % dsNum))
    chainFSRN = TChain("skimTree")
    [chainFSRN.Add(f) for f in fileFSRN]
    hFSRN = wl.H1D(chainFSRN,nBins,eLo,eHi,"trapENFCal","","trapENFCal","","hFSRN")
    getIntegralCounts(hFSRN, int2, int5, int10, int50, int250)

    fileFSWF = sorted(glob.glob("/global/homes/w/wisecg/project/cuts/fs_wf/fs_wf-DS%d-*.root" % dsNum))
    chainFSWF = TChain("skimTree")
    [chainFSWF.Add(f) for f in fileFSWF]
    hFSWF = wl.H1D(chainFSWF,nBins,eLo,eHi,"trapENFCal","","trapENFCal","","hFSWF")
    getIntegralCounts(hFSWF, int2, int5, int10, int50, int250)

    fileRNWF = sorted(glob.glob("/global/homes/w/wisecg/project/cuts/rn_wf/rn_wf-DS%d-*.root" % dsNum))
    chainRNWF = TChain("skimTree")
    [chainRNWF.Add(f) for f in fileRNWF]
    hRNWF = wl.H1D(chainRNWF,nBins,eLo,eHi,"trapENFCal","","trapENFCal","","hRNWF")
    getIntegralCounts(hRNWF, int2, int5, int10, int50, int250)

    fileFSRNWF = sorted(glob.glob("/global/homes/w/wisecg/project/cuts/fs_rn_wf/fs_rn_wf-DS%d-*.root" % dsNum))
    chainFSRNWF = TChain("skimTree")
    [chainFSRNWF.Add(f) for f in fileFSRNWF]
    hFSRNWF = wl.H1D(chainFSRNWF,nBins,eLo,eHi,"trapENFCal","","trapENFCal","","hFSRNWF")
    getIntegralCounts(hFSRNWF, int2, int5, int10, int50, int250)

    # stackoverflow trick to make list printing prettier
    class prettyfloat(float):
        def __repr__(self): return "%.2f" % self

    print "Int2:", map(prettyfloat, int2)
    print "Int5:", map(prettyfloat, int5)
    print "Int10:", map(prettyfloat, int10)
    print "Int50:", map(prettyfloat, int50)
    print "Int250:", map(prettyfloat, int250)


def makeCutSpectra(dsNum):
    """ ./plots2.py -cut [dsNum]
    Order to plot cut reductions:
        raw, fs, fs+wf, fs+rn+wf.
    """
    print "Generating cut spectra for DS-%d ..." % (dsNum)

    kpb = 0.1
    eLo, eHi = 0., 15.
    nBins = int((eHi-eLo)/kpb)

    fIn = TFile("../data/latRaw_DS%d.root" % dsNum)
    hRaw = fIn.Get("hSum")
    rebinFactor = int(kpb/0.01)
    hRaw = hRaw.Rebin(rebinFactor)
    hRaw.GetXaxis().SetRangeUser(eLo, eHi)
    hRaw.SetTitle("")
    hRaw.SetMinimum(0.1)

    fileFS = sorted(glob.glob("/global/homes/w/wisecg/project/cuts/fs/fitSlo-DS%d-*.root" % dsNum))
    chainFS = TChain("skimTree")
    [chainFS.Add(f) for f in fileFS]
    hFS = wl.H1D(chainFS,nBins,eLo,eHi,"trapENFCal","","trapENFCal","","hFS")

    fileFSWF = sorted(glob.glob("/global/homes/w/wisecg/project/cuts/fs_wf/fs_wf-DS%d-*.root" % dsNum))
    chainFSWF = TChain("skimTree")
    [chainFSWF.Add(f) for f in fileFSWF]
    hFSWF = wl.H1D(chainFSWF,nBins,eLo,eHi,"trapENFCal","","trapENFCal","","hFSWF")

    fileFSRNWF = sorted(glob.glob("/global/homes/w/wisecg/project/cuts/fs_rn_wf/fs_rn_wf-DS%d-*.root" % dsNum))
    chainFSRNWF = TChain("skimTree")
    [chainFSRNWF.Add(f) for f in fileFSRNWF]
    hFSRNWF = wl.H1D(chainFSRNWF,nBins,eLo,eHi,"trapENFCal","","trapENFCal","","hFSRNWF")

    c = TCanvas("c","Bob Ross's Canvas",800,600)
    c.SetLogy(1)
    hRaw.SetLineColor(ROOT.kBlack)
    hRaw.Draw("hist")
    hFS.SetLineColor(ROOT.kRed)
    hFS.Draw("hist same")
    hFSWF.SetLineColor(ROOT.kMagenta)
    hFSWF.Draw("hist same")
    hFSRNWF.SetLineColor(ROOT.kBlue)
    hFSRNWF.Draw("hist same")

    leg = TLegend(0.83,0.6,0.97,0.9)
    leg.AddEntry(hRaw,"raw","l")
    leg.AddEntry(hFS,"fitSlo","l")
    leg.AddEntry(hFSWF,"+wfStd","l")
    leg.AddEntry(hFSRNWF,"+riseNoise","l")
    leg.Draw("same")

    c.Print("../plots/cutSpecDS%d.pdf" % dsNum)


def applyThresholdCut(dsNum):
    """ ./plots2.py -thr [dsNum] """

    kpb = 0.1
    eLo, eHi = 0., 15.
    nBins = int((eHi-eLo)/kpb)

    chList = ds.GetGoodChanList(dsNum)
    chList = [ch for ch in chList if ch!=692 and ch!=1232]
    # if dsNum==5 and modNum==1:
    #     chList = [ch for ch in chList if ch < 1000 and ch!=692]
    # if dsNum==5 and modNum==2:
    #     chList = [ch for ch in chList if ch > 1000 and ch!=1232]

    # build the cut file list sequentially by bkgIdx and channel
    fList, fMissing = [], []
    dsPath = "/global/homes/w/wisecg/project/cuts/fs_rn_wf"
    for bkgIdx in range(ds.dsMap[dsNum]+1):
        tmp = glob.glob("%s/fs_rn_wf-DS%d-%d-*.root" % (dsPath, dsNum, bkgIdx))
        nFiles = len(tmp)
        for chan in chList:
            fName = "%s/fs_rn_wf-DS%d-%d-ch%d.root" % (dsPath, dsNum, bkgIdx, chan)
            if os.path.isfile(fName) == True:
                fList.append(fName)
            else:
                fMissing.append(fName)

    print "nFound: %d  nMissing: %d" % (len(fList),len(fMissing))

    # save ds+bkgIdx+ch as a key to do a quick threshold lookup?  Maybe fList is a dict?

    # fileFSRNWF = sorted(glob.glob("/global/homes/w/wisecg/project/cuts/fs_rn_wf/fs_rn_wf-DS%d-*.root" % dsNum))

    # chainFSRNWF = TChain("skimTree")
    # [chainFSRNWF.Add(f) for f in fileFSRNWF]
    # hFSRNWF = wl.H1D(chainFSRNWF,nBins,eLo,eHi,"trapENFCal","","trapENFCal","","hFSRNWF")

    # c = TCanvas("c","Shan is pretty",800,600)




if __name__=="__main__":
    main(sys.argv[1:])