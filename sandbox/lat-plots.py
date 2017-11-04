#!/usr/bin/env python
import sys, imp, glob, os
sys.argv.append("-b") # kill all interactive crap
ds = imp.load_source('DataSetInfo','../DataSetInfo.py')
wl = imp.load_source('waveLibs','../waveLibs.py')
import ROOT
from ROOT import gROOT, gStyle, gPad
from ROOT import TFile, TTree, TChain, TCanvas, TH1D, TH2D, TF1, TLegend, TLine, TGraph

from matplotlib import pyplot as plt
import seaborn as sns
sns.set(style='whitegrid', context='talk')

homePath = os.path.expanduser('~')
bgDir = homePath + "/project/bg-lat"
calDir = homePath + "/project/cal-lat"

c = TCanvas("c","Bob Ross's Canvas",800,600)

def fitMu():

    # load cal files the way job-panda does
    dsNum = 3
    # nIdx = wl.getNCalIdxs(dsNum, module=1)
    calIdx = 15
    fileList = wl.getCalFiles(dsNum, calIdx)
    lat = TChain("skimTree")
    for f in fileList: lat.Add(f)
    print "Found",lat.GetEntries(),"entries."

    h1 = wl.H2D(lat,50,0,50,500,-1000,25000,"fitMu:trapENFCalC","","Energy (keV)","fitMu","")
    h1.Draw("colz")
    c.SetLogz(1)
    c.Print("../plots/fitMu.pdf")


def wfStd():

    dsNum = 1
    calIdx = 10
    fileList = wl.getCalFiles(dsNum, calIdx)
    lat = TChain("skimTree")
    for f in fileList:
        print f
        lat.Add(f)
    cutFile = TFile(fileList[0])
    calCut = cutFile.Get("theCut").GetTitle()
    print "Found",lat.GetEntries(),"entries.  Using cut:",calCut

    # the wfstd draw segfaulted for one of the calIdx's.  ds 1, calIdx 10
    b1 = lat.GetBranch("wfstd")
    b2 = lat.GetBranch("trapENFCalC")
    print type(b1), type(b2)
    print b1.GetEntries(), b2.GetEntries()
    # lat.Draw("wfstd")

    # theCut = calCut + "&& gain==0"
    # h1 = wl.H2D(lat,50,0,5,50,0,5,"wfstd:trapENFCalC",theCut,"Energy (keV)","wfStd","")
    # h1.Draw("colz")
    # c.SetLogz(1)
    # c.Print("../plots/wfStd.pdf")


def bkgHists():
    """ brian wrote this one """

    dsNum, modNum = 1, 1

    chList = ds.GetGoodChanList(dsNum)
    if dsNum==5 and modNum == 1: # remove 692 and 1232 (both beges, so who cares)
        chList = [584, 592, 598, 608, 610, 614, 624, 626, 628, 632, 640, 648, 658, 660, 662, 672, 678, 680, 688, 690, 694]
    if dsNum==5 and modNum == 2:
        chList = [1106, 1110, 1120, 1124, 1128, 1170, 1172, 1174, 1176, 1204, 1208, 1298, 1302, 1330, 1332]

    bins,lower,upper = 2500,0,250

    skimTree = ROOT.TChain("skimTree")
    skimTree.Add("/projecta/projectdirs/majorana/users/wisecg/bg-lat/latSkimDS%d_*.root"%(dsNum))

    cuts = "gain==0 && trapENFCal>0.7 && mHL==1 && isGood && !muVeto && !(C==1&&isLNFill1) && !(C==2&&isLNFill2) && C!=0&&P!=0&&D!=0"

    skimCut = ROOT.TChain("skimTree")
    hList = []
    hCutList = []
    outFile = ROOT.TFile("./data/BkgHisto_DS%d.root"%(dsNum), "RECREATE")
    outFile.cd()
    hFullTotal = ROOT.TH1D("hFullTotal", "", bins,lower,upper)
    hfitSloTotal = ROOT.TH1D("hfitSloTotal", "", bins,lower,upper)
    for idx, ch in enumerate(chList):
        print "Drawing histograms for channel", ch
        hList.append(ROOT.TH1D())
        hCutList.append(ROOT.TH1D())
        skimCut.Reset()
        skimCut.Add("/projecta/projectdirs/majorana/users/wisecg/cuts/fs/fitSlo-DS%d-*-ch%d.root"%(dsNum, ch))
        hCutList[idx] = wl.H1D(skimCut,bins,lower,upper, "trapENFCalC", cuts+"&&channel==%d"%(ch), Title="hfitSlo_ch%d"%(ch))
        hList[idx] = wl.H1D(skimTree,bins,lower,upper,"trapENFCalC", cuts+"&&channel==%d"%(ch),Title="hFull_ch%d"%(ch))
        hfitSloTotal.Add(hCutList[idx])
        hFullTotal.Add(hList[idx])
        hCutList[idx].Write()
        hList[idx].Write()

    hfitSloTotal.Write()
    hFullTotal.Write()
    outFile.Close()


def fitSlo():
    """ make a raw spectrum + fitslo cut for each ds.
        rn we have DS1234.
    """

    dsNum = 1
    fsList = glob.glob("/global/homes/w/wisecg/project/cuts/fs/fitSlo-DS%d*" % dsNum)
    latList = glob.glob("/globa/homes/w/wisecg/project/bg-lat/latSkimDS%d" % dsNum)

    # chList = ds.GetGoodChanList(dsNum)
    # outFile = "/global/homes/w/wisecg/project/cuts/fs/fitSlo-DS%d-%d-ch%d.root" % (dsNum, subNum, ch)



if __name__=="__main__":
    gStyle.SetOptStat(0)
    gROOT.ProcessLine("gErrorIgnoreLevel = 3001;")

    # fitMu()
    # wfStd()
    fitSlo()