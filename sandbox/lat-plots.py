#!/usr/bin/env python
import sys, imp, glob, os
sys.argv.append("-b") # kill all interactive crap
ds = imp.load_source('DataSetInfo','../DataSetInfo.py')
wl = imp.load_source('waveLibs','../waveLibs.py')
import ROOT
from ROOT import gROOT, gStyle, gPad
from ROOT import TFile, TTree, TChain, TCanvas, TH1D, TH2D, TF1, TLegend, TLine, TGraph
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
    calIdx = 7
    fileList = wl.getCalFiles(dsNum, calIdx)
    lat = TChain("skimTree")
    for f in fileList: lat.Add(f)
    cutFile = TFile(fileList[0])
    calCut = cutFile.Get("theCut").GetTitle()
    print "Found",lat.GetEntries(),"entries.  Using cut:",calCut

    b1 = lat.GetBranch("wfstd")
    b2 = lat.GetBranch("trapENFCalC")
    print type(b1), type(b2)
    print b1.GetEntries(), b2.GetEntries()

    # lat.Draw("wfstd")


    h1 = wl.H2D(lat,500,0,50,500,0,50,"wfstd:trapENFCalC",calCut,"Energy (keV)","wfStd","")
    h1.Draw("colz")
    c.SetLogz(1)
    c.Print("./plots/wfStd.pdf")




if __name__=="__main__":
    gStyle.SetOptStat(0)
    gROOT.ProcessLine("gErrorIgnoreLevel = 3001;")

    # fitMu()
    wfStd()