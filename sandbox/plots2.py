#!/usr/bin/env python
import sys, imp, glob, os
sys.argv.append("-b") # kill interactivity before loading ROOT (PDSF)
ds = imp.load_source('DataSetInfo','../DataSetInfo.py')
wl = imp.load_source('waveLibs','../waveLibs.py')
from ROOT import gStyle, gROOT
from ROOT import TFile, TChain, TTree, TNamed, TCanvas

def main(argv):

    gStyle.SetOptStat(0)
    gROOT.ProcessLine("gErrorIgnoreLevel = 3001;")

    # genRawHists()
    readRawHists()


def genRawHists():
    """ Create "raw" LAT histograms for each dataset.
        Use a super fine binning and huge energy range for versatility.
    """
    kpb = 0.01
    eLo, eHi = 0., 15000.
    nBins = int((eHi-eLo)/kpb)

    for dsNum in [0]:

        # PLAN:
        # sum histogram
        # CPD:energy
        # CPD:run
        # could also do channel-by-channel, bkgIdx by bkgIdx, but hopefully CPD:run will be enough for diagnostics.

        fList = sorted(glob.glob("/global/homes/w/wisecg/project/bg-lat/latSkimDS%d*.root" % dsNum))
        fCut = TFile(fList[0])
        dcCut = fCut.Get("theCut").GetTitle()

        ch = TChain("skimTree")
        [ch.Add(f) for f in fList[:3]]

        fOut = TFile("../data/latRaw_DS%d.root" % dsNum,"RECREATE")

        h1 = wl.H1D(ch,nBins,eLo,eHi,"trapENFCal",dcCut,"trapENFCal","Counts","raw sum LAT spectrum","hSum")
        h1.Write()

        # def H1D(tree,bins,xlo,xhi,drawStr,cutStr,xTitle="",yTitle="",Title=None, Name=None):
        # def H2D(tree,xbins,xlo,xhi,ybins,ylo,yhi,drawStr,cutStr,xTitle="",yTitle="",Title=None, Name=None):

        fOut.Close()


def readRawHists():

    dsNum = 0
    fIn = TFile("../data/latRaw_DS%d.root" % dsNum)
    fIn.Print()

    h1 = fIn.Get("hSum")
    print h1.GetEntries()

    c = TCanvas("c","Bob Ross's Canvas",800,600)
    h1.Draw("hist")
    c.Print("../plots/latRawDS%d.pdf" % dsNum)



if __name__=="__main__":
    main(sys.argv[1:])