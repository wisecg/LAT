#!/usr/bin/env python
import sys, imp
ds = imp.load_source('DataSetInfo','../DataSetInfo.py')
wl = imp.load_source('waveLibs','../waveLibs.py')
sys.argv += [ '-b' ]  # force ROOT to be loaded in batch mode
import ROOT
from ROOT import gROOT, gStyle, gPad
from ROOT import TFile, TTree, TChain, TCanvas, TH1D, TH2D, TF1, TLegend, TLine, TGraph
"""
Control function is at the bottom.
h1 = wl.H1D(tree,bins,xlo,xhi,drawStr,cutStr,xTitle="",yTitle="",Title=None, Name=None)
h2 = wl.H2D(tree,xbins,xlo,xhi,ybins,ylo,yhi,drawStr,cutStr,xTitle="",yTitle="",Title=None, Name=None)
"""

def fitMu():

    c = TCanvas()






if __name__=="__main__":
    gStyle.SetOptStat(0)
    gROOT.ProcessLine("gErrorIgnoreLevel = 3001;")
    fitMu()