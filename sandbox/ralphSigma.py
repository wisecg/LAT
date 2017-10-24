#!/usr/bin/env python
from ROOT import TFile, TTree, TCanvas, TH1D
import numpy as np

def main():

    f = TFile("./data/latSkimDS1_run10247_0.root","READ")
    t = f.Get("skimTree")
    nEnt = t.GetEntries()

    theCut = f.Get("theCut").GetTitle()
    # theCut += " && trapENFCal  1.1"

    print "Using cut:\n",theCut,"\n"
    waveTree.Draw(">>elist", theCut, "entrylist")
    elist = gDirectory.Get("elist")
    waveTree.SetEntryList(elist)
    nList = elist.GetN()
    print "Found",nList,"entries passing cuts."

    fig = plt.figure(figsize=(10,7), facecolor='w')
    p1 = plt.subplot(211)
    p2 = plt.subplot(212)

    






if __name__=="__main__":
    main()