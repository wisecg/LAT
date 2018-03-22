#!/usr/bin/env python3
import sys
# import waveLibs as wl
from ROOT import TFile, TTree, gROOT
gROOT.ProcessLine("gErrorIgnoreLevel = 3001;")

def main():

    f1 = TFile("../skimDS1_run9456_low.root")
    f2 = TFile("../skimDS1_run9456_low_s.root")

    t1 = f1.Get("skimTree")
    t2 = f2.Get("skimTree")

    print("low:",f1.GetSize()/1e6, "low_s",f2.GetSize()/1e6)

    n1 = t1.Draw("trapENFCal","","goff")
    n2 = t2.Draw("trapENFCal","","goff")
    print("hits: low %d  low_s %d" % (n1, n2))

    n1 = t1.Draw("trapENFCal","gain==0","goff")
    n2 = t2.Draw("trapENFCal","gain==0","goff")
    print("HG hits: low %d  low_s %d" % (n1, n2))

    n1 = t1.Draw("trapENFCal","gain==1","goff")
    n2 = t2.Draw("trapENFCal","gain==1","goff")
    print("LG hits: low %d  low_s %d" % (n1, n2))

    n1 = t1.Draw("sumEH","","goff")
    n2 = t2.Draw("sumEH","","goff")
    print("sumEH hits: low %d  low_s %d" % (n1, n2))


if __name__=="__main__":
    main()