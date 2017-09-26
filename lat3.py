#!/usr/local/bin/python
"""
===================== LAT3.py =====================

PLACEHOLDER CODE.  THIS SHOULD BE A LAT SKIMMER,
WHICH READS DATASETINFO.PY AND PRODUCES CHANNEL-SPECIFIC
FILES WITH RUN CUTS APPLIED.

or should this be where "tuneCuts" is done?

v1: 07 Aug 2017

========= C. Wiseman (USC), B. Zhu (LANL) =========
"""
import sys, time
import numpy as np

def main(argv):

    print "========================================"
    print "LAT3 started:",time.strftime('%X %x %Z')
    startT = time.clock()

    dsNum, subNum = -1, -1
    bAll, bbcMax, b = False, False
    bCustom = False
    customPar = ""

    if len(argv) == 0:
        return
    for i, opt in enumerate(argv):
        if opt == "-all":
            bAll = True
            print "Tuning all cuts"
        if opt == "-bcMax":
            bbcMax = True
            print "Tuning bcMax"
        if opt == "-Custom":
            bCustom = True
            customPar = str(argv[i+1])
            print "Tuning custom cut parameter: ", customPar
        if opt == "-s":
            dsNum, subNum = int(argv[i+1]), int(argv[i+2])



    stopT = time.clock()
    print "Stopped:",time.strftime('%X %x %Z'),"\nProcess time (min):",(stopT - startT)/60


def TuneCut(dsNum, subNum, cal, chList, fastMode, cutName):
    c = TCanvas("c","c",1600,600)
    c.Divide(3,1,0.00001,0.00001)
    for ch in chList:

        eb, elo, ehi = 10,0,30
        e1dCut = 5.
        vb, vlo, vhi = 1000,100,2000
        d2Draw = "%s:trapENFCal>>b(%d,%d,%d,%d,%d,%d)"%(cutName,eb,elo,ehi,vb,vlo,vhi)
        d2Cut = theCut + " && channel==%d" % ch
        d1Cut = theCut + " && trapENFCal > %d && trapENFCal < %d && channel==%d" % (e1dCut,ehi,ch)
        outPlot = "./plots/%s/%s_ds%d_ch%d.png" % (cutName, cutName, dsNum,ch)
        cut99,cut95,cut01,cut05,cut90 = MakeCutPlot(c,cal,var,eb,elo,ehi,vb,vlo,vhi,d2Draw,d2Cut,d1Cut,outPlot,fastMode)
        
        # print "DS-%d  Ch %d  bcMax99 %.2f  bcMax95 %.2f  bcMax90 %.2f" % (dsNum,ch,cut99,cut95,cut90)


def MakeCutPlot(c,cal,var,eb,elo,ehi,vb,vlo,vhi,d2Draw,d2Cut,d1Cut,outPlot,fastMode):
    """ Repeated code is the DEVIL.  Even if you have to pass in 1,000,000 arguments. """

    # Calculate cut vals (assumes plot range is correct)
    h1 = wl.H1D(cal,vb,vlo,vhi,var,d1Cut)
    h1Sum = h1.Integral()
    if h1Sum == 0:
        return 0,0,0,0,0
    h1.Scale(1/h1Sum)
    cut99,cut95,cut01,cut05,cut90 = wl.GetIntegralPoints(h1)
    if fastMode:
        return cut99,cut95,cut01,cut05,cut90

    # Generate the plot for inspection.
    c.cd(2)
    gPad.SetLogy(0)
    h1.Draw("hist")

    c.cd(1)
    gPad.SetLogy(0)
    cal.Draw(d2Draw,d2Cut)

    l1, l2 = TLine(), TLine()
    l1.SetLineColor(ROOT.kGreen)
    l2.SetLineColor(ROOT.kRed)

    l1.DrawLine(elo, cut99, ehi, cut99)
    l2.DrawLine(elo, cut95, ehi, cut95)
    l2.DrawLine(elo, cut05, ehi, cut05)
    l1.DrawLine(elo, cut01, ehi, cut01)

    c.cd(3)
    x_h1, y_h1 = wl.npTH1D(h1)
    int_h1 = wl.integFunc(y_h1)
    g2 = TGraph(len(x_h1), x_h1, int_h1)
    g2.Draw("ACP")
    l1.DrawLine(cut99, 0, cut99, 1)
    l2.DrawLine(cut95, 0, cut95, 1)
    l1.DrawLine(cut01, 0, cut01, 1)
    l2.DrawLine(cut05, 0, cut05, 1)

    c.Print(outPlot)
    return cut99,cut95,cut01,cut05,cut90


if __name__ == "__main__":
    main(sys.argv[1:])
