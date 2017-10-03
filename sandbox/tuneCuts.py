#!/usr/bin/env python
import sys, shlex, glob, os
sys.argv += [ '-b' ] # force ROOT to be loaded in batch mode.
import ROOT
from ROOT import gROOT, gStyle, gPad
from ROOT import TFile, TTree, TChain, TCanvas, TH1D, TH2D, TF1, TLegend, TLine, TGraph
import numpy as np
import matplotlib.pyplot as plt
import waveLibs as wl
import DataSetInfo as ds

def main(argv):
    """ Calculates channel-by-channel cut values for LAT parameters.
        Cuts are calculated independent of other cut values.
    """
    # -- ROOT options --
    gROOT.ProcessLine(".x ~/env/MJDClintPlotStyle.C")
    gROOT.ProcessLine("gErrorIgnoreLevel = 3001;") # suppress ROOT messages

    # -- Get user args --
    dsNum = int(argv[0])
    fastMode = False
    f = dict.fromkeys(shlex.split('a b c d e'),False) # make a bunch of bools
    for i,opt in enumerate(argv):
        if opt == "-f":           fastMode = True
        if opt == "-bcMax":       f['a'] = True
        if opt == "-noiseWeight": f['b'] = True
        if opt == "-bcTime":      f['c'] = True
        if opt == "-tailSlope":   f['d'] = True
        if opt == "-fitSlo":      f['e'] = True
        if opt == "-burst":       burstCut(dsNum); return

    # -- Load chains for this DS --
    homePath = os.path.expanduser('~') # glob doesn't know how to expand this
    inPath = homePath + "/project/lat/latSkimDS%d*.root" % dsNum
    fileList = glob.glob(inPath)
    cutFile = TFile(fileList[0])
    theCut = cutFile.Get("theCut").GetTitle()

    cal = TChain("skimTree")
    cal.Add("~/project/cal-lat/latSkimDS%d*.root" % dsNum)

    # -- Load channel list --
    chList = ds.GetGoodChanList(dsNum)
    if dsNum==5: # remove 692 and 1232
        chList = [584, 592, 598, 608, 610, 614, 624, 626, 628, 632, 640, 648, 658, 660, 662, 672, 678, 680, 688, 690, 694, 1106, 1110, 1120, 1124, 1128, 1170, 1172, 1174, 1176, 1204, 1208, 1298, 1302, 1330, 1332]
    # chList = [578]
    print chList

    # -- Run routines --
    if f['a']: bcMax(dsNum,theCut,cal,chList,fastMode)
    if f['b']: noiseWeight(dsNum,theCut,cal,chList,fastMode)
    if f['c']: bcTime(dsNum,theCut,cal,chList,fastMode)
    if f['d']: tailSlope(dsNum,theCut,cal,chList,fastMode)
    if f['e']: fitSlo(dsNum,theCut,cal,chList,fastMode)


def bcMax(dsNum,theCut,cal,chList,fastMode):
    """ ./tuneCuts.py [dsNum] -bcMax
    Tune the LAT high-frequency parameter 'bcMax' by keeping 95%
    of the events between 5-30 kev in calibration data.
    Independent of other cuts.
    """
    c = TCanvas("c","c",1600,600)
    c.Divide(3,1,0.00001,0.00001)
    for ch in chList:

        var = "bcMax"
        eb, elo, ehi = 10,0,30
        e1dCut = 5.
        vb, vlo, vhi = 1000,100,2000
        d2Draw = "bcMax:trapENFCal>>b(%d,%d,%d,%d,%d,%d)"%(eb,elo,ehi,vb,vlo,vhi)
        d2Cut = theCut + " && channel==%d" % ch
        d1Cut = theCut + " && trapENFCal > %d && trapENFCal < %d && channel==%d" % (e1dCut,ehi,ch)
        outPlot = "./plots/bcMax/bcMax_ds%d_ch%d.png" % (dsNum,ch)
        cut99,cut95,cut01,cut05,cut90 = MakeCutPlot(c,cal,var,eb,elo,ehi,vb,vlo,vhi,d2Draw,d2Cut,d1Cut,outPlot,fastMode)
        print "DS-%d  Ch %d  bcMax99 %.2f  bcMax95 %.2f  bcMax90 %.2f" % (dsNum,ch,cut99,cut95,cut90)


def noiseWeight(dsNum,theCut,cal,chList,fastMode):
    """ ./tuneCuts.py [dsNum] -noiseWeight
    Tune the LAT (combination) parameter 'noiseWeight' by keeping the interval
    1% - 99% of the events between 5-30 kev in calibration data.
    Independent of other cuts.
    """
    c = TCanvas("c","c",1600,600)
    c.Divide(3,1,0.00001,0.00001)
    for ch in chList:

        var = "(waveS4-waveS1)/bcMax/trapENFCal"
        eb, elo, ehi = 10,0,30
        e1dCut = 5.
        vb, vlo, vhi = 1000,0,2
        d2Draw = "(waveS4-waveS1)/bcMax/trapENFCal:trapENFCal>>b(%d,%d,%d,%d,%d,%d)"%(eb,elo,ehi,vb,vlo,vhi)
        d2Cut = theCut + " && channel==%d" % ch
        d1Cut = theCut + " && trapENFCal > %d && trapENFCal < %d && channel==%d" % (e1dCut,ehi,ch)
        outPlot = "./plots/noiseWeight/noiseWeight_ds%d_ch%d.png" % (dsNum,ch)

        cut99,cut95,cut01,cut05,cut90 = MakeCutPlot(c,cal,var,eb,elo,ehi,vb,vlo,vhi,d2Draw,d2Cut,d1Cut,outPlot,fastMode)
        print "DS-%d  Ch %d  noiseWt99 %.2f  noiseWt01 %.2f" % (dsNum,ch,cut99,cut01)


def tailSlope(dsNum,theCut,cal,chList,fastMode):
    """ ./tuneCuts.py [dsNum] -tailSlope
    Tune the LAT 'tailSlope' parameters by keeping the interval
    1% - 99% of the events between 5-30 kev in calibration data.
    Do this for "pol2" and "pol3".
    Independent of other cuts.
    """
    c = TCanvas("c","c",1600,600)
    c.Divide(3,1,0.00001,0.00001)
    for ch in chList:

        # pol2
        var = "pol2"
        eb, elo, ehi = 10,0,30
        e1dCut = 5.
        vb, vlo, vhi = 1000,-6e-6,6e-6
        # d2Draw = "pol2:trapENFCal>>b(%d,%d,%d,%d,%d,%d)"%(eb,elo,ehi,vb,vlo,vhi) # ROOT doesn't like this
        d2Draw = "pol2:trapENFCal>>b(10,0,30,1000,-6e-6,6e-6)"
        d2Cut = theCut + " && channel==%d" % ch
        d1Cut = theCut + " && trapENFCal > %d && trapENFCal < %d && channel==%d" % (e1dCut,ehi,ch)
        outPlot = "./plots/tailSlope/pol2_ds%d_ch%d.png" % (dsNum,ch)

        cut99,cut95,cut01,cut05,cut90 = MakeCutPlot(c,cal,var,eb,elo,ehi,vb,vlo,vhi,d2Draw,d2Cut,d1Cut,outPlot,fastMode)
        print "DS-%d  Ch %d  pol2-99 %.2e  pol2-01 %.2e" % (dsNum,ch,cut99,cut01)

        # pol3
        var = "pol3"
        eb, elo, ehi = 10,0,30
        e1dCut = 5.
        vb, vlo, vhi = 1000,-0.1e-9,0.1e-9
        d2Draw = "pol3:trapENFCal>>b(10,0,30,1000,-0.1e-9,0.1e-9)"
        d2Cut = theCut + " && channel==%d" % ch
        d1Cut = theCut + " && trapENFCal > %d && trapENFCal < %d && channel==%d" % (e1dCut,ehi,ch)
        outPlot = "./plots/tailSlope/pol3_ds%d_ch%d.png" % (dsNum,ch)

        cut99,cut95,cut01,cut05,cut90 = MakeCutPlot(c,cal,var,eb,elo,ehi,vb,vlo,vhi,d2Draw,d2Cut,d1Cut,outPlot,fastMode)
        print "DS-%d  Ch %d  pol3-99 %.2e  pol3-01 %.2e" % (dsNum,ch,cut99,cut01)


def bcTime(dsNum,theCut,cal,chList,fastMode):
    """ ./tuneCuts.py [dsNum] -bcTime
    Tune the LAT 'bcTime' parameters by keeping events ABOVE 1%
    of the events between 5-30 kev in calibration data.

    NOTE: The bandTime parameter in LATv1 is wrong by an offset.
    After a re-skimming, the offset will need to be removed in these commands.
    """
    c = TCanvas("c","c",1600,600)
    c.Divide(3,1,0.00001,0.00001)
    for ch in chList:

        var = "(bandTime-tOffset-1100)/(matchTime-tOffset)"
        eb, elo, ehi = 10,0,5
        e1dCut = 0.
        vb, vlo, vhi = 1000,-0.5,3
        d2Draw = "(bandTime-tOffset-1100)/(matchTime-tOffset):trapENFCal>>b(%d,%d,%d,%d,%d,%d)"%(eb,elo,ehi,vb,vlo,vhi)
        d2Cut = theCut + " && channel==%d" % ch
        d1Cut = theCut + " && trapENFCal > %d && trapENFCal < %d && channel==%d" % (e1dCut,ehi,ch)
        outPlot = "./plots/bcTime/bctime_ds%d_ch%d.png" % (dsNum,ch)

        cut99,cut95,cut01,cut05,cut90 = MakeCutPlot(c,cal,var,eb,elo,ehi,vb,vlo,vhi,d2Draw,d2Cut,d1Cut,outPlot,fastMode)
        print "DS-%d  Ch %d  bcTime01 %.2e" % (dsNum,ch,cut01)


def fitSlo(dsNum,theCut,cal,chList,fastMode):
    """ ./tuneCuts.py [dsNum] -fitSlo
    Tune the LAT slowness parameter 'fitSlo' by keeping 95%
    of the events in the 238 peak of calibration data.

    NOTE: Brian and I think this (physics-based) cut should be applied
    AFTER the noise-rejection cuts - DEPENDENT on them.
    This makes the cut WAAAY tighter at 95%, so I changed to 99%.

    NOTE: the slowness model might not be good at the 238 peak (S/N too high?)
    David proposed testing this on a set of dummy waveforms, taken from each detector's DEP,
    and scaled down to random energies.
    """
    c = TCanvas("c","c",1600,600)
    c.Divide(3,1,0.00001,0.00001)
    for ch in chList:

        avseCut = " && avse > -1"

        bandTimeCut = " && bandTime-tOffset-1100 < 11000"

        channelCut = " && channel==%d" % ch

        bcMaxCut = " && bcMax < %.2e" % ds.bcMax[dsNum][ch][1] # 95% value

        noiseWt = "(waveS4-waveS1)/bcMax/trapENFCal"
        noiseWtCut = " && %s < %.2e && %s > %.2e" % (noiseWt,ds.noiseWt[dsNum][ch][0],noiseWt,ds.noiseWt[dsNum][ch][1])

        tailSlopeCut = " && pol2 > %.2e && pol2 < %.2e && pol3 > %.2e && pol3 < %.2e" % (ds.pol2[dsNum][ch][1],ds.pol2[dsNum][ch][0],ds.pol3[dsNum][ch][1],ds.pol3[dsNum][ch][0])

        bcTimeCut = " && (bandTime-tOffset-1100)/(matchTime-tOffset) > %.2e" % ds.bcTime[dsNum][ch]

        megaCut = theCut + avseCut + bandTimeCut + channelCut + bcMaxCut + noiseWtCut + tailSlopeCut + bcTimeCut

        var = "fitSlo"
        eb, elo, ehi = 10,237,240
        vb, vlo, vhi = 200,0,30
        d2Draw = "fitSlo:trapENFCal>>b(%d,%d,%d,%d,%d,%d)"%(eb,elo,ehi,vb,vlo,vhi)
        d2Cut = megaCut
        d1Cut = megaCut + " && trapENFCal > %d && trapENFCal < %d && channel==%d" % (elo,ehi,ch)

        outPlot = "./plots/fitSlo/fitSlo_ds%d_ch%d.png" % (dsNum,ch)

        cut99,cut95,cut01,cut05,cut90 = MakeCutPlot(c,cal,var,eb,elo,ehi,vb,vlo,vhi,d2Draw,d2Cut,d1Cut,outPlot,fastMode)
        print "DS-%d  Ch %d  fitSlo99 %.2f  fitSlo95 %.2f  fitSlo90 %.2f" % (dsNum,ch,cut99,cut95,cut90)


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


def burstCut(dsNum):
    """ ./tuneCuts.py [dsNum] -burst """

    # rates = {0:(30,5), 1:(20,5), 3:(50,5), 4:(20,3), 5:(150.,5.)} # v1 - before fitSLo
    rates = {0:(20,5), 1:(20,5), 3:(20,5), 4:(20,5), 5:(40.,5)} # v2 - after fitSlo

    chList = ds.GetGoodChanList(dsNum)
    nDets = len(chList)
    maxRate = rates[dsNum][0]
    maxChanRate = rates[dsNum][0] * rates[dsNum][1] / float(nDets)

    print "maxRate %d  nDets %d  factor %d  maxChanRate %.2f" % (rates[dsNum][0],nDets,rates[dsNum][1],maxChanRate)
    energyCut = "trapENFCal >= 1"

    ignoreList = {0:[656], 3:[592,692], 4:[1332], 5:[692,1232,1124]}
    bkg = ROOT.TChain("skimTree")
    for ch in chList:
        if ch not in ignoreList[dsNum]:
            f = "~/project/latskim/latSkimDS%d_ch%d.root" % (dsNum,ch)
            print "Added",f
            bkg.Add(f)

    # append to the text file (input to ds_livetime.cc)
    ds_livetimeList = open("burstCut_v1.txt", 'a')
    for key in ignoreList:
        killChs = ""
        for val in ignoreList[key]: killChs += " %d " % val
        ds_livetimeList.write("%s %s \n" % (key, killChs) )

    c0 = ROOT.TCanvas("c0","c0",1000,600)

    rlo, rhi = ds.dsRanges[dsNum][0], ds.dsRanges[dsNum][1]

    clo, chi = 570, 700
    if dsNum == 4: clo, chi = 1100, 1400
    if dsNum == 5: clo, chi = 570, 1400

    h0 = wl.H2D(bkg,rhi-rlo,rlo,rhi,chi-clo,clo,chi,"channel:run",energyCut,"Run Number","Channel")
    h0.Draw("colz")
    c0.SetLogz(1)

    c0.Print("./plots/burst/channelRateDS%d.pdf" % dsNum)

    c1 = ROOT.TCanvas("c","c",1600,600)
    c1.Divide(2,1,0)

    c1.cd(1)
    ROOT.gPad.SetLogy(1)

    h1 = wl.H1D(bkg,rhi-rlo,rlo,rhi,"run",energyCut,"Run Number","Counts")
    h1.SetLineColor(ROOT.kRed)
    h1.Draw()

    # Run & channel-based burst cut.

    runs,rates = wl.npTH1D(h1,"i")
    idx = np.where(rates > maxRate)
    print "Noisy runs:", runs[idx]

    burstCut, invBurstCut = energyCut + " && ", energyCut + " && ("

    print "maxChanRate:",maxChanRate
    for i,run in enumerate(runs[idx]):

        runCut = " && run==%d" % run
        # h2.append(ROOT.TH1D())

        h2 = wl.H1D(bkg,chi-clo,clo,chi,"channel",energyCut + runCut,"channel","Counts")

        if h2.GetEntries() == 0:
            continue

        chans,chRates = wl.npTH1D(h2,"i")
        idx2 = np.where(chRates > maxChanRate)
        print "run",int(run),"noisy chans",chans[idx2],"rates",chRates[idx2],"total entries",h2.GetEntries(),"runCut:",energyCut+runCut

        # Write run + channel groups to the file (input to ds_livetime.cc)
        noisyChans = ""
        for ch in chans[idx2]: noisyChans += "%d " % ch
        ds_livetimeList.write("%d %s \n" % (run,noisyChans))

        # Make the TCut
        runBurst = "&& !(run==%d && (" % run
        invBurst = "|| (run==%d && (" % run
        if i==0:
            runBurst = runBurst[3:]
            invBurst = invBurst[3:]

        for ch in chans[idx2]:
            runBurst += "channel==%d|| " % ch
            invBurst += "channel==%d|| " % ch
        runBurst = runBurst[:-3] + ")) "
        invBurst = invBurst[:-3] + ")) "
        burstCut += runBurst
        invBurstCut += invBurst

    invBurstCut += ")"

    ds_livetimeList.close()

    if len(runs[idx])==0:
        burstCut, invBurstCut = energyCut, energyCut

    # add the dead channels back in
    for ch in ignoreList[dsNum]:
        burstCut += " && channel!=%d" % ch

    print "\nBURST CUT:"
    print burstCut

    print "\nINVERSE BURST CUT:"
    print invBurstCut
    print ""

    h1a = wl.H1D(bkg,rhi-rlo,rlo,rhi,"run",burstCut)
    h1a.SetLineColor(ROOT.kBlue)
    h1a.Draw("same")

    c1.cd(2)
    ROOT.gPad.SetLogy(1)

    h1.Draw("hist")

    h1b = wl.H1D(bkg,rhi-rlo,rlo,rhi,"run",invBurstCut)
    h1b.SetLineColor(ROOT.kBlue)
    h1b.Draw("hist same")

    c1.Print("./plots/burst/burstCutDS%d.pdf" % dsNum)

    c2 = TCanvas("c2","c2",1000,600)
    c2.SetLogy(1)

    eb,elo,ehi = 150,0,30 # 5 bins/kev

    h2 = wl.H1D(bkg,eb,elo,ehi,"trapENFCal","","Energy","Counts")
    h2.SetLineColor(ROOT.kRed)

    h3 = wl.H1D(bkg,eb,elo,ehi,"trapENFCal",burstCut)
    h3.SetLineColor(ROOT.kBlue)

    h2.Draw("hist")
    h3.Draw("hist same")

    c2.Print("./plots/burst/energySpecDS%d.pdf" % dsNum)


if __name__ == "__main__":
    main(sys.argv[1:])
