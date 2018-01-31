#!/usr/bin/env python3
import sys, os, imp, glob
sys.argv.append("-b")
ds = imp.load_source('DataSetInfo',os.environ['LATDIR']+'/DataSetInfo.py')
wl = imp.load_source('waveLibs',os.environ['LATDIR']+'/waveLibs.py')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np

def main():

    # runByRun()
    # runTuneCut()
    getThreshold()


def runByRun():
    """ Directly confirm settings of ext pulser scripts. """
    import time
    from ROOT import TChain, GATDataSet

    calInfo = ds.CalInfo()

    runList = calInfo.GetSpecialRuns("extPulser",19)
    for run in runList:
        if run != 7234:
            continue

        fileList = []
        subFiles = glob.glob("%s/lat/latSkimDS%d_run%d_*.root" % (ds.specialDir, ds.GetDSNum(run), run))
        for idx in range(len(subFiles)):
            thisFile = "%s/lat/latSkimDS%d_run%d_%d.root" % (ds.specialDir, ds.GetDSNum(run), run, idx)
            if not os.path.isfile(thisFile):
                print("File doesn't exist: ",thisFile)
            else:
                fileList.append(thisFile)
        latChain = TChain("skimTree")
        for f in fileList: latChain.Add(f)

        detPos = wl.getDetPos(run)
        syncChan = wl.getChan(0,10,0) # 672
        extPChan = 688

        # theCut = "channel==%d || channel==%d" % (syncChan,extPChan)
        # theCut = "mH==2 && Entry$ < 50"
        # theCut = "mH==2 && Entry$ < 50 && !muVeto"
        theCut = "Entry$ < 100"
        # theCut = "Entry$ < 50 && trapENFCal > 1"
        # theCut = "mH==4 && (channel==640 || channel==646)"
        # theCut = "Entry$ < 100 && (channel==%d || channel==%d) && trapENFCal > 2" % (extPChan,syncChan)

        start = time.time()
        tNames = ["Entry$","run","channel","mH","trapENFCal","den90","den10","fitSlo","localTime_s","tOffset"]
        tVals = wl.GetVX(latChain,tNames,theCut)
        print("took",time.time()-start)

        for idx in range(tVals["run"].size):
            ent    = tVals["Entry$"][idx]
            run    = tVals["run"][idx]
            chan   = tVals["channel"][idx]
            mH     = tVals["mH"][idx]
            enf    = tVals["trapENFCal"][idx]
            d90    = tVals["den90"][idx]
            d10    = tVals["den10"][idx]
            fitSlo = tVals["fitSlo"][idx]
            gt     = tVals["localTime_s"][idx]
            tOff   = tVals["tOffset"][idx]*1e-9
            hitTime = gt+tOff
            print("%d  e%d  m%d  t%.8f  c%-4d  %-9.2f  %-8.2f  %.2f" % (run,ent,mH,hitTime,chan,enf,d90-d10,fitSlo))


def runTuneCut():
    """ Run TuneCut just like it would be in LAT3, but for ext pulser data. """
    global ROOT, TChain, TCanvas, gPad, TGraph, TLine
    from ROOT import TChain, TCanvas, gPad, TGraph, TLine
    import ROOT

    dsNum, subNum, tName, eLo, eHi = 0, 19, "extPulser", 0, 250
    chList = [624]
    par, parName = "fitSlo", "fitSlo"
    theCut = "channel==624"
    fastMode = False
    tRange = [eLo, eHi]

    cal = ds.CalInfo()
    # runList = cal.GetSpecialRuns("extPulser",subNum)
    runList = [7244]
    fList = ds.getLATRunList(runList, ds.specialDir+"/lat")
    skimTree = TChain("skimTree")
    for f in fList: skimTree.Add("%s/%s" % (ds.specialDir+"/lat", f))

    # returns the dict we would fill in the DB with.
    # cutDict = TuneCut(dsNum, subNum, tRange[0], tRange[1], tName, skimTree, chList, par, parName, theCut, fastMode)

    # make a plot to check the TuneCut plot.
    tNames = ["run","Entry$","channel","trapENFCal","fitSlo"]
    tVals = wl.GetVX(skimTree, tNames, "(channel==624 || channel==672) && Entry$ < 200")
    for idx in range(tVals["trapENFCal"].size):
        run = tVals["run"][idx]
        ent = tVals["Entry$"][idx]
        chan = tVals["channel"][idx]
        enf = tVals["trapENFCal"][idx]
        slo = tVals["fitSlo"][idx]
        print("r%d  i%d  c%d  e%-8.3f  s%.3f" % (run,ent,chan,enf,slo))
    return

    fig = plt.figure(figsize=(10,5),facecolor='w')
    xArr, yArr = tVals["trapENFCal"], tVals["fitSlo"]
    xLo, xHi, yLo, yHi = 0, 50, -20, 300
    bpY, bpX = 1., 0.2
    nBinsY, nBinsX = int((yHi-yLo)/bpY), int((xHi-xLo)/bpX)
    plt.hist2d(xArr, yArr, bins=[nBinsX,nBinsY], range=[[xLo,xHi],[yLo,yHi]], norm=LogNorm())
    plt.colorbar()
    plt.xlabel("trapENFCal (keV)",horizontalalignment='right',x=1.0)
    plt.ylabel("fitSlo",horizontalalignment='right',y=1.0)
    plt.savefig("../plots/extPulser_idx%d.pdf" % subNum)

    nOverflow = len([val for val in tVals["fitSlo"] if val > yHi or val < yLo])
    print("Found %d overflows of %d total (%.2f%%)" % (nOverflow, len(tVals["fitSlo"]), 100.*nOverflow/len(tVals["fitSlo"])))


def TuneCut(dsNum, subNum, tMin, tMax, tName, cal, chList, par, parName, theCut, fastMode):
    c = TCanvas("%s"%(parName),"%s"%(parName),1600,600)
    c.Divide(3,1,0.00001,0.00001)
    cutDict = {}
    for ch in chList:
        cutDict[ch] = [0,0,0,0,0]
        eb, elo, ehi = (tMax-tMin),tMin,tMax
        d1Cut = theCut + " && trapENFCal > %d && trapENFCal < %d && channel==%d" % (elo,ehi,ch)
        d2Cut = theCut + " && channel==%d" % ch
        nPass = cal.Draw("trapENFCal:%s"%(par), d1Cut, "goff")
        nEnergy = cal.GetV1()
        nCut = cal.GetV2()
        nCutList = list(float(nCut[n]) for n in range(nPass))
        nEnergyList = list(float(nEnergy[n]) for n in range(nPass))

        # Error and warning messages
        if len(nCutList) == 0 or len(nEnergyList) == 0:
            print("Error: Channel %d has no entries, cut cannot be set properly, setting to [0,0,0,0,0,0,0]"%(ch))
            cutDict[ch] = [0,0,0,0,0]
            continue
        if len(nCutList) <= 1000 or len(nEnergyList) <= 1000:
            print("Warning: Channel %d has less than 1000 entries, cut values may not be accurate"%(ch))

        vb, v5, v95 = 100000, np.percentile(nCutList, 5), np.percentile(nCutList,95)
        vlo, vhi = v5-5*abs(v5), v95+5*abs(v95)
        nCutListReduced = [x for x in nCutList if x > v5 and x < v95]
        outPlot = "../plots/%s_ds%d_idx%d_%s_ch%d.png" % (parName,dsNum,subNum,tName,ch)
        cut99,cut95,cut01,cut05,cut90 = MakeCutPlot(c,cal,par,eb,elo,ehi,vb,vlo,vhi,d2Cut,d1Cut,outPlot,fastMode)
        cutDict[ch] = [cut01,cut05,cut90,cut95,cut99]
    return cutDict


def MakeCutPlot(c,cal,var,eb,elo,ehi,vb,vlo,vhi,d2Cut,d1Cut,outPlot,fastMode):
    """ Creates a channel-specific energy calibration plot. """

    # Calculate cut vals (assumes plot range is correct)
    h1 = wl.H1D(cal,vb,vlo,vhi,var,d1Cut)
    h1Sum = h1.Integral()
    if h1Sum == 0:
        print("Error: Failed %s, histogram sum is 0 so cannot normalize, setting to [0,0,0,0,0]"%(var))
        return 0,0,0,0,0
    h1.Scale(1/h1Sum)
    try:
        cut99,cut95,cut01,cut05,cut90 = wl.GetIntegralPoints(h1)
    except:
        print("Error: Failed %s using cut %s, setting to [0,0,0,0,0]"%(var,d1Cut))
        return 0,0,0,0,0
    if fastMode:
        print("Returning fastMode output: ", cut99,cut95,cut01,cut05,cut90)
        return cut99,cut95,cut01,cut05,cut90

    # Generate the plot for inspection.
    c.cd(2)
    gPad.SetLogy(0)
    h1.GetXaxis().SetRangeUser(cut01-abs(0.25*cut01), cut99 + abs(0.25*cut99) )
    h1.SetTitle("")
    h1.GetXaxis().SetTitle(var)
    h1.Draw("hist")

    c.cd(1)
    gPad.SetLogy(0)
    cal.Draw("%s:trapENFCal>>b(%d,%d,%d,%d,%.3E,%.3E)"%(var,eb+10,elo-5,ehi+5,vb,cut01-abs(0.25*cut01),cut99+abs(0.25*cut99)) ,d2Cut)

    l1, l2, l3 = TLine(), TLine(), TLine()
    l1.SetLineColor(ROOT.kGreen)
    l2.SetLineColor(ROOT.kRed)
    l3.SetLineColor(ROOT.kMagenta)

    l1.DrawLine(elo-5, cut99, ehi+5, cut99)
    l2.DrawLine(elo-5, cut95, ehi+5, cut95)
    l2.DrawLine(elo-5, cut05, ehi+5, cut05)
    l1.DrawLine(elo-5, cut01, ehi+5, cut01)

    c.cd(3)
    x_h1, y_h1 = wl.npTH1D(h1)
    int_h1 = wl.integFunc(y_h1)
    g2 = TGraph(len(x_h1), x_h1, int_h1)
    g2.GetXaxis().SetRangeUser(cut01-abs(0.3*cut01), cut99 + abs(0.3*cut99) )
    g2.SetTitle("")
    g2.GetXaxis().SetTitle(var)
    g2.GetYaxis().SetTitle("Percentile")
    g2.Draw("ACP")
    l1.DrawLine(cut99, 0, cut99, 1)
    l2.DrawLine(cut95, 0, cut95, 1)
    l1.DrawLine(cut01, 0, cut01, 1)
    l2.DrawLine(cut05, 0, cut05, 1)

    c.Print(outPlot)
    return cut99,cut95,cut01,cut05,cut90


def getThreshold():
    """ get a threshold value for a particular run/channel.
    threshold - uses bkgIdx
        keys: thresh_ds[i]_bkgidx[j]
        vals: {[chan]:[50 pct mean, sigma]}
    """
    run, chan = 7244, 624

    # calDB = db.TinyDB('../calDB.json')
    # pars = db.Query()

    # this isn't gonna work for ext pulser runs, they're not in a bkgIdx ...
    # thD = ds.getDBRecord("thresh_ds%d_bkgidx%d" % (dsNum, bkgIdx), False, calDB, pars)
    # chThreshList = (ch for ch in chList if ch in thD.keys())
    # for ch in chThreshList:
    #     mu, sig = thD[ch][0], thD[ch][1]
    #     if mu > 5 or mu==0:
    #         print "Threshold for ch%d bkgidx%d zero or too high, skipping" % (ch, bkgIdx)
    #         continue
    #     yThresh = 0.5 * (1 + sp.erf( (xThresh - mu) / (2**(.5) * abs(sig))))
    #     yThresh *= expDict[ch][bkgIdx]
    #     yThreshTot += yThresh


if __name__=="__main__":
    main()