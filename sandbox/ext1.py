#!/usr/bin/env python3
import sys, os, imp, glob
sys.argv.append("-b")
ds = imp.load_source('DataSetInfo',os.environ['LATDIR']+'/DataSetInfo.py')
wl = imp.load_source('waveLibs',os.environ['LATDIR']+'/waveLibs.py')
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np
import tinydb as db

extPulserInfo = {
        # TODO: relabel special runs list and DataSetInfo to match this list
        # Test 1 - attenuation, --donotbuild was used
        6:  [[5942, 5945], [10,14,18,22], 200, 674], # att, rt, chan
        7:  [[5947, 5960], [10,14,18,22,26,30,34,38,42,46,50,54,58,62], 190, 674],
        8:  [[5964, 5977], [10,14,18,22,26,30,34,38,42,46,50,54,58,62], 190, 624],
        9:  [[5979, 5992], [10,14,18,22,26,30,34,38,42,46,50,54,58,62], 190, 688],
        10: [[6191, 6204], [10,14,18,22,26,30,34,38,42,46,50,54,58,62], 190, 662],
        11: [[6206, 6219], [10,14,18,22,26,30,34,38,42,46,50,54,58,62], 190, 608],
        # Test 2 - rise time
        12: [[6934, 6944], [140,145,150,155,160,165,170,175,180,185,190], 18, 674], # rt, att, chan
        13: [[6964, 6970], [4354,1257,1296,654,1278,1278,1278],0,[674,624,688,662,608,608,608]], # adc,att,chan
        14: [[6971, 6976], [140,150,160,170,180,190], 18, 614],
        15: [[6977, 6982], [140,150,160,170,180,190], 18, 624],
        16: [[7002, 7007], [140,150,160,170,180,190], 18, 688],
        17: [[7008, 7013], [140,150,160,170,180,190], 18, 662],
        # Test 3 - attenuation
        18: [[7219, 7233], [14,18,22,26,30,999,30,34,38,42,46,50,54,58,62], 155, 674], # att, rt, chan
        19: [[7234, 7246], [14,18,22,26,30,34,38,42,46,50,54,58,62], 164, 624],
        20: [[7247, 7259], [14,18,22,26,30,34,38,42,46,50,54,58,62], 146, 688],
        21: [[7260, 7272], [14,18,22,26,30,34,38,42,46,50,54,58,62], 138, 662],
        22: [[13168, 13181], [14,18,22,26,30,34,38,42,46,50,54,58,62], 999, 690]
    }

originalList = {
    # 0: [4547, 4547], # ignore
    # 1: [4549, 4572], # ignore, whole BG range, full of regular pulsers, etc
    # 2: [4573, 4831], # ignore, whole BG range, full of regular pulsers, etc
    # 3: [5525, 5534], # ignore, "setup system"
    # 4: [5535, 5554], # ignore, "numerous problems"
    # 5: [5555, 5850], # ignore, whole BG range, full of regular pulsers, etc
    # 6: [5872, 5877], # ignore, short runs
    7: [5940, 5963],
    8: [5964, 5978],
    9: [5979, 5992],
    10: [6191, 6205],
    11: [6206, 6219],
    12: [6934, 6944],
    13: [6964, 6970],
    14: [6971, 6976],
    15: [6977, 6982],
    16: [7002, 7007],
    17: [7008, 7013],
    18: [7219, 7233],
    19: [7234, 7246],
    20: [7247, 7259],
    21: [7260, 7272],
    22: [13168, 13181]
    }


def main():

    # runTuneCut()
    runByRun()
    # fitSloEfficiency1()
    # fitSloEfficiency2()
    # riseTimeStudy()
    # combineData()
    # riseNoiseEfficiency()
    # TEEfficiency()
    # checkBadRuns()
    # createBadRunList()
    # checkFiles()


def hist2dExample():
    xArr, yArr = sloArr, rtArr
    xLo, xHi, bpX, yLo, yHi, bpY = 50, 100, 0.2, 335, 375, 0.2
    nBY, nBX = int((yHi-yLo)/bpY), int((xHi-xLo)/bpY)
    plt.hist2d(xArr, yArr, bins=[nBX,nBY], range=[[xLo,xHi],[yLo,yHi]], norm=LogNorm())
    plt.colorbar()
    plt.xlim(50, 100)
    plt.title("pIdx %d, channel %d" % (pIdx, extChan))
    plt.xlabel("fitSlo",horizontalalignment='right',x=1.0)
    plt.ylabel("riseTime (ns)",horizontalalignment='right',y=1.0)
    plt.legend(loc="best")
    plt.savefig("../plots/rtStudy_idx%d.pdf" % (pIdx))
    # plt.show()


def checkFiles():

    calInfo = ds.CalInfo()
    extPulserInfo = calInfo.GetSpecialList()["extPulserInfo"]

    for pIdx in range(7,23+1):
        runList = calInfo.GetSpecialRuns("extPulser",pIdx)
        attList = extPulserInfo[pIdx][0]
        extChan = extPulserInfo[pIdx][-1]

        for i, run in enumerate(runList):
            fileList = ds.getLATRunList([run],"%s/lat" % (ds.specialDir))

            if len(fileList)==0:
                attRun = extPulserInfo[pIdx][0][i]
                print("No files:",run,"Att:",attRun)

    # result:
    noFiles = [6936,6937,6940,6942,6944,6965,6968,6969,6974,6977,7224,7267,7268,13168]


def createBadRunList():

    origRunList = []
    for key in originalList:
        runLo, runHi = originalList[key][0], originalList[key][1]
        origRunList.extend([run for run in range(runLo, runHi+1)])

    newRunList = []
    for key in extPulserInfo:
        runLo, runHi = extPulserInfo[key][0][0], extPulserInfo[key][0][1]
        newRunList.extend([run for run in range(runLo, runHi+1)])

    print(len(origRunList), len(newRunList))

    for run in origRunList:
        if run not in newRunList:
            print(run)


def checkBadRuns():
    from ROOT import TChain

    cal = ds.CalInfo()
    for pIdx in range(14,22+1):

        extChan = extPulserInfo[pIdx][-1]
        syncChan = wl.getChan(0,10,0) # 672

        runList = cal.GetSpecialRuns("extPulser",pIdx)
        fList = ds.getLATRunList(runList, ds.specialDir+"/lat")

        for f in fList:
            print(f)
            latChain = TChain("skimTree")
            latChain.Add("%s/%s" % (ds.specialDir+"/lat", f))

            tNames = ["Entry$","mH","channel","trapENFCal","fitSlo","den90","den10"]
            theCut = "(channel==%d || channel==%d) && trapENFCal > 1" % (syncChan, extChan)
            tVals = wl.GetVX(latChain,tNames)
            nPass = len(tVals["Entry$"])
            if nPass == 0: continue

            for idx in range(nPass):
                chan = tVals["channel"][idx]
                enf = tVals["trapENFCal"][idx]
                rt = tVals["den90"][idx]-tVals["den10"][idx]

            print(pIdx,latChain.GetEntries(),nPass)


    skipList = []


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


def runByRun():
    """ Directly confirm settings of ext pulser scripts. """
    import time
    from ROOT import TFile, TChain, GATDataSet, gROOT
    gROOT.ProcessLine("gErrorIgnoreLevel = 3001;") # suppress ROOT error messages

    extPDict = {7:674, 8:624, 9:688, 10:662, 11:608, 12:674, 14:608, 15:624, 16:688, 17:662, 18:674, 19:624, 20:688, 21:662, 22:690}
    syncChan = wl.getChan(0,10,0) # 672

    syncRateNominal = 20 # Hz

    calInfo = ds.CalInfo()

    fig = plt.figure(figsize=(10,6),facecolor='w')

    # pIdxs = [12,14,15,16,17] # test 2 - rise time
    # pIdxs = [18,19,20,21] # test 3 - attenuation
    pIdxs = [22]
    for pIdx in pIdxs:
        runList = calInfo.GetSpecialRuns("extPulser",pIdx)
        print("Range",pIdx)

        extChan = extPDict[pIdx]
        xArr, yArr = [], [] # we're gonna plot these

        # runList = [7234]
        for run in runList:
            # if run in [6936,6937,6940,6942,6944, 6974, 6977]: continue # test 2
            # if run in [7224] or run > 7266: continue # test 3

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

            tNames = ["Entry$","run","channel","mH","trapENFCal","den90","den10","fitSlo","localTime_s","tOffset","fitAmp"]
            # theCut = "(channel==%d || channel==%d) && mH==2" % (syncChan,extChan)
            # theCut += " && Entry$ < 100"
            theCut = "Entry$ < 200 && gain==0"
            tVals = wl.GetVX(latChain,tNames,theCut)

            # don't delete this
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

            continue

            # make sure we only have hits from syncChan and extChan
            # for entry in set(tVals["Entry$"]):
            #     idxs = [idx for idx in range(len(tVals["Entry$"])) if tVals["Entry$"][idx]==entry]
            #     chans = [tVals["channel"][idx] for idx in idxs]
            #     if not set([extChan,syncChan]).issubset(set(chans)):
            #         print("NOPE:",chans)

            gds = GATDataSet(int(run))
            runTime = gds.GetRunTime()/1e9
            if len(tVals["Entry$"])==0:
                print("Run %d, %.2f sec. Found no cts." % (run,runTime))
                continue

            syncRate = len(set(tVals["Entry$"]))/runTime
            expectedCts = runTime * syncRateNominal
            extPCts = len([ch for ch in tVals["channel"] if ch==extChan])
            syncCts = len([ch for ch in tVals["channel"] if ch==syncChan])
            extPRate = extPCts/runTime
            syncRate = syncCts/runTime

            syncAmp = [tVals["fitAmp"][i] for i in range(len(tVals["fitAmp"])) if tVals["channel"][i]==syncChan]
            syncAmp = np.asarray(syncAmp)
            muS, sigS = 0, 0
            if len(syncAmp)>0:
                muS, sigS = np.mean(syncAmp), np.std(syncAmp)

            extENF = [tVals["trapENFCal"][i] for i in range(len(tVals["trapENFCal"])) if tVals["channel"][i]==extChan]
            extENF = np.asarray(extENF)
            muE, sigE = 0, 0
            if len(extENF)>0:
                muE, sigE = np.mean(extENF), np.std(extENF)

            print("Run %d, %.2f sec.  #Expect %d  #Sync %d  (%.2f Hz)  #extP %d (%.2f Hz)  muE %.2f  sigE %.2f  muS %.2f  sigS %.2f" % (run,runTime,expectedCts,syncCts,syncRate,extPCts,extPRate,muE,sigE,muS,sigS))

            # fill the plot arrays
            xArr.extend([tVals["trapENFCal"][i] for i in range(len(tVals["trapENFCal"])) if tVals["channel"][i]==extChan])
            yArr.extend([tVals["fitSlo"][i] for i in range(len(tVals["fitSlo"])) if tVals["channel"][i]==extChan])

        return

        # make a plot for this range
        fig.clear()
        xLo, xHi, yLo, yHi = 0, 10, -20, 300 # test 3
        # xLo, xHi, yLo, yHi = 50, 100, -20, 200 # test 2
        bpY, bpX = 2, 0.1
        nBinsY, nBinsX = int((yHi-yLo)/bpY), int((xHi-xLo)/bpX)
        try:
            plt.hist2d(xArr, yArr, bins=[nBinsX,nBinsY], range=[[xLo,xHi],[yLo,yHi]], norm=LogNorm())
            plt.colorbar()
            plt.xlabel("trapENFCal (keV)",horizontalalignment='right',x=1.0)
            plt.ylabel("fitSlo",horizontalalignment='right',y=1.0)
            plt.title("Range %d, Channel %d" % (pIdx, extChan))
            plt.tight_layout()
            plt.savefig("../plots/extPulser_idx%d.pdf" % pIdx)
        except ValueError:
            pass


def fitSloEfficiency1():
    """ Plot the efficiency vs. energy """
    from ROOT import TChain, gROOT
    gROOT.ProcessLine("gErrorIgnoreLevel = 3001;")

    # 20: [7247, 7259], (run) P42661C P1D3 688
    pIdx, extChan, syncChan = 20, 688, 672
    calInfo = ds.CalInfo()
    runList = calInfo.GetSpecialRuns("extPulser",pIdx)

    # get fitSlo value from the DB for this channel.
    # calIdx's
    # 33: [[6904,6925],6904,7274],
    # 34: [[7275,7279,7281,7295],7275,7614]
    dsNum, modNum, calIdx = 0, 1, 33
    calDB = db.TinyDB('../calDB.json')
    pars = db.Query()
    fsD = ds.getDBRecord("fitSlo_ds%d_idx%d_m%d_Peak" % (dsNum, calIdx, modNum), False, calDB, pars)
    fsCut = fsD[extChan][2] # 90% value (used in LAT3)

    xArr, yArr = [], []
    for run in runList:
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
        tNames = ["Entry$","run","channel","mH","trapENFCal","fitSlo"]
        theCut = "(channel==%d || channel==%d) && mH==2" % (syncChan,extChan)
        tVals = wl.GetVX(latChain,tNames,theCut)
        nPass = len(tVals["Entry$"])
        if nPass == 0: continue

        xArr.extend([tVals["trapENFCal"][i] for i in range(nPass) if tVals["channel"][i]==extChan])
        yArr.extend([tVals["fitSlo"][i] for i in range(nPass) if tVals["channel"][i]==extChan])

    # calculate efficiency
    xLo, xHi, bpX = 0, 10, 0.2
    xAcc = np.arange(xLo, xHi, bpX)
    xArr, yArr = np.asarray(xArr), np.asarray(yArr)
    yEff = []
    errBars = []
    for i in range(int((xHi-xLo)/bpX),0,-1):
        eHi, eLo = i*bpX, (i-1)*bpX
        idx = np.where((xArr > eLo) & (xArr < eHi))
        ySlice = yArr[idx]
        idx2 = np.where(ySlice < fsCut)
        nTot, nPass = len(ySlice), len(ySlice[idx2])

        if len(yEff)==0 and nTot==0:
            eff = 100
        elif nTot == 0:
            eff = yEff[-1]
        else:
            eff = 100 * nPass/nTot
            # err = sm.stats.proportion.proportion_confint(nPass, nTot, alpha=0.05, method='beta')
        print("%d  %.1f-%.1f keV  %.2f%%" % (i, eLo, eHi, eff))
        yEff.append(eff)
    yEff.reverse()

    # make a plot
    fig = plt.figure(figsize=(10,6),facecolor='w')
    p1 = plt.subplot(111)

    xLo, xHi, bpX = 0, 10, 0.2
    yLo, yHi, bpY = -20, 1000, 2.
    nBinsY, nBinsX = int((yHi-yLo)/bpY), int((xHi-xLo)/bpX)
    h = p1.hist2d(xArr, yArr, bins=[nBinsX,nBinsY], range=[[xLo,xHi],[yLo,yHi]], norm=LogNorm())

    p1.axhline(fsCut, color='black', linewidth=3)
    p1.set_xlabel("trapENFCal (keV)",horizontalalignment='right',x=1.0)
    p1.set_ylabel("fitSlo",horizontalalignment='right',y=1.0)
    p1.set_title("Range %d, Channel %d.  fitSlo cut > %.2f" % (pIdx, extChan, fsCut))
    # plt.colorbar(h[3], ax=p1)

    p2 = p1.twinx()
    p2.plot(xAcc, yEff, ls='steps-post',color='red',linewidth=2.)

    p2.set_ylim(0,110)
    p2.set_ylabel('% Efficiency', color='r', horizontalalignment='right',y=1.0)
    p2.tick_params('y', colors='r')

    plt.tight_layout()
    plt.savefig("../plots/efficiency_idx%d.pdf" % pIdx)


def fitSloEfficiency2():
    """ Plot efficiency vs energy for the 3 lowest ext pulser sets. """

    from ROOT import TChain, gROOT
    gROOT.ProcessLine("gErrorIgnoreLevel = 3001;")

    pIdxs = [18,19,20]

    for pIdx in pIdxs:

        avgE = {}
        effic = {}

        # pIdx 20 settings, taken from ORCA logs (run database)
        pIdx, extChan, syncChan = 20, 688, 672
        riseTime = 146 # ns
        attens = {7247:14, 7248:18, 7249:22, 7250:26, 7251:30, 7252:34, 7253:38, 7254:42, 7255:46, 7256:50, 7257:54, 7258:58, 7259:62}

        # get fitSlo value from the DB for this channel.
        dsNum, modNum, calIdx = 0, 1, 33
        calDB = db.TinyDB('../calDB.json')
        pars = db.Query()
        fsD = ds.getDBRecord("fitSlo_ds%d_idx%d_m%d_Peak" % (dsNum, calIdx, modNum), False, calDB, pars)
        fsCut = fsD[extChan][2] # 90% value (used in LAT3)

        fig = plt.figure(figsize=(8,8),facecolor='w')
        p1 = plt.subplot(111)

        calInfo = ds.CalInfo()
        runList = calInfo.GetSpecialRuns("extPulser",pIdx)
        for run in runList:

            fileList = ds.getLATRunList([run],"%s/lat" % (ds.specialDir))
            latChain = TChain("skimTree")
            for f in fileList: latChain.Add("%s/lat/%s" % (ds.specialDir,f))
            tNames = ["Entry$","channel","trapENFCal","fitSlo"]
            theCut = "(channel==%d || channel==%d) && mH==2" % (syncChan, extChan)
            tVals = wl.GetVX(latChain,tNames,theCut)
            nPass = len(tVals["Entry$"])
            if nPass == 0: continue
            enfArr = [tVals["trapENFCal"][i] for i in range(nPass) if tVals["channel"][i]==extChan]
            sloArr = [tVals["fitSlo"][i] for i in range(nPass) if tVals["channel"][i]==extChan]
            enfArr, sloArr = np.asarray(enfArr), np.asarray(sloArr)

            muE, stdE = np.mean(enfArr), np.std(enfArr)
            muF, stdF = np.mean(sloArr), np.std(sloArr)
            nTot = len(sloArr)
            nAcc = len([fs for fs in sloArr if fs < fsCut])
            eff = 100. * (nAcc/nTot)

            print("Run %d  Att %d  muE %-5.2f  stdE %-5.2f  muF %-5.2f  stdF %-5.2f  nTot %-5d  nAcc %-5d  eff %.8f" % (run,attens[run],muE,stdE,muF,stdF,nTot,nAcc,eff))

            avgE[run] = muE
            effic[run] = eff

            p1.cla()
            xLo, xHi, bpX = muE - stdE*3, muE + stdE*3, 0.2
            yLo, yHi, bpY = -10, 300, 0.1
            nBinsY, nBinsX = int((yHi-yLo)/bpY), int((xHi-xLo)/bpX)
            p1.hist2d(enfArr, sloArr, bins=[nBinsX,nBinsY], range=[[xLo,xHi],[yLo,yHi]], norm=LogNorm())
            p1.axhline(fsCut, color='black', linewidth=3)
            p1.set_xlabel("trapENFCal (keV)",horizontalalignment='right',x=1.0)
            p1.set_ylabel("fitSlo",horizontalalignment='right',y=1.0)
            p1.set_title("run %d  ch %d  fsCut %.2f  eff %.2f" % (run,extChan,fsCut,eff))
            plt.tight_layout()
            plt.savefig("../plots/efficiency_idx%d_run%d.pdf" % (pIdx, run))

    # make arrays for plotting
    yAtt, yEff, xEne = [], [], []
    for run in sorted(attens):
        if run in avgE.keys() and run in effic.keys():
            yAtt.append(attens[run])
            yEff.append(effic[run])
            xEne.append(avgE[run])
            print(run,attens[run],avgE[run],effic[run])

    # make plot
    fig.clear()
    p1 = plt.subplot(211)

    p1.plot(xEne, yAtt, ".")
    p1.set_title("Att. vs Energy, m %.2f  b %.2f" % (m,b))
    p1.set_xlabel("Energy (keV)")
    p1.set_ylabel("Atten (db)")

    p2 = plt.subplot(212)
    p2.plot(xEne, yEff, ".")
    p2.set_title("Efficiency vs Energy")
    p2.set_xlabel("Energy (keV)")
    p2.set_ylabel("Efficiency")
    p2.set_ylim(0,110 )

    plt.tight_layout()
    plt.savefig("../plots/attVsE.pdf")


def riseTimeStudy():
    """ How much does fitSlo jump around? """
    calInfo = ds.CalInfo()
    runList = calInfo.GetSpecialRuns("extPulser",pIdx)
    for run in runList:
        fileList = ds.getLATRunList([run],"%s/lat" % (ds.specialDir))



if __name__=="__main__":
    main()