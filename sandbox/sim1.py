#!/usr/bin/env python3
import os, math, glob, imp, sys
import numpy as np

import matplotlib as mpl
mpl.use('Agg')
sys.argv.append("-b")
import matplotlib.pyplot as plt
# plt.style.use('../pltReports.mplstyle')
# import seaborn as sns

# load LAT libraries
ds = imp.load_source('DataSetInfo',os.environ['LATDIR']+'/sandbox/DataSetInfo.py')
wl = imp.load_source('waveLibs',os.environ['LATDIR']+'/waveLibs.py')
cal = ds.CalInfo()

def main():

    # getSimData()
    # plotSimData()
    # getSimData2()
    # getThresh()

    # compare1DHistos()
    # plotDLSimSpec()
    plotResiduals()
    # testSimData()
    # testNewSimData()


def getSimData():
    from ROOT import TChain, TFile

    # inDir = "/global/projecta/projectdirs/majorana/users/mbuuck/sim/MJDWithGrahamTDL_byDetESmearing/5.0_TDL/0.75_transition_point/0.50_transition_level/MJDemonstrator/linesource/M1CalSource/A224_Z88"
    inDir = os.environ['SLURM_TMP']

    fileList = sorted(glob.glob("%s/*.root" % inDir))
    # fileList = fileList[:10]
    fileList = fileList[:]

    evtTotal, evtBulk, evtTrans = [], [], []

    for i, f in enumerate(fileList):
        print("%d/%d %s" % (i,len(fileList),f))
        tf = TFile(f)
        simTree = tf.Get("simTree")

        theCut = "fNWaveforms==2 && fTotalEnergy/fActiveness > 0.237 && fTotalEnergy/fActiveness < 0.24"

        n1 = simTree.Draw('fEnergy*1000/fActiveness', theCut, 'goff')
        eTot = simTree.GetV1()
        evtTotal.extend([eTot[i] for i in range(n1)])

        n2 = simTree.Draw('fEnergy*1000/fActiveness', theCut + '&& fActiveness == 1', 'goff')
        eBulk = simTree.GetV1()
        evtBulk.extend([eBulk[i] for i in range(n2)])

        n3 = simTree.Draw('fEnergy*1000/fActiveness', theCut + '&& fActiveness < 1', 'goff')
        eTrans = simTree.GetV1()
        evtTrans.extend([eTrans[i] for i in range(n3)])

        tf.Close()

    np.savez("../data/sim1-evtTrans.npz",evtTotal,evtBulk,evtTrans)


def plotSimData():

    f = np.load("../data/sim1-evtTrans.npz")
    evtTotal, evtBulk, evtTrans = f['arr_0'], f['arr_1'], f['arr_2']

    fig = plt.figure()

    xLo, xHi, xpb = 0, 3000, 1

    x, hTotal = wl.GetHisto(evtTotal, xLo, xHi, xpb)
    x, hBulk = wl.GetHisto(evtBulk, xLo, xHi, xpb)
    x, hTrans = wl.GetHisto(evtTrans, xLo, xHi, xpb)

    plt.semilogy(x, hTotal, ls='steps', c='r', lw=2., label='total')
    plt.semilogy(x, hBulk, ls='steps', c='g', lw=2., label='bulk')
    plt.semilogy(x, hTrans, ls='steps', c='b', lw=2., label='transition')

    plt.ylim(1., 8 * max(hTotal))
    plt.xlabel("fEnergy (keV)", ha='right', x=1.)
    plt.ylabel("Simulated Counts", ha='right', y=1.)
    plt.legend(loc=1)
    plt.savefig("../plots/sim1-transSpec.png")

    plt.cla()

    transPct = 100 * (np.divide(hTrans, hTotal, dtype=float))
    sns.regplot(x=x, y = transPct, scatter_kws={'s':20})

    plt.xlabel("fEnergy (keV)", ha='right', x=1.)
    plt.ylabel("% of Transition Layer Events", ha='right', y=1.)

    plt.savefig("../plots/sim1-transPct.png")


def getSimData2():
    from ROOT import TFile, TTree

    # config, module = "DS5", "M1"
    # basePath = "/global/projecta/projectdirs/majorana/sim/MJDG41003GAT/"
    # sourceType, partClass, segment = "linesource", "%sCalSource" % module, "A224_Z88"
    # simPath = "%s/MJDemonstrator/%s/%s/%s" % (basePath,sourceType,partClass,segment)
    # simInfo = ds.SimInfo(config)
    # detList = simInfo.GetActiveDets(config, module)

    # inDir = "/global/projecta/projectdirs/majorana/users/mbuuck/sim/MJDWithGrahamTDL_byDetESmearing/5.0_TDL/0.75_transition_point/0.50_transition_level/MJDemonstrator/linesource/M1CalSource/A224_Z88"

    simFile = "~/project/sims/MJDemonstrator_linesource_A224_Z88_from_A224_Z88_to_A208_Z81_in_M1CalSource_500000_-1089829.root"

    procFile = "/global/projecta/projectdirs/majorana/users/mbuuck/sim/MJDWithGrahamTDL_byDetESmearing/5.0_TDL/0.75_transition_point/0.50_transition_level/MJDemonstrator/linesource/M1CalSource/A224_Z88/DS3processed_MJDemonstrator_linesource_A224_Z88_from_A224_Z88_to_A208_Z81_in_M1CalSource_500000_-1089829.root"

    sf = TFile(simFile)

    simTree = sf.Get('fTree')

    pf = TFile(procFile)
    procTree = pf.Get('simTree')

    print(simTree.GetEntries(), procTree.GetEntries())

    # branch: eventSteps
    # fSteps.fProcessName
    # fSteps.fPhysVolName
    # fSteps.fEdep
    # fSteps.fKineticEnergy
    # fSteps.fX, fY, fZ
    # fSteps.fLocalX
    # fSteps.fStepNumber

    # branch: eventPrimaries


def getThresh():

    import tinydb as db
    runList = cal.GetSpecialRuns("longCal",2) # ds3 long cal run, 17183 - 17302
    dsNum, bkgIdx = 3, 3 # 3: [[17183,17203],17183,17422],

    calDB = db.TinyDB('../calDB.json')
    pars = db.Query()
    thD = ds.getDBRecord("thresh_ds%d_bkgidx%d" % (dsNum, bkgIdx), False, calDB, pars)

    for chan in thD.keys():
        mu, sig = thD[chan][0], thD[chan][1]
        print(chan,mu,sig,mu+3*sig)


def compare1DHistos():
    """ Result: method 2 is preferable.
    Call wl.npTH1D with the 'ctr' option to get the bins centered
    Then if you want to plot a line, you have to shift it down
    by half the bin width.
    """
    from ROOT import TFile, TH1D, TCanvas

    f1 = TFile("../data/hists_1.0_TDL_0.80_transition_point_0.20_transition_level.root")
    sim = ds.SimInfo('DS3')
    detList = sim.GetActiveDets('DS3','M1')
    for det in detList:
        tSim = f1.Get("hSim_%s" % det)
        # tData = f1.Get("hData_%s" % det)
        # tResid = f1.Get("hResid_%s" % det)
        # tSE = f1.Get("hSE_%s" % det)

        # method 1
        nb = tSim.GetNbinsX()
        xLo, xHi = tSim.GetXaxis().GetXmin(), tSim.GetXaxis().GetXmax()
        xpb = (xHi-xLo)/nb # kev per bin
        hSim = np.asarray([tSim.GetBinContent(i) for i in range(nb)])
        x = np.arange(xLo, xHi, xpb)

        # method 2
        x1, y1, xpb1 = wl.npTH1D(tSim,'ctr')

        # x, hData, xpb = wl.npTH1D(tData,'ctr')
        # if (len(x)!=len(x1)): print("lengths don't match")
        # if (x[0]!=x1[0] or x[-1]!=x1[-1]): print("limits don't match")

        break

    idx = np.where((x>2612) & (x<2618))
    xMax = np.argmax(hSim[idx])
    eMax = x[idx][xMax]

    idx = np.where((x1>2612) & (x1<2618))
    xMax = np.argmax(y1[idx])
    eMax2 = x1[idx][xMax] - xpb1/2.

    print("np:",eMax, eMax2)

    fig = plt.figure()
    # plt.plot(x, hSim, c='b', ls='steps')
    # plt.plot(x1+xpb1/2., y1, c='r', ls='steps')
    plt.plot(x1, y1, c='r', ls='steps')
    plt.xlim(2612,2618)
    plt.ylim(0,0.0005)
    # plt.axvline(eMax,c='b')
    plt.axvline(eMax2,c='r')
    plt.savefig("../plots/sim1-deadSpec.png")

    c = TCanvas("c","c",800,600)
    tSim.Draw("hist")
    tSim.GetXaxis().SetRangeUser(2612,2618)
    print("root:",tSim.GetXaxis().GetBinCenter(tSim.GetMaximumBin()))
    c.Print("../plots/sim1-deadTh1.png")


def plotDLSimSpec():
    from ROOT import TFile, TH1D, TCanvas

    # f1 = TFile("../data/hists_1.0_TDL_0.80_transition_point_0.20_transition_level.root")
    # f1 = TFile("../data/hists_3.0_TDL_0.95_transition_point_0.40_transition_level.root")
    f1 = TFile("../data/hists_5.0_TDL_0.90_transition_point_0.25_transition_level.root")

    tSim = f1.Get("hSim_All")
    tData = f1.Get("hData_All")
    x, hSim, xpb = wl.npTH1D(tSim, 'ctr')
    x, hData, xpb = wl.npTH1D(tData, 'ctr')
    hRes = hSim - hData

    fig = plt.figure(figsize=(9,8))
    p1 = plt.subplot2grid((3,1),(0,0), rowspan=2)
    p2 = plt.subplot2grid((3,1),(2,0), sharex=p1)

    p1.semilogy(x, hSim, ls='steps', lw=1., c='r', label='sim')
    p1.semilogy(x, hData, ls='steps', lw=1., c='b', label='data')
    p1.legend(loc=1)

    p2.plot(x, hRes, ls='steps', lw=1, c='b', label='residual, sim-data')
    p2.legend(loc=1)
    p2.set_xlabel("Energy (keV)", ha='right', x=1.)

    plt.savefig("../plots/sim1-deadSpec.png")


def plotResiduals():
    """ sim 2 is the really good fit """
    from ROOT import TFile

    # Most likely case (best fit with 1.0x DL thickness): /global/projecta/projectdirs/majorana/users/mbuuck/sim/sim2data/201802/hists_1.0_TDL_0.90_transition_point_0.30_transition_level.root (edited)
    # potential best case (best fit with floating DL thickness): /global/projecta/projectdirs/majorana/users/mbuuck/sim/sim2data/201802/hists_3.0_TDL_0.95_transition_point_0.40_transition_level.root
    # worst fit (restricting to 1.0x DL thickness): /global/projecta/projectdirs/majorana/users/mbuuck/sim/sim2data/201802/hists_1.0_TDL_0.80_transition_point_0.20_transition_level.root
    # worst fit (floating DL thickness): /global/projecta/projectdirs/majorana/users/mbuuck/sim/sim2data/201802/hists_5.0_TDL_0.90_transition_point_0.25_transition_level.root

    f1 = TFile("../data/hists_1.0_TDL_0.80_transition_point_0.20_transition_level.root")
    tSim, tData = f1.Get("hSim_All"), f1.Get("hData_All")
    x,hSim1,_ = wl.npTH1D(tSim,'ctr')
    h,hData,_ = wl.npTH1D(tData,'ctr')
    hRes1 = hSim1 - hData

    # this is the good one
    f2 = TFile("../data/hists_3.0_TDL_0.95_transition_point_0.40_transition_level.root")
    tSim, tData = f2.Get("hSim_All"), f2.Get("hData_All")
    x,hSim2,_ = wl.npTH1D(tSim,'ctr')
    hRes2 = hSim2 - hData

    # brian asked for nCounts from 0-5 kev between sim1 and sim2

    idx = np.where((x > 7) & (x < 10))
    nCts1 = sum(hSim1[idx])
    nCts2 = sum(hSim2[idx])
    nCtsD = sum(hData[idx])
    print('nCts1',nCts1,'nCts2',nCts2,'nData',nCtsD)

    f3 = TFile("../data/hists_5.0_TDL_0.90_transition_point_0.25_transition_level.root")
    tSim, tData = f3.Get("hSim_All"), f3.Get("hData_All")
    x,hSim3,_ = wl.npTH1D(tSim,'ctr')
    hRes3 = hSim3 - hData

    fig = plt.figure(figsize=(9,8))
    p1 = plt.subplot2grid((3,1),(0,0), rowspan=2)
    p2 = plt.subplot2grid((3,1),(2,0), sharex=p1)

    p1.semilogy(x, hSim3, ls='steps', lw=1, c='m', label='sim 3')
    p1.semilogy(x, hSim1, ls='steps', lw=1, c='r', label='sim 1')
    p1.semilogy(x, hSim2, ls='steps', lw=1, c='g', label='sim 2')
    p1.semilogy(x, hData, ls='steps', lw=5, c='b', alpha=0.5, label='data')
    p1.set_xlabel("Energy (keV)", ha='right', x=1.)
    p1.set_xlim(0,250)
    p1.set_ylim(5e-5,5e-4)
    p1.legend()

    # p2.plot(x, hRes1, ls='steps', lw=1, c='r', label='sim 1')
    p2.plot(x, hRes2, ls='steps', lw=1, c='g', label='sim 2')
    # p2.plot(x, hRes3, ls='steps', lw=1, c='m', label='sim 3')
    p2.set_xlabel("Energy (keV)", ha='right', x=1.)
    p2.set_xlim(0,250)
    p2.legend()

    plt.tight_layout()
    plt.savefig("../plots/sim1-residuals.png")


def testSimData():
    from ROOT import TChain, TFile
    inDir = os.environ['SLURM_TMP']
    fileList = sorted(glob.glob("%s/DS3processed*.root" % inDir))
    fileList = fileList[:40]
    # fileList = fileList[:]

    hits = {i:[] for i in range(7)}
    actList = []

    totESpec = []
    sumHitESpec = []

    nPk, nCont = 0,0

    for i, f in enumerate(fileList):
        # print("%d/%d %s" % (i,len(fileList),f))
        print("%d/%d" % (i,len(fileList)))
        tf = TFile(f)
        simTree = tf.Get("simTree")

        for iEnt in range(simTree.GetEntries()):
            simTree.GetEntry(iEnt)

            ae = simTree.fAnalysisEvent
            sumE = ae.GetTotalEnergy()

            mH = ae.GetNElements()

            if 0.2 < sumE < 0.25 and mH==2:
                totESpec.append(sumE*1000)

                sumHitE = 0
                for iH in range(mH):
                    ele = ae.GetElement(iH)
                    ene = ele.GetEnergy()*1000
                    act = ele.GetActiveness()
                    if ene < 0.7 or np.isnan(ene) or np.isnan(act):
                        continue
                    sumHitE += ene
                sumHitESpec.append(sumHitE)

            if 0.23775 < sumE < 0.23924 and mH==2 : # sims window
            # if 0.23728 < sumE < 0.23946 and mH==2 : # data window

                hitE = []
                for iH in range(mH):
                    ele = ae.GetElement(iH)
                    ene = ele.GetEnergy()*1000
                    act = ele.GetActiveness()
                    # if ene < 0.7 or np.isnan(ene) or np.isnan(act):
                        # continue
                    actList.append(act)
                    hits[mH].append(ene/act)
                    hitE.append(ene)

                hitE = np.asarray(hitE)
                idx = np.where((hitE > 237.5) & (hitE < 239.2))
                if len(idx[0]) > 0:
                    nPk += 1
                    print("%-8d  pk   sumE %.2f  hits:[" % (iEnt, sumE*1000), "  ".join("%.2f" % e for e in hitE),"]")

                idx2 = np.where((hitE > 233.6) & (hitE < 235.55))
                if len(idx2[0]) > 0:
                    nCont += 1
                    print("%-8d  con  sumE %.2f  hits:[" % (iEnt, sumE*1000), "  ".join("%.2f" % e for e in hitE),"]")

    print("Peak cts: %d  Cont cts: %d" % (nPk,nCont))

    np.savez("../data/mult4-simTest.npz", hits, actList, totESpec, sumHitESpec)





if __name__=="__main__":
    main()

