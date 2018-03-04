#!/usr/bin/env python3
import os, math, glob, imp
import numpy as np
import matplotlib.pyplot as plt
plt.style.use('../pltReports.mplstyle')
import seaborn as sns

# load LAT libraries
ds = imp.load_source('DataSetInfo',os.environ['LATDIR']+'/DataSetInfo.py')
wl = imp.load_source('waveLibs',os.environ['LATDIR']+'/waveLibs.py')
calInfo = ds.CalInfo()

def main():

    # getSimData()
    plotSimData()


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

    np.savez("../plots/sim1-evtTrans.npz",evtTotal,evtBulk,evtTrans)


def plotSimData():

    f = np.load("../plots/sim1-evtTrans.npz")
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


if __name__=="__main__":
    main()

