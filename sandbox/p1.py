#!/usr/bin/env python3
import sys, os
import numpy as np

import matplotlib as mpl
mpl.use('Agg')
sys.argv.append("-b")
import matplotlib.pyplot as plt
plt.style.use('../pltReports.mplstyle')
# from matplotlib.colors import LogNorm, Normalize
# from matplotlib import gridspec

import waveLibs as wl
import dsi
det = dsi.DetInfo()
bkg = dsi.BkgInfo()

def main(argv):

    loadCutFiles()

def loadCutFiles():
    from ROOT import TChain

    dsList = [0,1,2,3,4,"5B","5C"]
    # dsList = [1]

    fList = []
    for ds in dsList:
        dsNum = int(ds[0]) if isinstance(ds, str) else int(ds)

        chList = det.getGoodChanList(dsNum)
        bkgRanges = bkg.getRanges(ds)
        for bIdx in bkgRanges:
            for ch in chList:
                cutFile = "%s/fs/fs_ds%d_%d_ch%d.root" % (dsi.cutDir,dsNum,bIdx,ch)
                # print(cutFile)
                if not os.path.isfile(cutFile):
                    # print("Couldn't find file:",cutFile)
                    continue
                fList.append(cutFile)

    tt = TChain("skimTree")
    for f in fList: tt.Add(f)

    fig = plt.figure()

    # n = tt.Draw("trapENFCal","","goff")
    # hitE = tt.GetV1()
    # hitE = np.asarray([hitE[i] for i in range(n)])
    # xLo, xHi, xpb = 0, 30, 0.2
    # x, y = wl.GetHisto(hitE, xLo, xHi, xpb)
    # plt.semilogy(x, y, ls='steps', c='b', label='bkg')
    # plt.xlabel("Energy (keV)", ha='right', x=1.)
    # plt.legend()
    # plt.savefig("../plots/p1-spec.png")

    n = tt.Draw("trapENFCal:fitSlo:riseNoise","trapENFCal < 30","goff")
    hitE, fSlo, rise = tt.GetV1(), tt.GetV2(), tt.GetV3()
    hitE = np.asarray([hitE[i] for i in range(n)])
    fSlo = np.asarray([fSlo[i] for i in range(n)])
    rise = np.asarray([rise[i] for i in range(n)])

    # plt.cla()
    # plt.plot(hitE, fSlo, ".", ms=2)
    # plt.xlabel("Energy (keV)", ha='right', x=1)
    # plt.ylabel("fitSlo", ha='right',y=1)
    # plt.savefig("../plots/p1-fSlo.png")

    plt.cla()
    plt.plot(hitE, rise, ".", ms=2)
    plt.xlabel("Energy (keV)", ha='right', x=1)
    plt.ylabel("riseNoise", ha='right', y=1)
    plt.ylim(0,5)
    plt.savefig("../plots/p1-rise.png")



if __name__=="__main__":
    main(sys.argv[1:])