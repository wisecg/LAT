#!/usr/bin/env python3
import sys
import glob
import numpy as np
import matplotlib.pyplot as plt
plt.style.use('./pltReports.mplstyle')

import dsi
import waveLibs as wl

def main(argv):

    # like a genius, i used a ton of different output file names
    # lat :     ~/project/bkg/lat/latSkimDS5_40_16.root
    # th :      ~/project/bkg/cut/th/th_ds0_10_ch592.root
    # rn :      none
    # fs :      ~/project/bkg/cut/fs/fs_ds0_0_ch594.root
    # fr95:     ~/project/bkg/cut/fr95/fr_ds0_0_ch576.root
    # frb90     ~/project/bkg/cut/frb90/frb90_ds0_0_ch576.root
    # frb95:    ~/project/bkg/cut/frb95/frb95_ds0_0_ch576.root
    # final90:  ~/project/bkg/cut/final90/final90_DS1.root
    # final95t: ~/project/bkg/cut/final95t/final95t_DS1.root  ***

    # "lat", "th", "rn", "fs", "fr90", "fr95", "frb90", "frb95", "final90", "final95"
    # type = "final95"

    # "fitSlo","riseNoise","tOffset"
    # par = "fitSlo"

    # which ds's to combine
    # dsList = [0,1,2,3,4,"5A","5B","5C"]

    for i, opt in enumerate(argv):
        if opt=="-g":
            getLAT()
        if opt=="-c":
            getCut(argv[i+1])
        if opt=="-f":
            getFinal()

    plotSpec()


def getLAT():

    from ROOT import TFile, TTree

    dsList = [0,1,2,3,4,"5A","5B","5C"]

    # histo output for each ds
    xLo, xHi, xpb = 0, 100, 0.1
    nbx = int((xHi-xLo)/xpb)
    xTot = np.arange(xLo, xHi, xpb)
    hTot = {ds:np.zeros(len(xTot)+1) for ds in dsList} # these 100% match w/ wl.GetHisto's output

    for ds in dsList:
        print("Scanning DS-%s" % (str(ds)))

        dsNum = int(ds[0]) if isinstance(ds,str) else ds

        fList = []
        if ds not in ["5A","5B","5C"]:
            fList = glob.glob("%s/latSkimDS%d*" % (dsi.latDir, ds))
        else:
            if ds=="5A": bLo, bHi = 0, 79
            if ds=="5B": bLo, bHi = 80, 112
            if ds=="5C": bLo, bHi = 113, 121
            for bIdx in range(bLo, bHi+1):
                fList.extend(glob.glob("%s/latSkimDS%d_%d_*" % (dsi.latDir, dsNum, bIdx)))

        for idx, f in enumerate(fList):

            fr = np.fabs(100*idx/len(fList) % 10)
            if fr < 0.5:
                print("%d/%d (file %s) %.1f%% done." % (idx, len(fList), f, 100*idx/len(fList)))

            tf = TFile(f)
            tt = tf.Get("skimTree")

            n = tt.Draw("trapENFCal","","goff")
            hitE = tt.GetV1()
            hitE = [hitE[i] for i in range(n)]

            xF, hF = wl.GetHisto(hitE, xLo, xHi, xpb)
            hTot[ds] = np.add(hTot[ds], hF)

            tf.Close()

    np.savez("./data/lat-cutspec.npz", xTot, hTot)


def getCut(cutType):
    """ ./cut-spec.py -c [cutType] """

    fPrefix = {
        "th"  : "%s/th/th" % dsi.cutDir,          # thresholds only.  ok
        # "fs"  : "%s/fs/fs" % dsi.cutDir,        # can't use this one, it was an old version of the 90% cut
        "fr95": "%s/fr95/fr" % dsi.cutDir,        # this is sort of PSA amalgalm, ok
        "frb95"  : "%s/frb95/frb95" % dsi.cutDir  # this is with burst cut, ok
    }

    from ROOT import TFile, TTree

    dsList = [0,1,2,3,4,"5A","5B","5C"]
    # dsList = ["5B"]

    # histo output for each ds
    xLo, xHi, xpb = 0, 100, 0.1
    nbx = int((xHi-xLo)/xpb)
    xTot = np.arange(xLo, xHi, xpb)
    hTot = {ds:np.zeros(len(xTot)+1) for ds in dsList} # these 100% match w/ wl.GetHisto's output

    for ds in dsList:
        print("Scanning DS-%s" % (str(ds)))

        dsNum = int(ds[0]) if isinstance(ds,str) else ds

        fList = []
        if ds not in ["5A","5B","5C"]:
            fList = glob.glob("%s_ds%d_*" % (fPrefix[cutType], dsNum))
        else:
            if ds=="5A": bLo, bHi = 0, 79
            if ds=="5B": bLo, bHi = 80, 112
            if ds=="5C": bLo, bHi = 113, 121
            for bIdx in range(bLo, bHi+1):
                fList.extend(glob.glob("%s_ds%d_%d_*" % (fPrefix[cutType], dsNum, bIdx)))

        for idx, f in enumerate(fList):

            fr = np.fabs(100*idx/len(fList) % 10)
            if fr < 0.5:
                print("%d/%d (file %s) %.1f%% done." % (idx, len(fList), f, 100*idx/len(fList)))

            tf = TFile(f)
            tt = tf.Get("skimTree")

            n = tt.Draw("trapENFCal","","goff")
            hitE = tt.GetV1()
            hitE = [hitE[i] for i in range(n)]

            xF, hF = wl.GetHisto(hitE, xLo, xHi, xpb)
            hTot[ds] = np.add(hTot[ds], hF)

            tf.Close()

    np.savez("./data/lat-cut-%s-spec.npz" % cutType, xTot, hTot)


def getFinal():
    """ ./cut-spec.py -f """

    from ROOT import TFile, TTree

    dsList = [0,1,2,3,4,"5A","5B","5C"]
    # dsList = ["5B"]

    # cutType = "final90"
    cutType = "final95t" # <-- use this one

    # histo output for each ds
    xLo, xHi, xpb = 0, 100, 0.1
    nbx = int((xHi-xLo)/xpb)
    xTot = np.arange(xLo, xHi, xpb)
    hTot = {ds:np.zeros(len(xTot)+1) for ds in dsList} # these 100% match w/ wl.GetHisto's output

    for ds in dsList:

        tf = TFile("%s/%s/%s_DS%s.root" % (dsi.cutDir, cutType, cutType, str(ds)))

        tt = tf.Get("skimTree")
        n = tt.Draw("trapENFCal","","goff")
        hitE = tt.GetV1()
        hitE = [hitE[i] for i in range(n)]

        xF, hF = wl.GetHisto(hitE, xLo, xHi, xpb)
        hTot[ds] = np.add(hTot[ds], hF)

        tf.Close()

    np.savez("./data/lat-cut-%s-spec.npz" % cutType, xTot, hTot)


def plotSpec():

    # initial lat
    f_lat = np.load("./data/lat-cutspec.npz")
    xTot, hTotDS_lat = f_lat['arr_0'], f_lat['arr_1'].item()
    hTot_lat = np.zeros(len(hTotDS_lat[0]))
    for ds in hTotDS_lat:
        hTot_lat = np.add(hTot_lat, hTotDS_lat[ds])

    # th
    f_th = np.load("./data/lat-cut-th-spec.npz")
    xTot, hTotDS_th = f_th['arr_0'], f_th['arr_1'].item()
    hTot_th = np.zeros(len(hTotDS_th[0]))
    for ds in hTotDS_th:
        hTot_th = np.add(hTot_th, hTotDS_th[ds])

    # fr95
    f_fr = np.load("./data/lat-cut-fr95-spec.npz")
    xTot, hTotDS_fr = f_fr['arr_0'], f_fr['arr_1'].item()
    hTot_fr = np.zeros(len(hTotDS_fr[0]))
    for ds in hTotDS_fr:
        hTot_fr = np.add(hTot_fr, hTotDS_fr[ds])

    # frb95
    f_frb = np.load("./data/lat-cut-frb95-spec.npz")
    xTot, hTotDS_frb = f_frb['arr_0'], f_frb['arr_1'].item()
    hTot_frb = np.zeros(len(hTotDS_frb[0]))
    for ds in hTotDS_frb:
        hTot_frb = np.add(hTot_frb, hTotDS_frb[ds])

    # final
    f_fin = np.load("./data/lat-cut-final95t-spec.npz")
    xTot, hTotDS_fin = f_fin['arr_0'], f_fin['arr_1'].item()
    hTot_fin = np.zeros(len(hTotDS_fin[0]))
    for ds in hTotDS_fin:
        hTot_fin = np.add(hTot_fin, hTotDS_fin[ds])

    xTot = np.insert(xTot, 0, 0, axis=0)

    # limit energy range
    idx = np.where(xTot < 50)

    # cmap = plt.cm.get_cmap('jet',6)

    plt.semilogy(xTot[idx], hTot_lat[idx], c='k', ls='steps', lw=2, label="LAT Skims")
    plt.semilogy(xTot[idx], hTot_th[idx], c='b', ls='steps', lw=2, label="Threshold Cut")
    plt.semilogy(xTot[idx], hTot_fr[idx], c='g', ls='steps', lw=2, label="PSA Cuts")
    plt.semilogy(xTot[idx], hTot_frb[idx], c='orange', ls='steps', lw=2, label="Burst Cut")
    plt.semilogy(xTot[idx], hTot_fin[idx], c='r', ls='steps', lw=2, label="tOffset Cut, Final Spectrum")

    plt.xlabel("Energy (keV)", ha='right', x=1)
    plt.ylabel("Counts", ha='right', y=1)
    plt.legend()
    plt.tight_layout()
    # plt.show()
    plt.savefig("./plots/lat-specReduction.pdf")



if __name__=="__main__":
    main(sys.argv[1:])