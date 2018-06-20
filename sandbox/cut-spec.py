#!/usr/bin/env python3

import glob
import dsi

def main():

    # specify a bkg cut folder
    # "lat", "th", "rn", "fs", "fr90", "fr95", "frb90", "frb95", "final90", "final95"
    cut = "final95"

    # specify parameters you want 2D histograms vs energy for
    # "fitSlo","riseNoise","tOffset"
    par = "fitSlo"

    # specify binning of each parameter (assume bins[0] = energy, xLo, xHi, xpb)
    bins = [[0, 50, 0.1], [-50, 200, 1]]

    # specify which ds's to combine
    dsList = [0,1,2,3,4,"5A","5B","5C"]

    getSpectrumData(cutType)


def getSpectrumData(cut, par, bins):

    xLo, xHi, xpb = bins[0]
    yLo, yHi, ypb = bins[1]

    # like a genius, i used a ton of different output file names
    lat : latSkimDS5_40_16.root
    final90: final90_DS1.root
    final95: final95_DS1.root
    fr, fr90, fr95: fr_ds0_0_ch576.root
    frb90 : frb90_ds0_0_ch576.root
    frb95: frb95_ds0_0_ch576.root
    fs: fs_ds0_0_ch594.root
    th: th_ds0_10_ch592.root
    *no rn

    bLo, bHi =
    if ds=="5A": bLo, bHi = 0, 79
    if ds=="5B": bLo, bHi = 80, 112
    if ds=="5C": bLo, bHi = 113, 121


    dataDir = "%s/%s/" % (dsi.cutDir, par)
    if par == "lat": dataDir = "%s/" % dsi.latDir

    fList = []
    for ds in dsList:
        if cut in ["lat", "":
            fList.extend(glob.glob("%s/*DS%d*" % (dataDir,ds))
        else:
            fList.extend(glob.glob("%s/*ds%d*.root" % (dataDir,ds)))

        final95_DS0.root
        fr95: fr_ds5_9_ch678.root










if __name__=="__main__":
    main()