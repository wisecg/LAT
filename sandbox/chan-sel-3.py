#!/usr/bin/env python3
import os
import numpy as np
import matplotlib.pyplot as plt
plt.style.use('../pltReports.mplstyle')

import dsi
bkg = dsi.BkgInfo()
det = dsi.DetInfo()

def main():

    checkLTOutput()
    # checkRates()


def checkLTOutput():
    from ROOT import TFile, TTree

    tf = TFile("../ds_2_output.root")
    tt = tf.Get("dsTree")
    for i in range(tt.GetEntries()):
        tt.GetEntry(i)
        for j in range(tt.channel.size()):
            print(tt.run, tt.channel[j], tt.livetime[j])


def checkRates():

    from ROOT import TFile, TTree

    ds = 1
    nBkg = bkg.dsMap()[ds]
    chList = det.getGoodChanList(ds)

    for bIdx in range(nBkg+1):

        print("DS-%d  bIdx %d" % (ds, bIdx))

        for ch in chList:

            fName = "%s/bkg/cut/fr/fr_ds%d_%d_ch%d.root" % (dsi.dataDir, ds, bIdx, ch)

            # skip nonexistent files.  other parts of chsel should tell us why these aren't here.
            if not os.path.isfile(fName):
                print("nope",fName)
                continue

            tf = TFile(fName)
            tt = tf.Get("skimTree")

            print(bIdx, ch, tt.GetEntries())

        return




if __name__=="__main__":
    main()