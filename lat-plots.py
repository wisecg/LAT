#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
plt.style.use('./pltReports.mplstyle')

import dsi
import waveLibs as wl

def main():

    spec1()

def spec1():
    from ROOT import TFile, TChain, TTree

    ds = 0
    inFile = "%s/bkg/cut/final/final_DS%s.root" % (dsi.dataDir, ds)
    tf = TFile(inFile)
    tt = tf.Get("skimTree")
    enrExp = float(tf.Get("enrExp (kg-d)").GetTitle())
    natExp = float(tf.Get("natExp (kg-d)").GetTitle())

    tCut = "!isEnr"

    xLo, xHi, xpb = 0, 20, 0.1

    n = tt.Draw("trapENFCal",tCut,"goff")
    hitE = tt.GetV1()
    hitE = [hitE[i] for i in range(n)]
    x, h1 = wl.GetHisto(hitE, xLo, xHi, xpb, shift=False)

    # h1 = np.divide(h1, enrExp*xpb) # scale by exposure and binning

    plt.plot(x, h1, 'b', ls='steps', lw=2, label='DS-%s' % ds)

    plt.xlabel("Energy (keV)", ha='right', x=1)

    # TODO: how to do the per kev scaling again?
    plt.ylabel("Counts/keV-kg-d", ha='right', y=1)

    plt.legend(loc=1)

    plt.tight_layout()

    plt.show()




if __name__=="__main__":
    main()