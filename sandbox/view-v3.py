#!/usr/bin/env python
import sys
import numpy as np
import waveLibs as wl
from scipy.signal import butter, lfilter
from ROOT import TChain, TTree

import matplotlib.pyplot as plt
plt.style.use('../clint.mpl')

def main(argv):

    ds = 5
    tt = TChain("skimTree")
    tt.Add("~/project/cal/lat/*.root")

    # i'm just pasting the output of view-v2 here.

    # thesis plot, fast vs slow
    # trapENF 16.0, iE 20958, iH 0 # fast
    # trapENF 16.0, iE 33039, iH 0 # slow
    evtList = [[20958,0],[33039,0]]
    cols = ['k','r','g','m']
    labels = ['Fast, 16.0 keV','Slow, 16.0 keV']
    alphas = [1, 0.7]

    for i, (iE, iH) in enumerate(evtList):
        tt.GetEntry(iE)

        run = tt.run
        chan = tt.channel.at(iH)
        hitE = tt.trapENFCal.at(iH)
        wf = tt.MGTWaveforms.at(iH)
        truncLo, truncHi = 0, 2
        if ds==6 or ds==2: truncLo = 4
        signal = wl.processWaveform(wf,truncLo,truncHi)
        waveBLSub = signal.GetWaveBLSub()
        waveTS = signal.GetTS()
        print("%d / %d  Run %d  chan %d  trapENF %.1f" % (i,len(evtList),run,chan,hitE))

        # plt.plot(waveTS, waveBLSub, "-", lw=1, c=cols[i], label="%.1f keV" % hitE)

        plt.plot(waveTS, waveBLSub, "-", lw=2, c=cols[i], alpha=alphas[i], label=labels[i])

        # # standard energy trapezoid
        # eTrap = wl.trapFilter(waveBLSub,400,250,-7200)
        #
        # nPad = len(waveBLSub)-len(eTrap)
        # eTrap = np.pad(eTrap, (nPad,0), 'constant')
        # eTrapTS = np.arange(0, len(eTrap)*10., 10)
        # ePickoff = eTrapTS[nPad + 400 + 200]
        # # plt.axvline(ePickoff, c='m')
        #
        # # trigger trap
        # tTrap = wl.trapFilter(waveBLSub,100,150,-7200)
        # tTrap = np.pad(tTrap, (nPad,0), 'constant')
        # tTrapTS = np.arange(0, len(tTrap)*10., 10)
        #
        # # low pass filter
        # B, A = butter(2,1e6/(1e8/2), btype='lowpass')
        # waveLP = lfilter(B, A, waveBLSub)
        #
        # plt.cla()
        # plt.plot(waveTS, waveBLSub, 'b', label='Raw WF, %.2f keV' % (hitE))
        # # plt.plot(waveTS, waveLP, 'r', alpha=0.7, label='Low-pass filter')

    plt.ylim(ymin=-10)

    plt.xlabel("Time (ns)", ha='right', x=1)
    plt.ylabel("Voltage (ADC)", ha='right',y=1)
    plt.legend(loc=4)

    # plt.pause(0.0001) # scan speed, does plt.show(block=False) automatically
    # plt.show()

    plt.savefig('../plots/fast-and-slow.pdf')


if __name__ == "__main__":
    main(sys.argv[1:])
