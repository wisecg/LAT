#!/usr/bin/env python
import sys
import numpy as np
import waveLibs as wl
from scipy.signal import butter, lfilter
from ROOT import TChain, TTree

import matplotlib.pyplot as plt
plt.style.use('../pltReports.mplstyle')

def main(argv):

    ds = 5
    tt = TChain("skimTree")
    tt.Add("~/project/cal/lat/*.root")

    # trapENF 1.1, iE 9895, iH 0
    # trapENF 2.0, iE 61372, iH 0
    # trapENF 3.1, iE 35219, iH 0
    # trapENF 4.2, iE 33507, iH 0
    # evtList = [[9895,0],[61372,0],[35219,0],[33507,0]]
    evtList = [[33507,0],[35219,0],[61372,0],[9895,0]]
    cols = ['r','g','m','b']

    fig = plt.figure(figsize=(8,7))

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

        plt.plot(waveTS, waveBLSub, "-", lw=1, c=cols[i], label="%.1f keV" % hitE)

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

    plt.xlabel("Time (ns)", ha='right', x=1)
    plt.ylabel("Voltage (ADC)", ha='right',y=1)
    plt.legend(loc=4)

    # plt.pause(0.0001) # scan speed, does plt.show(block=False) automatically
    plt.show()


if __name__ == "__main__":
    main(sys.argv[1:])
