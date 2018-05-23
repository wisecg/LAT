#!/usr/bin/env python
import sys
import numpy as np
import waveLibs as wl
from ROOT import TChain, TTree

import matplotlib.pyplot as plt
plt.style.use('../pltReports.mplstyle')

def main(argv):

    quickDraw = False
    for i, opt in enumerate(argv):
        if opt == "-q":
            print("Quick draw mode selected.")
            quickDraw = True

    ds = 5
    tt = TChain("skimTree")
    tt.Add("~/project/cal/lat/*.root")

    tCut = "trapENFCal > 238 && trapENFCal < 239"

    n = tt.Draw("Entry$:Iteration$",tCut,"goff")
    evt, itr = tt.GetV1(), tt.GetV2()
    evtList = [[int(evt[i]),int(itr[i])] for i in range(n)]

    nEnt = len(set([evt[i] for i in range(n)]))
    print("Found %d total entries, %d passing cut: %s" % (tt.GetEntries(), nEnt, tCut))

    i, pEvt = -1, -1
    while(True):
        i += 1
        if not quickDraw and i!=0:
            val = input()
            if val == "q": break
            if val == "p": i -= 2
            if val.isdigit() : i = int(val)
            if val == "s":
                plt.savefig("../plots/wf-%d.pdf" % i)
        if i >= len(evtList): break
        iE, iH = evtList[i]

        if iE != pEvt:
            tt.GetEntry(iE)
        pEvt = iE

        run = tt.run
        chan = tt.channel.at(iH)
        hitE = tt.trapENFCal.at(iH)

        wf = tt.MGTWaveforms.at(iH)
        truncLo, truncHi = 0, 2
        if ds==6 or ds==2: truncLo = 4
        signal = wl.processWaveform(wf,truncLo,truncHi)

        # waveform
        waveBLSub = signal.GetWaveBLSub()
        waveTS = signal.GetTS()
        print("%d / %d  Run %d  chan %d  trapENF %.1f" % (i,len(evtList),run,chan,hitE))

        # standard energy trapezoid
        eTrap = wl.trapFilter(waveBLSub,400,250,-7200)

        nPad = len(waveBLSub)-len(eTrap)
        eTrap = np.pad(eTrap, (nPad,0), 'constant')
        eTrapTS = np.arange(0, len(eTrap)*10., 10)
        ePickoff = eTrapTS[nPad + 400 + 200]
        # plt.axvline(ePickoff, c='m')

        # trigger trap
        tTrap = wl.trapFilter(waveBLSub,100,150,-7200)
        tTrap = np.pad(tTrap, (nPad,0), 'constant')
        tTrapTS = np.arange(0, len(tTrap)*10., 10)

        plt.cla()
        plt.plot(waveTS, waveBLSub, 'b', lw=1.5, label='Raw WF, %.2f keV' % (hitE))
        plt.xlabel("Time (ns)", ha='right', x=1)
        plt.ylabel("Voltage (ADC)", ha='right',y=1)
        plt.legend(loc=4)

        plt.pause(0.0001) # scan speed, does plt.show(block=False) automatically


if __name__ == "__main__":
    main(sys.argv[1:])
