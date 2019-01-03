#!/usr/bin/env python
import sys
import numpy as np
import dsi
import waveLibs as wl


import matplotlib.pyplot as plt
plt.style.use('../clint.mpl')

def main():

    # saveWFs()
    plotWFs()


def saveWFs():
    from ROOT import TChain

    # 1 P3KJR 9472 5254 C1P3D2 610 2115.519087 0.831066 -0.000310 1
    # 1 P3KJR 9544 14210 C1P1D2 582 2083.442207 4.322209 -0.000296 2
    # 1 P3KJR 10614 8717 C1P1D2 582 2167.199298 1.929810 -0.000396 19
    # 1 P3KJR 13345 31547 C1P1D3 580 2013.263961 14.821248 -0.000213 46
    tt = TChain("skimTree")
    tt.Add("%s/latSkimDS%d*.root" % (dsi.latDir,1))

    chan = 582
    ds = 1
    truncLo, truncHi = 0, 2
    if ds==6 or ds==2: truncLo = 4

    # wf 1 - passing DCR
    tCut1 = "trapENFCal > 2000 && trapENFCal < 2200 && channel==%d && avse>-1 && dcr99<0" % chan
    n = tt.Draw("Entry$:Iteration$",tCut1,"goff")
    evt, itr = tt.GetV1(), tt.GetV2()
    evtList = [[int(evt[i]),int(itr[i])] for i in range(n)]
    nEnt = len(set([evt[i] for i in range(n)]))
    passWF = []
    for i, (iE, iH) in enumerate(evtList):
        tt.GetEntry(iE)
        run = tt.run
        hitE = tt.trapENFCal.at(iH)
        dcr99 = tt.dcr99.at(iH)
        wf = tt.MGTWaveforms.at(iH)
        signal = wl.processWaveform(wf,truncLo,truncHi)
        wf1 = signal.GetWaveBLSub()
        ts1 = signal.GetTS()
        passWF.append((ts1,wf1,hitE,dcr99))
        print("%d / %d  Run %d  chan %d  trapENF %.1f" % (i+1,len(evtList),run,chan,hitE))

    # wf 2 - failing DCR
    tCut2 = "trapENFCal > 1500 && channel==%d && dcr99>0.001" % chan
    n = tt.Draw("Entry$:Iteration$",tCut2,"goff")
    evt, itr = tt.GetV1(), tt.GetV2()
    evtList = [[int(evt[i]),int(itr[i])] for i in range(n)]
    nEnt = len(set([evt[i] for i in range(n)]))

    failWF = []
    for i, (iE, iH) in enumerate(evtList):
        tt.GetEntry(iE)
        run = tt.run
        hitE = tt.trapENFCal.at(iH)
        dcr99 = tt.dcr99.at(iH)
        wf = tt.MGTWaveforms.at(iH)
        signal = wl.processWaveform(wf,truncLo,truncHi)
        wf2 = signal.GetWaveBLSub()
        ts2 = signal.GetTS()
        failWF.append((ts2,wf2,hitE,dcr99))
        print("%d / %d  Run %d  chan %d  trapENF %.1f" % (i+1,len(evtList),run,chan,hitE))

    np.savez("../data/dcr-example.npz",passWF,failWF)


def plotWFs():

    f = np.load("../data/dcr-example.npz")
    passWF, failWF = f['arr_0'], f['arr_1']

    # passWF
    # 1 / 2  Run 10614  chan 582  trapENF 2167.1
    # 2 / 2  Run 9544  chan 582  trapENF 2083.4

    # failWF
    # 1 / 12  Run 10131  chan 582  trapENF 1772.3
    # 2 / 12  Run 10149  chan 582  trapENF 2136.9
    # 3 / 12  Run 10555  chan 582  trapENF 1805.4
    # 4 / 12  Run 9484  chan 582  trapENF 2885.4
    # 5 / 12  Run 11094  chan 582  trapENF 2970.3
    # 6 / 12  Run 11100  chan 582  trapENF 2734.8
    # 7 / 12  Run 11408  chan 582  trapENF 2087.0
    # 8 / 12  Run 11426  chan 582  trapENF 1805.4
    # 9 / 12  Run 9555  chan 582  trapENF 2379.7
    # 10 / 12  Run 13050  chan 582  trapENF 2114.2
    # 11 / 12  Run 9667  chan 582  trapENF 3255.4
    # 12 / 12  Run 10028  chan 582  trapENF 2931.0

    ts1, wf1, hit1, dcr1 = passWF[0]
    ts2, wf2, hit2, dcr2 = failWF[7]

    # for i, fail in enumerate(failWF):
        # print(i,fail[2],fail[3])
    # return

    max1 = np.max(wf1)
    plt.plot(ts1, wf1/max1, 'b', lw=2, label='Passing DCR, %.1f keV' % hit1)

    max2 = np.max(wf2)
    plt.plot(ts2, wf2/max2, 'r', lw=2, label='Failing DCR, %.1f keV' % hit2)
    plt.xlabel("Time (ns)", ha='right', x=1)
    plt.ylabel("Amplitude (normalized)", ha='right',y=1)
    plt.legend(loc=4)

    plt.xlim(9500,20000)
    plt.ylim(0.8,1.05)

    # plt.show()
    plt.savefig('../plots/dcr-plot.pdf')


if __name__=="__main__":
    main()