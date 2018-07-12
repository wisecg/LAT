#!/usr/bin/env python3
"""
    % Email chain: "Threshold info for skim files"
    % File: ~/Documents/mjd/MJD\ Reports/david_threshold_erfns.pdf

    % Also David's main talk: http://mjcalendar.npl.washington.edu/indico/event/2351/material/slides/0.pdf
    # using positive and negative trigger data?
    2016-9-15-P3LQG_Run60001538 C2P1D3 P427488

    Questions:
    - Have you worked out how to take into account the effect of the trap-max nature of the on-board threshold?

    - How different are the threshold functions from the standard erf shape?

      -- david:  We could evaluate that by fitting the trap-max E=0 noise peak with a Gaussian and looking at the residuals.
                 I think it's pretty darn close; but how much discrepancy would we want to allow?
      -- jason:  It's an issue of where the analysis threshold will be.  To first order, the lower it is, the more sensitive we
                 are to non-Gaussian.  It also matters if the non-gaussianity is similar or different on all the channels.
                 If it's different, a lot of the differences will average away in a joint fit over many channels.
                 Another option besides fitting gaussians to the E=0 peaks is to fit your threshold curves to erfs.
                 It should be mathematically identical but the latter puts the residuals on the scale that matters and might
                 make it easier to evaluate.
      -- david: (does study)
                The residual on the threshold curve is less than about 1%.  Do we think that's good enough?
      -- jason: Nice work.  You can try adding a skew gaussian, but i would think that 1\% uncertainty is good enough for now,
                especially since it wiggles about the efficiency curve, so when you sum up over many slight left-right shifts
                in the erf position they will largely cancel.

    David's study:
        Take one detector, two different trap filters (trap max and fixed-time pickoff) w/ same rise and flat times.
        Fit with Gaussians
        Integrate to make erfn-type threshold curves
        Look at residuals
        Erfn residuals are less than 1\%.
        Call it good.

    Clint's take:
        For my thesis we are staying 3 sigma above the trigger threshold (about 99%.)
        So we don't need to adjust the erf shape,
        and we've learned from David's study that the onboard trap-max trigger doesn't bias the results
        (we get the same result as if we used the fixed-time pickoff.)
        Also, we have great stats on this measurement for every background run, so we know we're staying at 99% efficiency.
        Ok, done.  Now let's make some plots.
"""
import numpy as np
from scipy.special import erf
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from matplotlib import cm
plt.style.use('../pltReports.mplstyle')
import waveLibs as wl

def main():

    # getTrap()
    plotTrap()
    # plotTrigEff()


def getTrap():

    ds, run, chan = 1, 11085, 580 # C1P1D3
    # trigger eff vals from old study: mu = 1.448, sig = 0.305

    # 2016-9-15-P3LQG_Run60001538 C2P1D3 P42748B
    # ds, run, chan = 4, 60001538, 1110

    from ROOT import GATDataSet
    gds = GATDataSet(run)
    tt = gds.GetGatifiedChain()
    bt = gds.GetBuiltChain()
    tt.AddFriend(bt)

    tCut = "channel==%d" % chan
    n = tt.Draw("Entry$:Iteration$",tCut,"goff")
    evt, itr = tt.GetV1(), tt.GetV2()
    evtList = [[int(evt[i]),int(itr[i])] for i in range(n)]

    trapOutput = []

    pEvt = -1
    for i, (iE, iH) in enumerate(evtList):
        if iE != pEvt:
            tt.GetEntry(iE)
        pEvt = iE

        hitE = tt.trapENFCal.at(iH)
        wf = bt.event.GetWaveform(iH)
        sig = wl.processWaveform(wf)
        waveRaw = sig.GetWaveBLSub()
        waveTS = sig.GetTS()

        # mimic onboard energy trapezoid
        eTrap = wl.trapFilter(waveRaw,400,180,0)
        nPad = len(waveRaw)-len(eTrap)
        eTrap = np.pad(eTrap, (nPad,0), 'constant')
        eTrapTS = np.arange(0, len(eTrap)*10., 10)

        # trap timestamps shouldn't change, don't recalculate them
        if i==0:
            pTS = eTrapTS
        if not np.array_equal(eTrapTS, pTS):
            print("Warning, the timestamps changed.  Uggggg")
            return

        # print(i, iE, iH, hitE)
        trapOutput.append(eTrap)

    # np.savez("../data/trapOutput.npz",trapOutput,eTrapTS) # original ds4 run
    np.savez("../data/trapOutput-2.npz",trapOutput,eTrapTS) # new ds1 run


def plotTrap():
    """ Make 2 plots, one w/ full trap output on top and one zoomed in. """

    f = np.load("../data/trapOutput-2.npz")
    trapOutput, eTrapTS = f['arr_0'], f['arr_1']

    eMax = 0
    for eTrap in trapOutput:
        tmpE = np.log(np.amax(eTrap))
        if tmpE > eMax: eMax = tmpE

    fig = plt.figure(figsize=(8,8))
    p1 = plt.subplot(211)
    p2 = plt.subplot(212)

    cmap = cm.get_cmap('seismic')
    for i, eTrap in enumerate(trapOutput):
        col = np.log(np.amax(eTrap))/eMax
        # col = np.amax(eTrap)/eMax
        p1.plot(eTrapTS, eTrap, ls='steps', lw=2, c=cmap(col))
        p2.plot(eTrapTS, eTrap, ls='steps', lw=2, c=cmap(col))
        if i > 50: break

    p1.set_xlim(9750, 20000)
    p1.set_ylim(ymin=-20)
    # p1.set_xlabel("Time (ns)", ha='right', x=1)
    p1.set_ylabel("Trap Output (ADC)", ha='right', y=1)

    p2.set_xlim(9750, 10100)
    p2.set_ylim(-1, 8)
    p2.set_xlabel("Time (ns)", ha='right', x=1)
    # p2.set_ylabel("Trap Output (ADC)", ha='right', y=1)

    plt.tight_layout()
    # plt.show()
    plt.savefig("../plots/trap-crossing-2.pdf")


def effFunc(x, mu, sig):
    return (1/2) * (1 + erf( (x-mu) / (np.sqrt(2) * sig) ) )


def plotTrigEff():
    """ Plot the noise gaussian, trigger crossing gaussian, and the efficiency function. """

    f = np.load("../data/trapOutput.npz")
    trapOutput, eTrapTS = f['arr_0'], f['arr_1']

    sampNoise, sampTrig = [], []
    for eTrap in trapOutput:
        idx = np.where(eTrap != 0)
        if eTrapTS[idx][0] != 9800:
            print("wtffff, return")
            return

        hitADC = np.amax(eTrap)
        if -1 < hitADC < 10 :
            sampTrig.append(eTrap[idx][8])
        if hitADC > 50:
            sampNoise.append(eTrap[idx][0])

    fig = plt.figure(figsize=(9,5))
    p1 = plt.subplot(121)
    p2 = plt.subplot(122)

    x, hNoise = wl.GetHisto(sampNoise, -1, 2, 0.01)
    x, hTrig = wl.GetHisto(sampTrig, -1, 2, 0.01)
    p1.plot(x, hNoise, 'b', ls='steps', label="Noise Sample")
    p1.plot(x, hTrig, 'r', ls='steps', label="Trigger Sample")

    xF = np.arange(-1, 2, 0.01)

    pNoise,_ = curve_fit(wl.gauss_function, x, hNoise, p0=[1,np.mean(sampNoise),np.std(sampNoise)])
    fitSig = pNoise[2]
    p1.plot(xF, wl.gauss_function(xF,*pNoise), '-m', alpha=0.7, label='Noise, %s = %.2f' % (r'$\sigma$', fitSig))

    pTrig,_ = curve_fit(wl.gauss_function, x, hTrig, p0=[1,np.mean(sampTrig),np.std(sampTrig)])
    fitMu = pTrig[1]
    p1.plot(xF, wl.gauss_function(xF, *pTrig), '-k', alpha=0.5, label='Trigger, %s = %.2f' % (r'$\mu$', fitMu))

    p1.set_xlabel("Trap Output (ADC)", ha='right', x=1)
    p1.set_ylabel("Counts", ha='right', y=1)
    p1.legend(loc=2)

    xE = np.arange(0, 4, 0.01)
    p2.plot(xE, effFunc(xE, pTrig[1], pNoise[2]), '-b')
    p2.axvline(fitMu, color="k", alpha=0.5)

    p2.set_xlabel("Energy (ADC, uncalibrated)", ha='right', x=1)
    p2.set_ylabel("Efficiency", ha='right', y=1)

    # plt.show()
    plt.savefig("../plots/trig-effic.pdf")


if __name__=="__main__":
    main()
