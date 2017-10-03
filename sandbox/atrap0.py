#!/usr/bin/env python
import sys, time, imp
import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as interpolate
wl = imp.load_source('waveLibs', '../waveLibs.py')

def main():
    """ Calculate some new trapezoid parameters: latEM, latEF, latEAF, latEFC, latEAFC """

    # load a waveform
    file5 = np.load("./data/ds5exampleWaveform5.npz")
    dataTS, data, dataENM = file5['arr_0'], file5['arr_1'], file5['arr_3']

    # Short Trapezoid
    shortTrap = wl.trapFilter(data, rampTime=100, flatTime=150, decayTime=7200.)

    # Asymmetric Trapezoid
    asymTrap = wl.asymTrapFilter(data, ramp=200, flat=100, fall=40, padAfter=True)

    # Standard Energy Trapezoid
    longTrap = wl.trapFilter(data, rampTime=400, flatTime=250, decayTime=7200.)
    longTrapTS = np.linspace(0, len(longTrap)*10, len(longTrap))

    # Standard Energy Trapezoid with Baseline padded waveform
    padTrap = wl.trapFilter(np.pad(data, (400, 0), mode='symmetric'), rampTime=400, flatTime=250, decayTime=7200.)
    padTrapTS = np.linspace(0, len(padTrap)*10, len(padTrap))

    longTrapInterp = interpolate.interp1d(longTrapTS, np.asarray(longTrap).squeeze())
    padTrapInterp = interpolate.interp1d(padTrapTS, np.asarray(padTrap).squeeze())

    # Limit the range from 0 to 1000 samples (0 to 10 us) -- using 2.0 threshold like data for now...
    # Returned start times are all in units of ns!
    t0_SLE,_ = wl.walkBackT0(shortTrap, thresh=2.0, rmin=0, rmax=1000) # Leading-Edge on short trapezoid
    t0_ALE,_ = wl.walkBackT0(asymTrap, thresh=2.0, rmin=0, rmax=1000) # Leading-Edge on asymmetric trapezoid

    # Amplitude Evaluation -- Standard
    latE = np.amax(longTrap) # Basically trapENM

    # If fixed pickoff time is < 0, fix to 0. for interpolation
    latEF = longTrapInterp(np.amax([t0_SLE-7000+4000+2000, 0.])) # This should be ~trapENM
    latEAF = longTrapInterp(np.amax([t0_ALE-7000+4000+2000, 0.]))

    # Amplitude Evaluation -- Functional correction for t0 triggerWalk and then use padTrap
    # Function is still in development...
    # Function is exp(p0 + p1*E) where p0 ~ 8.5 and p1 ~ -0.1
    # Functional walk back distance is the minimum of the function value or ~5800 (standard value)

    t0_corr = -7000+8000+2000 - np.amin([np.exp(8.5 - 0.1*latE),5500.])
    latEFC = padTrapInterp( np.amax([t0_SLE + t0_corr, 0.]) )
    latEAFC = padTrapInterp( np.amax([t0_ALE + t0_corr, 0.]) )


    print "trapENM %.2f || latEM %.2f  ef %.2f  eaf %.2f  efc %.2f  eafc %.2f" % (dataENM,latE,latEF,latEAF,latEFC,latEAFC)

    return

    # make some plots
    fig = plt.figure(figsize=(10,6),facecolor='w')

    plt.plot(dataTS, data, color='blue', label='data')
    plt.plot(pTrapTS, data_pad, color='blue', alpha=0.5)

    plt.plot(sTrapTS, sTrap, color='red', label='sTrap')
    plt.axvline(t0_SLE, color='red')

    plt.plot(aTrapTS, aTrap, color='orange', label='aTrap')
    plt.axvline(t0_ALE, color='orange')

    plt.plot(eTrapTS, eTrap, color='green', label='eTrap', linewidth=3)
    plt.axhline(latE,color='green',linewidth=3)

    plt.plot(pTrapTS, pTrap, color='magenta', label='pTrap')
    plt.axhline(latEAFC,color='magenta')

    plt.title("trapENM %.2f || latEM %.2f  ef %.2f  eaf %.2f  efc %.2f  eafc %.2f" % (dataENM,latE,latEF,latEAF,latEFC,latEAFC))
    plt.legend(loc='best')
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()