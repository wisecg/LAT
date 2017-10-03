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

    # ======= calculate trapezoids =======

    # standard trapezoid - prone to walking, less sensitive to noise.  use to find energy
    eTrap = wl.trapFilter(data, 400, 250, 7200.)
    eTrapTS = np.arange(dataTS[0], len(eTrap)*10., 10)

    # short trapezoid - triggers more quickly, sensitive to noise.  use to find t0
    sTrap = wl.trapFilter(data, 100, 150, 7200.)
    sTrapTS = np.arange(dataTS[0], len(sTrap)*10., 10)

    # asymmetric trapezoid - used to find the t0 only
    aTrap = wl.asymTrapFilter(data, 200, 100, 40, True)
    aTrapTS = np.arange(dataTS[0], len(aTrap)*10., 10)


    # ======= find leading edges (t0 times) =======

    # limit the range from 0 to 14us, and use an ADC threshold of 2.0 (like the data) for now ...
    t0_SLE,_ = wl.walkBackT0(sTrap, 2., 0, 1000) # (in ns) finds leading edge from short trap
    t0_ALE,_ = wl.walkBackT0(aTrap, 2., 0, 1000) # (in ns) finds leading edge from asymmetric trap

    # standard energy trapezoid w/ a baseline padded waveform
    data_pad = np.pad(data,(400,0),'symmetric')
    pTrap = wl.trapFilter(data_pad, 400, 250, 7200.)
    pTrapTS = np.linspace(0, len(pTrap)*10, len(pTrap))


    # ======= calculate energy parameters =======

    # standard amplitude.  basically trapEM, but w/o NL correction if the input WF doesn't have it.
    latE = np.amax(eTrap)

    # standard amplitude with t0 from the shorter traps
    # If either fixed pickoff time (t0) is < 0, use the first sample as the amplitude (energy).

    eTrapInterp = interpolate.interp1d(eTrapTS, eTrap)
    latEF = eTrapInterp( np.amax([t0_SLE-7000+4000+2000, 0.]) ) # This should be ~trapEF
    latEAF = eTrapInterp( np.amax([t0_ALE-7000+4000+2000, 0.]) )

    # amplitude from padded trapezoid, with t0 from short traps and a correction function

    # function is under development.  currently: f() = exp(p0 + p1*E), p0 ~ 8.5, p1 ~ -0.1
    # functional walk back distance is *either* the minimum of the function value, or 5500 (standard value)

    pTrapInterp = interpolate.interp1d(pTrapTS, pTrap)
    t0_corr = -7000+8000+2000 - np.amin([np.exp(8.5 - 0.1*latE),5500.])
    latEFC = pTrapInterp( np.amax([t0_SLE + t0_corr, 0.]) )
    latEAFC = pTrapInterp( np.amax([t0_ALE + t0_corr, 0.]) )

    print "trapENM %.2f || latEM %.2f  ef %.2f  eaf %.2f  efc %.2f  eafc %.2f" % (dataENM,latE,latEF,latEAF,latEFC,latEAFC)

    # ======= make some plots =======

    fig = plt.figure(figsize=(10,6),facecolor='w')

    plt.plot(dataTS, data, color='blue', label='data')

    plt.plot(sTrapTS, sTrap, color='red', label='sTrap')
    plt.axvline(t0_SLE, color='red')

    plt.plot(aTrapTS, aTrap, color='orange', label='aTrap')
    plt.axvline(t0_ALE, color='orange')

    plt.plot(eTrapTS, eTrap, color='green', label='eTrap')
    plt.axhline(latE,color='green')

    plt.plot(pTrapTS, pTrap, color='magenta', label='pTrap')
    plt.axhline(latEAFC, color='magenta')

    plt.title("trapENM %.2f || latEM %.2f  ef %.2f  eaf %.2f  efc %.2f  eafc %.2f" % (dataENM,latE,latEF,latEAF,latEFC,latEAFC))
    plt.legend(loc='best')
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()