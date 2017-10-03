#!/usr/bin/env python
import pywt
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import waveLibs as wl
import waveModel as wm

def main():
    # Load a data and a temp waveform
    npzfile = np.load("./data/tailSlopeInputs.npz")
    dl, tl = npzfile['arr_0'], npzfile['arr_1']
    data, dataTS, dataE, dataRT = dl[0], dl[1], dl[2], dl[3]
    temp, tempTS, tempE, tempRT = tl[0], tl[1], tl[2], tl[3]

    print len(dataTS)

    # Denoise, just for fun
    wp = pywt.WaveletPacket(data=data, wavelet='haar', mode='symmetric',maxlevel=3)
    new_wp = pywt.WaveletPacket(data=None, wavelet='haar', mode='symmetric')
    new_wp['aaa'] = wp['aaa'].data
    waveDenoised = new_wp.reconstruct(update=False)
    # ix = len(newWave)
    # waveDenoised = np.delete(newWave, [ix-1,ix-2,ix-3])

    # Start the tail slope fit ranged based on the calculated start time
    # Plot the template "guess" to illustrate what we would expect to see
    guessFull, guessFullTS = wm.MakeModel(dl, tl, [dataRT,dataE,1.], opt="full")
    idxMax = np.where(guessFull == guessFull.max()) # returns an array/tuple
    idxMax = idxMax[0][0] # "cast" to int
    idxMax += 100 # move it ahead in case this is a slow pulse
    tail, tailTS = data[idxMax:], dataTS[idxMax:]

    # Exponential fit
    a, b = dataE, 72000
    popt, pcov = curve_fit(wl.tailModelExp, tailTS, tail, p0=[a,b])
    a, b = popt[0], popt[1]

    # Polynomial curve fit
    popt2, pcov2 = curve_fit(wl.tailModelPol, tailTS, tail)
    d, e, f, g = popt2[0], popt2[1], popt2[2], popt2[3]

    # Plots
    fig = plt.figure(figsize=(8,5), facecolor='w')
    a1 = plt.subplot(111)
    a1.set_ylabel("ADC",y=0.95, ha='right')
    a1.set_xlabel("Time (ns)",x=0.95, ha='right')
    a1.set_title("Energy %.1f  StartTime %.1f  %.1e * exp(-t/%.0f[us]) " % (dataE, dataRT, a, b/1000))
    a1.plot(dataTS,data,color='blue',alpha=0.7,label='Data')
    a1.plot(dataTS,waveDenoised,color='red',alpha=0.8,label='Lvl3 Haar Denoised')
    a1.axvline(x=guessFullTS[np.where(guessFull==guessFull.max())],color='black',label='Guess Max ADC')
    a1.plot(guessFullTS,guessFull,color='orange',label='Scaled Template WF')
    a1.plot(tailTS, wl.tailModelExp(tailTS, *popt), color='cyan',linewidth=2,label='Tail Exp. Fit')
    a1.plot(tailTS, wl.tailModelPol(tailTS, *popt2), color='green',linewidth=2,label='Tail Curve Fit')
    a1.legend(loc=4)
    plt.show()
    # plt.savefig("./plots/tailSlope.pdf")

if __name__ == "__main__":
    main()
