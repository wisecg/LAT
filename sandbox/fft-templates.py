#!/usr/bin/env python
import sys, random
from ROOT import TFile,gDirectory
import matplotlib.pyplot as plt
import waveLibs as wl
import numpy as np

""" This is an adaptation of the LIGO "Find an Inspiral" tutorial,
    applied to MJD data. https://losc.ligo.org/tutorial06/
"""
def main(argv):

    # Make a figure
    fig = plt.figure(figsize=(8,5), facecolor='w')
    p1 = plt.subplot(111)
    plt.show(block=False)

    # Sampling frequency - 100MHz
    fs = 1e8

    total_fft = np.zeros(1008)
    total_x_fft = np.zeros(1008)

    for i in range(100):

        # Make a template
        # samp, r, z, tempE, tempStart, smooth = 5000, 0, 15, 10, 2500, 100  # huge template
        # samp, r, z, tempE, tempStart, smooth = 2016, 0, 15, 10, 1000, 100 # regular size template
        samp, r, z, tempE, tempStart, smooth = 2016, 30, 30, 10, 1000, 100 # regular size temp, slower rise
        # samp, r, z, tempE, tempStart, smooth = 500, 0, 15, 10, 100, 100  # small template

        samp = 2016
        r = random.uniform(1, 30)
        z = random.uniform(1, 30)
        tempE = 10
        tempStart = 1000
        smooth = random.uniform(1,200)

        temp, tempTS = wl.MakeSiggenWaveform(samp,r,z,tempE,tempStart,smooth)

        x_fft, _, temp_fft = wl.fourierTransform(temp) # returns a power spectrum
        total_x_fft = x_fft
        total_fft += temp_fft

        p1.cla()
        idx = np.where(x_fft > 1000000)
        p1.loglog(x_fft[idx],total_fft[idx],color='blue')

        plt.pause(0.000001)


    plt.show()


if __name__ == "__main__":
    main(sys.argv[1:])
