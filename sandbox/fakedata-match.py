#!/usr/local/bin/python
import sys
from ROOT import TFile,TTree,TChain,TEntryList,gDirectory,gROOT,MGTWaveform
import ROOT
import waveLibs as wl
import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from scipy.signal import correlate
import pywt
import time

def main(argv):
    scanSpeed = 1.0
    opt1, opt2 = "", ""
    intMode, printWF = False, False
    if (len(argv) >= 1): opt1 = argv[0]
    if (len(argv) >= 2): opt2 = argv[1]
    if "-i" in (opt1, opt2):
        intMode = True
        print "Interactive mode selected."
    if "-s" in (opt1, opt2):
        printWF = True
        print "Saving WF plots to current directory."

    # Set input file
    inputFile = TFile("./data/SimWF_Baseline.root") # simulated noise pulses from Brian
    waveTree = inputFile.Get("skimTree")
    theCut = ""

    # Make a figure
    fig = plt.figure(figsize=(8,5), facecolor='w')
    # fig, (ax_orig, ax_noise, ax_corr) = plt.subplots(3, 1, sharex=True)
    # p1 = plt.subplot2grid((6,7), (0,0), colspan=4, rowspan=2) # original
    p1 = plt.subplot(111)
    plt.show(block=False)

    # Make a template waveform
    samp, r, z, ene, t0, smooth = 5000, 0, 15, 10, 2500, 100 # huge floating waveform
    tempOrig, tempOrigTS = wl.MakeSiggenWaveform(samp,r,z,ene,t0,smooth)

    # Loop over events
    iEnt = -1
    while(True):
        iEnt += 1
        if intMode==True and iEnt!=0:
            value = raw_input()
            if value=='q': break
            if value=='p': iEnt -= 2  # previous
            if (value.isdigit()): iEnt = int(value) # go to entry
        if iEnt >= waveTree.GetEntriesFast(): break
        waveTree.GetEntry(iEnt)

        # Make a fake "data" pulse
        signal = wl.processWaveform(waveTree.MGTWaveforms.at(0)) # brian only has 1 pulse per event
        noise = signal.GetWaveBLSub()
        dataTS = signal.GetTS() * 10.
        samp, r, z, ene, t0, smooth = 2016, 0, 15, 10, 1000, 100 # realistic waveform
        sigg, siggTS = wl.MakeSiggenWaveform(samp,r,z,ene,t0,smooth)
        data = sigg + noise

        # how long does this take?
        start = time.time()
        wp = pywt.WaveletPacket(data=data, wavelet='haar', mode='symmetric',maxlevel=3)
        new_wp = pywt.WaveletPacket(data=None, wavelet='haar', mode='symmetric')
        new_wp['aaa'] = wp['aaa'].data
        newWave = new_wp.reconstruct(update=False)
        stop = time.time()
        print "denoising took",stop-start

        # what about this part?
        start = time.time()
        deriv = wl.wfDerivative(newWave)
        stop = time.time()
        print "derivative took",stop-start

        p1.cla()
        p1.set_ylabel("ADC")
        p1.plot(dataTS,data,color='blue')
        p1.plot(siggTS,sigg,color='red')

        # Scale the template down and window it to match the fake data
        # cmap = cm.get_cmap('Spectral')
        # maxCorr = 0
        # bestShift = 0
        # for i in range(10):
        #     temp = tempOrig
        #     tempTS = tempOrigTS
        #     timeShift = -10000 - 1000*i
        #     tempTS = tempTS + timeShift
        #
        #     # May need to manually window if array matching becomes important
        #     idx = np.where((tempTS >= dataTS[0]) & (tempTS <= dataTS[-1]+5))
        #     tempTS = tempTS[idx]
        #     temp = temp[idx]
        #
        #     # Matched filter thing - flip it and reverse it
        #     # temp = -1.0 * temp
        #     # temp = np.flip(temp,0)
        #
        #     # apply correlation with template
        #     corr = correlate(temp,data,mode='same')
        #
        #     localMax = np.amax(corr)
        #     if localMax > maxCorr:
        #         maxCorr = localMax
        #         bestShift = timeShift
        #     print "timeShift %i  localMaxCorr %.1f  len(tempTS) %i  len(corr) %i" % (timeShift,localMax,len(tempTS),len(corr))
        #
        #     p1.plot(tempTS,temp,color=cmap(i/10.))
        #     p1.plot(tempTS,corr/10000,color=cmap(i/10.))
        #
        # print "final corr value:",maxCorr," bestShift ",bestShift
        #
        # plt.pause(scanSpeed)



if __name__ == "__main__":
    main(sys.argv[1:])
