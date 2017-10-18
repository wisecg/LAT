#!/usr/bin/env python
import sys
from ROOT import TFile,TTree,TChain,TEntryList,gDirectory,gROOT,MGTWaveform
import ROOT
import waveLibs as wl
import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from scipy.signal import correlate
from scipy.ndimage.filters import gaussian_filter
import pywt
import time
import matplotlib.cm as cm

def main(argv):
    scanSpeed = 0.0001
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
    inputFile = TFile("~/project/match-skim/waveletSkimDS5_run23920.root")
    waveTree = inputFile.Get("skimTree")

    # Set cut
    theCut = inputFile.Get("cutUsedHere").GetTitle()
    theCut += " && waveS5/trapENFCal < 1200 && trapENFCal > 2 && trapENFCal < 5"
    waveTree.Draw(">>elist", theCut, "entrylist")
    elist = gDirectory.Get("elist")
    waveTree.SetEntryList(elist)
    nList = elist.GetN()
    print "Using cut:\n",theCut,"\n"
    print "Found",nList,"entries passing cuts."

    # Make a figure
    fig = plt.figure(figsize=(8,7), facecolor='w')
    p1 = plt.subplot(211)
    p2 = plt.subplot(212)
    plt.show(block=False)

    # Make a huge floating template waveform
    samp, r, z, tempE, tempStart, smooth = 5000, 0, 15, 10, 2500, 100
    tempOrig, tempOrigTS = wl.MakeSiggenWaveform(samp,r,z,tempE,tempStart,smooth)

    # Loop over events
    iList = -1
    while(True):
        iList += 1
        if intMode==True and iList!=0:
            value = raw_input()
            if value=='q': break
            if value=='p': iList -= 2  # previous
            if (value.isdigit()): iList = int(value) # go to entry
        if iList >= elist.GetN(): break

        entry = waveTree.GetEntryNumber(iList);
        waveTree.LoadTree(entry)
        waveTree.GetEntry(entry)
        nChans = waveTree.channel.size()
        nWaves = waveTree.MGTWaveforms.size()
        numPass = waveTree.Draw("channel",theCut,"GOFF",1,iList)
        chans = waveTree.GetV1()
        chanList = list(set(int(chans[n]) for n in xrange(numPass)))

        # Loop over hits passing cuts
        hitList = (iH for iH in xrange(nChans) if waveTree.channel.at(iH) in chanList)  # a 'generator expression'
        for iH in hitList:

            # Load waveform data
            run = waveTree.run
            chan = waveTree.channel.at(iH)
            dataE = waveTree.trapENFCal.at(iH)
            dataENM = waveTree.trapENM.at(iH)
            dataMaxTime = waveTree.trapENMSample.at(iH)*10. - 4000
            # dataMaxTime = waveTree.blrwfFMR50.at(iH)
            signal = wl.processWaveform(waveTree.MGTWaveforms.at(iH))
            data = signal.GetWaveBLSub()
            dataTS = signal.GetTS()
            print "%d / %d  Run %d  nCh %d  chan %d  trapENF %.1f" % (iList,nList,run,nChans,chan,dataE)

            # Draw a trapezoid filter, append with 1000 samples of zeros
            # trap = np.zeros(1000)
            # trap = np.append(trap,wl.trapezoidalFilter(data))

            # 1. Signal-template convolution. Time domain, technically not quite a full matched filter?
            temp, tempTS = flipAndAlign(tempOrig,tempOrigTS,dataTS,dataMaxTime,dataENM,tempE)
            smf = gaussian_filter(temp * data,sigma=float( 10 ))

            # 2. crazy stuff that doesn't work
            deriv = wl.wfDerivative(data)
            wp = pywt.WaveletPacket(data=deriv, wavelet='haar', mode='symmetric',maxlevel=3)
            new_wp = pywt.WaveletPacket(data=None, wavelet='haar', mode='symmetric')
            new_wp['aaa'] = wp['aaa'].data
            denoised = new_wp.reconstruct(update=False)

            # 3. LIGO-inspired matched filter


            # cmap = cm.get_cmap('Spectral') # then in plot: color=cmap(fraction)
            p1.cla()
            p1.set_ylabel("ADC")
            p1.plot(dataTS,data,color='blue')
            p1.plot(tempTS,temp,color='red')
            # p1.plot(tempTS,trap,color='green')
            p1.plot(tempTS,smf,color='red')

            p2.cla()
            p2.plot(dataTS,denoised,color='blue')






            plt.pause(scanSpeed)


def flipAndAlign(tempOrig,tempOrigTS,dataTS,dataMaxTime,dataENM,tempE):
    temp, tempTS = tempOrig, tempOrigTS
    temp = temp * (dataENM / tempE)
    temp = np.flip(temp,0)
    tempMaxTime = np.argmax(temp)*10 # find max after flipping
    tempTS = tempTS - (tempMaxTime - dataMaxTime)
    idx = np.where((tempTS >= dataTS[0]-5) & (tempTS <= dataTS[-1]+5))
    temp, tempTS = temp[idx], tempTS[idx]
    return temp, tempTS


if __name__ == "__main__":
    main(sys.argv[1:])
