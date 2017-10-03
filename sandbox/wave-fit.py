#!/usr/bin/env python
import time
startT = time.clock()

import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
from matplotlib import gridspec
from ROOT import TFile,TTree,TChain,TEntryList,gDirectory,gROOT,MGTWaveform
import waveLibs as wl
import waveModel as wm
import pymc
from scipy.ndimage.filters import gaussian_filter

stopT = time.clock()
print "Import time (s): ", (stopT-startT)

def main(argv):
    """ Interactive-fit or 'rapid'-fit waveforms that pass a given TCut.
        BUG: Doesn't always work with a TChain.  Add input files together
        with hadd and use a single TFile.
    """
    scanSpeed = 0.2
    iList = -1
    opt1, opt2 = "", "-"
    intMode = False
    if (len(argv) >= 1): opt1 = argv[0]
    if (len(argv) >= 2): opt2 = argv[1]
    if "-i" in (opt1, opt2):
        intMode = True
        print "Interactive mode selected."
    # if opt2.isdigit:
        # iList = int(opt2) - 1

    startT = time.clock()

    # Load template waveforms
    templateFile = TFile("./data/wave-m1cal.root")
    waveTemplates = templateFile.Get("waveTree")
    calDict = {640:28, 672:18, 610:16, 580:19, 582:34, 648:38, 600:21, 578:39, 592:27, 664:55, 626:62, 692:8, 598:22, 690:52, 632:9, 608:7}
    # PlotTemplateWaveforms(waveTemplates,calDict)

    # Set input data and cuts
    # waveTree = TChain("skimTree") # BUG: Can't always recognize MGTWaveforms branch
    # waveTree.Add("~/project/wavelet-skim/waveletSkimDS3*")
    # inputFile = TFile("~/project/wavelet-skim/hadd/waveletSkimDS3.root")
    inputFile = TFile("~/project/wavelet-skim/waveletSkimDS3_13.root")
    waveTree = inputFile.Get("skimTree")
    print "Found",waveTree.GetEntries(),"input entries."

    # theCut = inputFile.Get("cutUsedHere").GetTitle()
    # theCut = "trapENFCal > 0.8 && gain==0 && mH==1 && isGood && !muVeto && !wfDCBits && !isLNFill1 && !isLNFill2 && trapETailMin < 0 && channel!=596 && channel!=676 && channel!=676 && channel!=612 && channel!=1104 && channel!=1200 && channel!=1334 && channel!=1336"

    # DS3 big cut
    theCut = "trapENFCal > 0.8 && gain==0 && mHClean==1 && isGood && !muVeto && !wfDCBits && !isLNFill1 && !isLNFill2 && trapETailMin < 0 && channel!=596 && channel!=676 && channel!=676 && channel!=612 && channel!=1104 && channel!=1200 && channel!=1334 && channel!=1336 && tOffset < 10 && waveS5/TMath::Power(trapENFCal,1/4) < 1200 && (waveS3-waveS2)/trapENFCal > 100 && (waveS3-waveS2)/trapENFCal < 300 && !(channel==692 && (run==16974 || run==16975 || run==16976 || run==16977 || run==16978 || run==16979))"

    # theCut += " && trapENFCal > 2 && trapENFCal < 2.1 && channel==578 && run==16868"  # template fast WF (2kev)
    theCut += " && trapENFCal > 3 && trapENFCal < 4"
    # theCut += " && trapENFCal < 10 && trapENFCal > 1"

    # Print cut and events passing cut
    print "Using cut:\n",theCut,"\n"
    waveTree.Draw(">>elist", theCut, "entrylist")
    elist = gDirectory.Get("elist")
    waveTree.SetEntryList(elist)
    nList = elist.GetN()
    print "Found",nList,"entries passing cuts."

    stopT=time.clock()
    print "Data loading time (s): ", (stopT-startT)

    # Make a figure
    fig = plt.figure(figsize=(8,7), facecolor='w')
    gs = gridspec.GridSpec(3, 1, height_ratios=[2, 3, 1])
    a1 = plt.subplot(gs[0])
    a2 = plt.subplot(gs[1])
    a3 = plt.subplot(gs[2])
    plt.show(block=False)


    # Loop over events
    while(True):
        iList += 1
        if intMode==True and iList != 0:
            value = raw_input()
            if value=='q': break
            if value=='p': iList -= 2  # previous
            if (value.isdigit()): iList = int(value) # go to entry
        if iList >= elist.GetN(): break

        entry = waveTree.GetEntryNumber(iList);
        waveTree.LoadTree(entry)
        waveTree.GetEntry(entry)
        nChans = waveTree.channel.size()
        numPass = waveTree.Draw("channel",theCut,"GOFF",1,iList)
        chans = waveTree.GetV1()
        chanList = list(set(int(chans[n]) for n in xrange(numPass)))

        # Loop over hits passing cuts
        hitList = (iH for iH in xrange(nChans) if waveTree.channel.at(iH) in chanList)  # a 'generator expression'
        for iH in hitList:

            # Load waveform for this hit
            run = waveTree.run
            chan = waveTree.channel.at(iH)
            energy = waveTree.trapENFCal.at(iH)
            print "%d / %d  Run %d  nCh %d  chan %d  trapENF %.1f" % (iList,nList,run,nChans,chan,energy)
            signal = wl.processWaveform(waveTree.MGTWaveforms.at(iH),opt='full')
            waveBLSub = signal.GetWaveBLSub()
            waveFilt = signal.GetWaveFilt()
            waveTS = signal.GetTS()
            baseAvg, noiseAvg = signal.GetBaseNoise()
            riseTime = waveTree.butterTime.at(iH) # replace with blrwfFMR50?

            # Load channel template waveform
            try: tempEntry = calDict[chan]
            except:
                print "ain't got channel",chan
                continue
            waveTemplates.GetEntry(calDict[chan])
            template = wl.processWaveform(waveTemplates.event.GetWaveform(waveTemplates.itr),opt='full')
            templateENF = waveTemplates.trapENFCal
            templateT50 = waveTemplates.blrwfFMR50
            tempTS = template.GetTS()
            tempBLSub = template.GetWaveBLSub()

            # Guess scaling parameters & fill arrays
            scaleGuess = energy/templateENF
            tDiffGuess = riseTime - templateT50
            tempGuess = tempBLSub * scaleGuess
            tempTSGuess = tempTS + tDiffGuess
            tempGuessRT = templateT50 + tDiffGuess

            # Fit around rising edge - rise time calculator method
            lo, hi = 300, 700  # samples before and after the start time
            loWin = riseTime/10 - lo
            hiWin = riseTime/10 + hi
            waveEndsRT = np.concatenate((np.arange(0,loWin), np.arange(hiWin, waveBLSub.size)), axis=0)
            waveEdgeRT = np.delete(waveBLSub, waveEndsRT)
            tsEdgeRT = np.delete(waveTS, waveEndsRT)

            # Fit with MCMC
            numSteps, burnin = 10000, 4000
            # numSteps, burnin = 1000, 400
            waveModel = pymc.Model( wm.WaveformModel(waveEdgeRT, tsEdgeRT, tempBLSub, tempTS, tDiffGuess, scaleGuess, noiseAvg) )
            M = pymc.MCMC( waveModel )
            M.use_step_method(pymc.Metropolis, M.riseTime, proposal_sd=4., proposal_distribution='Normal')
            M.use_step_method(pymc.Metropolis, M.slowness, proposal_sd=0.1, proposal_distribution='Normal')
            M.use_step_method(pymc.Metropolis, M.scale, proposal_sd=0.1, proposal_distribution='Normal')
            M.sample(iter=numSteps, verbose=0) # do the fit
            rt = np.median(M.trace('riseTime')[burnin:])
            slowness = np.median(M.trace('slowness')[burnin:])
            scale = np.median(M.trace('scale')[burnin:])

            # Recreate the guess, the best-fit signal, and the residual
            templateTSRT = tempTS + rt
            lo, hi = tsEdgeRT[0], tsEdgeRT[-1]+10
            idx = np.where((templateTSRT >= lo) & (templateTSRT <= hi))
            scaled = tempBLSub[idx] * scale
            smoothed = gaussian_filter(scaled, sigma=float(slowness))
            smoothedTS = templateTSRT[idx]
            guess = tempBLSub[idx] * scaleGuess
            guessTS = templateTSRT[idx]
            residual1 = smoothed - guess
            residual2 = smoothed - waveEdgeRT
            residual3 = smoothed - waveFilt[idx]

            # Calculate the Chi2 of the fit


            # Fill the figure
            a1.cla()
            a1.set_ylabel("ADC")
            a1.set_title("Run %d  Channel %d  Entry %d  trapENFCal %.1f" % (run,chan,iList,energy))
            a1.plot(waveTS,waveBLSub,color='blue')
            a1.plot(waveTS,waveFilt,color='green',alpha=0.6)
            a1.plot(tempTSGuess,tempGuess,color='orange')
            a1.axvline(x = riseTime, color='red',linewidth=3)
            a1.axvline(x = tempGuessRT, color='green',alpha=0.8)
            a1.axvline(x = loWin*10, color='cyan')
            a1.axvline(x = hiWin*10, color='cyan')

            a2.cla()
            a2.plot(tsEdgeRT,waveEdgeRT,color='blue')
            a2.plot(smoothedTS,smoothed,color='red',linewidth=4)
            a2.plot(guessTS,guess,color='orange')

            a3.cla()
            a3.set_xlabel("Time [ns]", x=0.95, ha='right')
            a3.set_ylabel("Residual [ADC]")
            a3.plot(smoothedTS,residual1, color='red')
            # a3.plot(smoothedTS,residual2, color='cyan')
            # a3.plot(smoothedTS,residual3, color='cyan')
            a3.axhline(y = 0, color='blue', alpha=0.3)

            plt.tight_layout()
            plt.subplots_adjust(hspace=0.15)
            plt.pause(scanSpeed)

# ===========================================================================

def PlotTemplateWaveforms(waveTree,calDict):
    """ Check the wf's grabbed from an M1 calibration run. """
    # plt.style.use('mjWaveforms')  # see ~/.matplotlib/stylelib
    # optional: redefine calDict with only one or two channels for debugging
    # calDict = {672:18}

    # Make a figure
    fig = plt.figure(figsize=(9, 7), facecolor='w')
    ax = plt.subplot(111)
    ax.set_xlabel("time (ns)")
    ax.set_ylabel("ADC")
    jet = plt.get_cmap('jet')
    cNorm  = colors.Normalize(vmin=578, vmax=672)
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)

    # Plot multiple waveforms
    for chan in calDict:
        waveTree.GetEntry(calDict[chan])
        waveMGT = waveTree.event.GetWaveform(waveTree.itr)
        signal = wl.processWaveform(waveMGT)
        waveRaw = signal.GetWaveRaw()
        waveOff = signal.GetOffset()/1e9 # seconds
        waveTS = signal.GetTS()
        colorVal = scalarMap.to_rgba(chan)
        label = ('%d' % chan)
        ax.plot(waveTS, waveRaw, color=colorVal, label=label)

    # make a legend and show the plot.
    handles,labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels, loc='upper right')
    plt.show()

if __name__ == "__main__":
    main(sys.argv[1:])
