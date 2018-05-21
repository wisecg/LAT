#!/usr/bin/env python
#
# wave-movie.py
# Clint Wiseman, USC
# 12 July 2017
#
# Modified to take in a chain, apply a cut onto a dummy tree and then load the dummy tree
# this allows loading entire datasets as a TChain

import sys, imp, glob, os
sys.argv.append("-b") # kill all interactive crap

import matplotlib as mpl
mpl.use('Agg')
sys.argv.append("-b")
import matplotlib.pyplot as plt
from matplotlib import animation

import waveLibs as wl
from ROOT import TFile,TTree,TChain,TEntryList,gDirectory,gROOT,MGTWaveform
import ROOT
import numpy as np


def main(argv):
    """
    Matplotlib animation tutorial:
    http://jakevdp.github.io/blog/2012/08/18/matplotlib-animation-tutorial/
    Requires ffmpeg.  (brew install ffmpeg)
    """
    chans = [1176, 672, 648, 660]
    chan = chans[0]

    dummyTree = ROOT.TChain("skimTree")
    # dummyTree.Add("/global/homes/w/wisecg/project/cal/lat/latSkimDS5_run2549*.root")
    # dummyTree.Add("/global/homes/w/wisecg/project/cal/lat/latSkimDS5_run2550*.root")
    dummyTree.Add("/global/homes/w/wisecg/project/bkg/waves/waveSkimDS6_*.root")
    print("Found %d entries." % (dummyTree.GetEntries()))

    theCut = "trapENFCal > 2000 && trapENFCal < 2200 && channel==%d && avse>-1" % chan

    gatTree = TTree()
    gatTree = dummyTree.CopyTree(theCut)
    outFile = "../plots/movie-5c-%d.mp4" % chan

    # Print cut and events passing cut
    print("Using cut:\n",theCut,"\n")
    gatTree.Draw(">>elist", theCut, "entrylist")
    elist = gDirectory.Get("elist")
    gatTree.SetEntryList(elist)
    nList = elist.GetN()
    print("Found",nList,"entries passing cuts.")

    nWFLimit = 1000 if nList > 1000 else nList

    # First set up the figure, the axis, and the plot element we want to animate
    fig = plt.figure(figsize=(20,10), facecolor='w')
    fig.set_size_inches(20,10, True)
    a1 = plt.subplot(111)
    a1.set_xlabel("Time (ns)")
    a1.set_ylabel("ADC")
    p1, = a1.plot(np.ones(1), np.ones(1), color='blue')

    # initialization function: plot the background of each frame
    def init():
        p1.set_data([],[])
        return p1,

    # animation function.  This is called sequentially (it's the loop over events.)
    def animate(iList):

        entry = gatTree.GetEntryNumber(iList);
        gatTree.LoadTree(entry)
        gatTree.GetEntry(entry)
        nChans = gatTree.channel.size()
        numPass = gatTree.Draw("channel",theCut,"GOFF",1,iList)
        chans = gatTree.GetV1()
        chanList = list(set(int(chans[n]) for n in range(numPass)))

        # Loop over hits passing cuts
        hitList = (iH for iH in range(nChans) if gatTree.channel.at(iH) in chanList)  # a 'generator expression'
        for iH in hitList:
            wf = gatTree.MGTWaveforms.at(iH)
            iEvent = entry
            run = gatTree.run
            chan = gatTree.channel.at(iH)
            energy = gatTree.trapENFCal.at(iH)
            dcr = gatTree.dcr99.at(iH)
            # fs = gatTree.fitSlo.at(iH)
            # rn = gatTree.riseNoise.at(iH)
            signal = wl.processWaveform(wf)
            waveRaw = signal.GetWaveRaw()
            waveTS = signal.GetTS()

            baseline = np.sum(waveRaw[:50])/50

            # fill the figure
            p1.set_ydata(waveRaw)
            p1.set_xdata(waveTS)
            # plt.title("Run %d  Ch %d  Entry %d  trapENFCal %.1f  fitSlo %.2f  riseNoise %.2f" % (run,chan,iList,energy, fs, rn))
            plt.title("Run %d  Ch %d  Entry %d  trapENFCal %.1f  dcr99 %.8f" % (run,chan,iList,energy,dcr))

            # dynamically scale the axes
            xmin, xmax = np.amin(waveTS), np.amax(waveTS)
            # a1.set_xlim([xmin,xmax])
            a1.set_xlim([9000,xmax])
            # if energy > 1000:
                # a1.set_ylim(0,1000)
            # else:
            ymin, ymax = np.amin(waveRaw), np.amax(waveRaw)
            # a1.set_ylim([ymin-abs(0.1*ymin),ymax+abs(0.1*ymax)])
            a1.set_ylim([ymax-abs(0.35*ymax),ymax+abs(0.2*ymax)])
            # a1.set_ylim([5000,6000])

            # print(progress)
            # print("%d / %d  Run %d  nCh %d  chan %d  trapE %.1f  samp %d" % (iList,nList,run,nChans,chan,energy,wf.GetLength()))
            if iList%500 == 0 and iList!=0:
                print("%d / %d entries saved (%.2f %% done)." % (iList,nList,100*(float(iList)/nList)))
            return p1,

    # call the animator.  blit=True means only re-draw the parts that have changed.
    # anim = animation.FuncAnimation(fig, animate, init_func=init, frames=elist.GetN(), interval=0, blit=True)
    anim = animation.FuncAnimation(fig, animate, init_func=init, frames=nWFLimit, interval=0, blit=True)

    # save the animation as an mp4.  This requires ffmpeg or mencoder to be
    # installed.  The extra_args ensure that the x264 codec is used, so that
    # the video can be embedded in html5.  You may need to adjust this for
    # your system: for more information, see
    # http://matplotlib.sourceforge.net/api/animation_api.html
    anim.save(outFile, fps=20)#, extra_args=['-vcodec', 'libx264'])


def trapFilt(data,ramp=400,flat=200,decay=72,trapThresh=2.,padAfter=False):

    # decay is unused
    # no charge trapping
    # no fixed time pickoff
    # no pole zero correction

    trap = np.zeros(len(data))
    trigger, triggerTS = False, 0

    # compute a moving average
    for i in range(len(data)-1000):
        w1 = ramp
        w2 = ramp+flat
        w3 = ramp+flat+ramp

        r1 = np.sum(data[i:w1+i])/ramp
        r2 = np.sum(data[w2+i:w3+i])/ramp

        if not padAfter:
            trap[i+1000] = r2 - r1
            if trap[i+1000] > trapThresh and not trigger:
                triggerTS = (i + 1000) * 10
                trigger = True
        else:
            trap[i] = r2 - r1
            if trap[i] > trapThresh and not trigger:
                triggerTS = i * 10
                trigger = True

    # if trigger: print("triggered!",triggerTS)
    return trap, trigger, triggerTS


if __name__ == "__main__":
    main(sys.argv[1:])
