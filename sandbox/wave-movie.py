#!/usr/bin/env python
#
# wave-movie.py
# Clint Wiseman, USC
# 12 July 2017

import sys, imp, glob, os
sys.argv.append("-b") # kill all interactive crap
ds = imp.load_source('DataSetInfo','../DataSetInfo.py')
wl = imp.load_source('waveLibs','../waveLibs.py')
from ROOT import TFile,TTree,TChain,TEntryList,gDirectory,gROOT,MGTWaveform,MJTMSWaveform,GATDataSet
import ROOT
import numpy as np
from scipy.signal import butter, lfilter
import matplotlib.pyplot as plt
from matplotlib import animation

def main(argv):
    """
    Matplotlib animation tutorial:
    http://jakevdp.github.io/blog/2012/08/18/matplotlib-animation-tutorial/
    Requires ffmpeg.  (brew install ffmpeg)
    """

    # Set input file and cuts
    # gatFile = TFile("~/project/mjddatadir/gatified/mjd_run27012.root") # generate w/ process_mjd_data_p1
    # bltFile = TFile("~/project/mjddatadir/built/OR_run27012.root")
    # gatTree = gatFile.Get("mjdTree")
    # bltTree = bltFile.Get("MGTree")
    run = 23725
    ds = GATDataSet(run)
    gatTree = ds.GetGatifiedChain()
    bltTree = ds.GetBuiltChain()
    gatTree.AddFriend(bltTree)
    print "Found",gatTree.GetEntries(),"input entries."

    # Get first timestamp
    gatTree.GetEntry(0)
    # print type(gatTree.localTime.localTime)


    theCut = "Entry$<50"
    # theCut = "channel == 674" # C1P7D3
    # theCut = "trapE < 6000 && channel==674 && timestamp/1e8 > 148359"
    # theCut = "channel==632 && timestamp/1e8 > 148359"

    outFile = "../plots/movie_run%d.mp4" % run

    # Print cut and events passing cut
    print "Using cut:\n",theCut,"\n"
    gatTree.Draw(">>elist", theCut, "entrylist")
    elist = gDirectory.Get("elist")
    gatTree.SetEntryList(elist)
    nList = elist.GetN()
    print "Found",nList,"entries passing cuts."



    # First set up the figure, the axis, and the plot element we want to animate
    fig = plt.figure(figsize=(11,5), facecolor='w')
    a1 = plt.subplot(111)
    # a2 = plt.subplot(111)
    # a1.set_xlim(0,50000)
    # a1.set_ylim(0,1000)
    a1.set_xlabel("time (ns)")
    a1.set_ylabel("ADC")
    p1, = a1.plot(np.ones(1), np.ones(1), color='blue')
    p2, = a1.plot(np.ones(1), np.ones(1), color='red')

    # initialization function: plot the background of each frame
    def init():
        p1.set_data([],[])
        p2.set_data([],[])
        return p1,p2,

    # animation function.  This is called sequentially (it's the loop over events.)
    def animate(iList):

        entry = gatTree.GetEntryNumber(iList);
        gatTree.LoadTree(entry)
        gatTree.GetEntry(entry)
        nChans = gatTree.channel.size()
        numPass = gatTree.Draw("channel",theCut,"GOFF",1,iList)
        chans = gatTree.GetV1()
        chanList = list(set(int(chans[n]) for n in xrange(numPass)))
        event = bltTree.event

        locTime = gatTree.localTime.localTime # wow, srsly?

        # Loop over hits passing cuts
        hitList = (iH for iH in xrange(nChans) if gatTree.channel.at(iH) in chanList)  # a 'generator expression'
        for iH in hitList:

            # wf = MGTWaveform()
            iEvent = 0
            wf_downsampled = event.GetWaveform(iH)
            wf_regular = event.GetAuxWaveform(iH)

            wf = MJTMSWaveform(wf_downsampled,wf_regular)


            iEvent = entry
            run = gatTree.run
            chan = gatTree.channel.at(iH)
            energy = gatTree.trapE.at(iH)
            # timestamp = gatTree.timestamp.at(iH) / 1e8
            # locTime = timestamp - firstTime


            signal = wl.processWaveform(wf)
            waveRaw = signal.GetWaveRaw()
            waveBLSub = signal.GetWaveBLSub()
            waveTS = signal.GetTS()

            baseline = np.sum(waveRaw[:50])/50

            # B, A = butter(2,1e6/(1e8/2), btype='lowpass')
            # data_lPass = lfilter(B, A, waveRaw)

            # simple trap filter (broken)
            trapTS = waveTS
            ADCThresh = 1.
            trap, trigger, triggerTS = trapFilt(waveBLSub,trapThresh=ADCThresh)
            trap = trap + baseline

            # fill the figure
            p1.set_ydata(waveRaw)
            p1.set_xdata(waveTS)
            p2.set_ydata(trap)
            p2.set_xdata(trapTS)
            plt.title("Run %d  Channel %d  Entry %d  trapE %.1f  locTime %.1f s" % (run,chan,iList,energy,locTime))

            # dynamically scale the axes
            xmin, xmax = np.amin(waveTS), np.amax(waveTS)
            a1.set_xlim([xmin,xmax])
            if energy < 1000:
                a1.set_ylim(0,1000)
            else:
                ymin, ymax = np.amin(waveRaw), np.amax(waveRaw)
                a1.set_ylim([ymin-abs(0.1*ymin),ymax+abs(0.1*ymax)])

            # print progress
            # print "%d / %d  Run %d  nCh %d  chan %d  trapE %.1f  samp %d" % (iList,nList,run,nChans,chan,energy,wf.GetLength())
            if iList%500 == 0 and iList!=0:
                print "%d / %d entries saved (%.2f %% done)." % (iList,nList,100*(float(iList)/nList))

            return p1,p2,

    # call the animator.  blit=True means only re-draw the parts that have changed.
    anim = animation.FuncAnimation(fig, animate, init_func=init, frames=elist.GetN(), interval=0, blit=True)

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

    # if trigger: print "triggered!",triggerTS
    return trap, trigger, triggerTS


if __name__ == "__main__":
    main(sys.argv[1:])
