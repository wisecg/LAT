#!/usr/bin/env python
import sys
import numpy as np
import matplotlib.pyplot as plt
from ROOT import TFile,TTree,TCanvas,TEntryList,gDirectory,MGTWaveform
import waveLibs as wl

def main(argv):
    scanSpeed = 0.0001
    opt1, opt2 = "", ""
    intMode, batMode = False, False
    if (len(argv) >= 1): opt1 = argv[0]
    if (len(argv) >= 2): opt2 = argv[1]
    if "-i" in (opt1, opt2):
        intMode = True
        print "Interactive mode selected."
    if "-b" in (opt1,opt2):
        batMode = True
        print "Batch mode selected."

    # Input: forced acquisition runs
    gatFile = TFile("~/project/mjddatadir/gatified/mjd_run13071.root")
    bltFile = TFile("~/project/mjddatadir/built/OR_run13071.root")
    waveTree = gatFile.Get("mjdTree")
    bltTree = bltFile.Get("MGTree")
    waveTree.AddFriend(bltTree)
    print "Found",waveTree.GetEntries(),"input entries."

    # Set cuts
    theCut = "rawWFMax < 200 && rawWFMin > -20"
    theCut += " && Entry$ < 10000"
    waveTree.Draw(">>elist", theCut, "entrylist")
    elist = gDirectory.Get("elist")
    waveTree.SetEntryList(elist)
    nList = elist.GetN()
    print "Using cut:\n",theCut,"\n"
    print "Found",nList,"entries passing cuts."

    # Loop figure (diagnostic)
    fig1 = plt.figure(figsize=(8,5), facecolor='w')
    p1 = plt.subplot(211)
    p2 = plt.subplot(212)
    if not batMode: plt.show(block=False)

    # Quantities to write to file
    avgPsd, xPsd = 0, 0 # pyplot.psd
    avgPwrSpec, xPwrSpec = 0, 0 # wl.fourierTransform
    data_forceAcq = 0 # single baseline
    data_fft = 0 # fft of single baseline

    nWF = 0
    ents = elist.GetN()

    # Loop over events
    iList = -1
    while(True):
        iList += 1
        if intMode==True and iList !=0:
            value = raw_input()
            if value=='q': break
            if value=='p': iList -= 2  # previous
            if (value.isdigit()): iList = int(value) # go to entry
        if iList >= ents: break

        # print progress
        if (np.fabs((float(iList)/ents)%0.1) < 0.001):
            print "%i done" % (100*iList/ents), np.fabs(float(iList)/ents)%0.1

        entry = waveTree.GetEntryNumber(iList)
        waveTree.LoadTree(entry)
        waveTree.GetEntry(entry)
        nChans = waveTree.channel.size()
        numPass = waveTree.Draw("channel",theCut,"GOFF",1,iList)
        chans = waveTree.GetV1()
        chanList = list(set(int(chans[n]) for n in xrange(numPass)))
        event = bltTree.event  # grab waveforms from built data

        # Loop over hits passing cuts
        hitList = (iH for iH in xrange(nChans) if waveTree.channel.at(iH) in chanList)  # a 'generator expression'
        for iH in hitList:

            # Get signal
            signal = wl.processWaveform(event.GetWaveform(iH)) # use when we have gat+blt data
            data = signal.GetWaveBLSub()
            dataTS = signal.GetTS() * 10.

            # Get power spectrum
            xPwrSpec,_,ypwr = wl.fourierTransform(data)
            avgPwrSpec += ypwr
            nWF += 1

            # Get full FFT (real, cplx)
            data_fft = np.fft.fft(data)
            data_forceAcq = data

            # power spectral density
            psd, xPsd = plt.psd(data, Fs=1e8, NFFT=2048, pad_to=2048, visible=True)
            # avgPsd += np.sqrt(psd)
            avgPsd += psd

            # Fill the figure
            if not batMode:
                p1.cla()
                p1.plot(dataTS,data,color='blue')
                p2.cla()
                plt.yscale('log')
                # p2.plot(xPwrSpec,avgPwrSpec,color='blue')
                p2.loglog(xPsd,avgPsd,color='blue')
                plt.pause(scanSpeed)


    # Plot average power spectrum
    fig2 = plt.figure(figsize=(8,5), facecolor='w')
    p3 = plt.subplot(111)
    plt.yscale('log')
    # plt.xscale('log')
    p3.set_xlabel("frequency")
    avgPwrSpec /= nWF
    idx = np.where(xPwrSpec > 10000) # low frequency cutoff
    p3.plot(xPwrSpec[idx],avgPwrSpec[idx],color='blue')
    fig2.savefig("./plots/powerSpectrum_run13071_linscale.pdf")

    # Get noise power density
    print "noise power density N_0: ", np.sum(avgPwrSpec[idx])

    # Plot power spectral density
    avgPsd /= nWF
    p3.cla()
    p3.loglog(xPsd,avgPsd,color='blue')
    fig2.savefig("./plots/psd_run13071.pdf")

    # Save stuff to file for optimal matched filter (ligo-mjd.py)
    np.savez("./data/fft_forcedAcqDS1.npz",avgPsd,xPsd,avgPwrSpec,xPwrSpec,data_forceAcq,data_fft)


if __name__ == "__main__":
    main(sys.argv[1:])
