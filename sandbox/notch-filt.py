#!/usr/local/bin/python
import sys, imp
import numpy as np
from scipy.signal import periodogram
import matplotlib.pyplot as plt
from matplotlib import gridspec
import ROOT
from ROOT import TFile,TTree,TChain,TEntryList,gDirectory,gROOT,MGTWaveform
ds = imp.load_source('DataSetInfo','../DataSetInfo.py')
wl = imp.load_source('waveLibs', '../waveLibs.py')

def main(argv):

    scanSpeed = 0.2
    opt1, opt2 = "", ""
    intMode = False
    if (len(argv) >= 1): opt1 = argv[0]
    if (len(argv) >= 2): opt2 = argv[1]
    if "-i" in (opt1, opt2):
        intMode = True
        print "Interactive mode selected."

    # Set input file and cuts
    # inputFile = TFile("~/project/cal-waves/waveSkimDS1_run10851.root")
    # inputFile = TFile("~/project/cal-waves/waveSkimDS1_run10852.root")
    # inputFile = TFile("~/project/cal-waves/waveSkimDS2_run14857_WithNLC.root")
    inputFile = TFile("~/project/cal-waves/waveSkimDS2_run14857.root")
    # inputFile = TFile("~/project/cal-waves/waveSkimDS5_run23895.root")
    waveTree = inputFile.Get("skimTree")
    print "Found",waveTree.GetEntries(),"input entries."

    theCut = inputFile.Get("theCut").GetTitle()
    # theCut += " && trapENFCal  1.1"
    # theCut = " trapENFCal > 0.5 && trapENFCal < 10"

    # Print cut and events passing cut
    print "Using cut:\n",theCut,"\n"
    waveTree.Draw(">>elist", theCut, "entrylist")
    elist = gDirectory.Get("elist")
    waveTree.SetEntryList(elist)
    nList = elist.GetN()
    print "Found",nList,"entries passing cuts."

    # Make a figure
    fig = plt.figure(figsize=(12,7), facecolor='w')
    p0 = plt.subplot(111)
    # gs = gridspec.GridSpec(2, 2, height_ratios=[2, 3])
    # p0 = plt.subplot(gs[:,0]) # power spec
    # p1 = plt.subplot(gs[0,1]) # baseline
    # p2 = plt.subplot(gs[1,1]) # full waveform
    plt.show(block=False)

    # Loop over events
    iList = -1
    while(True):
        iList += 1
        if intMode==True and iList !=0:
            value = raw_input()
            if value=='q': break
            if value=='p': iList -= 2  # previous
            if (value.isdigit()): iList = int(value) # go to entry
        if iList >= elist.GetN(): break

        entry = waveTree.GetEntryNumber(iList);
        waveTree.LoadTree(entry)
        waveTree.GetEntry(entry)
        nChans = waveTree.channel.size()
        nWFs = waveTree.MGTWaveforms.size()

        if (nWFs==0):
            print "Error - nWFs:",nWFs,"nChans",nChans
            continue

        numPass = waveTree.Draw("channel",theCut,"GOFF",1,iList)
        chans = waveTree.GetV1()
        chanList = list(set(int(chans[n]) for n in xrange(numPass)))

        # Loop over hits passing cuts
        hitList = (iH for iH in xrange(nChans) if waveTree.channel.at(iH) in chanList)  # a 'generator expression'
        for iH in hitList:
            run = waveTree.run
            chan = waveTree.channel.at(iH)
            dataENFCal = waveTree.trapENFCal.at(iH)
            wf = waveTree.MGTWaveforms.at(iH)
            signal = wl.processWaveform(wf)
            data = signal.GetWaveRaw()
            data_blSub = signal.GetWaveBLSub()
            dataTS = signal.GetTS()
            print "%d / %d  Run %d  nCh %d  chan %d  trapENF %.1f" % (iList,nList,run,nChans,chan,dataENFCal)

            # bl = data[0]
            # idx = np.where(data > 100000)
            # data[idx] = bl

            """
            # -- apply notch filter --
            # This algorithm was adapted from code by Ben Shanks, who was nice enough to share.

            # only look at power spectrum before the rising edge
            bl = data[:800]
            xf,power = periodogram(bl, fs=1E8)

            # find the preferred notch frequency: highest freq peak after some min freq
            min_freq = 0.25E7
            xf_idx_min = np.argmax(xf > min_freq)

            notch_freq = xf[ np.argmax(power[xf_idx_min:])+xf_idx_min ]

            # express that frequency as a fraction of the nyquist frequency (half the sampling rate)
            nyquist = 0.5*1E8
            w0 = notch_freq/nyquist

            # "quality factor" which determines width of notch
            Q = 10

            num, den = signal.iirnotch(w0, Q)

            ax1.plot(xf, power)
            ax1.set_ylabel('Noise Amplitude', color='b')
            ax1.tick_params('y', colors='b')

            ax_twin = ax1.twinx()
            w, h = signal.freqz(num, den)
            freq = w*1E8/(2*np.pi)
            ax_twin.plot(freq, 20*np.log10(abs(h)), color='red', ls=":", label="Notch")
            ax_twin.set_ylabel('Filter Amplitude [dB]', color='r')
            ax_twin.tick_params('y', colors='r')
            ax_twin.set_title(label="Q=%0.1f"%Q)
            """


            # -- plotting --

            # -- power spec --
            # p0.cla()
            # p0.axvline(notch_freq, color="green", ls="--")


            # -- baseline --
            p0.cla()
            # p0.margins(x=0)
            p0.plot(dataTS,data,color='blue',label='data')
            p0.set_title("Run %d  Entry %d  Channel %d  Energy %.2f" % (run,iList,chan,dataENFCal))
            p0.legend(loc='best')

            # -- full waveform --
            # p2.cla()

            plt.pause(scanSpeed)


if __name__ == "__main__":
    main(sys.argv[1:])
