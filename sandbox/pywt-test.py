#!/usr/bin/env python3
import sys, os, imp, pywt
sys.argv.append("-b")
wl = imp.load_source('waveLibs',os.environ['LATDIR']+'/waveLibs.py')
from ROOT import TChain,TFile,TTree,TEntryList,gDirectory,gROOT,std,TNamed,TObject
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal as sg

def main(argv):

    print("Hi Clint")

    # User args
    intMode, batMode = True, False
    for i,opt in enumerate(argv):
        if opt == "-b":
            intMode, batMode = False, True

    # Set file I/O and cuts
    inFile = "~/project/special/waves/waveSkimDS0_run4549.root"
    inputFile = TFile(inFile)
    waveTree = inputFile.Get("skimTree")
    # waveTree = TChain("skimTree")
    # waveTree.Add("~/project/wavelet-skim/hadd/waveletSkimDS3.root")
    print("Found",waveTree.GetEntries(),"input entries.")

    theCut = "trapENFCal > 0.8"
    print("Using the cut:\n",theCut)

    waveTree.Draw(">>elist", theCut, "GOFF entrylist")
    elist = gDirectory.Get("elist")
    waveTree.SetEntryList(elist)
    nList = elist.GetN()
    print("Found",nList,"entries passing cuts.")

    # Begin loop over events
    iList = -1
    while True:
        iList += 1
        if intMode and iList != 0:
            value = raw_input()
            if value=='q': break        # quit
            if value=='p': iList -= 2   # go to previous
            if (value.isdigit()):
                iList = int(value)      # go to entry number
        elif not intMode and batMode:
            plt.pause(0.001)            # rapid-draw mode
        if iList >= nList: break        # bail out, goose!

        entry = waveTree.GetEntryNumber(iList);
        waveTree.LoadTree(entry)
        waveTree.GetEntry(entry)
        nChans = waveTree.channel.size()

        # Loop over hits passing cuts
        numPass = waveTree.Draw("channel",theCut,"GOFF",1,iList)
        chans = waveTree.GetV1()
        chanList = list(set(int(chans[n]) for n in range(numPass)))
        hitList = (iH for iH in range(nChans) if waveTree.channel.at(iH) in chanList)  # a 'generator expression'
        for iH in hitList:
            run = waveTree.run
            iEvent = waveTree.iEvent
            chan = waveTree.channel.at(iH)
            energy = waveTree.trapENFCal.at(iH)
            wf = waveTree.MGTWaveforms.at(iH)

            if wf.GetID() != chan:
                print("Warning!!  Vector matching failed!  iList %d  run %d  iEvent %d" % (iList,run,iEvent))
                break

            # Let's start the show - grab a waveform.
            # Remove first 4 samples when we have multisampling
            # Remove last 2 samples to get rid of the ADC spike at the end of all wf's.
            truncLo, truncHi = 0, 2
            signal = wl.processWaveform(wf,truncLo,truncHi)
            data = signal.GetWaveRaw()
            data_blSub = signal.GetWaveBLSub()
            dataTS = signal.GetTS()
            # tOffset[iH] = signal.GetOffset()
            dataBL,dataNoise = signal.GetBaseNoise()

            # wavelet packet transform
            wp = pywt.WaveletPacket(data_blSub, 'db2', 'symmetric', maxlevel=4)
            nodes = wp.get_level(4, order='freq')
            wpCoeff = np.array([n.data for n in nodes],'d')
            wpCoeff = abs(wpCoeff)

            # wavelet parameters
            # First get length of wavelet on the time axis, the scale axis will
            # always be the same due to the number of levels in the wavelet
            wpLength = len(wpCoeff[1,:])
            waveS1 = np.sum(wpCoeff[0:1,1:wpLength//4+1]) # python3 : floor division (//) returns an int

            print("%d  run %d  iEvt %d  iH %d  chan %d  ene %.2f  waveS1 %.2f" % (iList,run,iEvent,iH,chan,energy,waveS1))

            if batMode: continue

            # Make a figure (only setting data in the loop is faster)
            plt.close("all")
            fig = plt.figure(figsize=(7,7), facecolor='w')
            a1 = plt.subplot(211)
            a1.set_xlabel("Time (ns)")
            a1.set_ylabel("ADC")
            a1.plot(waveTS,waveBLSub,color='blue')
            a1.plot(waveTS,waveFilt,color='red',linewidth=2.,alpha=0.5)

            # a2 = plt.subplot(212)
            # a2.plot(waveTS,waveFiltDeriv,color='red')
            # a2.axvline(x=butterMaxTS,color='green',alpha=0.5)

            a1.set_title("Run %d  Ch %d  Energy %.2f  \n S5 %.2f  (S3-S2)/E %.2f  (S3-S4)/E %.2f" % (run,chan,energy,s5,((s3-s2)/energy),((s3-s4)/energy)))

            # Wavelet plot
            a3 = plt.subplot(212)
            a3.cla()
            a3.imshow(waveletYTrans, interpolation='nearest',
                aspect="auto", origin="lower",extent=[0, 1, 0, len(waveletYTrans)],cmap='jet')

            plt.tight_layout()
            plt.subplots_adjust(hspace=.2,top=0.92)
            plt.pause(0.00001)
            # pdf.savefig()

if __name__ == "__main__":
    main(sys.argv[1:])
