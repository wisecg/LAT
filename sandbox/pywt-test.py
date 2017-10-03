#!/usr/bin/env python
import sys, time, os
from ROOT import TChain,TFile,TTree,TEntryList,gDirectory,gROOT,std,TNamed,TObject
import numpy as np
import waveLibs as wl
import matplotlib.pyplot as plt
from scipy import signal as sg

def main(argv):

    # User args
    intMode, batMode = True, False
    for i,opt in enumerate(argv):
        if opt == "-b":
            intMode, batMode = False, True

    # Set file I/O and cuts
    # inFile = "~/project/wavelet-skim/waveletSkimDS4_0.root"
    # inputFile = TFile(inFile)
    # waveTree = inputFile.Get("skimTree")
    waveTree = TChain("skimTree")
    waveTree.Add("~/project/wavelet-skim/hadd/waveletSkimDS3.root")
    print "Found",waveTree.GetEntries(),"input entries."

    # theCut = inputFile.Get("theCut").GetTitle()
    # theCut = "trapENFCal > 0.8 && gain==0 && mH==1 && isGood && !muVeto && !wfDCBits && !isLNFill1 && !isLNFill2 && trapETailMin < 0 && channel!=596 && channel!=676 && channel!=676 && channel!=612 && channel!=1104 && channel!=1200 && channel!=1334 && channel!=1336"
    # theCut += " && trapENFCal > 1.2 && trapENFCal < 4"

    theCut = "trapENFCal > 0.8 && gain==0 && mH==1 && isGood && !muVeto && !wfDCBits && !isLNFill1 && !isLNFill2 && trapETailMin < 0 && channel!=596 && channel!=676 && channel!=676 && channel!=612 && channel!=1104 && channel!=1200 && channel!=1334 && channel!=1336 && tOffset < 10 && waveS5/TMath::Power(trapENFCal,1/4) < 1200 && (waveS3-waveS2)/trapENFCal > 100 && (waveS3-waveS2)/trapENFCal < 300 && !(channel==692 && (run==16974 || run==16975 || run==16976 || run==16977 || run==16978 || run==16979)) && trapENFCal > 2000"

    print "Using the cut:\n",theCut

    waveTree.Draw(">>elist", theCut, "GOFF entrylist")
    elist = gDirectory.Get("elist")
    waveTree.SetEntryList(elist)
    nList = elist.GetN()
    print "Found",nList,"entries passing cuts."

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
        chanList = list(set(int(chans[n]) for n in xrange(numPass)))
        hitList = (iH for iH in xrange(nChans) if waveTree.channel.at(iH) in chanList)  # a 'generator expression'
        for iH in hitList:
            run = waveTree.run
            iEvent = waveTree.iEvent
            chan = waveTree.channel.at(iH)
            energy = waveTree.trapENFCal.at(iH)
            wf = waveTree.MGTWaveforms.at(iH)

            if wf.GetID() != chan:
                print "Warning!!  Vector matching failed!  iList %d  run %d  iEvent %d" % (iList,run,iEvent)
                break

            signal = wl.processWaveform(wf,opt='full',N=2, Wn=0.01) # N=2, Wn=0.01 is doing well
            waveBLSub = signal.GetWaveBLSub()
            waveFilt = signal.GetWaveFilt()
            waveTS = signal.GetTS()
            # waveletYOrig, waveletYTrans = signal.GetWavelet()
            _,waveletYTrans = wl.waveletTransform(waveBLSub,level=2)

            # waveFTX, waveFTY, waveFTPow = signal.GetFourier()
            # waveTrap = signal.GetTrapezoid()
            # waveTrapF = signal.GetFiltTrapezoid()
            # waveTrapX = np.linspace(0, len(waveTrap), num=len(waveTrap))
            # _,waveletFilt = wl.waveletTransform(waveFilt,level=3)  # testing other levels on the filtered WF

            # Wavelet cut parameters (only valid for numLevels=4)
            s0 = np.sum(waveletYTrans[0:1,1:-1])
            s1 = np.sum(waveletYTrans[0:1,1:33])
            s2 = np.sum(waveletYTrans[0:1,33:65])
            s3 = np.sum(waveletYTrans[0:1,65:97])
            s4 = np.sum(waveletYTrans[0:1,97:-1])
            s5 = np.sum(waveletYTrans[2:7,1:-1])
            wpar4 = np.amax(waveletYTrans[0:1,1:-1])

            # Smoothed derivative params
            waveFiltDeriv = wl.wfDerivative(waveFilt)
            butterMaxVal = np.amax(waveFiltDeriv)
            butterMaxTS = np.argmax(waveFiltDeriv[100:-10])*10 + signal.GetOffset() + 1000
            sigOffset = signal.GetOffset()
            print "max %.3f  ts %.0f ns  offset %.0f ns" % (butterMaxVal, butterMaxTS, sigOffset)

            if (batMode==False):
                print "i %d  run %d  iEvent(i) %d  iH(j+1) %d/%d  chan %d  energy %.2f  bt %.0f  b2 %.0f" % (iList,run,iEvent,iH+1,nChans,chan,energy,waveTree.butterTime.at(iH),butterMaxTS)


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
