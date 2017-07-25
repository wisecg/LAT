#!/usr/local/bin/python
import sys
from ROOT import TFile,TTree,TChain,TEntryList,gDirectory,gROOT,MGTWaveform
import ROOT
import waveLibs as wl
import numpy as np

def main(argv):
    """Interactive-draw or rapid-draw waveforms that pass a given TCut.

       Bug: Doesn't always work with a TChain.  Add input files together
       with hadd and use a single TFile.
    """

    scanSpeed = 0.2
    opt1, opt2 = "", ""
    intMode, printWF, warpMode = False, False, False
    if (len(argv) >= 1): opt1 = argv[0]
    if (len(argv) >= 2): opt2 = argv[1]
    if "-i" in (opt1, opt2):
        intMode = True
        print "Interactive mode selected."
    if "-s" in (opt1, opt2):
        printWF = True
        print "Saving WF plots to current directory."
    if "-w" in (opt1, opt2):
        warpMode = True
        import matplotlib
        matplotlib.use('TkAgg')
        print "Warp-speed drawing, don't trust figure axes."
    import matplotlib.pyplot as plt

    # Set input file and cuts
    # waveTree = TChain("skimTree") # NOTE: Can't always recognize MGTWaveforms branch
    # waveTree.Add("~/project/wavelet-skim/waveletSkimDS3*")
    # inputFile = TFile("~/project/wavelet-skim/hadd/waveletSkimDS3.root")
    inputFile = TFile("~/project/wavelet-skim/waveletSkimDS3_1.root")
    waveTree = inputFile.Get("skimTree")
    print "Found",waveTree.GetEntries(),"input entries."

    # theCut = inputFile.Get("cutUsedHere").GetTitle()
    theCut = "trapENFCal > 0.8 && gain==0 && mH==1 && isGood && !muVeto && !wfDCBits && !isLNFill1 && !isLNFill2 && trapETailMin < 0 && channel!=596 && channel!=676 && channel!=676 && channel!=612 && channel!=1104 && channel!=1200 && channel!=1334 && channel!=1336"

    # S5
    # theCut += " && waveS5 > 4500 && waveS5 < 6000 && trapENFCal > 30 && trapENFCal < 50"  # DS4 blob #1
    # theCut += " && waveS5 > 1500 && waveS5 < 2500 && trapENFCal > 20 && trapENFCal < 30"  # DS4 blob #2
    # theCut += " && waveS5 > 1500 && waveS5 < 2500 && trapENFCal > 30 && trapENFCal < 40"  # DS4 blob #3
    # theCut += " && waveS5 > 1500 && waveS5 < 7000 && trapENFCal > 0 && trapENFCal < 5"    # DS4 blob #4
    # theCut += " && waveS5 < 1500 && trapENFCal > 1 && trapENFCal < 2"
    # theCut += " && (waveS5/TMath::Power(trapENFCal,1/4) < 1200 && tOffset < 10)"

    # S32/E
    # theCut += " && (waveS3-waveS2)/trapENFCal > 400 && (waveS3-waveS2)/trapENFCal < 600 && trapENFCal > 20 && trapENFCal < 30"
    # theCut += " && (waveS3-waveS2)/trapENFCal > 150 && (waveS3-waveS2)/trapENFCal < 320 && trapENFCal > 1.5 && trapENFCal < 5"

    # S34/E
    # theCut += " && (waveS3-waveS4)/trapENFCal > 25 && (waveS3-waveS4)/trapENFCal < 60 && trapENFCal > 0.8 && trapENFCal < 10"
    # theCut += " && (waveS3-waveS4)/trapENFCal < 0 && trapENFCal > 0.8 && trapENFCal < 5"

    # DS3 noisy runs - channel 692 is the culprit
    # theCut += " && tOffset < 10 && waveS5/TMath::Power(trapENFCal,1/4) < 1200 && (waveS3-waveS2)/trapENFCal > 100 && (waveS3-waveS2)/trapENFCal < 300 && ((run==16974 || run==16975 || run==16976 || run==16977 || run==16979) && channel!=692)"

    # theCut = "trapENFCal > 0.8 && gain==0 && mH==1 && isGood && !muVeto && !wfDCBits && !isLNFill1 && !isLNFill2 && trapETailMin < 0 && channel!=596 && channel!=676 && channel!=676 && channel!=612 && channel!=1104 && channel!=1200 && channel!=1334 && channel!=1336 && tOffset < 10 && waveS5/TMath::Power(trapENFCal,1/4) < 1200 && (waveS3-waveS2)/trapENFCal > 100 && (waveS3-waveS2)/trapENFCal < 300 && !(channel==692 && (run==16974 || run==16975 || run==16976 || run==16977 || run==16978 || run==16979)) && butterTime > 19000"

    # Print cut and events passing cut
    print "Using cut:\n",theCut,"\n"
    waveTree.Draw(">>elist", theCut, "entrylist")
    elist = gDirectory.Get("elist")
    waveTree.SetEntryList(elist)
    nList = elist.GetN()
    print "Found",nList,"entries passing cuts."


    # Make a figure (only setting data in the loop is faster)
    fig = plt.figure(figsize=(8,5), facecolor='w')
    a1 = plt.subplot(111)
    a2 = plt.subplot(111)
    a1.set_xlabel("time (ns)")
    a1.set_ylabel("ADC")
    p1, = a1.plot(np.ones(1), np.ones(1), color='blue')
    p2, = a2.plot(np.ones(1), np.ones(1), color='red', linewidth=2., alpha=0.5)
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
        nWaves = waveTree.MGTWaveforms.size()
        numPass = waveTree.Draw("channel",theCut,"GOFF",1,iList)
        chans = waveTree.GetV1()
        chanList = list(set(int(chans[n]) for n in xrange(numPass)))

        # Loop over hits passing cuts
        hitList = (iH for iH in xrange(nChans) if waveTree.channel.at(iH) in chanList)  # a 'generator expression'
        for iH in hitList:
            run = waveTree.run
            chan = waveTree.channel.at(iH)
            energy = waveTree.trapENFCal.at(iH)
            signal = wl.processWaveform(waveTree.MGTWaveforms.at(iH))
            waveBLSub = signal.GetWaveBLSub()
            waveFilt = signal.GetWaveFilt()
            waveTS = signal.GetTS()
            print "%d / %d  Run %d  nCh %d  chan %d  trapENF %.1f" % (iList,nList,run,nChans,chan,energy)

            # fill the figure
            p1.set_ydata(waveBLSub)
            p1.set_xdata(waveTS)
            p2.set_ydata(waveFilt)
            p2.set_xdata(waveTS)
            xmin, xmax = np.amin(waveTS), np.amax(waveTS)
            ymin, ymax = np.amin(waveBLSub), np.amax(waveBLSub)
            a1.set_xlim([xmin,xmax])
            a1.set_ylim([ymin-abs(0.1*ymin),ymax+abs(0.1*ymax)])
            plt.title("Run %d  Channel %d  Entry %d  trapENFCal %.1f" % (run,chan,iList,energy))

            if warpMode:
                # faster, requires TkAgg backend, doesn't update axes
                a1.draw_artist(a1.patch)
                a1.draw_artist(p1)
                # a1.draw_artist(a1.yaxis) # doesn't work
                fig.canvas.blit(a1.bbox)
                fig.canvas.flush_events()
            else:
                plt.pause(scanSpeed)

            if (printWF): plt.savefig("./plots/wave-%d-%d-%d.pdf" % (run,iList,chan))


if __name__ == "__main__":
    main(sys.argv[1:])
