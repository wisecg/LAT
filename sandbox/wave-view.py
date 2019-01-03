#!/usr/bin/env python
import sys, pywt
from ROOT import TFile,TTree,TChain,TEntryList,gDirectory,gROOT,MGTWaveform
import ROOT
import dsi
import waveLibs as wl
import numpy as np
from scipy.signal import butter, lfilter

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
        print("Interactive mode selected.")
    if "-s" in (opt1, opt2):
        printWF = True
        print("Saving WF plots to current directory.")
    if "-w" in (opt1, opt2):
        warpMode = True
        import matplotlib
        matplotlib.use('TkAgg')
        print("Warp-speed drawing, don't trust figure axes.")
    import matplotlib.pyplot as plt
    plt.style.use('../clint.mpl')
    from matplotlib.colors import LogNorm, Normalize

    # Set input file and cuts
    dsNum = 5
    # waveTree = TChain("skimTree") # NOTE: Can't always recognize MGTWaveforms branch
    # waveTree.Add("~/project/lat/latSkimDS1*")
    inputFile = TFile("/Users/wisecg/project/cal/lat/latSkimDS5_run20044_0.root")
    # inputFile = TFile("~/project/cal-waves/waveSkimDS2_run14857.root")
    # inputFile = TFile("~/project/cal-waves/waveSkimDS5_run23895.root")

    waveTree = inputFile.Get("skimTree")
    print("Found",waveTree.GetEntries(),"input entries.")

    theCut = inputFile.Get("theCut").GetTitle()
    # theCut += " && trapENFCal  1.1"
    # theCut = " trapENFCal > 0.5 && trapENFCal < 10"

    # Print cut and events passing cut
    print("Using cut:\n",theCut,"\n")
    waveTree.Draw(">>elist", theCut, "entrylist")
    elist = gDirectory.Get("elist")
    waveTree.SetEntryList(elist)
    nList = elist.GetN()
    print("Found",nList,"entries passing cuts.")


    # Make a figure (only setting data in the loop is faster)
    fig = plt.figure()
    a1 = plt.subplot(111)
    a1.set_xlabel("Time (ns)", ha='right', x=1)
    a1.set_ylabel("Voltage (ADC)", ha='right', y=1)
    p1, = a1.plot(np.ones(1), np.ones(1), color='blue', alpha=0.8, label='data')
    p2, = a1.plot(np.ones(1), np.ones(1), color='magenta', linewidth=1.5, label='denoised')
    p3, = a1.plot(np.ones(1), np.ones(1), color='red', linewidth=3, label='lowpass')
    # a1.legend(loc=4)
    plt.show(block=False)

    # Loop over events
    iList = -1
    while(True):
        iList += 1
        if intMode==True and iList !=0:
            value = input()
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
            print("Error - nWFs:",nWFs,"nChans",nChans)
            continue

        numPass = waveTree.Draw("channel",theCut,"GOFF",1,iList)
        chans = waveTree.GetV1()
        chanList = list(set(int(chans[n]) for n in range(numPass)))

        # Loop over hits passing cuts
        hitList = (iH for iH in range(nChans) if waveTree.channel.at(iH) in chanList)
        for iH in hitList:
            run = waveTree.run
            chan = waveTree.channel.at(iH)
            energy = waveTree.trapENFCal.at(iH)
            wf = waveTree.MGTWaveforms.at(iH)

            truncLo, truncHi = 0, 2
            if dsNum==6 or dsNum==2: truncLo = 4
            signal = wl.processWaveform(wf,truncLo,truncHi)

            waveBLSub = signal.GetWaveBLSub()
            # waveRaw = signal.GetWaveRaw()
            waveTS = signal.GetTS()
            print("%d / %d  Run %d  nCh %d  chan %d  trapENF %.1f" % (iList,nList,run,nChans,chan,energy))


            # standard energy trapezoid
            eTrap = wl.trapFilter(waveBLSub,400,250,-7200)

            nPad = len(waveBLSub)-len(eTrap)
            eTrap = np.pad(eTrap, (nPad,0), 'constant')
            eTrapTS = np.arange(0, len(eTrap)*10., 10)
            ePickoff = eTrapTS[nPad + 400 + 200]
            # plt.axvline(ePickoff, c='m')

            # trigger trap
            tTrap = wl.trapFilter(waveBLSub,100,150,-7200)
            tTrap = np.pad(tTrap, (nPad,0), 'constant')
            tTrapTS = np.arange(0, len(tTrap)*10., 10)

            # -- fill the figure --
            p1.set_ydata(waveBLSub)
            p1.set_xdata(waveTS)

            p2.set_ydata(tTrap)
            p2.set_ydata(tTrapTS)

            p3.set_ydata(eTrap)
            p3.set_xdata(eTrapTS)

            xmin, xmax = np.amin(waveTS), np.amax(waveTS)
            ymin, ymax = np.amin(waveBLSub), np.amax(waveBLSub)
            a1.set_xlim([xmin,xmax])
            a1.set_ylim([ymin-abs(0.1*ymin),ymax+abs(0.1*ymax)])
            # plt.title("Run %d  Channel %d  Entry %d  trapENFCal %.1f" % (run,chan,iList,energy))

            if warpMode:
                # superfast, requires TkAgg backend, doesn't update axes
                a1.draw_artist(a1.patch)
                a1.draw_artist(p1)
                fig.canvas.blit(a1.bbox)
                fig.canvas.flush_events()
            else:
                plt.pause(scanSpeed)

            if (printWF): plt.savefig("./plots/wave-%d-%d-%d.pdf" % (run,iList,chan))


if __name__ == "__main__":
    main(sys.argv[1:])
