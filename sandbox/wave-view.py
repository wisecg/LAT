#!/usr/local/bin/python
import sys, pywt, imp
from ROOT import TFile,TTree,TChain,TEntryList,gDirectory,gROOT,MGTWaveform
import ROOT
# import waveLibs as wl
ds = imp.load_source('DataSetInfo','../DataSetInfo.py')
wl = imp.load_source('waveLibs', '../waveLibs.py')
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
    # waveTree.Add("~/project/lat/latSkimDS1*")
    inputFile = TFile("~/project/cal-waves/waveSkimDS1_run10851.root")
    # inputFile = TFile("~/project/cal-waves/waveSkimDS2_run14857.root")
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


    # Make a figure (only setting data in the loop is faster)
    fig = plt.figure(figsize=(10,7), facecolor='w')
    a1 = plt.subplot(111)
    a1.set_xlabel("time (ns)")
    a1.set_ylabel("ADC")
    p1, = a1.plot(np.ones(1), np.ones(1), color='blue', alpha=0.8, label='data')
    p2, = a1.plot(np.ones(1), np.ones(1), color='orange', linewidth=1.5, label='denoised')
    p3, = a1.plot(np.ones(1), np.ones(1), color='red', linewidth=3, label='lowpass')
    a1.legend(loc=4)
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
            energy = waveTree.trapENFCal.at(iH)
            wf = waveTree.MGTWaveforms.at(iH)
            signal = wl.processWaveform(wf,0,0)
            # waveBLSub = signal.GetWaveBLSub()
            waveRaw = signal.GetWaveRaw()
            waveTS = signal.GetTS()
            print "%d / %d  Run %d  nCh %d  chan %d  trapENF %.1f" % (iList,nList,run,nChans,chan,energy)

            # # low pass filter
            # B, A = butter(2,1e6/(1e8/2), btype='lowpass')
            # data_lPass = lfilter(B, A, waveBLSub)
            #
            # # wavelet packet denoised waveform
            # wp = pywt.WaveletPacket(waveBLSub, 'db2', 'symmetric', maxlevel=4)
            # nodes = wp.get_level(4, order='freq')
            # waveletYTrans = np.array([n.data for n in nodes],'d')
            # waveletYTrans = abs(waveletYTrans)
            #
            # # reconstruct waveform w/ only lowest frequency.
            # new_wp = pywt.WaveletPacket(data=None, wavelet='db2', mode='symmetric')
            # new_wp['aaa'] = wp['aaa'].data
            # data_wlDenoised = new_wp.reconstruct(update=False)
            #
            # # resize in a smart way
            # diff = len(data_wlDenoised) - len(waveBLSub)
            # if diff > 0: data_wlDenoised = data_wlDenoised[diff:]

            # print signal.GetOffset()
            # print waveTS[:10]
            # print waveRaw[:10]

            # -- fill the figure --
            # p1.set_ydata(waveBLSub)
            p1.set_ydata(waveRaw)
            p1.set_xdata(waveTS)

            # p2.set_ydata(data_wlDenoised)
            # p2.set_xdata(waveTS)
            #
            # p3.set_ydata(data_lPass)
            # p3.set_xdata(waveTS)

            xmin, xmax = np.amin(waveTS), np.amax(waveTS)
            ymin, ymax = np.amin(waveRaw), np.amax(waveRaw)
            a1.set_xlim([xmin,xmax])
            a1.set_ylim([ymin-abs(0.1*ymin),ymax+abs(0.1*ymax)])
            plt.title("Run %d  Channel %d  Entry %d  trapENFCal %.1f" % (run,chan,iList,energy))

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
