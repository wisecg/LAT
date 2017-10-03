#!/usr/bin/env python
import sys
import numpy as np
import matplotlib.pyplot as plt
from ROOT import TFile,TTree,TCanvas,TEntryList,gDirectory,MGTWaveform
import waveLibs as wl

def main(argv):

def main(argv):
    """ Bug: Doesn't always work with a TChain.
    Add input files together with hadd and use a single TFile.
    """
    scanSpeed = 0.2
    opt1, opt2 = "", ""
    intMode, printWF = False, False
    if (len(argv) >= 1): opt1 = argv[0]
    if (len(argv) >= 2): opt2 = argv[1]
    if "-i" in (opt1, opt2):
        intMode = True
        print "Interactive mode selected."
    if "-s" in (opt1, opt2):
        printWF = True
        print "Saving WF plots to current directory."

    # Set input file
    # DS5 calibration run (lots of all pulses)
    # inputFile = TFile("~/project/match-skim/waveletSkimDS5_run23884.root")

    # Simulated noise pulses from Brian
    inputFile = TFile("./data/SimWF_Baseline.root")
    waveTree = inputFile.Get("skimTree")

    # ext. pulsing runs (P2D2, chan 674) https://majorana.npl.washington.edu/elog/Run+Elog/1093
    # gatFile = TFile("~/project/gat-blt-data/gatified/mjd_run7219.root")
    # bltFile = TFile("~/project/gat-blt-data/built/OR_run7219.root")

    # ext. trigger runs https://majorana.npl.washington.edu/elog/Run+Elog/1207
    # gatFile = TFile("~/project/gat-blt-data/gatified/mjd_run9913.root")
    # bltFile = TFile("~/project/gat-blt-data/built/OR_run9913.root")

    # waveTree = gatFile.Get("mjdTree")
    # bltTree = bltFile.Get("MGTree")
    # waveTree.AddFriend(bltTree)

    print "Found",waveTree.GetEntries(),"input entries."

    # Set cuts
    # theCut = inputFile.Get("cutUsedHere").GetTitle()
    # theCut = "trapENFCal > 0.8 && gain==0 && mH==1 && isGood && !muVeto && !wfDCBits && !isLNFill1 && !isLNFill2 && trapETailMin < 0 && channel!=596 && channel!=676 && channel!=676 && channel!=612 && channel!=1104 && channel!=1200 && channel!=1334 && channel!=1336"
    # theCut = "channel==674"
    theCut = ""
    print "Using cut:\n",theCut,"\n"
    waveTree.Draw(">>elist", theCut, "entrylist")
    elist = gDirectory.Get("elist")
    waveTree.SetEntryList(elist)
    nList = elist.GetN()
    print "Found",nList,"entries passing cuts."

    # Make a figure
    fig = plt.figure(figsize=(8,5), facecolor='w')
    # p1 = plt.subplot2grid((6,7), (0,0), colspan=4, rowspan=2) # original
    p1 = plt.subplot(111)
    plt.show(block=False)

    # Make a template waveform
    samp, r, z, ene, t0, smooth = 2016, 0, 15, 10, 1000, 100
    temp, tempTS = MakeSiggenWaveform(samp,r,z,ene,t0,smooth)

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
        # nChans = waveTree.channel.size()
        nChans = waveTree.MGTWaveforms.size()
        # numPass = waveTree.Draw("channel",theCut,"GOFF",1,iList)
        # chans = waveTree.GetV1()
        # chanList = list(set(int(chans[n]) for n in xrange(numPass)))
        # event = bltTree.event  # grab waveforms from built data (comment this out if we're using clint-skim)

        # Loop over hits passing cuts
        # hitList = (iH for iH in xrange(nChans) if waveTree.channel.at(iH) in chanList)  # a 'generator expression'
        # for iH in hitList:
        for iH in range(nChans):
            # run = waveTree.run
            # chan = waveTree.channel.at(iH)
            # dataE = waveTree.trapENFCal.at(iH)
            signal = wl.processWaveform(waveTree.MGTWaveforms.at(iH)) # use when we have skim
            # signal = wl.processWaveform(event.GetWaveform(iH)) # use when we have gat+blt data
            waveBLSub = signal.GetWaveBLSub()
            waveFilt = signal.GetWaveFilt()
            waveTS = signal.GetTS() * 10.
            # print "%d / %d  Run %d  nCh %d  chan %d  trapENF %.1f" % (iList,nList,run,nChans,chan,dataE)

            # Generate a siggen WF
            samp, r, z, ene, t0, smooth = 2016, 0, 15, 10, 1000, 100
            sigg, siggTS = MakeSiggenWaveform(samp,r,z,ene,t0,smooth)

            data = sigg + waveBLSub

            # Fill the figure
            p1.cla()
            p1.set_ylabel("ADC")
            # p1.set_title("Run %d  Channel %d  Entry %d  trapENFCal %.1f" % (run,chan,iList,dataE))
            p1.plot(waveTS,data,color='blue')
            p1.plot(tempTS,temp,color='red')

            plt.pause(scanSpeed)
            # if (printWF): plt.savefig("./plots/wave-%d-%d-%d.pdf" % (run,iList,chan))



def MakeSiggenWaveform(samp,r,z,ene,t0,smooth=1,phi=np.pi/8):
    """ Use pysiggen to generate a waveform w/ arb. ADC amplitude. Thanks Ben. """

    wf_length = samp # in tens of ns.  This sets how long a wf you will simulate
    # r ranges from 0 to detector.detector_radius
    # phi ranges from 0 to np.pi/4                (Clint: Limit is the pcrad, set below)
    # z ranges from 0 to detector.detector_length (Clint: Limit is the pclen, set below)
    # energy sets the waveform amplitude
    # t0 sets the start point of the waveform
    # smooth is a gaussian smoothing parameter
    # default: r, phi, z, energy, t0, smooth = (15, np.pi/8, 15, 80, 10,15)

    # Probably don't mess around with anything in this block
    fitSamples = 1000 # ask me if you think you need to change this (you almost definitely don't)
    timeStepSize = 1 # don't change this, you'll break everything
    # Create a detector model (don't change any of this)
    detName = "./data/conf/P42574A_grad%0.2f_pcrad%0.2f_pclen%0.2f.conf" % (0.05,2.5, 1.65)
    detector =  Detector(detName, timeStep=timeStepSize, numSteps=fitSamples*10./timeStepSize, maxWfOutputLength=5000)
    detector.LoadFieldsGrad("./data/fields_impgrad.npz",pcLen=1.6, pcRad=2.5)
    # Sets the impurity gradient.  Don't bother changing this
    detector.SetFieldsGradIdx(0)

    # First 3 params control the preamp shaping.
    # The 4th param is the RC decay, you can change that if you want.
    # params 5 and 6 are set not to do anything (theyre for the 2 rc constant decay model)
    rc_decay = 72.6 #us
    detector.SetTransferFunction(50, -0.814072377576, 0.82162729751, rc_decay, 1, 1)

    wf_notrap = np.copy(detector.MakeSimWaveform(r, phi, z, ene, t0, wf_length, h_smoothing=smooth))
    timesteps = np.arange(0, wf_length) * 10 # makes it so your plot is in ns

    # Add charge trapping
    # trap_constants = [8, 28,288] # these are in microseconds
    # trap_constant_colors = ["red", "blue", "purple"]
    # for (idx, trap_rc) in enumerate(trap_constants):
    #     detector.trapping_rc = trap_rc
    #     wf = np.copy(detector.MakeSimWaveform(r, phi, z, ene, t0, wf_length, h_smoothing=smooth))
    #     ax0.plot(timesteps, wf, color = trap_constant_colors[idx],  label = "%0.1f us trapping" % trap_rc )
    #     ax1.plot(timesteps, wf_notrap - wf, color = trap_constant_colors[idx])
    #     print "amplitude diff: %f" % ( (np.amax(wf_notrap) - np.amax(wf)) /  np.amax(wf_notrap) )

    return wf_notrap, timesteps



def makeNoisePlots():

    # External pulsing runs
    # key: https://majorana.npl.washington.edu/elog/Run+Elog/1093
    # P2D2 is the magic detector.
    # 9  1  10  3  2  11  0  9  3  1500  "P42664A"  2  6  674  675  "C1P2"  2  679
    f1 = TFile("~/project/noise/mjd_run7219.root")
    puls = f1.Get("mjdTree")
    print "puls:",puls.GetEntries()

    # Forced trigger runs
    f2 = TFile("~/project/noise/mjd_run9913.root")
    trig = f2.Get("mjdTree")
    print "trig:",trig.GetEntries()

    # Make some quick plots
    c = TCanvas("c","Bob Ross's Canvas",800,600)

    puls.Draw("trapENFCal")
    c.Print("./plots/puls_energy.pdf")

    trig.Draw("trapENFCal")
    c.Print("./plots/trig_energy.pdf")

    trig.Draw("trapENFCal:timestamp")
    c.Print("./plots/trig_vsTime.pdf")


if __name__ == "__main__":
    main(sys.argv[1:])
