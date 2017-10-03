#!/usr/bin/env python
import sys, time, pywt
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
from matplotlib import gridspec
from ROOT import TFile,TTree,TChain,TEntryList,gDirectory,gROOT,MGTWaveform
import waveLibs as wl
import waveModel as wm
import pymc
from scipy.ndimage.filters import gaussian_filter
from scipy import interpolate
from scipy.optimize import curve_fit

def main(argv):
    """ Interactive-fit or 'rapid'-fit waveforms that pass a given TCut.
        BUG: Doesn't always work with a TChain.  Add input files together
        with hadd and use a single TFile.
    """
    scanSpeed = 0.2
    iList = -1
    opt1, opt2 = "", "-"
    intMode = False
    if (len(argv) >= 1): opt1 = argv[0]
    if (len(argv) >= 2): opt2 = argv[1]
    if "-i" in (opt1, opt2):
        intMode = True
        print "Interactive mode selected."

    startT = time.clock()

    # Load template waveforms
    templateFile = TFile("./data/wave-m1cal.root")
    waveTemplates = templateFile.Get("waveTree")
    calDict = {640:28, 672:18, 610:16, 580:19, 582:34, 648:38, 600:21, 578:39, 592:27, 664:55, 626:62, 692:8, 598:22, 690:52, 632:9, 608:7}
    # PlotTemplateWaveforms(waveTemplates,calDict)
    # return

    # Set input data and cuts
    inputFile = TFile("~/project/wavelet-skim/hadd/waveletSkimDS3.root")
    # inputFile = TFile("~/project/wavelet-skim/waveletSkimDS3_13.root")
    waveTree = inputFile.Get("skimTree")
    print "Found",waveTree.GetEntries(),"input entries."

    # theCut = inputFile.Get("cutUsedHere").GetTitle()

    # DS3 big cut
    theCut = "trapENFCal > 0.8 && gain==0 && mHClean==1 && isGood && !muVeto && !wfDCBits && !isLNFill1 && !isLNFill2 && trapETailMin < 0 && channel!=596 && channel!=676 && channel!=676 && channel!=612 && channel!=1104 && channel!=1200 && channel!=1334 && channel!=1336 && tOffset < 10 && waveS5/TMath::Power(trapENFCal,1/4) < 1200 && (waveS3-waveS2)/trapENFCal > 100 && (waveS3-waveS2)/trapENFCal < 300 && !(channel==692 && (run==16974 || run==16975 || run==16976 || run==16977 || run==16978 || run==16979)) && butterTime < 11000"

    # theCut += " && trapENFCal > 2 && trapENFCal < 2.1 && channel==578 && run==16868"  # template fast WF (2kev)
    theCut += " && trapENFCal > 7 && trapENFCal < 10"
    # theCut += " && trapENFCal < 20 && trapENFCal > 2 && run > 17480"
    # theCut += " && kvorrT/trapENFCal > 2.2 && trapENFCal > 2 && trapENFCal < 10"

    # Print cut and events passing cut
    print "Using cut:\n",theCut,"\n"
    waveTree.Draw(">>elist", theCut, "entrylist")
    elist = gDirectory.Get("elist")
    waveTree.SetEntryList(elist)
    nList = elist.GetN()
    print "Found",nList,"entries passing cuts."

    stopT=time.clock()
    print "Data loading time (s): ", (stopT-startT)

    # Make a figure
    fig = plt.figure(figsize=(11,7), facecolor='w')
    p1 = plt.subplot2grid((6,7), (0,0), colspan=4, rowspan=2) # original
    p2 = plt.subplot2grid((6,7), (2,0), colspan=4, rowspan=3) # rising edge
    p3 = plt.subplot2grid((6,7), (5,0), colspan=4, rowspan=1) # residual
    p4 = plt.subplot2grid((6,7), (0,4), colspan=3, rowspan=2 )           # trace 1
    p5 = plt.subplot2grid((6,7), (2,4), colspan=3, rowspan=2, sharex=p4) # trace 2
    p6 = plt.subplot2grid((6,7), (4,4), colspan=3, rowspan=2, sharex=p4) # trace 3
    plt.show(block=False)

    # Loop over events
    while(True):
        saveMe = False
        iList += 1
        if intMode==True and iList != 0:
            value = raw_input()
            if value=='q': break
            if value=='s': saveMe=True
            if value=='p': iList -= 2  # previous
            if (value.isdigit()): iList = int(value) # go to entry
        if iList >= elist.GetN(): break

        entry = waveTree.GetEntryNumber(iList);
        waveTree.LoadTree(entry)
        waveTree.GetEntry(entry)
        nChans = waveTree.channel.size()
        numPass = waveTree.Draw("channel",theCut,"GOFF",1,iList)
        chans = waveTree.GetV1()
        chanList = list(set(int(chans[n]) for n in xrange(numPass)))

        # Loop over hits passing cuts
        hitList = (iH for iH in xrange(nChans) if waveTree.channel.at(iH) in chanList)  # a 'generator expression'
        for iH in hitList:

            # Load waveform for this hit
            run = waveTree.run
            chan = waveTree.channel.at(iH)
            dataE = waveTree.trapENFCal.at(iH)
            dataST = waveTree.butterTime.at(iH) # replace with blrwfFMR50?
            toe = waveTree.kvorrT.at(iH)/dataE
            print "%d / %d  Run %d  nCh %d  chan %d  trapENF %.1f  t/e %.1f" % (iList,nList,run,nChans,chan,dataE,toe)
            signal = wl.processWaveform(waveTree.MGTWaveforms.at(iH),opt='full')
            waveBLSub = signal.GetWaveBLSub()
            waveFilt = signal.GetWaveFilt()
            waveTS = signal.GetTS()
            baseAvg, dataNoise = signal.GetBaseNoise()

            # Denoise the data waveform (take only lowest-frequency components)
            wp = pywt.WaveletPacket(data=waveBLSub, wavelet='haar', mode='symmetric',maxlevel=3)
            new_wp = pywt.WaveletPacket(data=None, wavelet='haar', mode='symmetric')
            new_wp['aaa'] = wp['aaa'].data
            waveDenoised = new_wp.reconstruct(update=False)

            # Load channel template waveform
            try: tempEntry = calDict[chan]
            except:
                print "ain't got channel",chan
                continue
            waveTemplates.GetEntry(calDict[chan])
            template = wl.processWaveform(waveTemplates.event.GetWaveform(waveTemplates.itr),opt='full')
            temp = template.GetWaveBLSub()
            tempTS = template.GetTS()
            tempE = waveTemplates.trapENFCal
            tempST = waveTemplates.blrwfFMR50

            # Window the fit around rising edge - start time calculator method
            loWin, hiWin = dataST - 1000, dataST + 4000 # ns
            if loWin < waveTS[0] or hiWin > waveTS[-1]:
                print "Window out of range!  dataST: %.1f  loWin %.1f  hiWin %.1f" % (dataST,loWin,hiWin)
            idx = np.where((waveTS >= loWin) & (waveTS <= hiWin))
            data = waveBLSub[idx]
            # data = waveDenoised[idx]
            dataTS = waveTS[idx]

            # Pack into lists
            rawList = [waveBLSub, waveTS, dataE, dataST]
            dataList = [data, dataTS, dataE, dataST, loWin, hiWin]
            tempList = [temp, tempTS, tempE, tempST]

            # Optionally save to a file
            if saveMe:
                print "Saved entry",iList,iH
                np.savez("./data/optimumInputs.npz",rawList,tempList)
                # np.savez("./data/rawInputs.npz",rawList)

            # Recreate the guess and the guess's rising edge
            guessFull, guessFullTS = wm.MakeModel(dataList, tempList, [dataST,dataE,1.], opt="full")
            guess, guessTS = wm.MakeModel(dataList, tempList, [dataST,dataE,1.], opt="!fancy")

            # Make an "almost complete" guess - no MCMC
            # st, en, slo = dataST-100, dataE, 5
            # InterpFn = interpolate.interp1d(tempTS, temp, kind="linear", copy="False", assume_sorted="True")
            # model, modelTS = wm.MakeModel(dataList, tempList, [st,en,slo], fn=InterpFn)

            # Fit with MCMC and get best-fit parameters
            numSteps, burnIn = 2000, 1000 # default - 10000, 5000.  try 3000, 2000
            myModel = wm.TemplateModel( dataList, dataNoise, tempList )
            waveModel = pymc.Model( myModel )
            M = pymc.MCMC( waveModel )
            M.use_step_method(pymc.Metropolis, M.startTime, proposal_sd=100., proposal_distribution='Normal')
            M.use_step_method(pymc.Metropolis, M.energy, proposal_sd=1., proposal_distribution='Normal')
            M.use_step_method(pymc.Metropolis, M.slowness, proposal_sd=100., proposal_distribution='Normal')
            M.sample(iter=numSteps, verbose=0) # do the fit
            st = np.median(M.trace('startTime')[burnIn:])
            en = np.median(M.trace('energy')[burnIn:])
            slo = np.median(M.trace('slowness')[burnIn:])
            InterpFn = interpolate.interp1d(tempTS, temp, kind="linear", copy="False", assume_sorted="True")
            model, modelTS = wm.MakeModel(dataList, tempList, [st,en,slo], fn=InterpFn)

            # Calculate residual, Chi2/NDF
            residual = model - data
            frac = (np.power(data - model, 2)) / np.abs(model)
            chi2NDF = np.sum(frac) / len(model)
            print "chi2/NDF:",chi2NDF

            # ** Do a separate & simple fit of the tail slope **
            # TODO: Add this to process-waveforms.py
            idxMax = np.where(guessFull == guessFull.max()) # returns an array/tuple
            idxMax = idxMax[0][0] # "cast" to int
            tail, tailTS = waveDenoised[idxMax:], waveTS[idxMax:]
            popt,_ = curve_fit(wl.tailModelPol, tailTS, tail) # poly fit
            a, b = dataE, 72000
            popt2,_ = curve_fit(wl.tailModelExp, tailTS, tail, p0=[a,b]) # expo fit

            # Fill the figure
            p1.cla()
            p1.set_ylabel("ADC")
            p1.set_title("Run %d  Channel %d  Entry %d\ntrapENFCal %.1f  T/E %.1f  ST %.1f" % (run,chan,iList,dataE,toe,dataST))
            p1.plot(waveTS,waveBLSub,color='blue')
            p1.plot(waveTS,waveDenoised,color='red',alpha=0.8)
            p1.plot(guessFullTS,guessFull,color='orange',linewidth=2)
            p1.axvline(x = dataST, color='green',linewidth=2)
            p1.axvline(x = loWin, color='black')
            p1.axvline(x = hiWin, color='black')
            p1.plot(tailTS, wl.tailModelPol(tailTS, *popt), color='cyan', linewidth=1) # tail poly fit
            p1.plot(tailTS, wl.tailModelExp(tailTS, *popt2), color='cyan', linewidth=1) # tail expo fit

            p2.cla()
            p2.plot(dataTS,data,color='blue',label='Data')
            p2.plot(guessTS,guess,color='orange',label='Guess')
            p2.plot(modelTS,model,color='red',linewidth=3,alpha=0.7,label='Fit')
            p2.legend(loc=4)

            p3.cla()
            p3.set_xlabel("Time [ns]", x=0.95, ha='right')
            p3.set_ylabel("Residual [ADC]")
            p3.plot(modelTS,residual, color='red')
            p3.axhline(y = 0, color='blue', alpha=0.3)
            p3.axhline(y = dataNoise, color='blue', alpha=0.3)
            p3.axhline(y = -1.0*dataNoise, color='blue', alpha=0.3)

            p4.cla()
            minST = tempST - tempTS[-1] + hiWin
            maxST = tempST - tempTS[0] + loWin
            p4.set_title("startTime %.1f  Energy %.2f\nSlow %.1f  chi2/ndf %.1f  Min %d  Max %d" % (st,en,slo,chi2NDF,minST,maxST))
            p4.plot(M.trace('startTime')[:])
            p4.set_ylabel('startTime')

            p5.cla()
            p5.plot(M.trace('energy')[:])
            p5.set_ylabel('energy')

            p6.cla()
            p6.plot(M.trace('slowness')[:])
            p6.set_ylabel('slowness')

            plt.tight_layout()
            plt.subplots_adjust(hspace=0.35)
            plt.pause(scanSpeed)

# ===========================================================================

def PlotTemplateWaveforms(waveTree,calDict):
    """ Check the wf's grabbed from an M1 calibration run. """
    # plt.style.use('mjWaveforms')  # see ~/.matplotlib/stylelib
    # optional: redefine calDict with only one or two channels for debugging
    # calDict = {672:18}
    calDict = {640: 28, 672: 18, 610: 16, 580: 19, 648: 38, 632: 9, 578: 39, 592: 27, 664: 55, 626: 62, 692: 8, 598: 22, 690: 52, 600: 21, 608: 7}

    # Make a figure
    fig = plt.figure(figsize=(8, 5), facecolor='w')
    ax = plt.subplot(111)
    ax.set_xlabel("Time (ns)",x=0.95)
    ax.set_ylabel("ADC",y=0.95)
    jet = plt.get_cmap('jet')
    cNorm  = colors.Normalize(vmin=578, vmax=672)
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)

    # Plot multiple waveforms
    for chan in calDict:
        waveTree.GetEntry(calDict[chan])
        waveMGT = waveTree.event.GetWaveform(waveTree.itr)
        signal = wl.processWaveform(waveMGT)
        waveRaw = signal.GetWaveRaw()
        waveOff = signal.GetOffset()/1e9 # seconds
        waveTS = signal.GetTS()
        colorVal = scalarMap.to_rgba(chan)
        label = ('%d' % chan)
        ax.plot(waveTS, waveRaw, color=colorVal, label=label)

        print np.amax(waveRaw)

    # make a legend and show the plot.
    handles,labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels, loc='upper right')
    plt.show()
    # plt.savefig("./plots/template-wfs.pdf")

if __name__ == "__main__":
    main(sys.argv[1:])
