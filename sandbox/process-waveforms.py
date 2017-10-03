#!/usr/bin/env python
""" process-waveforms.py
    Takes an MJD skim file (augmented with an MGTWaveform branch)
    and calculates wavelet parameters for events passing cuts.
    Usage:
    ./process-waveforms.py [-i interactive mode]
                           [-r [dsNum] [subNum] set range & sub-range]
                           [-b batch mode, create new file]
                           [-f [dsNum] [runNum] file mode]
    C. Wiseman, B. O'Zhu
    v1. 3/20/2017
"""
import sys, time, os
from ROOT import TFile,TTree,TEntryList,gDirectory,gROOT,std,TNamed,TObject, MGWFTimePointCalculator
from scipy.optimize import curve_fit
import numpy as np
import waveLibs as wl
import pywt

def main(argv):
    # gROOT.ProcessLine("gErrorIgnoreLevel = 3001;") # suppress ROOT error messages
    startT = time.clock()
    print "Started:",time.strftime('%X %x %Z')
    intMode, batMode, rangeMode, fileMode, dsNum, subNum, runNum = False, False, False, False, -1, -1, -1
    for i,opt in enumerate(argv):
        if opt == "-i":
            intMode = True
            print "Interactive mode selected. Use \"p\" for previous and \"q\" to exit."
        if opt == "-r":
            rangeMode = True
            dsNum, subNum = int(argv[i+1]), int(argv[i+2])
            print "Scanning DS-%d range %d" % (dsNum, subNum)
        if opt == "-f":
            fileMode = True
            dsNum, runNum = int(argv[i+1]), int(argv[i+2])
        if opt == "-b":
            batMode = True
            import matplotlib
            if os.environ.get('DISPLAY','') == '':
                print('No display found. Using non-interactive Agg backend')
                matplotlib.use('Agg')
            print "Batch mode selected.  A new file will be created."
    import matplotlib.pyplot as plt


    # Set file I/O and cuts
    inFile = "./data/waveSkimDS4_test.root"
    outFile = "./data/waveletSkimDS4_test.root"
    if (rangeMode):
        inFile = "~/project/v2-waveskim/waveSkimDS%d_%d.root" % (dsNum,subNum)
        outFile = "~/project/v2-processwfs/waveletSkimDS%d_%d.root" % (dsNum,subNum)
    if (fileMode):
        # inFile = "~/project/cal-wave-skim/waveSkimDS%d_run%d.root" % (dsNum,runNum)
        # outFile = "~/project/cal-wavelet-skim/waveletSkimDS%d_run%d.root" % (dsNum,runNum)
        inFile = "./waveSkimDS%d_run%d.root" % (dsNum, runNum)
        outFile = "./waveletSkimDS%d_run%d.root" % (dsNum, runNum)

    inputFile = TFile(inFile)
    waveTree = inputFile.Get("skimTree")
    print "Found",waveTree.GetEntries(),"input entries."

    theCut = inputFile.Get("theCut").GetTitle()
    # theCut = "trapENFCal > 0.8 && gain==0 && mHClean==1 && isGood && !muVeto && !isLNFill1 && !isLNFill2 && P!=0 && D!=0 && trapETailMin<0"
    # theCut += " && Entry$ < 2000"
    # theCut += " && trapENFCal > 20"
    print "Using cut:\n",theCut,"\n"

    print "Attempting entrylist draw ..."
    waveTree.Draw(">>elist", theCut, "GOFF entrylist")
    elist = gDirectory.Get("elist")
    waveTree.SetEntryList(elist)
    nList = elist.GetN()
    print "Found",nList,"entries passing cuts."


    # In batch mode ONLY, create an output file+tree & append new branches
    outputFile = TFile()
    outTree = TTree()
    if (batMode==True):
        outputFile = TFile(outFile, "RECREATE")
        print "Attempting tree copy to",outFile
        outTree = waveTree.CopyTree("")
        outTree.Write()
        print "Wrote",outTree.GetEntries(),"entries."
        cutUsed = TNamed("cutUsedHere",theCut)
        cutUsed.Write()

    waveS0 = std.vector("double")()
    waveS1 = std.vector("double")()
    waveS2 = std.vector("double")()
    waveS3 = std.vector("double")()
    waveS4 = std.vector("double")()
    waveS5 = std.vector("double")()
    wpar4 = std.vector("double")()
    waveEnt = std.vector("double")()
    butterMax = std.vector("double")()
    butterTime = std.vector("double")()
    tOffset = std.vector("double")()
    lastZeroTime = std.vector("double")()
    pol0 = std.vector("double")()
    pol1 = std.vector("double")()
    pol2 = std.vector("double")()
    pol3 = std.vector("double")()
    exp0 = std.vector("double")()
    exp1 = std.vector("double")()
    rt10 = std.vector("double")()
    rt20 = std.vector("double")()
    rt50 = std.vector("double")()
    rt90 = std.vector("double")()

    bWaveS0 = outTree.Branch("waveS0", waveS0)
    bWaveS1 = outTree.Branch("waveS1", waveS1)
    bWaveS2 = outTree.Branch("waveS2", waveS2)
    bWaveS3 = outTree.Branch("waveS3", waveS3)
    bWaveS4 = outTree.Branch("waveS4", waveS4)
    bWaveS5 = outTree.Branch("waveS5", waveS5)
    bWPar4 = outTree.Branch("wpar4", wpar4)
    bWaveEnt = outTree.Branch("waveEnt",waveEnt)
    bButterMax = outTree.Branch("butterMax",butterMax)
    bButterTime = outTree.Branch("butterTime",butterTime)
    bTOffset = outTree.Branch("tOffset",tOffset)
    bLastZeroTime = outTree.Branch("lastZeroTime",lastZeroTime)
    bPol0 = outTree.Branch("pol0",pol0)
    bPol1 = outTree.Branch("pol1",pol1)
    bPol2 = outTree.Branch("pol2",pol2)
    bPol3 = outTree.Branch("pol3",pol3)
    bExp0 = outTree.Branch("exp0",exp0)
    bExp1 = outTree.Branch("exp1",exp1)
    bRt10 = outTree.Branch("rt10",rt10)
    bRt20 = outTree.Branch("rt20",rt20)
    bRt50 = outTree.Branch("rt50",rt50)
    bRt90 = outTree.Branch("rt90",rt90)


    # Make a figure
    fig = plt.figure(figsize=(7,7), facecolor='w')
    p1 = plt.subplot(211)
    p2 = plt.subplot(212)
    if not batMode: plt.show(block=False)

    # Begin loop over events
    iList = -1
    while True:
        iList += 1
        if intMode==True and iList != 0:
            value = raw_input()
            if value=='q': break        # quit
            if value=='p': iList -= 2   # go to previous
            if (value.isdigit()):
                iList = int(value)      # go to entry number
        elif intMode==False and batMode==False:
            plt.pause(0.001)            # rapid-draw mode
        if iList >= nList: break        # bail out, goose!

        entry = waveTree.GetEntryNumber(iList);
        waveTree.LoadTree(entry)
        waveTree.GetEntry(entry)
        nChans = waveTree.channel.size()
        waveS1.assign(nChans,-999)
        waveS0.assign(nChans,-999)
        waveS1.assign(nChans,-999)
        waveS2.assign(nChans,-999)
        waveS3.assign(nChans,-999)
        waveS4.assign(nChans,-999)
        waveS5.assign(nChans,-999)
        wpar4.assign(nChans,-999)
        waveEnt.assign(nChans,-999)
        butterMax.assign(nChans,-999)
        butterTime.assign(nChans,-999)
        tOffset.assign(nChans,-999)
        lastZeroTime.assign(nChans,-999)
        pol0.assign(nChans,-999)
        pol1.assign(nChans,-999)
        pol2.assign(nChans,-999)
        pol3.assign(nChans,-999)
        exp0.assign(nChans,-999)
        exp1.assign(nChans,-999)
        rt10.assign(nChans,-999)
        rt20.assign(nChans,-999)
        rt50.assign(nChans,-999)
        rt90.assign(nChans,-999)

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
            startTime = waveTree.triggerTrapt0.at(iH)
            blrwfFMR50 = waveTree.blrwfFMR50.at(iH)

            if wf.GetID() != chan:
                print "Warning!!  Vector matching failed!  iList %d  run %d  iEvent %d" % (iList,run,iEvent)
                break

            signal = wl.processWaveform(wf,opt='full')
            waveBLSub = signal.GetWaveBLSub()
            waveFilt = signal.GetWaveFilt()
            waveTS = signal.GetTS()
            waveletYOrig, waveletYTrans = signal.GetWavelet()

            # waveFTX, waveFTY, waveFTPow = signal.GetFourier()
            # waveTrap = signal.GetTrapezoid()
            # waveTrapF = signal.GetFiltTrapezoid()
            # waveTrapX = np.linspace(0, len(waveTrap), num=len(waveTrap))
            # _,waveletFilt = wl.waveletTransform(waveFilt,level=3)  # testing other levels on the filtered WF

            # Wavelet cut parameters
            waveS0[iH] = np.sum(waveletYTrans[0:1,1:-1])
            waveS1[iH] = np.sum(waveletYTrans[0:1,1:33])
            waveS2[iH] = np.sum(waveletYTrans[0:1,33:65])
            waveS3[iH] = np.sum(waveletYTrans[0:1,65:97])
            waveS4[iH] = np.sum(waveletYTrans[0:1,97:-1])
            waveS5[iH] = np.sum(waveletYTrans[2:-1,1:-1])
            wpar4[iH] = np.amax(waveletYTrans[0:1,1:-1])

            # Waveform entropy parameters
            d1 = 2. * np.multiply(waveletYTrans[0:1,1:65], np.log(waveletYTrans[0:1,1:65]/waveS0[iH]/2.0))
            d2 = np.multiply(waveletYTrans[0:1,65:-1], np.log(waveletYTrans[0:1,65:-1]/waveS0[iH]))
            waveEnt[iH] = np.abs(np.sum(d1)) + np.abs(np.sum(d2))

            # Smoothed derivative params
            waveFiltDeriv = wl.wfDerivative(waveFilt)
            butterMax[iH] = np.amax(waveFiltDeriv)
            butterTime[iH] = np.argmax(waveFiltDeriv[100:])*10 + signal.GetOffset() + 1000
            tOffset[iH] = signal.GetOffset()
            # print "max %.3f  ts %.0f ns  offset %.0f ns" % (butterMax[iH], butterTime[iH], tOffset[iH])

            # Make a super denoised waveform
            wp = pywt.WaveletPacket(data=waveBLSub, wavelet='haar', mode='symmetric',maxlevel=3)
            new_wp = pywt.WaveletPacket(data=None, wavelet='haar', mode='symmetric')
            new_wp['aaa'] = wp['aaa'].data
            superDenoisedWF = new_wp.reconstruct(update=False)
            mgtSDWF = wl.MGTWFFromNpArray(superDenoisedWF)

            # Try to get start time by finding the super denoised last zero crossing
            lastZero = 0
            zeros = np.asarray(np.where(superDenoisedWF < 0.1))
            if (zeros.size > 0): lastZero = zeros[0,-1]
            lastZeroTime[iH] = waveTS[lastZero]

            # Time point calculator
            timePointCalc = MGWFTimePointCalculator();
            timePointCalc.AddPoint(.1)
            timePointCalc.AddPoint(.2)
            timePointCalc.AddPoint(.5)
            timePointCalc.AddPoint(.9)
            timePointCalc.FindTimePoints(mgtSDWF)
            rt10[iH] = timePointCalc.GetFromStartRiseTime(0)*10
            rt20[iH] = timePointCalc.GetFromStartRiseTime(1)*10
            rt50[iH] = timePointCalc.GetFromStartRiseTime(2)*10
            rt90[iH] = timePointCalc.GetFromStartRiseTime(3)*10
            # print "lastZero %.2f  startTime %.2f  10 %.2f  20 %.2f  50 %.2f  90 %.2f" % (lastZeroTime[iH], startTime,rt10, rt20, rt50, rt90)

            # Fit tail slope (2 methods).  Guard against fit fails
            idxMax = np.where(waveTS > rt90[iH]) # returns an array/tuple
            idxMax = idxMax[0][0] # "cast" to int
            tail, tailTS = waveBLSub[idxMax:], waveTS[idxMax:]
            try:
                popt1,_ = curve_fit(wl.tailModelPol, tailTS, tail)
                pol0[iH], pol1[iH], pol2[iH], pol3[iH] = popt1[0], popt1[1], popt1[2], popt1[3]
            except: pass
            try:
                popt2,_ = curve_fit(wl.tailModelExp, tailTS, tail, p0=[energy,72000])
                exp0[iH], exp1[iH] = popt2[0], popt2[1]
            except:
                # print "curve_fit tailModelExp failed, run %i  event %i  channel %i" % (run, iList, chan)
                pass

            # Make a plot
            if not batMode:
                print "i %d  run %d  iEvent(i) %d  iH(j+1) %d/%d  chan %d  energy %.2f  bt %.2f" % (iList,run,iEvent,iH+1,nChans,chan,energy,butterTime[iH])

                # Waveform plot
                p1.cla()
                p1.plot(waveTS,waveBLSub,color='blue')
                p1.plot(waveTS,superDenoisedWF,color='red')
                p1.plot(tailTS, wl.tailModelPol(tailTS, *popt1), color='cyan')
                p1.plot(tailTS, wl.tailModelExp(tailTS, *popt2), color='orange')
                p1.set_xlabel("Time (ns)")
                p1.set_ylabel("ADC")
                p1.set_title("Run %d  Ch %d  Energy %.2f  \n S5 %.2f  (S3-S2)/E %.2f  (S3-S4)/E %.2f" % (run,chan,energy,waveS5[iH],((waveS3[iH]-waveS2[iH])/energy),((waveS3[iH]-waveS4[iH])/energy)))

                # Wavelet plot
                p2.cla()
                p2.imshow(waveletYTrans, interpolation='nearest', aspect="auto", origin="lower",extent=[0, 1, 0, len(waveletYTrans)],cmap='jet')

                plt.tight_layout()
                plt.subplots_adjust(hspace=.2,top=0.92)
                plt.pause(0.00001)


        # End loop over hits
        if (batMode == True):
            bWaveS0.Fill()
            bWaveS1.Fill()
            bWaveS2.Fill()
            bWaveS3.Fill()
            bWaveS4.Fill()
            bWaveS5.Fill()
            bWPar4.Fill()
            bWaveEnt.Fill()
            bButterMax.Fill()
            bButterTime.Fill()
            bTOffset.Fill()
            bLastZeroTime.Fill()
            bPol0.Fill()
            bPol1.Fill()
            bPol2.Fill()
            bPol3.Fill()
            bExp0.Fill()
            bExp1.Fill()
            bRt10.Fill()
            bRt20.Fill()
            bRt50.Fill()
            bRt90.Fill()
            if iList%5000 == 0:
                outTree.Write("",TObject.kOverwrite)
                print "%d / %d entries saved (%.2f %% done)." % (iList,nList,100*(float(iList)/nList))

    # End loop over events
    if (batMode == True):
        outTree.Write("",TObject.kOverwrite)
        print "Wrote",outTree.GetBranch("channel").GetEntries(),"entries in the copied tree,"
        print "and wrote",bWaveS0.GetEntries(),"entries in the new branches."

    stopT = time.clock()
    print "Process time (min): ", (stopT-startT)/60


if __name__ == "__main__":
    main(sys.argv[1:])
