#!/usr/local/bin/python
#!/usr/common/usg/software/python/2.7.9/bin/python
#!/bin/bash
""":" `which python` "$0" "$@"; exit 1 ":""" # dammit, this works on PDSF but not my laptop
"""
============== lat.py: (L)ANL (A)nalysis (T)oolkit ==============

Secondary waveform processing for the Majorana Demonstrator.

Takes a single skim file (augmented with an MGTWaveform branch),
or built/gatified data.  Calculates various waveform parameters
for HITS passing cuts.

Does not handle TChains - pyroot can't always recognize
the vector<MGTWaveform*> branch in the skim tree. Damn you, ROOT.

Usage:
./lat.py [-r [dsNum] [subNum] use skim range & sub-range]
         [-p [inFile] [outFile] path mode: manual file locations]
         [-d [inPath] [outPath] set input and output directories]
         [-f [dsNum] [runNum] use single skim file]
         [-g [dsNum] [runNum] use single gat/blt file]
         [-s -- use a custom file, must combine w/ path mode]
         [-i [plotNum] interactive mode]
         [-b batch mode -- creates new file]

v1: 27 May 2017

================ C. Wiseman (USC), B. Zhu (LANL) ================
"""
import sys, time, os, pywt
from ROOT import TFile, TTree, TEntryList, gDirectory, TNamed, std, TObject
from ROOT import GATDataSet, MGTEvent, MGTWaveform, MGWFTimePointCalculator
import numpy as np
from scipy.signal import butter, lfilter, filtfilt
import scipy.optimize as op
from scipy.ndimage.filters import gaussian_filter
from scipy import interpolate
import waveLibs as wl
import waveModel as wm

def main(argv):
    print "======================================="
    print "LAT started:",time.strftime('%X %x %Z')
    startT = time.clock()
    # gROOT.ProcessLine("gErrorIgnoreLevel = 3001;") # suppress ROOT error messages

    # let's get some (m)args
    intMode, batMode, rangeMode, fileMode, gatMode, singleMode, pathMode = False, False, False, False, False, False, False
    dsNum, subNum, runNum, plotNum = -1, -1, -1, 1
    pathToInput, pathToOutput, manualInput, manualOutput = ".", ".", "", ""

    if len(argv)==0: return
    for i,opt in enumerate(argv):
        if opt == "-r":
            rangeMode, dsNum, subNum = True, int(argv[i+1]), int(argv[i+2])
            print "Scanning DS-%d sub-range %d" % (dsNum, subNum)
        if opt == "-p":
            pathMode, manualInput, manualOutput = True, argv[i+1], argv[i+2]
            print "Manually set input/output files:\nInput: %s\nOutput: %s" % (manualInput, manualOutput)
        if opt == "-d":
            pathToInput, pathToOutput = argv[i+1], argv[i+2]
            print "Custom paths: Input %s,  Output %s" % (pathToInput,pathToOutput)
        if opt == "-f":
            fileMode, dsNum, runNum = True, int(argv[i+1]), int(argv[i+2])
            print "Scanning DS-%d, run %d" % (dsNum, runNum)
        if opt == "-g":
            gatMode, runNum = True, int(argv[i+1])
            print "GATDataSet mode.  Scanning run %d" % (runNum)
        if opt == "-s":
            singleMode, runNum, pathToInput = True, int(argv[i+1]), argv[i+2]
            print "Single file mode.  Scanning run %d" % (runNum)
        if opt == "-i":
            intMode, plotNum = True, int(argv[i+1])
            print "Interactive mode selected. Use \"p\" for previous and \"q\" to exit."
        if opt == "-b":
            batMode = True
            import matplotlib
            if os.environ.get('DISPLAY','') == '':
                print('No display found. Using non-interactive Agg backend')
                matplotlib.use('Agg')
            print "Batch mode selected.  A new file will be created."
    import matplotlib.pyplot as plt
    from matplotlib import gridspec

    # File I/O
    inFile, outFile, bltFile = TFile(), TFile(), TFile()
    gatTree, bltTree, oTree = TTree(), TTree(), TTree()
    theCut, inPath, outPath = "", "", ""

    # Set input and output files
    if rangeMode:
        inPath = "%s/waveSkimDS%d_%d.root" % (pathToInput, dsNum, subNum)
        outPath = "%s/latSkimDS%d_%d.root" % (pathToOutput, dsNum, subNum)
    if fileMode:
        inFull = "%s/waveSkimDS%d_run%d.root" % (pathToInput, dsNum, runNum)
        outPath = "%s/latSkimDS%d_run%d.root" % (pathToOutput, dsNum, runNum)
    if pathMode:
        inPath, outPath = manualInput, manualOutput
    if gatMode:
        ds = GATDataSet()
        gatPath = ds.GetPathToRun(runNum,GATDataSet.kGatified)
        bltPath = ds.GetPathToRun(runNum,GATDataSet.kBuilt)
        outPath = "%s/lat_run%d.root" % (pathToOutput, runNum)
    if pathMode and gatMode:
        outPath = manualOutput

    # Initialize trees
    if rangeMode or fileMode or pathMode:
        inFile = TFile(inPath)
        gatTree = inFile.Get("skimTree")
        theCut = inFile.Get("theCut").GetTitle()
    elif gatMode:
        inFile = TFile(gatPath)
        bltFile = TFile(bltPath)
        gatTree = inFile.Get("mjdTree")
        bltTree = bltFile.Get("MGTree")
        gatTree.AddFriend(bltTree)
    if singleMode:
        gatTree = inFile.Get("skimTree")

    # ============================================================
    # Set cuts (maybe someday they should be read in ...)
    # theCut += " && trapENFCal < 10 && trapENFCal > 2 && kvorrT/trapENFCal < 1"
    # theCut += " && Entry$ < 50"
    # theCut = "trapENFCal < 10 && fitSlo > 30 && trapENFCal > 2"
    theCut = "trapENFCal > 2 && trapENFCal < 10"
    print "WARNING: Custom cut in use!"
    # ============================================================
    gatTree.Draw(">>elist", theCut, "entrylist")
    elist = gDirectory.Get("elist")
    gatTree.SetEntryList(elist)
    nList = elist.GetN()
    print "Using cut:\n",theCut
    print "Found",gatTree.GetEntries(),"input entries."
    print "Found",nList,"entries passing cuts."

    # Output: In batch mode (-b) only, create an output file+tree & append new branches.
    if batMode:
        outFile = TFile(outPath, "RECREATE")
        print "Attempting tree copy to",outPath
        oTree = gatTree.CopyTree("")
        oTree.Write()
        print "Wrote",oTree.GetEntries(),"entries."
        cutUsed = TNamed("theCut",theCut)
        cutUsed.Write()

    waveS0, waveS1, waveS2 = std.vector("double")(), std.vector("double")(), std.vector("double")()
    waveS3, waveS4, waveS5 = std.vector("double")(), std.vector("double")(), std.vector("double")()
    wpar4, waveEnt, tOffset = std.vector("double")(), std.vector("double")(), std.vector("double")()
    bcMax, bcMin = std.vector("double")(), std.vector("double")()
    derivMax, derivTime = std.vector("double")(), std.vector("double")()
    bandMax, bandTime = std.vector("double")(), std.vector("double")()
    raw10, raw50, raw90 = std.vector("double")(), std.vector("double")(), std.vector("double")()
    den10, den50, den90 = std.vector("double")(), std.vector("double")(), std.vector("double")()
    oppie = std.vector("double")()
    fitMatch, fitE, fitSlo = std.vector("double")(), std.vector("double")(), std.vector("double")()
    matchMax, matchWidth, matchTime = std.vector("double")(), std.vector("double")(), std.vector("double")()
    pol0, pol1, pol2, pol3 = std.vector("double")(), std.vector("double")(), std.vector("double")(), std.vector("double")()
    exp0, exp1 = std.vector("double")(), std.vector("double")()
    fitMax = std.vector("double")()
    fails = std.vector("int")()

    # It's not possible to put the "oTree.Branch" call into a class initializer (waveLibs::latBranch). You suck, ROOT.
    b1, b2, b3 = oTree.Branch("waveS0",waveS0), oTree.Branch("waveS1",waveS1), oTree.Branch("waveS2",waveS2)
    b4, b5, b6 = oTree.Branch("waveS3",waveS3), oTree.Branch("waveS4",waveS4), oTree.Branch("waveS5",waveS5)
    b7, b8, b9 = oTree.Branch("wpar4",wpar4), oTree.Branch("waveEnt",waveEnt), oTree.Branch("tOffset",tOffset)
    b10, b11 = oTree.Branch("bcMax",bcMax), oTree.Branch("bcMin",bcMin)
    b12, b13 = oTree.Branch("derivMax",derivMax), oTree.Branch("derivTime",derivTime)
    b14, b15 = oTree.Branch("bandMax",bandMax), oTree.Branch("bandTime",bandTime)
    b16, b17, b18 = oTree.Branch("raw10",raw10), oTree.Branch("raw50",raw50), oTree.Branch("raw90",raw90)
    b19, b20, b21 = oTree.Branch("den10",den10), oTree.Branch("den50",den50), oTree.Branch("den90",den90)
    b22 = oTree.Branch("oppie",oppie)
    b23, b24, b25 = oTree.Branch("fitMatch", fitMatch), oTree.Branch("fitE", fitE), oTree.Branch("fitSlo", fitSlo)
    b26, b27, b28 = oTree.Branch("matchMax", matchMax), oTree.Branch("matchWidth", matchWidth), oTree.Branch("matchTime", matchTime)
    b29, b30, b31, b32 = oTree.Branch("pol0", pol0), oTree.Branch("pol1", pol1), oTree.Branch("pol2", pol2), oTree.Branch("pol3", pol3)
    b33, b34 = oTree.Branch("exp0", exp0), oTree.Branch("exp1", exp1)
    b35 = oTree.Branch("fitMax",fitMax)
    b36 = oTree.Branch("fails",fails)

    # make a dictionary that can be iterated over (avoids code repetition in the loop)
    brDict = {
    "waveS0":[waveS0, b1], "waveS1":[waveS1, b2], "waveS2":[waveS2, b3],
    "waveS3":[waveS3, b4], "waveS4":[waveS4, b5], "waveS5":[waveS5, b6],
    "wpar4":[wpar4, b7], "waveEnt":[waveEnt, b8], "tOffset":[tOffset, b9],
    "bcMax":[bcMax, b10], "bcMin":[bcMin, b11],
    "derivMax":[derivMax, b12], "derivTime":[derivTime, b13],
    "bandMax":[bandMax, b14], "bandTime":[bandTime, b15],
    "raw10":[raw10, b16], "raw50":[raw50, b17], "raw90":[raw90, b18],
    "den10":[den10, b19], "den50":[den50, b20], "den90":[den90, b21],
    "oppie":[oppie, b22],
    "fitMatch":[fitMatch, b23], "fitE":[fitE, b24], "fitSlo":[fitSlo, b25],
    "matchMax":[matchMax, b26], "matchWidth":[matchWidth, b27], "matchTime":[matchTime, b28],
    "pol0":[pol0, b29], "pol1":[pol1, b30], "pol2":[pol2, b31], "pol3":[pol3, b32],
    "exp0":[exp0, b33], "exp1":[exp1, b34],
    "fitMax":[fitMax,b35],
    "fails":[fails,b36]
    }


    # Make a figure (-i option: select different diagnostic plots)
    # NOTE: If you want to add a custom or one-off figure, PLEASE don't change plotNum 0-4.
    fig = plt.figure(figsize=(10,7), facecolor='w')
    p = []
    for i in range(1,8): p.append(plt.subplot())
    if plotNum==0:
        p[0] = plt.subplot(111)
    elif plotNum==1 or plotNum==2:
        p[0] = plt.subplot(211)
        p[1] = plt.subplot(212)
    elif plotNum==3:
        p[0] = plt.subplot2grid((2,5), (0,0), colspan=3)
        p[1] = plt.subplot2grid((2,5), (0,3), colspan=2)
        p[2] = plt.subplot2grid((2,5), (1,0), colspan=3)
    elif plotNum==4:
        p[0] = plt.subplot2grid((6,7), (0,0), colspan=4, rowspan=3) # data & fit
        p[1] = plt.subplot2grid((6,7), (3,0), colspan=4, rowspan=3) # matched filter
        p[2] = plt.subplot2grid((6,7), (0,4), colspan=3, rowspan=2) # trace 1
        p[3] = plt.subplot2grid((6,7), (2,4), colspan=3, rowspan=2) # trace 2
        p[4] = plt.subplot2grid((6,7), (4,4), colspan=3, rowspan=2) # trace 3
    elif plotNum==5:
        p[0] = plt.subplot(111) # bandTime plot
    elif plotNum==6:
        # waveform fit
        p[0] = plt.subplot2grid((6,6), (0,0), colspan=6, rowspan=4) # data & fit
        p[1] = plt.subplot2grid((6,6), (4,0), colspan=2, rowspan=2) # trace 1
        p[2] = plt.subplot2grid((6,6), (4,2), colspan=2, rowspan=2) # trace 2
        p[3] = plt.subplot2grid((6,6), (4,4), colspan=2, rowspan=2) # trace 3
    elif plotNum==7:
        p[0] = plt.subplot(111) # matched filter plot

    if not batMode: plt.show(block=False)


    # Make a signal template
    print "Generating signal template ..."
    tSamp, tR, tZ, tAmp, tST, tSlo = 5000, 0, 15, 100, 2500, 10
    # tOrig, tOrigTS = wl.MakeSiggenWaveform(tSamp,tR,tZ,tAmp,tST,tSlo) # Damn you to hell, PDSF
    templateFile = np.load("./data/lat_template.npz")
    if dsNum==2 or dsNum==6: templateFile = np.load("./data/lat_ds2template.npz")
    tOrig, tOrigTS = templateFile['arr_0'], templateFile['arr_1']


    # Load stuff from DS1 forced acq. runs
    npzfile = np.load("./data/fft_forcedAcqDS1.npz")
    noise_asd, noise_xFreq, avgPwrSpec, xPwrSpec, data_forceAcq, data_fft = npzfile['arr_0'],npzfile['arr_1'],npzfile['arr_2'],npzfile['arr_3'],npzfile['arr_4'],npzfile['arr_5']


    # Loop over events
    print "Starting event loop ..."
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
            plt.pause(0.00001)          # rapid-draw mode
        if iList >= nList: break        # bail out, goose!

        entry = gatTree.GetEntryNumber(iList);
        gatTree.LoadTree(entry)
        gatTree.GetEntry(entry)
        nChans = gatTree.channel.size()
        event = MGTEvent()
        if gatMode: event = bltTree.event

        # Reset all branch vectors
        # NOTE: The events sometimes contain 'straggler' hits that do not pass the
        # given TCut.  This line sets ALL the new parameters to -88888 by default.
        # If you see this value in a plot, then you must be including hits that
        # passed the cut in wave-skim but did not pass the (different?) cut in LAT.
        for key in brDict: brDict[key][0].assign(nChans,-88888)
        brDict["fails"][0].assign(nChans,0) # set error code to 'true' by default
        errorCode = [0,0,0,0]

        # Loop over hits passing cuts
        numPass = gatTree.Draw("channel",theCut,"GOFF",1,iList)
        chans = gatTree.GetV1()
        chanList = list(set(int(chans[n]) for n in xrange(numPass)))
        hitList = (iH for iH in xrange(nChans) if gatTree.channel.at(iH) in chanList)  # a 'generator expression'
        for iH in hitList:

            # ------------------------------------------------------------------------
            # Waveform processing

            # load data
            run = gatTree.run
            chan = gatTree.channel.at(iH)
            dataENF = gatTree.trapENFCal.at(iH)
            dataENM = gatTree.trapENM.at(iH)
            dataMax = gatTree.trapENMSample.at(iH)*10. - 4000
            wf = MGTWaveform()
            iEvent = 0
            if gatMode:
                wf = event.GetWaveform(iH)
                iEvent = entry
            else:
                wf = gatTree.MGTWaveforms.at(iH)
                iEvent = gatTree.iEvent

            # print "%d:  run %d  chan %d  trapENFCal %.2f" % (iList, run, chan, dataENF)

            # be absolutely sure you're matching the right waveform to this hit
            if wf.GetID() != chan:
                print "ERROR -- Vector matching failed.  iList %d  run %d  iEvent %d" % (iList,run,iEvent)
                return


            # Let's start the show
            removeNBeg, removeNEnd = 0, 2
            if dsNum==2 or dsNum==6: removeNBeg = 3
            signal = wl.processWaveform(wf,removeNBeg,removeNEnd)
            data = signal.GetWaveBLSub()
            dataTS = signal.GetTS()
            tOffset[iH] = signal.GetOffset()
            _,dataNoise = signal.GetBaseNoise()

            # wavelet packet transform
            wp = pywt.WaveletPacket(data, 'db2', 'symmetric', maxlevel=4)
            nodes = wp.get_level(4, order='freq')
            waveletYTrans = np.array([n.data for n in nodes],'d')
            waveletYTrans = abs(waveletYTrans)

            # wavelet parameters

            waveS0[iH] = np.sum(waveletYTrans[0:1,1:-1])
            waveS1[iH] = np.sum(waveletYTrans[0:1,1:33])
            waveS2[iH] = np.sum(waveletYTrans[0:1,33:65])
            waveS3[iH] = np.sum(waveletYTrans[0:1,65:97])
            waveS4[iH] = np.sum(waveletYTrans[0:1,97:-1])
            waveS5[iH] = np.sum(waveletYTrans[2:-1,1:-1])
            wpar4[iH] = np.amax(waveletYTrans[0:1,1:-1])

            S6 = np.sum(waveletYTrans[2:9,1:33])
            S7 = np.sum(waveletYTrans[2:9,33:65])
            S8 = np.sum(waveletYTrans[2:9,65:97])
            S9 = np.sum(waveletYTrans[2:9,97:-1])
            S10 = np.sum(waveletYTrans[9:,1:33])
            S11 = np.sum(waveletYTrans[9:,33:65])
            S12 = np.sum(waveletYTrans[9:,65:97])
            S13 = np.sum(waveletYTrans[9:,97:-1])
            sumList = [S6, S7, S8, S9, S10, S11, S12, S13]
            bcMax[iH] = np.max(sumList)
            bcMin[iH] = np.min(sumList)

            # reconstruct waveform w/ only lowest frequency.
            new_wp = pywt.WaveletPacket(data=None, wavelet='db2', mode='symmetric')
            new_wp['aaa'] = wp['aaa'].data
            data_wlDenoised = new_wp.reconstruct(update=False)
            # resize in a smart way
            diff = len(data_wlDenoised) - len(data)
            if diff > 0: data_wlDenoised = data_wlDenoised[diff:]

            # calculate wf entropy
            d1 = 2. * np.multiply(waveletYTrans[0:1,1:65], np.log(waveletYTrans[0:1,1:65]/waveS0[iH]/2.0))
            d2 = np.multiply(waveletYTrans[0:1,65:-1], np.log(waveletYTrans[0:1,65:-1]/waveS0[iH]))
            waveEnt[iH] = np.abs(np.sum(d1)) + np.abs(np.sum(d2))


            # waveform high/lowpass filters - parameters are a little arbitrary

            B1,A1 = butter(2, [1e5/(1e8/2),1e6/(1e8/2)], btype='bandpass')
            data_bPass = lfilter(B1, A1, data)

            B2, A2 = butter(1, 0.08)
            data_filt = filtfilt(B2, A2, data)
            data_filtDeriv = wl.wfDerivative(data_filt)
            filtAmp = np.amax(data_filtDeriv) # scale the max to match the amplitude
            data_filtDeriv = data_filtDeriv * (dataENM / filtAmp)

            B3, A3 = butter(2,1e6/(1e8/2), btype='lowpass')
            data_lPass = lfilter(B3, A3, data)

            idx = np.where((dataTS > dataTS[0]+100) & (dataTS < dataTS[-1]-100))
            windowingOffset = dataTS[idx][0] - dataTS[0]

            bandMax[iH] = np.amax(data_bPass[idx])
            bandTime[iH] = dataTS[ np.argmax(data_bPass[idx])] + windowingOffset

            derivMax[iH] = np.amax(data_filtDeriv[idx])
            derivTime[iH] = dataTS[ np.argmax(data_filtDeriv[idx])] + windowingOffset


            # timepoints of raw and low-pass waveforms
            unnecessaryWF = wl.MGTWFFromNpArray(data)
            tpc = MGWFTimePointCalculator();
            tpc.AddPoint(.2)
            tpc.AddPoint(.5)
            tpc.AddPoint(.9)
            # tpc.FindTimePoints(wf) # FIXME: this should work but timepoints are drawed weirdedly
            tpc.FindTimePoints(unnecessaryWF)
            raw10[iH] = tpc.GetFromStartRiseTime(0)*10
            raw50[iH] = tpc.GetFromStartRiseTime(1)*10
            raw90[iH] = tpc.GetFromStartRiseTime(2)*10

            mgtLowPass = wl.MGTWFFromNpArray(data_lPass)
            tpc.FindTimePoints(mgtLowPass)
            den10[iH] = tpc.GetFromStartRiseTime(0)*10
            den50[iH] = tpc.GetFromStartRiseTime(1)*10
            den90[iH] = tpc.GetFromStartRiseTime(2)*10


            # load and scale template waveform

            temp, tempTS = tOrig, tOrigTS
            temp = temp * (dataENM / tAmp)         # scale by amplitudes (not energies)
            tempMax = np.argmax(temp) * 10         # convert to ns
            tempTS = tempTS - (tempMax - dataMax)  # align @ max of rising edge
            tempMax = tempTS[np.argmax(temp)]      # get the new max TS after the shifting
            tempE = dataENF # set 'calibrated' energy of template equal to data's trapENFCal.

            # set window, guess parameters, pack into lists
            loWin, hiWin = dataTS[0], dataTS[-1]
            mt, en, slo = dataMax, dataENF, 20.  # guesses
            InterpFn = interpolate.interp1d(tempTS, temp, kind="linear", copy="False", assume_sorted="True")
            dataList = [data, dataTS, dataENF, dataMax, loWin, hiWin, dataNoise]
            tempList = [temp, tempTS, tempE, tempMax]
            datas = [dataList, tempList, InterpFn]
            floats = [mt,en,slo]

            # save the initial guess
            guess, guessTS = wm.MakeModel(dataList, tempList, floats, fn=InterpFn)

            # run waveform fitter
            # nelder-mead (~0.5 sec).
            MakeTracesGlobal()
            result = op.minimize(findLnLike, floats, args=datas, method="Nelder-Mead")#,options={"maxiter":10})
            if not result["success"]:
                print "fit 'fail': ", result["message"]
                errorCode[0] = 1
            mt, en, slo = result["x"]
            nelder, nelderTS = wm.MakeModel(dataList, tempList, [mt,en,slo], fn=InterpFn)

            # save parameters.  take absval for slowness only.
            slo = abs(slo)
            fitMatch[iH], fitE[iH], fitSlo[iH] = mt, en, slo


            # optimal matched filter (freq. domain)

            # we use the guess (not the fit result) to keep this independent of the wf fitter.
            data_fft = np.fft.fft(data) # can also try taking fft of the low-pass data
            temp_fft = np.fft.fft(guess)

            datafreq = np.fft.fftfreq(data.size) * 1e8
            power_vec = np.interp(datafreq, noise_xFreq, noise_asd) # load power spectra from file

            # Apply the filter
            optimal = data_fft * temp_fft.conjugate() / power_vec
            optimal_time = 2 * np.fft.ifft(optimal)

            # Normalize the output
            df = np.abs(datafreq[1] - datafreq[0]) # freq. bin size
            sigmasq = 2 * (temp_fft * temp_fft.conjugate() / power_vec).sum() * df
            sigma = np.sqrt(np.abs(sigmasq))
            SNR = abs(optimal_time) / (sigma)
            oppie[iH] = np.amax(SNR)


            # time domain match filter

            match, matchTS = wm.MakeModel(dataList, tempList, [mt,en,slo], opt="nowindow")
            match = match[::-1]

            # line up the max of the match with the max of the best-fit signal
            modelMaxTime = nelderTS[np.argmax(nelder)]
            matchMaxTime = matchTS[np.argmax(match)]
            matchTS = matchTS + (modelMaxTime - matchMaxTime)
            fitMax[iH] = modelMaxTime

            idx = np.where((matchTS >= dataTS[0]-5) & (matchTS <= dataTS[-1]+5))
            match, matchTS = match[idx], matchTS[idx]

            # make sure match is same length as data, then compute convolution and parameters
            matchMax[iH], matchWidth[iH], matchTime[iH] = -888, -888, -888
            smoothMF = np.zeros(len(dataTS))
            if len(match)!=len(data):
                # print "array mismatch: len match %i  len data %i " % (len(match),len(data))
                errorCode[1] = 1
            else:
                # I can't decide which of these is better.
                # smoothMF = gaussian_filter(match * data, sigma=float(5))
                smoothMF = gaussian_filter(match * data_lPass,sigma=float(5))

                # compute match filter parameters
                matchMax[iH] = np.amax(smoothMF)
                matchTime[iH] = matchTS[ np.argmax(smoothMF) ]
                idx = np.where(smoothMF > matchMax[iH]/2.)
                if len(matchTS[idx]>1): matchWidth[iH] = matchTS[idx][-1] - matchTS[idx][0]


            # Fit tail slope (2 methods).  Guard against fit fails

            idx = np.where(dataTS >= modelMaxTime)
            tail, tailTS = data[idx], dataTS[idx]
            popt1,popt2 = 0,0
            try:
                popt1,_ = op.curve_fit(wl.tailModelPol, tailTS, tail)
                pol0[iH], pol1[iH], pol2[iH], pol3[iH] = popt1[0], popt1[1], popt1[2], popt1[3]
            except:
                # print "curve_fit tailModelPol failed, run %i  event %i  channel %i" % (run, iList, chan)
                errorCode[2] = 1
                pass
            try:
                popt2,_ = op.curve_fit(wl.tailModelExp, tailTS, tail, p0=[dataENM,72000])
                exp0[iH], exp1[iH] = popt2[0], popt2[1]
            except:
                # print "curve_fit tailModelExp failed, run %i  event %i  channel %i" % (run, iList, chan)
                errorCode[3] = 1
                pass



            # ------------------------------------------------------------------------
            # End waveform processing.

            # Calculate error code
            fails[iH] = 0
            for i,j in enumerate(errorCode):
                if j==1: fails[iH] += int(j)<<i
            # print "fails:",fails[iH]

            # Make plots!
            if batMode: continue
            if plotNum==0: # raw data
                p[0].cla()
                p[0].plot(dataTS,data,'b')
                p[0].set_title("Run %d  Entry %d  Channel %d  ENFCal %.2f" % (run,iEvent,chan,dataENF))

            elif plotNum==1: # wavelet plot
                p[0].cla()
                p[0].plot(dataTS,data,'b')
                p[0].plot(dataTS,data_wlDenoised,'r')
                p[0].set_title("Run %d  Entry %d  Channel %d  ENFCal %.2f" % (run,iEvent,chan,dataENF))

                p[1].cla()
                p[1].imshow(waveletYTrans, interpolation='nearest', aspect="auto", origin="lower",extent=[0, 1, 0, len(waveletYTrans)],cmap='jet')
                p[1].set_title("S5/E %.2f  (S3-S1)/E %.2f  bcMax/bcMin %.2f" % (waveS5[iH]/dataENF, (waveS3[iH]-waveS1[iH])/dataENF, bcMax[iH]/bcMin[iH]))

            elif plotNum==2: # time points, bandpass filters, tail slope
                p[0].cla()
                p[0].plot(dataTS,data,color='blue',label='data')
                p[0].axvline(raw10[iH],color='red',label='rawTP')
                p[0].axvline(raw50[iH],color='red')
                p[0].axvline(raw90[iH],color='red')
                p[0].axvline(den10[iH],color='black',label='lpTP')
                p[0].axvline(den50[iH],color='black')
                p[0].axvline(den90[iH],color='black')
                p[0].plot(nelderTS,nelder,color='magenta',label='bestfit')
                if errorCode[2]!=1: p[0].plot(tailTS, wl.tailModelPol(tailTS, *popt1), color='orange',linewidth=2, label='tailPol')
                if errorCode[3]!=1: p[0].plot(tailTS, wl.tailModelExp(tailTS, *popt2), color='magenta',linewidth=2, label='tailExp')
                p[0].legend(loc=4)
                p[0].set_title("Run %d  Entry %d  Channel %d  ENFCal %.2f" % (run,iEvent,chan,dataENF))

                p[1].cla()
                p[1].plot(dataTS,data_lPass,color='blue',label='lowpass')
                p[1].plot(dataTS,data_filtDeriv,color='green',label='filtDeriv')
                p[1].plot(dataTS,data_filt,color='black',label='filtfilt')
                p[1].plot(dataTS,data_bPass,color='red',label='bpass')
                p[1].axvline(bandTime[iH],color='orange',label='bandTime')
                p[1].axvline(derivTime[iH],color='magenta',label='derivTime')
                p[1].legend(loc=4)

            elif plotNum==3: # matched filter - freq
                p[0].cla()
                p[0].plot(dataTS,data,'b')
                p[0].plot(guessTS,guess,'r')
                p[0].plot(nelderTS,nelder,color='cyan')
                p[0].set_title("Run %d  Entry %d  Channel %d  ENFCal %.2f" % (run,iEvent,chan,dataENF))

                data_asd, data_xFreq = plt.psd(data, Fs=1e8, NFFT=2048, pad_to=2048, visible=False)
                temp_asd, temp_xFreq = plt.psd(temp, Fs=1e8, NFFT=2048, pad_to=2048, visible=False)

                p[1].cla()
                p[1].loglog(data_xFreq, np.sqrt(data_asd), 'b')
                p[1].loglog(noise_xFreq, np.sqrt(noise_asd), 'g')
                p[1].loglog(temp_xFreq, np.sqrt(temp_asd), 'r')
                p[1].set_xlabel('Frequency (Hz)')
                p[1].set_ylabel('ASD')
                p[1].grid('on')

                p[2].cla()
                p[2].plot(dataTS, SNR)
                p[2].set_title('oppie %.1f' % (oppie[iH]))
                p[2].set_xlabel('Offset time (s)')
                p[2].set_ylabel('SNR')

            elif plotNum==4: # matched filter - time + fitter
                p[0].cla()
                p[0].plot(dataTS,data,'b',label='data')
                idx = np.where((tempTS >= dataTS[0]-5) & (tempTS <= dataTS[-1]+5))
                p[0].plot(tempTS[idx],temp[idx],'r',label='template')
                p[0].plot(nelderTS,trap,color='k',label='bestfit-trap')
                p[0].plot(nelderTS,nelder,color='cyan',label='bestfit')
                p[0].axvline(fitMatch[iH],color='magenta',label='fitMatch')
                p[0].axvline(fitMax[iH],color='orange',label='fitMax')
                p[0].set_xlabel('Time (s)')
                p[0].set_ylabel('Voltage (arb)')
                p[0].set_title("Run %d  Entry %d  Channel %d  ENFCal %.2f" % (run,iEvent,chan,dataENF))
                p[0].legend(loc=4)

                p[1].cla()
                p[1].plot(dataTS,data,'b',alpha=0.8,label='data')
                p[1].plot(matchTS,match,'k',label='match')
                p[1].plot(nelderTS,nelder,color='orange',label='bestFit')
                p[1].plot(dataTS,data_lPass,'r',label='lowpass')
                if len(match)==len(data):
                    p[1].plot(matchTS,smoothMF,'g',label='smoothMF')
                idx = np.where(smoothMF > matchMax[iH]/2.)
                p[1].axvline(matchTS[idx][-1],color='green',alpha=0.7)
                p[1].axvline(matchTS[idx][0],color='green',alpha=0.7)
                p[1].set_title("mMax %.2f  mWidth %.2f  mMax/E %.2f" % (matchMax[iH],matchWidth[iH],matchMax[iH]/dataENF))
                p[1].legend(loc=4)

                p[2].cla()
                p[2].set_title("maxTime %.1f  Energy %.2f  Slow %.1f" % (mt,en,slo))
                p[2].plot(mtTrace[1:])
                p[2].set_ylabel('maxTime')
                p[3].cla()
                p[3].plot(enTrace[1:])
                p[3].set_ylabel('energy')
                p[4].cla()
                p[4].plot(sloTrace[1:])
                p[4].set_ylabel('slowness')

            if plotNum==5: # bandTime plot
                p[0].cla()
                p[0].plot(dataTS,data,color='blue',label='data',alpha=0.7)
                p[0].plot(dataTS,data_lPass,color='magenta',label='lowpass',linewidth=4)
                p[0].plot(dataTS,data_bPass,color='red',label='bpass',linewidth=4)
                p[0].axvline(bandTime[iH],color='orange',label='bandTime',linewidth=4)
                p[0].legend(loc=4)
                p[0].set_xlabel('Time (ns)')
                p[0].set_ylabel('ADC (arb)')
                p[0].set_title("Run %d  Entry %d  Channel %d  ENFCal %.2f" % (run,iEvent,chan,dataENF))

            if plotNum==6: # waveform fit plot
                p[0].cla()
                p[0].plot(dataTS,data,color='blue',label='data',alpha=0.8)
                idx = np.where((tempTS >= dataTS[0]-5) & (tempTS <= dataTS[-1]+5))
                p[0].plot(tempTS[idx],temp[idx],color='orange',label='template')
                p[0].plot(nelderTS,nelder,color='red',label='bestfit',linewidth=3)
                p[0].axvline(fitMatch[iH],color='magenta',label='fitMatch',linewidth=4,alpha=0.5)
                p[0].set_title("Run %d  Entry %d  Channel %d  ENFCal %.2f  fitMatch %.1f  fitE %.2f  fitSlo %.1f" % (run,iEvent,chan,dataENF,mt,en,slo))
                p[0].legend(loc=4)
                p[1].cla()
                p[1].plot(mtTrace[1:],label='fitMatch',color='red')
                p[1].legend(loc=4)
                p[1].set_xlabel('Fit Steps')
                p[2].cla()
                p[2].plot(enTrace[1:],label='fitE',color='green')
                p[2].legend(loc=4)
                p[3].cla()
                p[3].plot(sloTrace[1:],label='fitSlo',color='blue')
                p[3].legend(loc=4)

            if plotNum==7: # match filter plot
                p[0].cla()
                p[0].plot(dataTS,data,color='blue',label='data',alpha=0.7)
                p[0].plot(nelderTS,nelder,color='red',label='bestfit',linewidth=3)
                p[0].axvline(matchTime[iH],color='orange',label='matchTime',linewidth=2)
                p[0].plot(matchTS,smoothMF,color='magenta',label='smoothMF',linewidth=4)
                p[0].plot(matchTS,match,color='cyan',label='match',linewidth=3)
                p[0].set_xlabel('Time (s)')
                p[0].set_ylabel('Voltage (arb)')
                p[0].legend(loc=4)
                p[0].set_title("Run %d  Entry %d  Channel %d  ENFCal %.2f  matchMax %.2f  matchTime %.2f  matchWidth %.2f" % (run,iEvent,chan,dataENF,matchMax[iH],matchTime[iH],matchWidth[iH]))

            plt.tight_layout()
            plt.pause(0.000001)
            # ------------------------------------------------------------------------

        # End loop over hits, fill branches
        if batMode:
            for key in brDict:
                brDict[key][1].Fill()
            if iList%5000 == 0 and iList!=0:
                oTree.Write("",TObject.kOverwrite)
                print "%d / %d entries saved (%.2f %% done), time: %s" % (iList,nList,100*(float(iList)/nList),time.strftime('%X %x %Z'))

    # End loop over events
    if batMode:
        oTree.Write("",TObject.kOverwrite)
        print "Wrote",oTree.GetBranch("channel").GetEntries(),"entries in the copied tree,"
        print "and wrote",b1.GetEntries(),"entries in the new branches."

    stopT = time.clock()
    print "Stopped:",time.strftime('%X %x %Z'),"\nProcess time (min):",(stopT - startT)/60
    print float(nList)/((stopT-startT)/60.),"entries per minute."


def findLnLike(floats, datas):
    # SciPy Minimizer.
    #     This could be moved to waveModel.py IF you decide you
    #     don't need to look at the trace arrays anymore.
    #     Also, this version tries to align the MAX time
    #     (not the start time) of the waveforms.
    global mtTrace, enTrace, sloTrace

    # Unpack parameters.
    dataList, tempList, InterpFn = datas
    data, dataTS, dataE, dataMax, loWin, hiWin, dataNoise = dataList
    temp, tempTS, tempE, tempMax = tempList

    # Can fix parameters here by setting them back to their input guess values
    mt, en, slo = 0, 0, 0
    if len(floats)==1: # basinhopping case
        slo, mt, en = floats[0], dataMax, dataE
    if len(floats)==3: # minimize case (no fixing)
        mt, en, slo = floats

    # Make a trace
    mtTrace = np.append(mtTrace,mt)
    enTrace = np.append(enTrace,en)
    sloTrace = np.append(sloTrace,slo)
    # print mt, en, slo

    # Build the model to compare with data
    model, modelTS = wm.MakeModel(dataList, tempList, [mt,en,slo], fn=InterpFn)

    # Find LL of data vs. model
    lnLike = 0.5 * np.sum ( np.power((data-model)/dataNoise, 2) - np.log( 1 / np.power(dataNoise,2) ) )
    return lnLike


def MakeTracesGlobal():
    # This is so 'findLnLike' can write to the trace arrays.
    # It has to remain in this file to work.
    tmp1 = tmp2 = tmp3 = np.empty([1,])
    global mtTrace, enTrace, sloTrace
    mtTrace, enTrace, sloTrace = tmp1, tmp2, tmp3


if __name__ == "__main__":
    main(sys.argv[1:])