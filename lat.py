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
v2: 04 Aug 2017 - improvements to wf fitting, handle multisampling, etc.

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
import scipy.special as sp
import waveLibs as wl
import waveModel as wm
limit = sys.float_info.max # equivalent to std::numeric_limits::max() in C++

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
    import matplotlib.ticker as mtick

    # File I/O
    inFile, outFile, bltFile = TFile(), TFile(), TFile()
    gatTree, bltTree, oTree = TTree(), TTree(), TTree()
    theCut, inPath, outPath = "", "", ""

    # Set input and output files
    if rangeMode:
        inPath = "%s/waveSkimDS%d_%d.root" % (pathToInput, dsNum, subNum)
        outPath = "%s/latSkimDS%d_%d.root" % (pathToOutput, dsNum, subNum)
    if fileMode:
        inPath = "%s/waveSkimDS%d_run%d.root" % (pathToInput, dsNum, runNum)
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
    theCut += " && trapENFCal < 6 && trapENFCal > 1"
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
    if batMode and not intMode:
        outFile = TFile(outPath, "RECREATE")
        print "Attempting tree copy to",outPath
        oTree = gatTree.CopyTree("")
        oTree.Write()
        print "Wrote",oTree.GetEntries(),"entries."
        cutUsed = TNamed("theCut",theCut)
        cutUsed.Write()

    waveS1, waveS2 = std.vector("double")(), std.vector("double")()
    waveS3, waveS4, waveS5 = std.vector("double")(), std.vector("double")(), std.vector("double")()
    tOffset = std.vector("double")()
    bcMax, bcMin = std.vector("double")(), std.vector("double")()
    bandMax, bandTime = std.vector("double")(), std.vector("double")()
    den10, den50, den90 = std.vector("double")(), std.vector("double")(), std.vector("double")()
    oppie = std.vector("double")()
    fitMu, fitAmp, fitSlo = std.vector("double")(), std.vector("double")(), std.vector("double")()
    fitTau, fitBL = std.vector("double")(), std.vector("double")()
    matchMax, matchWidth, matchTime = std.vector("double")(), std.vector("double")(), std.vector("double")()
    pol0, pol1, pol2, pol3 = std.vector("double")(), std.vector("double")(), std.vector("double")(), std.vector("double")()
    fails, fitChi2, fitLL = std.vector("int")(), std.vector("double")(), std.vector("double")()
    wpRiseNoise = std.vector("double")()
    t0_SLE, t0_ALE, lat, latF = std.vector("double")(), std.vector("double")(), std.vector("double")(), std.vector("double")()
    latAF, latFC, latAFC = std.vector("double")(), std.vector("double")(), std.vector("double")()

    # It's not possible to put the "oTree.Branch" call into a class initializer (waveLibs::latBranch). You suck, ROOT.
    b1, b2 = oTree.Branch("waveS1",waveS1), oTree.Branch("waveS2",waveS2)
    b3, b4, b5 = oTree.Branch("waveS3",waveS3), oTree.Branch("waveS4",waveS4), oTree.Branch("waveS5",waveS5)
    b6 = oTree.Branch("tOffset",tOffset)
    b7, b8 = oTree.Branch("bcMax",bcMax), oTree.Branch("bcMin",bcMin)
    b9, b10 = oTree.Branch("bandMax",bandMax), oTree.Branch("bandTime",bandTime)
    b11, b12, b13 = oTree.Branch("den10",den10), oTree.Branch("den50",den50), oTree.Branch("den90",den90)
    b14 = oTree.Branch("oppie",oppie)
    b15, b16, b17 = oTree.Branch("fitMu", fitMu), oTree.Branch("fitAmp", fitAmp), oTree.Branch("fitSlo", fitSlo)
    b18, b19 = oTree.Branch("fitTau",fitTau), oTree.Branch("fitBL",fitBL)
    b20, b21, b22 = oTree.Branch("matchMax", matchMax), oTree.Branch("matchWidth", matchWidth), oTree.Branch("matchTime", matchTime)
    b23, b24, b25, b26 = oTree.Branch("pol0", pol0), oTree.Branch("pol1", pol1), oTree.Branch("pol2", pol2), oTree.Branch("pol3", pol3)
    b27, b28, b29 = oTree.Branch("fails",fails), oTree.Branch("fitChi2",fitChi2), oTree.Branch("fitLL",fitLL)
    b30 = oTree.Branch("wpRiseNoise",wpRiseNoise)
    b31, b32, b33, b34 = oTree.Branch("t0_SLE",t0_SLE), oTree.Branch("t0_ALE",t0_ALE), oTree.Branch("lat",lat), oTree.Branch("latF",latF)
    b35, b36, b37 = oTree.Branch("latAF",latAF), oTree.Branch("latFC",latFC), oTree.Branch("latAFC",latAFC)

    # make a dictionary that can be iterated over (avoids code repetition in the loop)
    brDict = {
        "waveS1":[waveS1, b1], "waveS2":[waveS2, b2],
        "waveS3":[waveS3, b3], "waveS4":[waveS4, b4], "waveS5":[waveS5, b5],
        "tOffset":[tOffset, b6],
        "bcMax":[bcMax, b7], "bcMin":[bcMin, b8],
        "bandMax":[bandMax, b9], "bandTime":[bandTime, b10],
        "den10":[den10, b11], "den50":[den50, b12], "den90":[den90, b13],
        "oppie":[oppie, b14],
        "fitMu":[fitMu, b15], "fitAmp":[fitAmp, b16], "fitSlo":[fitSlo, b17],
        "fitTau":[fitTau, b18], "fitBL":[fitBL,b19],
        "matchMax":[matchMax, b20], "matchWidth":[matchWidth, b21], "matchTime":[matchTime, b22],
        "pol0":[pol0, b23], "pol1":[pol1, b24], "pol2":[pol2, b25], "pol3":[pol3, b26],
        "fails":[fails,b27], "fitChi2":[fitChi2,b28], "fitLL":[fitLL,b29],
        "wpRiseNoise":[wpRiseNoise,b30],
        "t0_SLE":[t0_SLE,b31], "t0_ALE":[t0_ALE,b32], "lat":[lat,b33], "latF":[latF,b34],
        "latAF":[latAF,b35], "latFC":[latFC,b36], "latAFC":[latAFC,b37]
    }


    # Make a figure (-i option: select different plots)
    fig = plt.figure(figsize=(12,7), facecolor='w')
    p = []
    for i in range(1,8): p.append(plt.subplot())
    if plotNum==0 or plotNum==7:
        p[0] = plt.subplot(111)  # 0-raw waveform, 7-new trap filters
    elif plotNum==1 or plotNum==2:
        p[0] = plt.subplot(211)  # 1-wavelet, 2-time points, bandpass filters, tail slope
        p[1] = plt.subplot(212)
    elif plotNum==3:
        p[0] = plt.subplot2grid((2,5), (0,0), colspan=3)  # oppie / freq-domain matched filter
        p[1] = plt.subplot2grid((2,5), (0,3), colspan=2)
        p[2] = plt.subplot2grid((2,5), (1,0), colspan=3)
    elif plotNum==4:
        p[0] = plt.subplot(111)  # time-domain matched filter
    elif plotNum==5:
        p[0] = plt.subplot(111)  # bandpass / bandTime
    elif plotNum==6:
        p[0] = plt.subplot2grid((6,10), (0,0), colspan=10, rowspan=3) # waveform fit
        p[1] = plt.subplot2grid((6,10), (3,0), colspan=10, rowspan=1) # residual
        p[2] = plt.subplot2grid((6,10), (4,0), colspan=2, rowspan=2) # traces
        p[3] = plt.subplot2grid((6,10), (4,2), colspan=2, rowspan=2)
        p[4] = plt.subplot2grid((6,10), (4,4), colspan=2, rowspan=2)
        p[5] = plt.subplot2grid((6,10), (4,6), colspan=2, rowspan=2)
        p[6] = plt.subplot2grid((6,10), (4,8), colspan=2, rowspan=2)
    if not batMode: plt.show(block=False)


    # Load a fast signal template - used w/ the freq-domain matched filter
    # print "Generating signal template ..."
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
            dataENFCal = gatTree.trapENFCal.at(iH)
            dataENM = gatTree.trapENM.at(iH)
            dataTSMax = gatTree.trapENMSample.at(iH)*10. - 4000
            wf = MGTWaveform()
            iEvent = 0
            if gatMode:
                wf = event.GetWaveform(iH)
                iEvent = entry
            else:
                wf = gatTree.MGTWaveforms.at(iH)
                iEvent = gatTree.iEvent

            # print "%d:  run %d  chan %d  trapENFCal %.2f" % (iList, run, chan, dataENFCal)

            # be absolutely sure you're matching the right waveform to this hit
            if wf.GetID() != chan:
                print "ERROR -- Vector matching failed.  iList %d  run %d  iEvent %d" % (iList,run,iEvent)
                return

            # Let's start the show
            # remove first 4 samples when we have multisampling
            # remove last 2 samples to get rid of the ADC spike at the end of all wf's.
            removeNBeg, removeNEnd = 0, 2
            if dsNum==6 or dsNum==2: removeNBeg = 4
            signal = wl.processWaveform(wf,removeNBeg,removeNEnd)
            data = signal.GetWaveRaw()
            data_blSub = signal.GetWaveBLSub()
            dataTS = signal.GetTS()
            tOffset[iH] = signal.GetOffset()
            dataBL,dataNoise = signal.GetBaseNoise()

            # wavelet packet transform
            wp = pywt.WaveletPacket(data_blSub, 'db2', 'symmetric', maxlevel=4)
            nodes = wp.get_level(4, order='freq')
            wpCoeff = np.array([n.data for n in nodes],'d')
            wpCoeff = abs(wpCoeff)

            # wavelet parameters
            waveS1[iH] = np.sum(wpCoeff[0:1,1:33])
            waveS2[iH] = np.sum(wpCoeff[0:1,33:65])
            waveS3[iH] = np.sum(wpCoeff[0:1,65:97])
            waveS4[iH] = np.sum(wpCoeff[0:1,97:-1])
            waveS5[iH] = np.sum(wpCoeff[2:-1,1:-1])
            S6 = np.sum(wpCoeff[2:9,1:33])
            S7 = np.sum(wpCoeff[2:9,33:65])
            S8 = np.sum(wpCoeff[2:9,65:97])
            S9 = np.sum(wpCoeff[2:9,97:-1])
            S10 = np.sum(wpCoeff[9:,1:33])
            S11 = np.sum(wpCoeff[9:,33:65])
            S12 = np.sum(wpCoeff[9:,65:97])
            S13 = np.sum(wpCoeff[9:,97:-1])
            sumList = [S6, S7, S8, S9, S10, S11, S12, S13]
            bcMax[iH] = np.max(sumList)
            bcMin[iH] = np.min(sumList)

            # reconstruct waveform w/ only lowest frequency.
            new_wp = pywt.WaveletPacket(data=None, wavelet='db2', mode='symmetric')
            new_wp['aaa'] = wp['aaa'].data
            data_wlDenoised = new_wp.reconstruct(update=False)
            # resize in a smart way
            diff = len(data_wlDenoised) - len(data_blSub)
            if diff > 0: data_wlDenoised = data_wlDenoised[diff:]

            # waveform high/lowpass filters - parameters are a little arbitrary

            B1,A1 = butter(2, [1e5/(1e8/2),1e6/(1e8/2)], btype='bandpass')
            data_bPass = lfilter(B1, A1, data_blSub)

            B2, A2 = butter(1, 0.08)
            data_filt = filtfilt(B2, A2, data_blSub)
            data_filtDeriv = wl.wfDerivative(data_filt)
            filtAmp = np.amax(data_filtDeriv) # scale the max to match the amplitude
            data_filtDeriv = data_filtDeriv * (dataENM / filtAmp)

            B3, A3 = butter(2,1e6/(1e8/2), btype='lowpass')
            data_lPass = lfilter(B3, A3, data_blSub)

            idx = np.where((dataTS > dataTS[0]+100) & (dataTS < dataTS[-1]-100))
            windowingOffset = dataTS[idx][0] - dataTS[0]

            bandMax[iH] = np.amax(data_bPass[idx])
            bandTime[iH] = dataTS[ np.argmax(data_bPass[idx])] - windowingOffset


            # timepoints of low-pass waveforms
            tpc = MGWFTimePointCalculator();
            tpc.AddPoint(.2)
            tpc.AddPoint(.5)
            tpc.AddPoint(.9)
            mgtLowPass = wl.MGTWFFromNpArray(data_lPass)
            tpc.FindTimePoints(mgtLowPass)
            den10[iH] = tpc.GetFromStartRiseTime(0)*10
            den50[iH] = tpc.GetFromStartRiseTime(1)*10
            den90[iH] = tpc.GetFromStartRiseTime(2)*10


            # ================ xgauss waveform fitting ================

            amp, mu, sig, tau, bl = dataENM, dataTSMax, 600., -72000., dataBL
            floats = np.asarray([amp, mu, sig, tau, bl])
            temp = xgModelWF(dataTS, floats)
            MakeTracesGlobal()

            # get the noise of the denoised wf
            denoisedNoise,_,_ = wl.baselineParameters(data_wlDenoised)

            # NOTE: fit is to wavelet-denoised data, BECAUSE there are no HF components in the model,
            # AND we'll still calculate fitChi2 w/r/t the data, not the denoised data.
            # datas = [dataTS, data, dataNoise] # fit data
            datas = [dataTS, data_wlDenoised + dataBL, denoisedNoise] # fit wavelet-denoised data w/ Bl added back in

            # L-BGFS-B
            bnd = ((None,None),(None,None),(None,None),(None,None),(None,None)) # A,mu,sig,tau.  Unbounded rn.
            # opts = {'disp': None,   # None, True to print convergence messages
            #         'maxls': 100,   # 20, max line search steps.  This option doesn't exist on PDSF (scipy 0.15)
            #         'iprint': -1,   # -1
            #         'gtol': 1e-08,  # 1e-05
            #         'eps': 1e-08,   # 1e-08
            #         'maxiter': 15000,   # 15000
            #         'ftol': 2.220446049250313e-09,
            #         'maxcor': 10,   # 10
            #         'maxfun': 15000}    # 15000

            # numerical gradient - seems trustworthy
            start = time.clock()
            result = op.minimize(lnLike, floats, args=datas, method="L-BFGS-B", options=None, bounds=None)
            fitSpeed = time.clock() - start

            if not result["success"]:
                # print "fit 'fail': ", result["message"]
                errorCode[0] = 1

            # save parameters
            amp, mu, sig, tau, bl = result["x"]
            fitMu[iH], fitAmp[iH], fitSlo[iH], fitTau[iH], fitBL[iH] = mu, amp, sig, tau, bl
            floats = [amp, mu, sig, tau, bl]
            fit = xgModelWF(dataTS, floats)

            # log-likelihood of this fit
            fitLL[iH] = result["fun"]

            # chi-square of this fit
            # Textbook is (observed - expected)^2 / expected,
            # but we'll follow MGWFCalculateChiSquare.cc and do (observed - expected)^2 / NDF.
            # NOTE: we're doing the chi2 against the DATA, though the FIT is to the DENOISED DATA.
            fitChi2[iH] = np.sum(np.square(data-fit)) / len(data)


            # get wavelet coeff's for rising edge only.  normalize to bcMin
            # view this w/ plot 1

            # find the window of rising edge
            fit_blSub = fit - bl
            fitMaxTime = dataTS[np.argmax(fit_blSub)]
            fitStartTime = dataTS[0]
            idx = np.where(fit_blSub < 0.1)
            if len(dataTS[idx] > 0): fitStartTime = dataTS[idx][-1]
            fitRiseTime50 = (fitMaxTime + fitStartTime)/2.

            # bcMin is 32 samples long in the x-direction.
            # if we make the window half as wide, it'll have the same # of coeff's as bcMin.
            # this is still 'cheating' since we're not summing over the same rows.
            numXRows = wpCoeff.shape[1]
            wpCtrRise = int((fitRiseTime50 - dataTS[0]) / (dataTS[-1] - dataTS[0]) * numXRows)
            wpLoRise = wpCtrRise - 8
            if wpLoRise < 0: wpLoRise = 0
            wpHiRise = wpCtrRise + 8
            if wpHiRise > numXRows: wpHiRise = numXRows

            # sum all HF wavelet components for this edge.
            wpRiseNoise[iH] = np.sum(wpCoeff[2:-1,wpLoRise:wpHiRise]) / (bcMin[iH])

            # =========================================================

            # optimal matched filter (freq. domain)
            # we use the pysiggen fast template (not the fit result) to keep this independent of the wf fitter.

            # pull in the template, shift it, and make sure it's the same length as the data
            guessTS = tOrigTS - 15000.
            idx = np.where((guessTS > -5) & (guessTS < dataTS[-1]))
            guessTS, guess = guessTS[idx], tOrig[idx]
            if len(guess)!=len(data):
                if len(guess)>len(data):
                    guess, guessTS = guess[0:len(data)], guessTS[0:len(data)]
                else:
                    guess = np.pad(guess, (0,len(data)-len(guess)), 'edge')
                    guessTS = np.pad(guessTS, (0,len(data)-len(guessTS)), 'edge')

            data_fft = np.fft.fft(data_blSub) # can also try taking fft of the low-pass data
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


            # time-domain matched filter.  use the baseline-subtracted wf as data, and fit_blSub too.

            # make a longer best-fit waveform s/t it can be shifted L/R.
            matchTS = np.append(dataTS, np.arange(dataTS[-1], dataTS[-1] + 20000, 10)) # add 2000 samples
            match = xgModelWF(matchTS, [amp, mu+10000., sig, tau, bl]) # shift mu accordingly
            match = match[::-1] - bl # time flip and subtract off bl

            # line up the max of the 'match' (flipped wf) with the max of the best-fit wf
            matchMaxTime = matchTS[np.argmax(match)]
            matchTS = matchTS + (fitMaxTime - matchMaxTime)

            # resize match wf and the window function to have same # samples as dataTS.
            if matchTS[0] <= dataTS[0] and matchTS[-1] >= dataTS[-1]:
                idx = np.where((matchTS >= dataTS[0]) & (matchTS <= dataTS[-1]-5))
                match, matchTS = match[idx], matchTS[idx]

            elif matchTS[-1] < dataTS[-1]: # too early
                idx = np.where(matchTS >= dataTS[0])
                tmpData = match[idx]
                match = np.pad(tmpData, (0,len(data)-len(tmpData)), 'edge')
                matchTS, windowTS = dataTS, dataTS

            elif matchTS[0] > dataTS[0]: # too late
                idx = np.where(matchTS <= dataTS[-1])
                tmpData = match[idx]
                match = np.pad(tmpData, (len(data)-len(tmpData),0), 'edge')
                matchTS, windowTS = dataTS, dataTS

            # make sure match is same length as data, then compute convolution and parameters
            matchMax[iH], matchWidth[iH], matchTime[iH] = -888, -888, -888
            smoothMF, windMF = np.zeros(len(dataTS)), np.zeros(len(dataTS))
            if len(match)!=len(data):
                print "array mismatch: len match %i  len data %i " % (len(match),len(data))
                # errorCode[1] = 1
            else:
                # compute match filter parameters
                smoothMF = gaussian_filter(match * data_blSub, sigma=5.)
                matchMax[iH] = np.amax(smoothMF)
                matchTime[iH] = matchTS[ np.argmax(smoothMF) ]
                idx = np.where(smoothMF > matchMax[iH]/2.)
                if len(matchTS[idx]>1): matchWidth[iH] = matchTS[idx][-1] - matchTS[idx][0]


            # Fit tail slope to polynomial.  Guard against fit fails

            idx = np.where(dataTS >= fitMaxTime)
            tail, tailTS = data[idx], dataTS[idx]
            popt1,popt2 = 0,0
            try:
                popt1,_ = op.curve_fit(wl.tailModelPol, tailTS, tail)
                pol0[iH], pol1[iH], pol2[iH], pol3[iH] = popt1[0], popt1[1], popt1[2], popt1[3]
            except:
                # print "curve_fit tailModelPol failed, run %i  event %i  channel %i" % (run, iList, chan)
                errorCode[2] = 1
                pass

            # =========================================================

            # new trap filters.
            # params: t0_SLE, t0_ALE, lat, latF, latAF, latFC, latAFC

            # calculate trapezoids

            # standard trapezoid - prone to walking, less sensitive to noise.  use to find energy
            eTrap = wl.trapFilter(data_blSub, 400, 250, 7200.)
            eTrapTS = np.arange(0, len(eTrap)*10., 10)
            eTrapInterp = interpolate.interp1d(eTrapTS, eTrap)

            # short trapezoid - triggers more quickly, sensitive to noise.  use to find t0
            sTrap = wl.trapFilter(data_blSub, 100, 150, 7200.)
            sTrapTS = np.arange(0, len(sTrap)*10., 10)

            # asymmetric trapezoid - used to find the t0 only
            aTrap = wl.asymTrapFilter(data_blSub, 200, 100, 40, True)
            aTrapTS = np.arange(0, len(aTrap)*10., 10)


            # find leading edges (t0 times)

            # limit the range from 0 to 14us, and use an ADC threshold of 2.0 (like the data) for now ...
            t0_SLE[iH],_ = wl.walkBackT0(sTrap, eTrapTS[-1]+7000-4000-2000, 2., 0, 1000) # (in ns) finds leading edge from short trap
            t0_ALE[iH],_ = wl.walkBackT0(aTrap, eTrapTS[-1]+7000-4000-2000, 2., 0, 1000) # (in ns) finds leading edge from asymmetric trap

            # standard energy trapezoid w/ a baseline padded waveform
            data_pad = np.pad(data_blSub,(400,0),'symmetric')
            pTrap = wl.trapFilter(data_pad, 400, 250, 7200.)
            pTrapTS = np.linspace(0, len(pTrap)*10, len(pTrap))
            pTrapInterp = interpolate.interp1d(pTrapTS, pTrap)

            # calculate energy parameters
            # standard amplitude.  basically trapEM, but w/o NL correction if the input WF doesn't have it.
            lat[iH] = np.amax(eTrap)

            # standard amplitude with t0 from the shorter traps
            # If either fixed pickoff time (t0) is < 0, use the first sample as the amplitude (energy).
            latF[iH] = eTrapInterp( np.amax([t0_SLE[iH]-7000+4000+2000, 0.]) ) # This should be ~trapEF
            latAF[iH] = eTrapInterp( np.amax([t0_ALE[iH]-7000+4000+2000, 0.]) )

            # amplitude from padded trapezoid, with t0 from short traps and a correction function
            # function is under development.  currently: f() = exp(p0 + p1*E), p0 ~ 7.8, p1 ~ -0.45 and -0.66
            # functional walk back distance is *either* the minimum of the function value, or 5500 (standard value)

            # t0_corr = -7000+8000+2000 # no correction
            t0_corr = -7000+8000+2000 - np.amin([np.exp(7.8 - 0.45*lat[iH]),1000.])
            t0A_corr = -7000+8000+2000 - np.amin([np.exp(7.8 - 0.66*lat[iH]),1000.])

            latFC[iH] = pTrapInterp( np.amax([t0_SLE[iH] + t0_corr, 0.]) )
            latAFC[iH] = pTrapInterp( np.amax([t0_ALE[iH] + t0A_corr, 0.]) )


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
                p[0].set_title("Run %d  Entry %d  Channel %d  ENFCal %.2f" % (run,iList,chan,dataENFCal))

            if plotNum==1: # wavelet plot
                p[0].cla()
                p[0].margins(x=0)
                p[0].plot(dataTS,data_blSub,color='blue',label='data')
                p[0].plot(dataTS,data_wlDenoised,color='cyan',label='denoised',alpha=0.7)
                p[0].axvline(fitRiseTime50,color='green',label='fit 50%',linewidth=2)
                p[0].plot(dataTS,fit_blSub,color='red',label='bestfit',linewidth=2)
                p[0].set_title("Run %d  Entry %d  Channel %d  ENFCal %.2f  flo %.0f  fhi %.0f  fhi-flo %.0f" % (run,iList,chan,dataENFCal,fitStartTime,fitMaxTime,fitMaxTime-fitStartTime))
                p[0].legend(loc='best')

                p[1].cla()
                p[1].imshow(wpCoeff, interpolation='nearest', aspect="auto", origin="lower",extent=[0, 1, 0, len(wpCoeff)],cmap='jet')
                p[1].axvline(float(wpLoRise)/numXRows,color='orange',linewidth=2)
                p[1].axvline(float(wpHiRise)/numXRows,color='orange',linewidth=2)
                p[1].set_title("waveS5 %.2f  bcMax %.2f  bcMin %.2f  wpRiseNoise %.2f" % (waveS5[iH], bcMax[iH], bcMin[iH], wpRiseNoise[iH]))


            if plotNum==2: # time points, bandpass filters, tail slope
                p[0].cla()
                p[0].plot(dataTS,data,color='blue',label='data')
                p[0].axvline(den10[iH],color='black',label='lpTP')
                p[0].axvline(den50[iH],color='black')
                p[0].axvline(den90[iH],color='black')
                p[0].plot(dataTS,fit,color='magenta',label='bestfit')
                if errorCode[2]!=1: p[0].plot(tailTS, wl.tailModelPol(tailTS, *popt1), color='orange',linewidth=2, label='tailPol')
                p[0].legend(loc='best')
                p[0].set_title("Run %d  Entry %d  Channel %d  ENFCal %.2f" % (run,iEvent,chan,dataENFCal))

                p[1].cla()
                p[1].plot(dataTS,data_lPass,color='blue',label='lowpass')
                p[1].plot(dataTS,data_filtDeriv,color='green',label='filtDeriv')
                p[1].plot(dataTS,data_filt,color='black',label='filtfilt')
                p[1].plot(dataTS,data_bPass,color='red',label='bpass')
                p[1].axvline(bandTime[iH],color='orange',label='bandTime')
                p[1].legend(loc='best')

            if plotNum==3: # freq-domain matched filter
                p[0].cla()
                p[0].plot(dataTS,data,'b')
                p[0].plot(dataTS,temp,'r')
                p[0].plot(dataTS,fit,color='cyan')
                p[0].set_title("Run %d  Entry %d  Channel %d  ENFCal %.2f" % (run,iEvent,chan,dataENFCal))

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

            if plotNum==4: # time domain match filter plot
                p[0].cla()
                p[0].plot(dataTS,data_blSub,color='blue',label='data',alpha=0.7)
                p[0].plot(dataTS,fit_blSub,color='red',label='bestfit',linewidth=3)
                p[0].axvline(matchTime[iH],color='orange',label='matchTime',linewidth=2)
                p[0].plot(matchTS,smoothMF,color='magenta',label='smoothMF',linewidth=3)
                p[0].plot(matchTS,match,color='cyan',label='match',linewidth=3)
                p[0].set_xlabel('Time (s)')
                p[0].set_ylabel('Voltage (arb)')
                p[0].legend(loc='best')
                p[0].set_title("Run %d  Entry %d  Channel %d  ENFCal %.2f  matchMax %.2f  matchTime %.2f  matchWidth %.2f" % (run,iEvent,chan,dataENFCal,matchMax[iH],matchTime[iH],matchWidth[iH]))

            if plotNum==5: # bandTime plot
                p[0].cla()
                p[0].plot(dataTS,data_blSub,color='blue',label='data',alpha=0.7)
                p[0].plot(dataTS,data_lPass,color='magenta',label='lowpass',linewidth=4)
                p[0].plot(dataTS,data_bPass,color='red',label='bpass',linewidth=4)
                p[0].axvline(bandTime[iH],color='orange',label='bandTime',linewidth=4)
                p[0].legend(loc='best')
                p[0].set_xlabel('Time (ns)')
                p[0].set_ylabel('ADC (arb)')
                p[0].set_title("Run %d  Entry %d  Channel %d  ENFCal %.2f" % (run,iEvent,chan,dataENFCal))

            if plotNum==6: # waveform fit plot
                p[0].cla()
                p[0].plot(dataTS,data,color='blue',label='data')
                # p[0].plot(dataTS,data_wlDenoised,color='cyan',label='wlDenoised',alpha=0.5)
                p[0].plot(dataTS,temp,color='orange',label='xgauss guess')
                p[0].plot(dataTS,fit,color='red',label='xgauss fit')
                p[0].set_title("Run %d  evt %d  chan %d  trapENFCal %.1f  trapENM %.1f  deltaBL %.1f\n  amp %.2f  mu %.2f  sig %.2f  tau %.2f  chi2 %.2f  spd %.3f" % (run,iList,chan,dataENFCal,dataENM,dataBL-bl,amp,mu,sig,tau,fitChi2[iH],fitSpeed))
                p[0].legend(loc='best')
                p[1].cla()
                p[1].plot(dataTS,data-fit,color='blue',label='residual')
                p[1].legend(loc='best')
                p[2].cla()
                p[2].plot(ampTr[1:],label='amp',color='red')
                p[2].legend(loc='best')
                p[3].cla()
                p[3].plot(muTr[1:],label='mu',color='green')
                p[3].legend(loc='best')
                p[4].cla()
                p[4].plot(sigTr[1:],label='sig',color='blue')
                p[4].yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))
                p[4].legend(loc='best')
                p[5].cla()
                p[5].plot(tauTr[1:],label='tau',color='black')
                p[5].legend(loc='best')
                p[6].cla()
                p[6].plot(blTr[1:],label='bl',color='magenta')
                p[6].legend(loc='best')

            if plotNum==7: # new traps plot
                p[0].cla()
                p[0].plot(dataTS, data_blSub, color='blue', label='data')
                p[0].plot(sTrapTS, sTrap, color='red', label='sTrap')
                p[0].axvline(t0_SLE[iH], color='red')
                p[0].plot(aTrapTS, aTrap, color='orange', label='aTrap')
                p[0].axvline(t0_ALE[iH], color='orange')
                p[0].plot(eTrapTS, eTrap, color='green', label='eTrap')
                p[0].axhline(lat[iH],color='green')
                p[0].plot(pTrapTS, pTrap, color='magenta', label='pTrap')
                p[0].axhline(latAFC[iH], color='magenta')
                p[0].set_title("trapENFCal %.2f  trapENM %.2f || latEM %.2f  latEF %.2f  latEAF %.2f  latEFC %.2f  latEAFC %.2f" % (dataENFCal,dataENM,lat[iH],latF[iH],latAF[iH],latFC[iH],latAFC[iH]))
                p[0].legend(loc='best')

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
    if batMode and not intMode:
        oTree.Write("",TObject.kOverwrite)
        print "Wrote",oTree.GetBranch("channel").GetEntries(),"entries in the copied tree,"
        print "and wrote",b1.GetEntries(),"entries in the new branches."

    stopT = time.clock()
    print "Stopped:",time.strftime('%X %x %Z'),"\nProcess time (min):",(stopT - startT)/60
    print float(nList)/((stopT-startT)/60.),"entries per minute."


def evalXGaus(x,mu,sig,tau):
    """ Ported from GAT/BaseClasses/GATPeakShapeUtil.cc
        Negative tau: Regular WF, high tail
        Positive tau: Backwards WF, low tail
    """
    tmp = (x-mu + sig**2./2./tau)/tau
    if all(tmp < limit):
        return np.exp(tmp)/2./np.fabs(tau) * sp.erfc((tau*(x-mu)/sig + sig)/np.sqrt(2.)/np.fabs(tau))
    else:
        print "Exceeded limit ..."
        # Here, exp returns NaN (in C++).  So use an approx. derived from the asymptotic expansion for erfc, listed on wikipedia.
        den = 1./(sig + tau*(x-mu)/sig)
        return sig * gaus(x,mu,sig) * den * (1.-tau**2. * den**2.)

def xgModelWF(dataTS, floats):
    """ Make a model waveform: Take a timestamp vector, generate an
        xGauss model, normalize to 1, then scale its max value to amp.
    """
    amp, mu, sig, tau, bl = floats
    model = evalXGaus(dataTS,mu,sig,tau)

    # pin max value of function to amp
    model = model * 1./np.sum(model)
    xMax = np.argmax(model)
    model = model * (amp / model[xMax])

    # float the baseline
    model = model + bl

    return model

def MakeTracesGlobal():
    """ This is so 'lnLike' can write to the trace arrays. Has to remain in this file to work. """
    tmp1 = tmp2 = tmp3 = tmp4 = tmp5 = np.empty([1,])
    global ampTr, muTr, sigTr, tauTr, blTr
    ampTr, muTr, sigTr, tauTr, blTr = tmp1, tmp2, tmp3, tmp4, tmp5

def lnLike(floats, *datas):
    """ log-likelihood function: L(A,mu,sig,tau)
    To make this work with op.minimize, 'datas' is passed in as a tuple (the asterisk),
    where the original list is the 1st element.
    """
    # Add to traces.
    global ampTr, muTr, sigTr, tauTr, blTr
    amp, mu, sig, tau, bl = floats
    ampTr = np.append(ampTr, amp)
    muTr = np.append(muTr, mu)
    sigTr = np.append(sigTr, sig)
    tauTr = np.append(tauTr, tau)
    blTr = np.append(blTr, bl)

    dataTS, data, dataNoise = datas[0][0], datas[0][1], datas[0][2]

    model = xgModelWF(dataTS, floats)
    lnLike = 0.5 * np.sum ( np.power((data-model)/dataNoise, 2) - np.log( 1 / np.power(dataNoise,2) ) )
    return lnLike


if __name__ == "__main__":
    main(sys.argv[1:])