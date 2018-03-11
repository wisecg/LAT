#!/usr/bin/env python3
"""
============== lat.py: (L)ANL (A)nalysis (T)oolkit ==============

Secondary waveform processing for the Majorana Demonstrator.

Takes a single skim file (augmented with an MGTWaveform branch),
or built/gatified data.  Calculates various waveform parameters
for HITS passing cuts.

Does not handle TChains - pyroot can't always recognize
the vector<MGTWaveform*> branch in the skim tree when a TEntryList
is in use.  Damn you, ROOT.

Usage:
./lat.py [-r [dsNum] [subNum] use skim range & sub-range]
         [-p [inFile] [outFile] path mode: manual file locations]
         [-d [inPath] [outPath] set input and output directories]
         [-f [dsNum] [runNum] use single skim file]
         [-g [dsNum] [runNum] use single gat/blt file]
         [-s [fileName] uses a custom file, use full file path]
         [-i [plotNum] interactive mode]
         [-c "custom cut" -- adds custom cut application]
         [-b batch mode -- creates new file]

v1: 27 May 2017
v2: 04 Aug 2017 - improvements to wf fitting, handle multisampling, etc.
v3: 18 Jan 2018 - update to python3
v4: 09 Mar 2018 - wf fitting error handling (scipy v1.0 improves convergence!)

================ C. Wiseman (USC), B. Zhu (LANL) ================
"""
import sys, time, os, pywt, imp
from ROOT import TFile, TTree, TEntryList, gDirectory, TNamed, std, TObject, gROOT
from ROOT import GATDataSet, MGTEvent, MGTWaveform, MGWFTimePointCalculator
import numpy as np
from scipy.signal import butter, lfilter, filtfilt
import scipy.optimize as op
from scipy.ndimage.filters import gaussian_filter
from scipy import interpolate
import scipy.special as sp
# import waveLibs as wl
wl = imp.load_source('waveLibs',os.environ['LATDIR']+'/waveLibs.py')

def main(argv):

    print("=======================================")
    print("LAT started:",time.strftime('%X %x %Z'))
    startT = time.clock()
    # gROOT.ProcessLine("gErrorIgnoreLevel = 3001;") # suppress ROOT error messages
    global batMode
    intMode, batMode, rangeMode, fileMode, gatMode, singleMode, pathMode, cutMode = False, False, False, False, False, False, False, False
    dontUseTCuts = False
    dsNum, subNum, runNum, plotNum = -1, -1, -1, 1
    pathToInput, pathToOutput, manualInput, manualOutput, customPar = ".", ".", "", "", ""

    if len(argv)==0: return
    for i,opt in enumerate(argv):
        if opt == "-r":
            rangeMode, dsNum, subNum = True, int(argv[i+1]), int(argv[i+2])
            print("Scanning DS-%d sub-range %d" % (dsNum, subNum))
        if opt == "-p":
            pathMode, manualInput, manualOutput = True, argv[i+1], argv[i+2]
            print("Manually set input/output files:\nInput: %s\nOutput: %s" % (manualInput, manualOutput))
        if opt == "-d":
            pathToInput, pathToOutput = argv[i+1], argv[i+2]
            print("Custom paths: Input %s,  Output %s" % (pathToInput,pathToOutput))
        if opt == "-f":
            fileMode, dsNum, runNum = True, int(argv[i+1]), int(argv[i+2])
            print("Scanning DS-%d, run %d" % (dsNum, runNum))
        if opt == "-g":
            gatMode, runNum = True, int(argv[i+1])
            print("GATDataSet mode.  Scanning run %d" % (runNum))
        if opt == "-s":
            singleMode, pathToInput = True, argv[i+1]
            print("Single file mode.  Scanning {}".format(pathToInput))
        if opt == "-i":
            intMode, plotNum = True, int(argv[i+1])
            print("Interactive mode selected. Use \"p\" for previous and \"q\" to exit.")
        if opt == "-x":
            dontUseTCuts = True
            print("DC TCuts deactivated.  Retaining all events ...")
        if opt == "-c":
            cutMode, customPar = True, str(argv[i+1])
            print("Using custom cut parameter: {}".format(customPar))
        if opt == "-b":
            batMode = True
            import matplotlib
            # if os.environ.get('DISPLAY','') == '':
                # print('No display found. Using non-interactive Agg backend')
            matplotlib.use('Agg')
            print("Batch mode selected.  A new file will be created.")
    import matplotlib.pyplot as plt
    from matplotlib import gridspec
    import matplotlib.ticker as mtick

    # File I/O
    inFile, outFile, bltFile = TFile(), TFile(), TFile()
    gatTree, bltTree, out = TTree(), TTree(), TTree()
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
        print(gatTree.GetEntries(),"entries in input tree.")
    elif gatMode:
        inFile = TFile(gatPath)
        bltFile = TFile(bltPath)
        gatTree = inFile.Get("mjdTree")
        bltTree = bltFile.Get("MGTree")
        gatTree.AddFriend(bltTree)
    if singleMode:
        inFile = TFile(pathToInput)
        gatTree = inFile.Get("skimTree")

    # apply cut to tree
    if (rangeMode or fileMode or pathMode) and not dontUseTCuts:
        try:
            theCut = inFile.Get("theCut").GetTitle()
        except ReferenceError:
            theCut = ""
    if cutMode:
        # theCut += customPar
        # theCut = "(channel==672 || channel==674) && mH==2" # sync chan: 672, extp chan: 674
        theCut = "channel==674 && mH==2"
        theCut += " && fitSlo < 10"
        print("WARNING: Custom cut in use! : ",theCut)

    gatTree.Draw(">>elist", theCut, "entrylist")
    elist = gDirectory.Get("elist")
    gatTree.SetEntryList(elist)
    nList = elist.GetN()
    print("Using cut:\n",theCut)
    print("Found",gatTree.GetEntries(),"input entries.")
    print("Found",nList,"entries passing cuts.")

    # Output: In batch mode (-b) only, create an output file+tree & append new branches.
    if batMode and not intMode:
        outFile = TFile(outPath, "RECREATE")
        print("Attempting tree copy to",outPath)
        out = gatTree.CopyTree("")
        out.Write()
        print("Wrote",out.GetEntries(),"entries.")
        cutUsed = TNamed("theCut",theCut)
        cutUsed.Write()

    waveS1, waveS2 = std.vector("double")(), std.vector("double")()
    waveS3, waveS4, waveS5 = std.vector("double")(), std.vector("double")(), std.vector("double")()
    bcMax, bcMin = std.vector("double")(), std.vector("double")()
    bandMax, bandTime = std.vector("double")(), std.vector("double")()
    den10, den50, den90 = std.vector("double")(), std.vector("double")(), std.vector("double")()
    oppie = std.vector("double")()
    fitMu, fitAmp, fitSlo = std.vector("double")(), std.vector("double")(), std.vector("double")()
    fitTau, fitBL = std.vector("double")(), std.vector("double")()
    matchMax, matchWidth, matchTime = std.vector("double")(), std.vector("double")(), std.vector("double")()
    pol0, pol1, pol2, pol3 = std.vector("double")(), std.vector("double")(), std.vector("double")(), std.vector("double")()
    fails, fitChi2, fitLL = std.vector("int")(), std.vector("double")(), std.vector("double")()
    riseNoise = std.vector("double")()
    t0_SLE, t0_ALE, lat, latF = std.vector("double")(), std.vector("double")(), std.vector("double")(), std.vector("double")()
    latAF, latFC, latAFC = std.vector("double")(), std.vector("double")(), std.vector("double")()
    nMS = std.vector("int")()
    tE50, latE50, wfStd = std.vector("double")(), std.vector("double")(), std.vector("double")()
    wfAvgBL, wfRMSBL = std.vector("double")(), std.vector("double")()
    fitErr = std.vector("int")()

    # It's not possible to put the "out.Branch" call into a class initializer (waveLibs::latBranch). You suck, ROOT.
    b1, b2 = out.Branch("waveS1",waveS1), out.Branch("waveS2",waveS2)
    b3, b4, b5 = out.Branch("waveS3",waveS3), out.Branch("waveS4",waveS4), out.Branch("waveS5",waveS5)
    b7, b8 = out.Branch("bcMax",bcMax), out.Branch("bcMin",bcMin)
    b9, b10 = out.Branch("bandMax",bandMax), out.Branch("bandTime",bandTime)
    b11, b12, b13 = out.Branch("den10",den10), out.Branch("den50",den50), out.Branch("den90",den90)
    b14 = out.Branch("oppie",oppie)
    b15, b16, b17 = out.Branch("fitMu", fitMu), out.Branch("fitAmp", fitAmp), out.Branch("fitSlo", fitSlo)
    b18, b19 = out.Branch("fitTau",fitTau), out.Branch("fitBL",fitBL)
    b20, b21, b22 = out.Branch("matchMax", matchMax), out.Branch("matchWidth", matchWidth), out.Branch("matchTime", matchTime)
    b23, b24, b25, b26 = out.Branch("pol0", pol0), out.Branch("pol1", pol1), out.Branch("pol2", pol2), out.Branch("pol3", pol3)
    b27, b28, b29 = out.Branch("fails",fails), out.Branch("fitChi2",fitChi2), out.Branch("fitLL",fitLL)
    b30 = out.Branch("riseNoise",riseNoise)
    b31, b32, b33, b34 = out.Branch("t0_SLE",t0_SLE), out.Branch("t0_ALE",t0_ALE), out.Branch("lat",lat), out.Branch("latF",latF)
    b35, b36, b37 = out.Branch("latAF",latAF), out.Branch("latFC",latFC), out.Branch("latAFC",latAFC)
    b38 = out.Branch("nMS",nMS)
    b39, b40, b41 = out.Branch("tE50", tE50), out.Branch("latE50", latE50), out.Branch("wfStd", wfStd)
    b42, b43 = out.Branch("wfAvgBL", wfAvgBL), out.Branch("wfRMSBL", wfRMSBL)
    b44 = out.Branch("fitErr",fitErr)

    # make a dictionary that can be iterated over (avoids code repetition in the loop)
    brDict = {
        "waveS1":[waveS1, b1], "waveS2":[waveS2, b2],
        "waveS3":[waveS3, b3], "waveS4":[waveS4, b4], "waveS5":[waveS5, b5],
        "bcMax":[bcMax, b7], "bcMin":[bcMin, b8],
        "bandMax":[bandMax, b9], "bandTime":[bandTime, b10],
        "den10":[den10, b11], "den50":[den50, b12], "den90":[den90, b13],
        "oppie":[oppie, b14],
        "fitMu":[fitMu, b15], "fitAmp":[fitAmp, b16], "fitSlo":[fitSlo, b17],
        "fitTau":[fitTau, b18], "fitBL":[fitBL,b19],
        "matchMax":[matchMax, b20], "matchWidth":[matchWidth, b21], "matchTime":[matchTime, b22],
        "pol0":[pol0, b23], "pol1":[pol1, b24], "pol2":[pol2, b25], "pol3":[pol3, b26],
        "fails":[fails,b27], "fitChi2":[fitChi2,b28], "fitLL":[fitLL,b29],
        "riseNoise":[riseNoise,b30],
        "t0_SLE":[t0_SLE,b31], "t0_ALE":[t0_ALE,b32], "lat":[lat,b33], "latF":[latF,b34],
        "latAF":[latAF,b35], "latFC":[latFC,b36], "latAFC":[latAFC,b37],
        "nMS":[nMS,b38], "tE50":[tE50,b39], "latE50":[latE50,b40], "wfStd":[wfStd,b41],
        "wfAvgBL":[wfAvgBL,b42], "wfRMSBL":[wfRMSBL,b43],
        "fitErr":[fitErr,b44]
    }

    # Make a figure (-i option: select different plots)
    fig = plt.figure(figsize=(12,7), facecolor='w')
    if plotNum==0 or plotNum==7 or plotNum==8:
        p0 = plt.subplot(111)  # 0-raw waveform, 7-new trap filters
    elif plotNum==1 or plotNum==2:
        p0 = plt.subplot(211)  # 1-wavelet, 2-time points, bandpass filters, tail slope
        p1 = plt.subplot(212)
    elif plotNum==3:
        p0 = plt.subplot2grid((2,5), (0,0), colspan=3)  # oppie / freq-domain matched filter
        p1 = plt.subplot2grid((2,5), (0,3), colspan=2)
        p2 = plt.subplot2grid((2,5), (1,0), colspan=3)
    elif plotNum==4:
        p0 = plt.subplot(111)  # time-domain matched filter
    elif plotNum==5:
        p0 = plt.subplot(111)  # bandpass / bandTime
    elif plotNum==6:
        p0 = plt.subplot2grid((6,10), (0,0), colspan=10, rowspan=3) # waveform fit
        p1 = plt.subplot2grid((6,10), (3,0), colspan=10, rowspan=1) # residual
        p2 = plt.subplot2grid((6,10), (4,0), colspan=2, rowspan=2) # traces
        p3 = plt.subplot2grid((6,10), (4,2), colspan=2, rowspan=2)
        p4 = plt.subplot2grid((6,10), (4,4), colspan=2, rowspan=2)
        p5 = plt.subplot2grid((6,10), (4,6), colspan=2, rowspan=2)
        p6 = plt.subplot2grid((6,10), (4,8), colspan=2, rowspan=2)
    elif plotNum==9:
        p0 = plt.subplot2grid((5,1), (0,0)) # 9- wpt on wf fit residual
        p1 = plt.subplot2grid((5,1), (1,0), rowspan=2)
        p2 = plt.subplot2grid((5,1), (3,0), rowspan=2)
    if not batMode: plt.show(block=False)


    # Load a fast signal template - used w/ the freq-domain matched filter
    # print("Generating signal template ...")
    tSamp, tR, tZ, tAmp, tST, tSlo = 5000, 0, 15, 100, 2500, 10
    # tOrig, tOrigTS = wl.MakeSiggenWaveform(tSamp,tR,tZ,tAmp,tST,tSlo) # Damn you to hell, PDSF
    templateFile = np.load("%s/data/lat_template.npz" % os.environ['LATDIR'])
    if dsNum==2 or dsNum==6: templateFile = np.load("%s/data/lat_ds2template.npz" % os.environ['LATDIR'])
    tOrig, tOrigTS = templateFile['arr_0'], templateFile['arr_1']


    # Load stuff from DS1 forced acq. runs
    npzfile = np.load("%s/data/fft_forcedAcqDS1.npz" % os.environ['LATDIR'])
    noise_asd, noise_xFreq, avgPwrSpec, xPwrSpec, data_forceAcq, data_fft = npzfile['arr_0'],npzfile['arr_1'],npzfile['arr_2'],npzfile['arr_3'],npzfile['arr_4'],npzfile['arr_5']


    # Loop over events
    print("Starting event loop ...")
    iList = -1
    while True:
        iList += 1
        if intMode==True and iList != 0:
            value = input()
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
        chanList = list(set(int(chans[n]) for n in range(numPass)))
        hitList = (iH for iH in range(nChans) if gatTree.channel.at(iH) in chanList)  # a 'generator expression'
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

            # print("%d:  run %d  chan %d  trapENFCal %.2f" % (iList, run, chan, dataENFCal))

            # be absolutely sure you're matching the right waveform to this hit
            if wf.GetID() != chan:
                print("ERROR -- Vector matching failed.  iList %d  run %d  iEvent %d" % (iList,run,iEvent))
                return

            # Let's start the show - grab a waveform.
            # Remove first 4 samples when we have multisampling
            # Remove last 2 samples to get rid of the ADC spike at the end of all wf's.
            truncLo, truncHi = 0, 2
            if dsNum==6 or dsNum==2: truncLo = 4
            signal = wl.processWaveform(wf,truncLo,truncHi)
            data = signal.GetWaveRaw()
            data_blSub = signal.GetWaveBLSub()
            dataTS = signal.GetTS()
            dataBL,dataNoise = signal.GetBaseNoise()

            # wavelet packet transform
            wp = pywt.WaveletPacket(data_blSub, 'db2', 'symmetric', maxlevel=4)
            nodes = wp.get_level(4, order='freq')
            wpCoeff = np.array([n.data for n in nodes],'d')
            wpCoeff = abs(wpCoeff)

            # wavelet parameters
            # First get length of wavelet on the time axis, the scale axis will always be the same
            # due to the number of levels in the wavelet
            wpLength = len(wpCoeff[1,:])
            waveS1[iH] = np.sum(wpCoeff[0:1,1:wpLength//4+1]) # python3 : floor division (//) returns an int
            waveS2[iH] = np.sum(wpCoeff[0:1,wpLength//4+1:wpLength//2+1])
            waveS3[iH] = np.sum(wpCoeff[0:1,wpLength//2+1:3*wpLength//4+1])
            waveS4[iH] = np.sum(wpCoeff[0:1,3*wpLength//4+1:-1])
            waveS5[iH] = np.sum(wpCoeff[2:-1,1:-1])
            S6 = np.sum(wpCoeff[2:9,1:wpLength//4+1])
            S7 = np.sum(wpCoeff[2:9,wpLength//4+1:wpLength//2+1])
            S8 = np.sum(wpCoeff[2:9,wpLength//2+1:3*wpLength//4+1])
            S9 = np.sum(wpCoeff[2:9,3*wpLength//4+1:-1])
            S10 = np.sum(wpCoeff[9:,1:wpLength//4+1])
            S11 = np.sum(wpCoeff[9:,wpLength//4+1:wpLength//2+1])
            S12 = np.sum(wpCoeff[9:,wpLength//2+1:3*wpLength//4+1])
            S13 = np.sum(wpCoeff[9:,3*wpLength//4+1:-1])
            sumList = [S6, S7, S8, S9, S10, S11, S12, S13]
            bcMax[iH] = np.max(sumList)
            bcMin[iH] = 1. if np.min(sumList) < 1 else np.min(sumList)

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

            # used in the multisite tagger
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
            if not batMode: MakeTracesGlobal()

            # get the noise of the denoised wf
            denoisedNoise,_,_ = wl.baselineParameters(data_wlDenoised)

            # NOTE: fit is to wavelet-denoised data, BECAUSE there are no HF components in the model,
            # AND we'll still calculate fitChi2 w/r/t the data, not the denoised data.
            # datas = [dataTS, data, dataNoise] # fit data
            datas = [dataTS, data_wlDenoised + dataBL, denoisedNoise] # fit wavelet-denoised data w/ Bl added back in

            # Set bounds - A,mu,sig,tau,bl.
            # bnd = ((None,None),(None,None),(None,None),(None,None),(None,None))   # often gets caught at sig=0
            bnd = ((None,None),(None,None),(2.,None),(-72001.,-71999.),(None,None)) # gets caught much less often.

            # L-BGFS-B with numerical gradient.
            start = time.clock()
            result = op.minimize(lnLike, floats, args=datas, method="L-BFGS-B", options=None, bounds=bnd)
            fitSpeed = time.clock() - start

            fitErr[iH] = 0
            if not result["success"]:
                # print("fit fail: ", result["message"])
                fitErr[iH] = 1
                errorCode[0] = 1

            amp, mu, sig, tau, bl = result["x"]

            # save parameters

            fitMu[iH], fitAmp[iH], fitSlo[iH], fitTau[iH], fitBL[iH] = mu, amp, sig, tau, bl
            floats = np.asarray([amp, mu, sig, tau, bl])
            fit = xgModelWF(dataTS, floats)

            # print("%d/%d iH %d  e %-10.2f  fs %-8.2f  f %d" % (iList, nList, iH, dataENFCal, fitSlo[iH], fitErr[iH]))

            # log-likelihood of this fit
            fitLL[iH] = result["fun"]

            # chi-square of this fit
            # Textbook is (observed - expected)^2 / expected,
            # but we'll follow MGWFCalculateChiSquare.cc and do (observed - expected)^2 / NDF.
            # NOTE: we're doing the chi2 against the DATA, though the FIT is to the DENOISED DATA.
            fitChi2[iH] = np.sum(np.square(data-fit)) / (len(data)-1)/dataNoise


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
            riseNoise[iH] = np.sum(wpCoeff[2:-1,wpLoRise:wpHiRise]) / bcMin[iH]

            # print("%d %d %d %d e %-5.2f  bmax %-6.2f  bmin %-6.2f  mu %-5.2f  a %-5.2f  s %-5.2f  bl %-5.2f  rn %.2f" % (run,iList,iH,chan,dataENFCal,bcMax[iH],bcMin[iH],fitMu[iH],fitAmp[iH],fitSlo[iH],fitBL[iH],riseNoise[iH]))

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
            # this kills the 1-1 matching between matchTS and dataTS (each TS has some offset)
            matchMaxTime = matchTS[np.argmax(match)]
            matchTS = matchTS + (fitMaxTime - matchMaxTime)

            # resize match, matchTS to have same # samples as data, dataTS.
            # this is the only case we really care about
            # ("too early" and "too late" also happen, but the shift is larger than the trigger walk, making it unphysical)
            if matchTS[0] <= dataTS[0] and matchTS[-1] >= dataTS[-1]:
                idx = np.where((matchTS >= dataTS[0]) & (matchTS <= dataTS[-1]))
                match, matchTS = match[idx], matchTS[idx]
                sizeDiff = len(dataTS)-len(matchTS)
                if sizeDiff < 0:
                    match, matchTS = match[:sizeDiff], matchTS[:sizeDiff]
                elif sizeDiff > 0:
                    match = np.hstack((match, np.zeros(sizeDiff)))
                    matchTS = np.hstack((matchTS, dataTS[-1*sizeDiff:]))
                if len(match) != len(data):
                    print("FIXME: match filter array manip is still broken.")

            # compute match filter parameters
            matchMax[iH], matchWidth[iH], matchTime[iH] = -888, -888, -888
            if len(match)==len(data):
                smoothMF = gaussian_filter(match * data_blSub, sigma=5.)
                matchMax[iH] = np.amax(smoothMF)
                matchTime[iH] = matchTS[ np.argmax(smoothMF) ]
                idx = np.where(smoothMF > matchMax[iH]/2.)
                if len(matchTS[idx]>1):
                    matchWidth[iH] = matchTS[idx][-1] - matchTS[idx][0]


            # Fit tail slope to polynomial.  Guard against fit fails

            idx = np.where(dataTS >= fitMaxTime)
            tail, tailTS = data[idx], dataTS[idx]
            popt1,popt2 = 0,0
            try:
                popt1,_ = op.curve_fit(wl.tailModelPol, tailTS, tail)
                pol0[iH], pol1[iH], pol2[iH], pol3[iH] = popt1[0], popt1[1], popt1[2], popt1[3]
            except:
                # print("curve_fit tailModelPol failed, run %i  event %i  channel %i" % (run, iList, chan))
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
            aTrap = wl.asymTrapFilter(data_blSub, 4, 10, 200, True) # (0.04us, 0.1us, 2.0us)
            aTrapTS = np.arange(0, len(aTrap)*10., 10)

            # find leading edges (t0 times)

            # limit the range from 0 to 10us, and use an ADC threshold of 1.0 as suggested by DCR
            t0_SLE[iH],_ = wl.walkBackT0(sTrap, eTrapTS[-1]+7000-4000-2000, 1., 0, 1000) # (in ns) finds leading edge from short trap
            t0_ALE[iH],_ = wl.walkBackT0(aTrap, eTrapTS[-1]+7000-4000-2000, 1., 0, 1000) # (in ns) finds leading edge from asymmetric trap

            # standard energy trapezoid w/ a baseline padded waveform
            data_pad = np.pad(data_blSub,(200,0),'symmetric')
            pTrap = wl.trapFilter(data_pad, 400, 250, 7200.)
            pTrapTS = np.linspace(0, len(pTrap)*10, len(pTrap))
            pTrapInterp = interpolate.interp1d(pTrapTS, pTrap)

            # calculate energy parameters
            # standard amplitude.  basically trapEM, but w/o NL correction if the input WF doesn't have it.
            lat[iH] = np.amax(eTrap)

            # Calculate DCR suggested amplitude, using the 50% to the left and right of the maximum point
            t0_F50,t0fail1 = wl.walkBackT0(pTrap, thresh=lat[iH]*0.5, rmin=0, rmax=len(pTrap)-1)
            t0_B50,t0fail2 = wl.walkBackT0(pTrap, thresh=lat[iH]*0.5, rmin=0, rmax=len(pTrap)-1, forward=True)
            t0_E50 = (t0_F50 + t0_B50)/2.0

            #TODO -- if it's necessary due to the trigger walk, we could potentially add a way to recursively increase the threshold until a timepoint is found, however it will still always fail for most noise events
            if not t0fail1 or not t0fail2:
                latE50[iH] = 0 # Set amplitude to 0 if one of the evaluations failed
            else:
                latE50[iH] = pTrapInterp(t0_E50) # Maybe I should call this latDCR50 to confuse people
            tE50[iH] = t0_B50 - t0_F50 # Save the difference between the middle points, can be used as a cut later

            # standard amplitude with t0 from the shorter traps
            # If either fixed pickoff time (t0) is < 0, use the first sample as the amplitude (energy).
            latF[iH] = eTrapInterp( np.amax([t0_SLE[iH]-7000+4000+2000, 0.]) ) # This should be ~trapEF
            latAF[iH] = eTrapInterp( np.amax([t0_ALE[iH]-7000+4000+2000, 0.]) )

            # amplitude from padded trapezoid, with t0 from short traps and a correction function
            # function is under development.  currently: f() = exp(p0 + p1*E), p0 ~ 7.8, p1 ~ -0.45 and -0.66
            # functional walk back distance is *either* the minimum of the function value, or 5500 (standard value)

            # t0_corr = -7000+6000+2000 # no correction
            t0_corr = -7000+6000+2000 - np.amin([np.exp(7.8 - 0.45*lat[iH]),1000.])
            t0A_corr = -7000+6000+2000 - np.amin([np.exp(7.8 - 0.66*lat[iH]),1000.])

            latFC[iH] = pTrapInterp( np.amax([t0_SLE[iH] + t0_corr, 0.]) )
            latAFC[iH] = pTrapInterp( np.amax([t0_ALE[iH] + t0A_corr, 0.]) )


            # =========================================================

            # the genius multisite event tagger - plot 8

            # decide a threshold
            dIdx = np.argmax(data_filtDeriv)
            dMax = data_filtDeriv[dIdx]
            dRMS,_,_ = wl.baselineParameters(data_filtDeriv)
            # msThresh = np.amax([dMax * .2, dRMS * 5.])
            # msThresh = dMax * .15
            msThresh = 50.  # I don't know.  this seems like a good value

            # run peak detect algorithm
            maxtab,_ = wl.peakdet(data_filtDeriv, msThresh)

            # profit
            msList = []
            for iMax in range(len(maxtab)):
                idx = int(maxtab[iMax][0])
                val = maxtab[iMax][1]
                msList.append(dataTS[idx])
                # print("%d  idx %d  TS %d  val %.2f  thresh %.2f" % (iList, idx, dataTS[idx], val, msThresh))
            nMS[iH] = len(maxtab)

            # =========================================================
            # wfStd analysis
            wfAvgBL[iH] = dataBL
            wfRMSBL[iH] = dataNoise
            wfStd[iH] = np.std(data[5:-5])

            # ------------------------------------------------------------------------
            # End waveform processing.

            # Calculate error code
            fails[iH] = 0
            for i,j in enumerate(errorCode):
                if j==1: fails[iH] += int(j)<<i
            # print("fails:",fails[iH])

            # Make plots!
            if batMode: continue
            if plotNum==0: # raw data
                p0.cla()
                p0.plot(dataTS,data,'b')
                p0.set_title("Run %d  Entry %d  Channel %d  ENFCal %.2f" % (run,iList,chan,dataENFCal))

            if plotNum==1: # wavelet plot
                p0.cla()
                p0.margins(x=0)
                p0.plot(dataTS,data_blSub,color='blue',label='data')
                p0.plot(dataTS,data_wlDenoised,color='cyan',label='denoised',alpha=0.7)
                p0.axvline(fitRiseTime50,color='green',label='fit 50%',linewidth=2)
                p0.plot(dataTS,fit_blSub,color='red',label='bestfit',linewidth=2)
                p0.set_title("Run %d  Entry %d  Channel %d  ENFCal %.2f  flo %.0f  fhi %.0f  fhi-flo %.0f" % (run,iList,chan,dataENFCal,fitStartTime,fitMaxTime,fitMaxTime-fitStartTime))
                p0.legend(loc='best')

                p1.cla()
                p1.imshow(wpCoeff, interpolation='nearest', aspect="auto", origin="lower",extent=[0, 1, 0, len(wpCoeff)],cmap='viridis')
                p1.axvline(float(wpLoRise)/numXRows,color='orange',linewidth=2)
                p1.axvline(float(wpHiRise)/numXRows,color='orange',linewidth=2)
                p1.set_title("waveS5 %.2f  bcMax %.2f  bcMin %.2f  riseNoise %.2f" % (waveS5[iH], bcMax[iH], bcMin[iH], riseNoise[iH]))

            if plotNum==2: # time points, bandpass filters, tail slope
                p0.cla()
                p0.plot(dataTS,data,color='blue',label='data')
                p0.axvline(den10[iH],color='black',label='lpTP')
                p0.axvline(den50[iH],color='black')
                p0.axvline(den90[iH],color='black')
                p0.plot(dataTS,fit,color='magenta',label='bestfit')
                if errorCode[2]!=1: p0.plot(tailTS, wl.tailModelPol(tailTS, *popt1), color='orange',linewidth=2, label='tailPol')
                p0.legend(loc='best')
                p0.set_title("Run %d  Entry %d  Channel %d  ENFCal %.2f" % (run,iEvent,chan,dataENFCal))

                p1.cla()
                p1.plot(dataTS,data_lPass,color='blue',label='lowpass')
                p1.plot(dataTS,data_filtDeriv,color='green',label='filtDeriv')
                p1.plot(dataTS,data_filt,color='black',label='filtfilt')
                p1.plot(dataTS,data_bPass,color='red',label='bpass')
                p1.axvline(bandTime[iH],color='orange',label='bandTime')
                p1.legend(loc='best')

            if plotNum==3: # freq-domain matched filter
                p0.cla()
                p0.plot(dataTS,data,'b')
                p0.plot(dataTS,temp,'r')
                p0.plot(dataTS,fit,color='cyan')
                p0.set_title("Run %d  Entry %d  Channel %d  ENFCal %.2f" % (run,iEvent,chan,dataENFCal))

                data_asd, data_xFreq = plt.psd(data, Fs=1e8, NFFT=2048, pad_to=2048, visible=False)
                temp_asd, temp_xFreq = plt.psd(temp, Fs=1e8, NFFT=2048, pad_to=2048, visible=False)

                p1.cla()
                p1.loglog(data_xFreq, np.sqrt(data_asd), 'b')
                p1.loglog(noise_xFreq, np.sqrt(noise_asd), 'g')
                p1.loglog(temp_xFreq, np.sqrt(temp_asd), 'r')
                p1.set_xlabel('Frequency (Hz)')
                p1.set_ylabel('ASD')
                p1.grid('on')

                p2.cla()
                p2.plot(dataTS, SNR)
                p2.set_title('oppie %.1f' % (oppie[iH]))
                p2.set_xlabel('Offset time (s)')
                p2.set_ylabel('SNR')

            if plotNum==4: # time domain match filter plot
                p0.cla()
                p0.plot(dataTS,data_blSub,color='blue',label='data',alpha=0.7)
                p0.plot(dataTS,fit_blSub,color='red',label='bestfit',linewidth=3)
                p0.axvline(matchTime[iH],color='orange',label='matchTime',linewidth=2)
                p0.plot(matchTS,smoothMF,color='magenta',label='smoothMF',linewidth=3)
                p0.plot(matchTS,match,color='cyan',label='match',linewidth=3)
                p0.set_xlabel('Time (s)')
                p0.set_ylabel('Voltage (arb)')
                p0.legend(loc='best')
                p0.set_title("Run %d  Entry %d  Channel %d  ENFCal %.2f  matchMax %.2f  matchTime %.2f  matchWidth %.2f" % (run,iEvent,chan,dataENFCal,matchMax[iH],matchTime[iH],matchWidth[iH]))

            if plotNum==5: # bandTime plot
                p0.cla()
                p0.plot(dataTS,data_blSub,color='blue',label='data',alpha=0.7)
                p0.plot(dataTS,data_lPass,color='magenta',label='lowpass',linewidth=4)
                p0.plot(dataTS,data_bPass,color='red',label='bpass',linewidth=4)
                p0.axvline(bandTime[iH],color='orange',label='bandTime',linewidth=4)
                p0.legend(loc='best')
                p0.set_xlabel('Time (ns)')
                p0.set_ylabel('ADC (arb)')
                p0.set_title("Run %d  Entry %d  Channel %d  ENFCal %.2f" % (run,iEvent,chan,dataENFCal))

            if plotNum==6: # waveform fit plot
                p0.cla()
                p0.plot(dataTS,data,color='blue',label='data')
                # p0.plot(dataTS,data_wlDenoised,color='cyan',label='wlDenoised',alpha=0.5)
                p0.plot(dataTS,temp,color='orange',label='xgauss guess')
                p0.plot(dataTS,fit,color='red',label='xgauss fit')
                p0.set_title("Run %d  evt %d  chan %d  trapENFCal %.1f  trapENM %.1f  deltaBL %.1f\n  amp %.2f  mu %.2f  sig %.2f  tau %.2f  chi2 %.2f  spd %.3f" % (run,iList,chan,dataENFCal,dataENM,dataBL-bl,amp,mu,sig,tau,fitChi2[iH],fitSpeed))
                p0.legend(loc='best')
                p1.cla()
                p1.plot(dataTS,data-fit,color='blue',label='residual')
                p1.legend(loc='best')
                p2.cla()
                p2.plot(ampTr[1:],label='amp',color='red')
                p2.legend(loc='best')
                p3.cla()
                p3.plot(muTr[1:],label='mu',color='green')
                p3.legend(loc='best')
                p4.cla()
                p4.plot(sigTr[1:],label='sig',color='blue')
                p4.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))
                p4.legend(loc='best')
                p5.cla()
                p5.plot(tauTr[1:],label='tau',color='black')
                p5.legend(loc='best')
                p6.cla()
                p6.plot(blTr[1:],label='bl',color='magenta')
                p6.legend(loc='best')

                print(gatTree.fitSlo.at(iH), sig)

            if plotNum==7: # new traps plot
                p0.cla()
                p0.plot(dataTS, data_blSub, color='blue', label='data')
                p0.plot(sTrapTS, sTrap, color='red', label='sTrap')
                p0.axvline(t0_SLE[iH], color='red')
                p0.plot(aTrapTS, aTrap, color='orange', label='aTrap')
                p0.axvline(t0_ALE[iH], color='orange')
                p0.plot(eTrapTS, eTrap, color='green', label='eTrap')
                p0.axhline(lat[iH],color='green')
                p0.plot(pTrapTS, pTrap, color='magenta', label='pTrap')
                p0.axhline(latAFC[iH], color='magenta')
                p0.axhline(latE50[iH], color='cyan')
                p0.set_title("trapENFCal %.2f  trapENM %.2f || latEM %.2f  latEF %.2f  latEAF %.2f  latEFC %.2f  latEAFC %.2f  latE50 %.2f" % (dataENFCal,dataENM,lat[iH],latF[iH],latAF[iH],latFC[iH],latAFC[iH], latE50[iH]))
                p0.legend(loc='best')

            if plotNum==8: # multisite tag plot
                p0.cla()
                p0.plot(dataTS, data_blSub, color='blue', label='data')
                p0.plot(dataTS, data_filtDeriv, color='red', label='filtDeriv')
                for mse in msList: p0.axvline(mse, color='green')
                p0.axhline(msThresh,color='red')
                p0.legend()

            if plotNum==9: # wavelet vs wf fit residual plot

                # wavelet packet transform on wf fit residual
                fitResid = data-fit
                wpRes = pywt.WaveletPacket(fitResid, 'db2', 'symmetric', maxlevel=4)
                nodesRes = wpRes.get_level(4, order='freq')
                wpCoeffRes = np.array([n.data for n in nodesRes], 'd')
                wpCoeffRes = abs(wpCoeffRes)
                R6 = np.sum(wpCoeffRes[2:9,1:wpLength//4+1])
                R7 = np.sum(wpCoeffRes[2:9,wpLength//4+1:wpLength//2+1])
                R8 = np.sum(wpCoeffRes[2:9,wpLength//2+1:3*wpLength//4+1])
                R9 = np.sum(wpCoeffRes[2:9,3*wpLength//4+1:-1])
                R10 = np.sum(wpCoeffRes[9:,1:wpLength//4+1])
                R11 = np.sum(wpCoeffRes[9:,wpLength//4+1:wpLength//2+1])
                R12 = np.sum(wpCoeffRes[9:,wpLength//2+1:3*wpLength//4+1])
                R13 = np.sum(wpCoeffRes[9:,3*wpLength//4+1:-1])
                RsumList = [R6, R7, R8, R9, R10, R11, R12, R13]
                bcMinRes = 1. if np.min(RsumList) < 1 else np.min(RsumList)
                riseNoiseRes = np.sum(wpCoeffRes[2:-1,wpLoRise:wpHiRise]) / bcMinRes
                rnCut = 1.1762 + 0.00116 * np.log(1 + np.exp((dataENFCal-7.312)/0.341))

                p0.cla()
                p0.margins(x=0)
                p0.plot(dataTS,data_blSub,color='blue',label='data')
                # p0.plot(dataTS,data_wlDenoised,color='cyan',label='denoised',alpha=0.7)
                # p0.axvline(fitRiseTime50,color='green',label='fit 50%',linewidth=2)
                p0.plot(dataTS,fit_blSub,color='red',label='bestfit',linewidth=2)
                # p0.set_title("Run %d  Entry %d  Channel %d  ENFCal %.2f  flo %.0f  fhi %.0f  fhi-flo %.0f" % (run,iList,chan,dataENFCal,fitStartTime,fitMaxTime,fitMaxTime-fitStartTime))
                # p0.legend(loc='best')
                p0.set_title("Run %d  Entry %d  Channel %d  ENFCal %.2f  flo %.0f  fhi %.0f  fhi-flo %.0f  approxFitE %.2f" % (run,iList,chan,dataENFCal,fitStartTime,fitMaxTime,fitMaxTime-fitStartTime,fitAmp[iH]*0.4))

                p1.cla()
                p1.plot(dataTS,fitResid,color='blue')

                p2.cla()
                p2.set_title("riseNoise %.2f  rnCut %.2f  riseNoiseRes %.2f  bcMinRes %.2f  bcMin %.2f  max %.2f" % (riseNoise[iH],rnCut,riseNoiseRes,bcMinRes,bcMin[iH],wpCoeffRes.max()))
                p2.imshow(wpCoeffRes, interpolation='nearest', aspect="auto", origin="lower",extent=[0, 1, 0, len(wpCoeff)],cmap='viridis')

            plt.tight_layout()
            plt.pause(0.000001)
            # ------------------------------------------------------------------------

        # End loop over hits, fill branches
        if batMode:
            for key in brDict:
                brDict[key][1].Fill()
            if iList%5000 == 0 and iList!=0:
                out.Write("",TObject.kOverwrite)
                print("%d / %d entries saved (%.2f %% done), time: %s" % (iList,nList,100*(float(iList)/nList),time.strftime('%X %x %Z')))

    # End loop over events
    if batMode and not intMode:
        out.Write("",TObject.kOverwrite)
        print("Wrote",out.GetBranch("channel").GetEntries(),"entries in the copied tree,")
        print("and wrote",b1.GetEntries(),"entries in the new branches.")

    stopT = time.clock()
    print("Stopped:",time.strftime('%X %x %Z'),"\nProcess time (min):",(stopT - startT)/60)
    print(float(nList)/((stopT-startT)/60.),"entries per minute.")


def evalGaus(x,mu,sig):
    return np.exp(-((x-mu)**2./2./sig**2.))


def evalXGaus(x,mu,sig,tau):
    """ Ported from GAT/BaseClasses/GATPeakShapeUtil.cc
        Negative tau: Regular WF, high tail
        Positive tau: Backwards WF, low tail
    """
    tmp = (x-mu + sig**2./2./tau)/tau

    # np.exp of this is 1.7964120280206387e+308, the largest python float value: sys.float_info.max
    fLimit = 709.782

    if all(tmp < fLimit):
        return np.exp(tmp)/2./np.fabs(tau) * sp.erfc((tau*(x-mu)/sig + sig)/np.sqrt(2.)/np.fabs(tau))
    else:
        # print("Exceeded limit ...")
        # Use an approx. derived from the asymptotic expansion for erfc, listed on wikipedia.
        den = 1./(sig + tau*(x-mu)/sig)
        return sig * evalGaus(x,mu,sig) * den * (1.-tau**2. * den**2.)


def xgModelWF(dataTS, floats):
    """ Make a model waveform: Take a timestamp vector, generate an
        xGauss model, normalize to 1, then scale its max value to amp.
    """
    amp, mu, sig, tau, bl = floats
    model = evalXGaus(dataTS,mu,sig,tau)

    if np.isnan(model).any():
        # print("d'oh!",model)
        return(np.zeros(len(model)))
    if np.sum(model)==0:
        # print("dooh!",model)
        return(np.zeros(len(model)))

    # pin max value of function to amp
    model = model * 1./np.sum(model)
    xMax = np.argmax(model)
    model = model * (amp / model[xMax])

    # float the baseline
    model = model + bl

    return model


def MakeTracesGlobal():
    """ This is so 'lnLike' can write to the trace arrays. Has to remain in this file to work. """
    tmp1, tmp2, tmp3, tmp4, tmp5 = [], [], [], [], []
    global ampTr, muTr, sigTr, tauTr, blTr
    ampTr, muTr, sigTr, tauTr, blTr = tmp1, tmp2, tmp3, tmp4, tmp5


def lnLike(floats, *datas):
    """ log-likelihood function: L(A,mu,sig,tau)
    To make this work with op.minimize, 'datas' is passed in as a tuple (the asterisk),
    where the original list is the 1st element.
    """
    amp, mu, sig, tau, bl = floats
    tau = -72000. # manually fix tau

    # Only fill traces when we're not in batch mode
    global batMode
    if not batMode:
        global ampTr, muTr, sigTr, tauTr, blTr
        ampTr.append(amp)
        muTr.append(mu)
        sigTr.append(sig)
        tauTr.append(tau)
        blTr.append(bl)

    dataTS, data, dataNoise = datas[0][0], datas[0][1], datas[0][2]

    model = xgModelWF(dataTS, [amp, mu, sig, tau, bl])
    lnLike = 0.5 * np.sum ( np.power((data-model)/dataNoise, 2) - np.log( 1 / np.power(dataNoise,2) ) )
    return lnLike


if __name__ == "__main__":
    main(sys.argv[1:])
