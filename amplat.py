#!/usr/local/bin/python
#!/usr/common/usg/software/python/2.7.9/bin/python
#!/bin/bash
""":" `which python` "$0" "$@"; exit 1 ":""" # dammit, this works on PDSF but not my laptop
"""
============== amplat.py: (A)uto-calibration (M)odified (P)eripheral to (L)ANL (A)nalysis (T)oolkit ==============

Stripped down version of LAT for calculating various energy estimators only!

Takes a single skim file (augmented with an MGTWaveform branch),
or built/gatified data.  Calculates various waveform parameters
for HITS passing cuts.

Does not handle TChains - pyroot can't always recognize
the vector<MGTWaveform*> branch in the skim tree. Damn you, ROOT.

Usage:
./amplat.py [-r [dsNum] [subNum] use skim range & sub-range]
         [-p [inFile] [outFile] path mode: manual file locations]
         [-d [inPath] [outPath] set input and output directories]
         [-f [dsNum] [runNum] use single skim file]
         [-g [dsNum] [runNum] use single gat/blt file]
         [-s -- use a custom file, must combine w/ path mode]
         [-b batch mode -- creates new file]

v1: 27 May 2017

================ C. Wiseman (USC), B. Zhu (LANL) ================
"""
import sys, time, os, pywt
from ROOT import TFile, TTree, TEntryList, gDirectory, TNamed, std, TObject
from ROOT import GATDataSet, MGTEvent, MGTWaveform
import numpy as np
from scipy.signal import butter, lfilter, filtfilt
import scipy.optimize as op
from scipy.ndimage.filters import gaussian_filter
from scipy import interpolate
import waveLibs as wl
import waveModel as wm

def main(argv):
    print "================================================="
    print "Auto Calibration Still Sucks... wait... started!",time.strftime('%X %x %Z')
    print "================================================="
    startT = time.clock()
    # gROOT.ProcessLine("gErrorIgnoreLevel = 3001;") # suppress ROOT error messages

    # let's get some (m)args
    batMode, rangeMode, fileMode, gatMode, singleMode, pathMode = False, False, False, False, False, False
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
        if opt == "-b":
            batMode = True
            import matplotlib
            if os.environ.get('DISPLAY','') == '':
                print('No display found. Using non-interactive Agg backend')
                matplotlib.use('Agg')
            print "Batch mode selected.  A new file will be created."

    # File I/O
    inFile, outFile, bltFile = TFile(), TFile(), TFile()
    gatTree, bltTree, oTree = TTree(), TTree(), TTree()
    theCut, inPath, outPath = "", "", ""

    # Set input and output files
    if rangeMode:
        inPath = "%s/waveSkimDS%d_%d.root" % (pathToInput, dsNum, subNum)
        outPath = "%s/amplatSkimDS%d_%d.root" % (pathToOutput, dsNum, subNum)
    if fileMode:
        inPath = "%s/waveSkimDS%d_run%d.root" % (pathToInput, dsNum, runNum)
        outPath = "%s/amplatSkimDS%d_run%d.root" % (pathToOutput, dsNum, runNum)
    if pathMode:
        inPath, outPath = manualInput, manualOutput
    if gatMode:
        ds = GATDataSet()
        gatPath = ds.GetPathToRun(runNum,GATDataSet.kGatified)
        bltPath = ds.GetPathToRun(runNum,GATDataSet.kBuilt)
        outPath = "%s/amplat_run%d.root" % (pathToOutput, runNum)
    if pathMode and gatMode:
        outPath = manualOutput

    # Initialize trees
    if rangeMode or fileMode or pathMode:
        print "Loading: ", inPath
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
    # theCut = "trapENFCal > 2 && trapENFCal < 10"
    # print "WARNING: Custom cut in use!"
    
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

    Max, ENF, ENAF = std.vector("double")(), std.vector("double")(), std.vector("double")()
    ENFC, ENAFC = std.vector("double")(), std.vector("double")()

    # It's not possible to put the "oTree.Branch" call into a class initializer (waveLibs::latBranch). You suck, ROOT.
    b1, b2, b3 = oTree.Branch("Max",Max), oTree.Branch("ENF",ENF), oTree.Branch("ENAF",ENAF)
    b4, b5 = oTree.Branch("ENFC",ENFC), oTree.Branch("ENAFC",ENAFC)

    # make a dictionary that can be iterated over (avoids code repetition in the loop)
    brDict = {
    "Max":[Max, b1], "ENF":[ENF, b2], "ENAF":[ENAF, b3],
    "ENFC":[ENFC, b4], "ENAFC":[ENAFC, b5],
    }

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

            # be absolutely sure you're matching the right waveform to this hit
            if wf.GetID() != chan:
                print "ERROR -- Vector matching failed.  iList %d  run %d  iEvent %d" % (iList,run,iEvent)
                return

            # Let's start the show
            removeNBeg, removeNEnd = 0, 2
            if dsNum==2 or dsNum==6: removeNBeg = 3
            signal = wl.processWaveform(wf,removeNBeg,removeNEnd)
            data = signal.GetWaveBLSub()

            # Short Trapezoid 
            shortTrap = wl.trapFilter(data, rampTime=100, flatTime=150, decayTime=7200.)
            # Asymmetric Trapezoid
            asymTrap = wl.asymTrapFilter(data, ramp=200, flat=100, fall=40, padAfter=True)
            # Standard Energy Trapezoid
            longTrap = wl.trapFilter(data, rampTime=400, flatTime=250, decayTime=7200.)
            longTrapTS = np.linspace(0, len(longTrap)*10, len(longTrap))
            # Standard Energy Trapezoid with Baseline padded waveform
            padTrap = wl.trapFilter(np.pad(data, (400, 0), mode='symmetric'), rampTime=400, flatTime=250, decayTime=7200.)
            padTrapTS = np.linspace(0, len(padTrap)*10, len(padTrap))

            longTrapInterp = interpolate.interp1d(longTrapTS, np.asarray(longTrap).squeeze())
            padTrapInterp = interpolate.interp1d(padTrapTS, np.asarray(padTrap).squeeze()) 

            # Limit the range from 0 to 1400 samples (0 to 14 us) -- using 2.0 threshold like data for now...
            # Returned start times are all in units of ns! 
            LEt0,_ = wl.walkBackT0(shortTrap, thresh=2.0, rmin=0, rmax=1000) # Leading-Edge on short trapezoid
            Asymt0,_ = wl.walkBackT0(asymTrap, thresh=2.0, rmin=0, rmax=1000) # Leading-Edge on asymmetric trapezoid

            # Amplitude Evaluation -- Standard
            Max[iH] = np.amax(longTrap) # Basically trapENM

            # If fixed pickoff time is < 0, fix to 0. for interpolation
            ENF[iH] = longTrapInterp(np.amax([LEt0-7000+4000+2000, 0.])) # This should be ~trapENF
            ENAF[iH] = longTrapInterp(np.amax([Asymt0-7000+4000+2000, 0.]))

            # Amplitude Evaluation -- Functional correction for t0 triggerWalk and then use padTrap
            # Function is still in development...
            # Function is exp(p0 + p1*E) where p0 ~ 7.8 and p1 ~ -0.45
            # Functional walk back distance is the minimum of the function value or ~5800 (standard value)
            ENFC[iH] = padTrapInterp( np.amax([LEt0-7000+8000+2000-np.amin([np.exp(7.8 - 0.45*Max[iH]),1000.]), 0.]) )
            ENAFC[iH] = padTrapInterp( np.amax([Asymt0-7000+8000+2000-np.amin([np.exp(7.8 - 0.66*Max[iH]),1000.]), 0.]) )

            # WF fitting amplitude


            # Oppie


            # Band Max


            # ------------------------------------------------------------------------
            # End waveform processing.
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