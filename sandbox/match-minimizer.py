#!/usr/bin/env python
import sys, time, pymc, pywt
import numpy as np
from scipy import interpolate
import scipy.optimize as op
import scipy.signal as sig
from scipy.ndimage.filters import gaussian_filter
import matplotlib.pyplot as plt
from ROOT import TFile,TTree,TChain,TEntryList,gDirectory,gROOT,MGTWaveform
import waveLibs as wl
import waveModel as wm

def main(argv):
    """ Get a sweet minimizer + time domain match filter working. Also graduate. """
    scanSpeed = 1.
    iList = -1
    opt1 = ""
    intMode = False
    if (len(argv) >= 1): opt1 = argv[0]
    if "-i" in (opt1): intMode = True

    # Set input data and cuts
    inputFile = TFile("~/project/match-skim/waveletSkimDS5_run23920.root") # calibration
    # inputFile = TFile("~/project/v2-processwfs/waveletSkimDS5_90.root") # noisy BG
    waveTree = inputFile.Get("skimTree")
    theCut = inputFile.Get("cutUsedHere").GetTitle()
    # theCut += " && waveS5/trapENFCal < 1200 && trapENFCal < 10"
    theCut += " && trapENFCal < 30"

    # Print cut and events passing cut
    waveTree.Draw(">>elist", theCut, "entrylist")
    elist = gDirectory.Get("elist")
    waveTree.SetEntryList(elist)
    nList = elist.GetN()
    print "Using cut:\n",theCut,"\n"
    print "Found",waveTree.GetEntries(),"input entries."
    print "Found",nList,"entries passing cuts."

    # Make a figure
    fig = plt.figure(figsize=(15,10), facecolor='w')
    p1 = plt.subplot2grid((6,7), (0,0), colspan=4, rowspan=3) # data & fit
    p2 = plt.subplot2grid((6,7), (3,0), colspan=4, rowspan=3) # matched filter
    p3 = plt.subplot2grid((6,7), (0,4), colspan=3, rowspan=2) # trace 1
    p4 = plt.subplot2grid((6,7), (2,4), colspan=3, rowspan=2) # trace 2
    p5 = plt.subplot2grid((6,7), (4,4), colspan=3, rowspan=2) # trace 3
    plt.show(block=False)

    # Make template(s)
    tSamp, tR, tZ, tAmp, tST, tSlo = 5000, 0, 15, 100, 2500, 10
    tOrig, tOrigTS = wl.MakeSiggenWaveform(tSamp,tR,tZ,tAmp,tST,tSlo)

    # Loop over events
    while True:
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

            # ------------------------------------------------------------------------

            # Load data
            run = waveTree.run
            chan = waveTree.channel.at(iH)
            dataE = waveTree.trapENFCal.at(iH)
            dataENM = waveTree.trapENM.at(iH)
            # dataMax = waveTree.trapENMSample.at(iH)*10. - 4000
            dataMax = waveTree.rt90.at(iH)
            signal = wl.processWaveform(waveTree.MGTWaveforms.at(iH))
            data = signal.GetWaveBLSub()
            dataTS = signal.GetTS()
            dataBL, dataNoise = signal.GetBaseNoise()
            print "%d / %d  Run %d  nCh %d  chan %d  trapENF %.1f  len %d" % (iList,nList,run,nChans,chan,dataE,len(data))

            # Load and scale template
            temp, tempTS = tOrig, tOrigTS
            temp = temp * (dataENM / tAmp)         # scale by amplitudes (not energies)
            tempMax = np.argmax(temp) * 10         # convert to ns
            tempTS = tempTS - (tempMax - dataMax)  # align @ max of rising edge
            tempMax = tempTS[np.argmax(temp)]      # get the new max TS after the shifting
            tempE = dataE # set 'calibrated' energy of template equal to data's trapENFCal.

            # Lowpass or bandpass filter the data
            # NOTE: This is tough.  There might be a feature due to rising edges, between 3e6 and 5e6 Hz.
            #       But that's only from data. Doing it w/ siggen templates yielded squat.
            # B,A = sig.butter(2, [3e6/(1e8/2),5e6/(1e8/2)], btype='bandpass')
            B,A = sig.butter(2,1e6/(1e8/2),btype='lowpass') # nice lowpass filter
            data_LowPass = sig.lfilter(B, A, data)
            B1,A1 = sig.butter(2, [1e5/(1e8/2),1e6/(1e8/2)], btype='bandpass') # attempt at bandpass
            data_BanPass = sig.lfilter(B1, A1, data)


            # -- Run minimizers

            # Set window and parameters
            loWin, hiWin = dataTS[0], dataTS[-1]
            mt, en, slo = dataMax-5, dataE+1, 20.  # guesses
            # print "guesses:        mt %.0f  en %.3f  slo %.2f" % (mt, en, slo)

            # pack into lists
            InterpFn = interpolate.interp1d(tempTS, temp, kind="linear", copy="False", assume_sorted="True")
            dataList = [data, dataTS, dataE, dataMax, loWin, hiWin, dataNoise]
            tempList = [temp, tempTS, tempE, tempMax]
            datas = [dataList, tempList, InterpFn]
            floats = [mt,en,slo]
            guess, guessTS = wm.MakeModel(dataList, tempList, floats, fn=InterpFn)

            start = time.clock()
            MakeTracesGlobal()

            fitFail = False # you should output this to ROOT

            # 0. brute - trying to make a quick initial guess for the slowness parameter
            # FIXME: this keeps returning the lowest allowed value for slowness.
            # bruteArgs = [en, datas]
            # ranges = (slice(dataMax-500,dataMax+300,100), slice(1,100,10))
            # resbrute = op.brute(bruteLnLike, ranges, args=(bruteArgs,), full_output=False, finish=None, disp=False)
            # mt, slo = resbrute[0], resbrute[1]
            # print "brute results:  mt %.0f  en %.3f  slo %.2f" % (mt, en, slo)
            # floats = [mt, en, slo]

            # 1. nelder-mead (~0.5 sec).  nobody does it better.  works almost like magic.
            result = op.minimize(findLnLike, floats, args=datas, method="Nelder-Mead")#,options={"maxiter":10})
            if not result["success"]:
                print "fit 'fail': ", result["message"]
                fitFail = True
            mt, en, slo = result["x"]
            # print "nelder results: mt %.0f  en %.3f  slo %.2f" % (mt, en, slo)
            nelder, nelderTS = wm.MakeModel(dataList, tempList, [mt,en,slo], fn=InterpFn)
            nelderStop = time.clock()

            # 2. powell "polish step" (~0.1 sec).  Sometimes screws up Nelder's good results!
            # floats = [mt,en,slo]
            # powell = op.minimize(findLnLike, floats, args=datas, method="Powell")
            # if not powell["success"]:
            #     fitFail= True
            #     print powell["message"]
            # mt, en, slo = powell["x"]
            # # print "powell results: mt %.0f  en %.3f  slo %.2f" % (mt, en, slo)
            # model, modelTS = wm.MakeModel(dataList, tempList, [mt,en,slo], fn=InterpFn)
            #
            # powellStop = time.clock()
            # print "Time: nelder %.3f  powell %.3f  total %.3f" % (nelderStop-start, powellStop-nelderStop, powellStop-start)

            # take absolute values for parameters
            mt, en, slo = abs(mt), abs(en), abs(slo)

            model, modelTS = wm.MakeModel(dataList, tempList, [mt,en,slo], fn=InterpFn)


            # compute a trap energy off the model (need pinghan's calibration constants)
            trap = np.zeros(1000)
            trap = np.append(trap,wl.trapezoidalFilter(model))
            trapTS = modelTS
            trapMax = np.amax(trap)



            # -- Time domain match filter

            match, matchTS = wm.MakeModel(dataList, tempList, [mt,en,slo], opt="nowindow")
            match = np.flip(match,0)

            # line up the max of the match with the max of the best-fit signal
            modelMaxTime = modelTS[np.argmax(model)]
            matchMaxTime = matchTS[np.argmax(match)]
            matchTS = matchTS + (modelMaxTime - matchMaxTime)

            idx = np.where((matchTS >= dataTS[0]-5) & (matchTS <= dataTS[-1]+5))
            match, matchTS = match[idx], matchTS[idx]

            # make sure match is same length as data, then compute convolution and parameters
            matchMax, matchWidth, matchArea = 0, 9999, 0
            smoothMF = np.zeros(len(dataTS))
            if len(match)!=len(data):
                print "array mismatch: len match %i  len data %i " % (len(match),len(data))
            else:
                # I can't decide which of these is better.
                # smoothMF = gaussian_filter(match * data, sigma=float(5))
                smoothMF = gaussian_filter(match * data_LowPass,sigma=float(5))

                # computer AWESOME match filter parameters
                matchMax = np.amax(smoothMF)
                idx = np.where(smoothMF > matchMax/2.)
                matchWidth = matchTS[idx][-1] - matchTS[idx][0]
                matchArea = np.sum(match[idx]) / float(len(match[idx]))


            # Plot data and template
            p1.cla()
            p1.plot(dataTS,data,'b',label='data')
            idx2 = np.where((tempTS >= dataTS[0]-5) & (tempTS <= dataTS[-1]+5))
            p1.plot(tempTS[idx2],temp[idx2],'r',label='template')
            # p1.plot(dataTS, data_LowPass, label='lowpass')
            # p1.plot(dataTS, data_BanPass, label='bandpass')
            p1.plot(trapTS,trap,color='k',label='bestfit-trap')
            p1.plot(nelderTS,nelder,color='cyan',label='bestfit')
            # p1.plot(modelTS,model,color='magenta',label='powells')
            p1.set_xlabel('Time (s)')
            p1.set_ylabel('Voltage (arb)')
            p1.set_title("Run %i  Ch %i  E %.2f  ENM %.2f  trapMax %.2f" % (run,chan,dataE,dataENM,trapMax))
            p1.legend(loc=4)

            p2.cla()
            p2.plot(dataTS,data,'b',alpha=0.8,label='data')
            # p2.plot(guessTS,guess,'y',alpha=0.8,label='guess')
            p2.plot(matchTS,match,'k',label='match')
            p2.plot(modelTS,model,color='orange',label='model')

            p2.plot(dataTS,data_LowPass,'r',label='lowpass')

            if len(match)==len(data): p2.plot(matchTS,smoothMF,'g',label='smoothMF')

            p2.axvline(matchTS[idx][-1],color='green',alpha=0.7)
            p2.axvline(matchTS[idx][0],color='green',alpha=0.7)

            p2.set_title("mMax %.2f  mWidth %.2f  mAreaNorm %.2f  mM/mW %.4f" % (matchMax,matchWidth,matchArea,matchMax/matchWidth))
            p2.legend(loc=4)


            # plot traces
            p3.cla()
            p3.set_title("maxTime %.1f  Energy %.2f  Slow %.1f" % (mt,en,slo))
            p3.plot(mtTrace[1:])
            p3.set_ylabel('maxTime')
            p4.cla()
            p4.plot(enTrace[1:])
            p4.set_ylabel('energy')
            p5.cla()
            p5.plot(sloTrace[1:])
            p5.set_ylabel('slowness')

            plt.tight_layout()
            plt.pause(scanSpeed)

            # ------------------------------------------------------------------------


def bruteLnLike(x0, *bruteArgs):
    mt, slo = x0
    en, datas = bruteArgs[0]
    floats = [mt, en, slo]
    # print mt, en, slo
    findLnLike(floats, datas)

def findLnLike(floats, datas):
    """ SciPy Minimizer.
        This could be moved to waveModel.py IF you decide you
        don't need to look at the trace arrays anymore.

        Also, the version in this file (match-minimizer.py)
        is different because it tries to align the MAX time
        (not the start time) of the waveforms.
    """
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
    """ This is so 'findLnLike' can write to the trace arrays.
        It has to remain in this file to work. """
    tmp1 = tmp2 = tmp3 = np.empty([1,])
    global mtTrace, enTrace, sloTrace
    mtTrace, enTrace, sloTrace = tmp1, tmp2, tmp3


if __name__ == "__main__":
    main(sys.argv[1:])
