#!/usr/local/bin/python
import sys, time, pymc, pywt
import numpy as np
from scipy import interpolate
import scipy.optimize as op
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from ROOT import TFile,TTree,TChain,TEntryList,gDirectory,gROOT,MGTWaveform
import waveLibs as wl
import waveModel as wm



def main(argv):
    """ Interactive-fit or 'rapid'-fit waveforms that pass a given TCut.
        BUG: Doesn't always work with a TChain.  Add input files together
        with hadd and use a single TFile.
    """
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
    theCut += " && waveS5/trapENFCal < 1200 && trapENFCal > 1.5 && trapENFCal < 3"

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
    p1 = plt.subplot2grid((6,7), (0,0), colspan=4, rowspan=2) # original
    p2 = plt.subplot2grid((6,7), (2,0), colspan=4, rowspan=3) # rising edge
    p3 = plt.subplot2grid((6,7), (5,0), colspan=4, rowspan=1) # residual
    p4 = plt.subplot2grid((6,7), (0,4), colspan=3, rowspan=2) # trace 1
    p5 = plt.subplot2grid((6,7), (2,4), colspan=3, rowspan=2) # trace 2
    p6 = plt.subplot2grid((6,7), (4,4), colspan=3, rowspan=2) # trace 3
    plt.show(block=False)

    # Make template(s)
    # npzfile = np.load("./data/genTemplateWF.npz") # gen-template.py
    # temp, tempTS, tempE, tempST = npzfile['arr_0']+1, npzfile['arr_1'], npzfile['arr_2'], npzfile['arr_3']*10

    samp, r, z, tempE, tempST, smooth = 2016, 0, 15, 3938, 1000, 100 # gen-template.py
    # samp, r, z, tempE, tempST, smooth = 5000, 0, 15, 10, 2500, 100  # huge template
    # samp, r, z, tempE, tempST, smooth = 2016, 0, 15, 10, 1000, 100 # regular size template
    # samp, r, z, tempE, tempST, smooth = 2016, 30, 30, 10, 1000, 100 # regular size temp, slower rise
    # samp, r, z, tempE, tempST, smooth = 500, 0, 15, 10, 100, 100  # small template
    temp, tempTS = wl.MakeSiggenWaveform(samp,r,z,tempE,tempST,smooth)
    tempST = tempST * 10 # convert to ns


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
            dataList = [data, dataTS, dataE, dataST, loWin, hiWin, dataNoise]
            tempList = [temp, tempTS, tempE, tempST]

            # Optionally save something to a file
            if saveMe: np.savez("./data/tailSlopeInputs.npz",rawList,tempList)

            # Recreate the guess and the guess's rising edge
            guessFull, guessFullTS = wm.MakeModel(dataList, tempList, [dataST,dataE,1.], opt="full")
            guess, guessTS = wm.MakeModel(dataList, tempList, [dataST,dataE,1.], opt="!fancy")

            # Make an "almost complete" guess - no fitting
            # st, en, slo = dataST-100, dataE, 5
            InterpFn = interpolate.interp1d(tempTS, temp, kind="linear", copy="False", assume_sorted="True")
            # model, modelTS = wm.MakeModel(dataList, tempList, [st,en,slo], fn=InterpFn)

            # Fit with MCMC and get best-fit parameters
            # numSteps, burnIn = 10000, 5000 # default - 10000, 5000.  try 3000, 2000.  long test: 20000,10000
            # myModel = TemplateModel( dataList, dataNoise, tempList )
            # waveModel = pymc.Model( myModel )
            # M = pymc.MCMC( waveModel )
            # M.use_step_method(pymc.Metropolis, M.startTime, proposal_sd=100., proposal_distribution='Normal')
            # M.use_step_method(pymc.Metropolis, M.energy, proposal_sd=1., proposal_distribution='Normal')
            # M.use_step_method(pymc.Metropolis, M.slowness, proposal_sd=100., proposal_distribution='Normal')
            # M.sample(iter=numSteps, verbose=0) # do the fit
            # st = np.median(M.trace('startTime')[burnIn:])
            # en = np.median(M.trace('energy')[burnIn:])
            # slo = np.median(M.trace('slowness')[burnIn:])
            # InterpFn = interpolate.interp1d(tempTS, temp, kind="linear", copy="False", assume_sorted="True")
            # model, modelTS = wm.MakeModel(dataList, tempList, [st,en,slo], fn=InterpFn)


            # Fit with SciPy minimizer(s) and get best-fit parameters

            # how long does it take?
            start = time.time()

            MakeTracesGlobal() # creates 3 global arrays: startTrace, enTrace, sloTrace
            datas = [dataList, tempList, InterpFn]

            # floats = [1.]
            # result = op.basinhopping(findLnLike, floats, T=1000, stepsize=15, niter_success=2, minimizer_kwargs={"args":datas})
            # print "slowness:",result["x"],"dataE",dataE,"dataST",dataST
            # floats = [dataST, dataE, result["x"]]
            # MakeTracesGlobal() # resets traces

            floats = [dataST, dataE, 20.]
            result = op.minimize(findLnLike, floats, args=datas, method="Nelder-Mead")
            if not result["success"]: print result["message"]
            st, en, slo = result["x"]

            model, modelTS = wm.MakeModel(dataList, tempList, [st,en,slo], fn=InterpFn)

            stop = time.time()
            print "fitting took",stop-start



            # Calculate residual, Chi2/NDF, likelihood, etc.
            residual = model - data
            frac = (np.power(data - model, 2)) / np.abs(model)
            chi2NDF = np.sum(frac) / len(model)
            inv_sigma2 = 1.0/(dataNoise**2)
            lnLike = -0.5*(np.sum((data - model)**2*inv_sigma2 - np.log(inv_sigma2)))

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

            p2.cla()
            p2.plot(dataTS,data,color='blue',label='Data')
            p2.plot(guessTS,guess,color='orange',label='Guess')
            p2.plot(modelTS,model,color='red',linewidth=3,label='Best Fit')
            p2.legend(loc=4)

            p3.cla()
            p3.set_xlabel("Time [ns]", x=0.95, ha='right')
            p3.set_ylabel("Residual [ADC]")
            p3.plot(guessTS,residual, color='red')
            p3.axhline(y = 0, color='blue', alpha=0.3)
            p3.axhline(y = dataNoise, color='blue', alpha=0.3)
            p3.axhline(y = -1.0*dataNoise, color='blue', alpha=0.3)

            p4.cla()
            p4.set_title("startTime %.1f  Energy %.2f  Slow %.1f" % (st,en,slo))

            p4.plot(stTrace[1:])
            p4.set_ylabel('startTime')

            p5.cla()
            p5.plot(enTrace[1:])
            p5.set_ylabel('energy')

            p6.cla()
            p6.plot(sloTrace[1:])
            p6.set_ylabel('slowness')

            plt.tight_layout()
            plt.subplots_adjust(hspace=0.35)
            plt.pause(scanSpeed)

            # ------------------------------------------------------------------------


def findLnLike(floats, datas):
    """ SciPy Minimizer.
        This could be moved to waveModel.py IF you decide you
        don't need to look at the trace arrays anymore. """
    global stTrace, enTrace, sloTrace

    # Yeezus.  At this point should you just write a damn class?
    # Unpack parameters.
    dataList, tempList, InterpFn = datas
    data, dataTS, dataE, dataST, loWin, hiWin, dataNoise = dataList
    temp, tempTS, tempE, tempST = tempList
    st = en = slo = 0
    if len(floats)==1: # basinhopping case
        slo, st, en = floats[0], dataST, dataE

    if len(floats)==3: # minimize case
        st, en, slo = floats

    # Make a trace
    stTrace = np.append(stTrace,st)
    enTrace = np.append(enTrace,en)
    sloTrace = np.append(sloTrace,slo)
    # print st, en, slo

    # Build the model to compare with data
    model, modelTS = wm.MakeModel(dataList, tempList, [st,en,slo], fn=InterpFn)

    # Find LL and return
    lnLike = 0.5 * np.sum ( np.power((data-model)/dataNoise, 2) - np.log( 1 / np.power(dataNoise,2) ) )
    return lnLike


def MakeTracesGlobal():
    """ This is so 'findLnLike' can write to the trace arrays.
        It has to remain in this file to work. """
    tmp1 = tmp2 = tmp3 = np.empty([1,])
    global stTrace, enTrace, sloTrace
    stTrace, enTrace, sloTrace = tmp1, tmp2, tmp3


if __name__ == "__main__":
    main(sys.argv[1:])
