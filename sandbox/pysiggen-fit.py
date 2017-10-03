#!/usr/bin/env python
import sys, time, pywt
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
from matplotlib import gridspec
from ROOT import TFile,TTree,TChain,TEntryList,TNamed,TObject,gDirectory,gROOT,MGTWaveform,std
import waveLibs as wl
import waveModel as wm
import pymc
from scipy.ndimage.filters import gaussian_filter
from scipy import interpolate
from scipy.optimize import curve_fit
from matplotlib.backends.backend_pdf import PdfPages
from multiprocessing import Process, Manager

def main(argv):
    """ Interactive-fit or 'rapid'-fit waveforms that pass a given TCut.
        BUG: Doesn't always work with a TChain.  Add input files together
        with hadd and use a single TFile.
    """
    scanSpeed = 0.2
    iList = -1
    opt1 = ""
    intMode, batMode = False, False
    if (len(argv) >= 1): opt1 = argv[0]
    if "-i" in (opt1):
        intMode = True
        print "Interactive mode selected."
    if "-b" in (opt1):
        batMode = True
        print "Batch mode selected.  A new file will be created."

    # Load template waveform (created by gen-template.py)
    npzfile = np.load("./data/genTemplateWF.npz")
    temp, tempTS, tempE, tempST = npzfile['arr_0']+1, npzfile['arr_1'], npzfile['arr_2'], npzfile['arr_3']*10

    # Set cuts
    # theCut = inputFile.Get("cutUsedHere").GetTitle()
    # DS3 "big DC" cut - PRELIMINARY
    theCut = "trapENFCal > 0.8 && gain==0 && mHClean==1 && isGood && !muVeto && !wfDCBits && !isLNFill1 && !isLNFill2 && trapETailMin < 0 && channel!=596 && channel!=676 && channel!=676 && channel!=612 && channel!=1104 && channel!=1200 && channel!=1334 && channel!=1336 && tOffset < 10 && waveS5/TMath::Power(trapENFCal,1/4) < 1200 && (waveS3-waveS2)/trapENFCal > 100 && (waveS3-waveS2)/trapENFCal < 300 && !(channel==692 && (run==16974 || run==16975 || run==16976 || run==16977 || run==16978 || run==16979)) && butterTime < 11000"
    # theCut += " && trapENFCal > 1.5 && trapENFCal < 2.1"
    # theCut += " && trapENFCal < 20 && trapENFCal > 2 && run > 17480"
    # theCut += " && kvorrT/trapENFCal > 2.2 && trapENFCal > 2 && trapENFCal < 10"

    # Set file I/O and create entry lists
    startT = time.clock()
    inFile = "~/project/wavelet-skim/waveletSkimDS3_1.root"
    # inFile = "~/project/wavelet-skim/hadd/waveletSkimDS3.root"
    outFile = "~/project/fit-skim/fitSkimDS3_1.root"
    inputFile = TFile(inFile)
    waveTree = inputFile.Get("skimTree")
    print "Found",waveTree.GetEntries(),"input entries.  Using cut:\n",theCut,"\n"
    waveTree.Draw(">>elist", theCut, "entrylist")
    elist = gDirectory.Get("elist")
    waveTree.SetEntryList(elist)
    nList = elist.GetN()
    print "Found",nList,"entries passing cuts."
    stopT=time.clock()
    print "Data loading time (s): ", (stopT-startT)

    # In batch mode ONLY, create an output file+tree & append new branches
    outputFile = TFile()
    outTree = TTree()
    if batMode:
        outputFile = TFile(outFile, "RECREATE")
        print "Attempting tree copy to",outFile
        outTree = waveTree.CopyTree("")
        outTree.Write()
        print "Wrote",outTree.GetEntries(),"entries."
        cutUsed = TNamed("cutUsedHere",theCut)
        cutUsed.Write()
    fitStart = std.vector("double")()
    fitE = std.vector("double")()
    fitSlo = std.vector("double")()
    fitStartSD = std.vector("double")()
    fitESD = std.vector("double")()
    fitSloSD = std.vector("double")()
    fitChi2NDF = std.vector("double")()
    fitLnLike = std.vector("double")()
    tExp1 = std.vector("double")()
    tExp2 = std.vector("double")()
    tPol0 = std.vector("double")()
    tPol1 = std.vector("double")()
    tPol2 = std.vector("double")()
    tPol3 = std.vector("double")()
    baseAvg = std.vector("double")()
    baseNoise = std.vector("double")()
    bFitStart = outTree.Branch("fitStart", fitStart)
    bFitE = outTree.Branch("fitE", fitE)
    bFitSlo = outTree.Branch("fitSlo", fitSlo)
    bFitStart_sd = outTree.Branch("fitStartSD", fitStartSD)
    bFitE_sd = outTree.Branch("fitESD", fitESD)
    bFitSlo_sd = outTree.Branch("fitSloSD", fitSloSD)
    bFitChi2NDF = outTree.Branch("fitChi2NDF", fitChi2NDF)
    bFitLnLike = outTree.Branch("fitLnLike", fitLnLike)
    bTExp1 = outTree.Branch("tExp1", tExp1)
    bTExp2 = outTree.Branch("tExp2", tExp2)
    bTPol0 = outTree.Branch("tPol0", tPol0)
    bTPol1 = outTree.Branch("tPol1", tPol1)
    bTPol2 = outTree.Branch("tPol2", tPol2)
    bTPol3 = outTree.Branch("tPol3", tPol3)
    bBaseAvg = outTree.Branch("baseAvg", baseAvg)
    bBaseNoise = outTree.Branch("baseNoise", baseNoise)

    # Make a figure
    # with PdfPages('multipage_pdf.pdf') as pdf:
    fig = plt.figure(figsize=(11,7), facecolor='w')
    p1 = plt.subplot2grid((6,7), (0,0), colspan=4, rowspan=2) # original
    p2 = plt.subplot2grid((6,7), (2,0), colspan=4, rowspan=3) # rising edge
    p3 = plt.subplot2grid((6,7), (5,0), colspan=4, rowspan=1) # residual
    p4 = plt.subplot2grid((6,7), (0,4), colspan=3, rowspan=2 )           # trace 1
    p5 = plt.subplot2grid((6,7), (2,4), colspan=3, rowspan=2, sharex=p4) # trace 2
    p6 = plt.subplot2grid((6,7), (4,4), colspan=3, rowspan=2, sharex=p4) # trace 3
    if not batMode: plt.show(block=False)

    # Setup multiprocessing
    manager = Manager()
    returnDict = manager.dict()

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

        fitStart.assign(nChans,-9999)
        fitE.assign(nChans,-9999)
        fitSlo.assign(nChans,-9999)
        fitStartSD.assign(nChans,-9999)
        fitESD.assign(nChans,-9999)
        fitSloSD.assign(nChans,-9999)
        fitChi2NDF.assign(nChans,-9999)
        fitLnLike.assign(nChans,-9999)
        tExp1.assign(nChans,-9999)
        tExp2.assign(nChans,-9999)
        tPol0.assign(nChans,-9999)
        tPol1.assign(nChans,-9999)
        tPol2.assign(nChans,-9999)
        tPol3.assign(nChans,-9999)
        baseAvg.assign(nChans,-9999)
        baseNoise.assign(nChans,-9999)

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
            dataBaseline, dataNoise = signal.GetBaseNoise()

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
            dataList = [data, dataTS, dataE, dataST, loWin, hiWin]
            tempList = [temp, tempTS, tempE, tempST]

            # Optionally save something to a file
            if saveMe:
                print "Saved entry",iList,iH
                np.savez("./data/tailSlopeInputs.npz",rawList,tempList)

            # Recreate the guess and the guess's rising edge
            guessFull, guessFullTS = wm.MakeModel(dataList, tempList, [dataST,dataE,1.], opt="full")
            guess, guessTS = wm.MakeModel(dataList, tempList, [dataST,dataE,1.], opt="!fancy")

            # Make an "almost complete" guess - no MCMC
            # st, en, slo = dataST-100, dataE, 5
            # InterpFn = interpolate.interp1d(tempTS, temp, kind="linear", copy="False", assume_sorted="True")
            # model, modelTS = wm.MakeModel(dataList, tempList, [st,en,slo], fn=InterpFn)

            # Fit with MCMC and get best-fit parameters
            numSteps, burnIn = 3000, 1800  # default: 10000, 5000.  fast: 3000, 1800  long test: 20000,10000
            wfModel = wm.TemplateModel( dataList, dataNoise, tempList )
            p = Process(target=RunMCMC, args=(wfModel, numSteps, burnIn, returnDict))
            p.start()
            p.join()
            startTimeTr = returnDict["startTimeTr"]
            energyTr = returnDict["energyTr"]
            slownessTr = returnDict["slownessTr"]
            st = np.median(startTimeTr[burnIn:])
            en = np.median(energyTr[burnIn:])
            slo = np.median(slownessTr[burnIn:])
            InterpFn = interpolate.interp1d(tempTS, temp, kind="linear", copy="False", assume_sorted="True")
            model, modelTS = wm.MakeModel(dataList, tempList, [st,en,slo], fn=InterpFn)

            # Save some extra parameters for the ROOT output
            # Calculate residual, Chi2/NDF, likelihood, etc.
            st_std = np.std(startTimeTr[burnIn:])
            en_std = np.std(energyTr[burnIn:])
            slo_std = np.std(slownessTr[burnIn:])
            residual = model - data
            frac = (np.power(data - model, 2)) / np.abs(model)
            chi2NDF = np.sum(frac) / len(model)
            inv_sigma2 = 1.0/(dataNoise**2)
            lnLike = -0.5*(np.sum((data-model)**2*inv_sigma2 - np.log(inv_sigma2)))

            # ** Do a separate & simple fit of the tail slope **
            # TODO: Add this to process-waveforms.py
            idxMax = np.where(guessFull == guessFull.max()) # returns an array/tuple
            idxMax = idxMax[0][0] # "cast" to int
            tail, tailTS = waveDenoised[idxMax:], waveTS[idxMax:]
            popt,_ = curve_fit(wl.tailModelPol, tailTS, tail) # poly fit
            pol0, pol1, pol2, pol3 = popt[0], popt[1], popt[2], popt[3]
            a, b = dataE, 72000
            popt2,_ = curve_fit(wl.tailModelExp, tailTS, tail, p0=[a,b]) # expo fit
            e1, e2 = popt2[0], popt2[1]

            # Assign values to output vectors and fill branches
            fitStart[iH], fitStartSD[iH] = st, st_std
            fitE[iH], fitESD[iH] = en, en_std
            fitSlo[iH], fitSloSD[iH] = slo, slo_std
            fitChi2NDF[iH] = chi2NDF
            fitLnLike[iH] = lnLike
            tExp1[iH], tExp2[iH] = e1, e2
            tPol0[iH], tPol1[iH], tPol2[iH], tPol3[iH] = pol0, pol1, pol2, pol3
            baseAvg[iH] = dataBaseline
            baseNoise[iH] = dataNoise
            if batMode:
                bFitStart.Fill()
                bFitE.Fill()
                bFitSlo.Fill()
                bFitStart_sd.Fill()
                bFitE_sd.Fill()
                bFitSlo_sd.Fill()
                bFitChi2NDF.Fill()
                bFitLnLike.Fill()
                bTExp1.Fill()
                bTExp2.Fill()
                bTPol0.Fill()
                bTPol1.Fill()
                bTPol2.Fill()
                bTPol3.Fill()
                bBaseAvg.Fill()
                bBaseNoise.Fill()
                if iList % 5000 == 0:
                    outTree.Write("",TObject.kOverwrite)
                    print "%d / %d entries saved (%.2f %% done)." % (iList,nList,100*(float(iList)/nList))

            # If not in batch mode, fill the figure
            if batMode: continue
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
            # special: plot the values of the trace after burn-in
            # to see how the model is covering the "money-zone"/rising edge after it's converged.
            # for i in range(burnIn,numSteps):
                # st_tr, en_tr, slo_tr = M.trace('startTime')[i], M.trace('energy')[i], M.trace('slowness')[i]
                # trace, traceTS = wm.MakeModel(dataList, tempList, [st_tr,en_tr,slo_tr], fn=InterpFn)
                # p2.plot(traceTS, trace, color='red',alpha=0.1,linewidth=2)
            p2.plot(modelTS,model,color='red',linewidth=3,alpha=0.7,label='Best Fit')
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
            p4.plot(startTimeTr[:])
            p4.set_ylabel('startTime')
            p4.axvline(x=burnIn,color='red',alpha=0.5)

            p5.cla()
            p5.plot(energyTr[:])
            p5.set_ylabel('energy')
            p5.axvline(x=burnIn,color='red',alpha=0.5)

            p6.cla()
            p6.plot(slownessTr[:])
            p6.set_ylabel('slowness')
            p6.axvline(x=burnIn,color='red',alpha=0.5)

            plt.tight_layout()
            plt.subplots_adjust(hspace=0.35)
            plt.pause(scanSpeed)
            # pdf.savefig()

    # End loop over events
    if batMode:
        outTree.Write("",TObject.kOverwrite)
        print "Wrote",outTree.GetBranch("channel").GetEntries(),"entries in the copied tree,"
        print "and wrote",bFitStart.GetEntries(),"entries in the new branches."


def RunMCMC(wfModel, numSteps, burnIn, returnDict):
    """ This is outside main() and called with 'multiprocessing' to avoid memory leaks. """
    M = pymc.MCMC( pymc.Model( wfModel ) )
    M.use_step_method(pymc.Metropolis, M.startTime, proposal_sd=100., proposal_distribution='Normal')
    M.use_step_method(pymc.Metropolis, M.energy, proposal_sd=1., proposal_distribution='Normal')
    M.use_step_method(pymc.Metropolis, M.slowness, proposal_sd=100., proposal_distribution='Normal')
    M.sample(iter=numSteps, verbose=0)
    returnDict["startTimeTr"] = M.trace("startTime")[:]
    returnDict["energyTr"] = M.trace("energy")[:]
    returnDict["slownessTr"] = M.trace("slowness")[:]

if __name__ == "__main__":
    main(sys.argv[1:])
