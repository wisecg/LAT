#!/usr/bin/env python
"""
    This code performs the Coherent Bragg-Primakov Axion analysis using pymc3

    0) This code requires a couple of inputs:
        a) Angle-averaged Axion PDF (provided by Wenqin)
            => In the current implementation, this PDF is loaded and wrangled everytime the code is run, which is admittedly very inefficient in terms of overhead.
        b) Final low energy spectra from 5 - 20 keV
            => In the current implementation, the ROOT file is converted to a HDF file (only once)
        c) Final low energy efficiencies (in a ROOT file) along with data-set specific run lists (in a text file)
            => In the current implementation of the code, these inputs are calculated in dsi.py

    Limits we want to beat (for Lambda):
    lambda = np.power(g_agg * 1e8, 4)
    CDMS: 0.00332
    EDELWEISS: 0.002137
    MJD: Hopefully somewhere around 0.001
    DAMA: 0.00084

    Energy Ranges (as percentage of total axion content):
    From 2.4 - 18: 0.988110 of total axion
    From 3.0 - 18: 0.958015 of total axion
    From 3.2 - 18: 0.943619 of total axion
    From 3.5 - 18: 0.914782 of total axion
    From 3.8 - 18: 0.877228 of total axion
    From 4.0 - 18: 0.848106 of total axion

    Wenqin's Current Fits (mean):
    5 -- 20:
    lambda: 175.3/206623.6666 = 0.00085
    4 -- 20:
    lambda: 181.7/258812.4304 = 0.0007

    Brian's results (95 Limit):
    New:
    5 -- 20: 159.16 => 0.00077

    ok, the upper limit was 143. and now 146.
    but, the expected number of axions was 213100.962 and now 194123.1755
    194123.1755/213100.96=91%
    so the change is 9%
    or 146/194123.1755 v.s. 143./213100.962, which is 0.00075 v.s. 0.00067
    so 12% more conservative

    Axion Integrals:
    194112.5744
    Central: 194123.1755 -> 159.62/194123.1755 = 0.000822
    Lo90: 192770.2718 -> 164.75/192770.2718 = 0.000855
    IS: 188168.4693 -> 162.33/188168.4693 = 0.000863
    Sideband: 202775.0857 -> 144.308/202775.0857 = 0.000712
    -- Not using these --
    Hi90: 195454.8897 ->
    Hi68: 194939.535 ->

    DS6a Central: 151.98 -> 151.98/194123.1691 = 0.00078

    Efficiency Studies:
    http://mjcalendar.npl.washington.edu/indico/event/2761/material/slides/0.pdf
    http://mjcalendar.npl.washington.edu/indico/event/2813/material/slides/0.pdf
    https://mjcalendar.npl.washington.edu/indico/event/2781/material/slides/0.pdf

"""

import sys, os, imp, time
import numpy as np
import pymc3 as pm
import matplotlib.pyplot as plt
from matplotlib import gridspec
import seaborn as sns
import theano.tensor as tt
from scipy.interpolate import interp1d
from scipy.stats import norm
import pandas as pd
sns.set(style='darkgrid', context='talk')

wl = imp.load_source('waveLibs', '{}/waveLibs.py'.format(os.environ['LATDIR']))
dsi = imp.load_source('dsi', '{}/dsi.py'.format(os.environ['LATDIR']))

bkg = dsi.BkgInfo()
cal = dsi.CalInfo()
det = dsi.DetInfo()

# Define some global parameters -- I can definitely organize this better but whatever
inDir = os.environ['LATDIR']+'/data/MCMC'
dsList = ['5B', '5C', '6']
# dsList = ['5B']

# The start and end time are for DS5b
startTime = 1485575988
# endTime = 1489762153 # DS5b

# The start and end time are for DS5c
# startTime = 1489773072
# endTime = 1494223870 # DS5c

# The start and end time are for DS5a
# startTime = 1494535909
endTime = 1523845062 # DS6a end day (April 16th 2018)

# This is for calculating stuff... that I haven't set yet
totExposure = 2621.72447852 # DS5b, DS5c, DS6

# Fake start times for testing
# startTime = 1483228800
# endTime = 1514678400 # End of 2017
# endTime = 1546300800 # End of 2018
yearInSec = 31536000 # Number of seconds in 1 year
nBurn = 1000 # Number of burn-in samples for the MCMC

# Energy range (also defines fitting range!)
# The tritium endpoint is 18 keV, the axion PDF stops at 12 keV.
# We want to extend the PDF to slightly beyond the tritium endpoint in order to capture the flat bkg
energyThresh, energyThreshMax, binSize = 5.0, 20.0, 0.1
energyBins = np.linspace(energyThresh,energyThreshMax, round((energyThreshMax-energyThresh)/binSize)+1)

def main(argv):
    # Some random booleans
    bDebug, bDrawPDF, bDiagnostic, bFinalSpectra, bBackend, bPPC = False, False, False, False, False, False
    # Default backendDir is my own directory!
    backendDir = '/Users/brianzhu/code/LAT/data/MCMC/AveragedAxion2_Enr_Central_DS5breso_5_20'
    bSample, bFull = False, True
    axionBinSize = 5 # This is in units of minutes
    pdfArrDict, pdfFlatDict = {}, {}
    # Set initial seed number here -- this is an option that can be changed
    dsNum, seedNum = '5b', 1
    effMode = None

    if len(argv)==0:
        help_options()
        return
    for i,opt in enumerate(argv):
        if opt == '-h' or opt == '--help':
            help_options()
            return
        if opt == '-reduce':
            # Save data into dataframe -- if this hasn't been done before -- exits immediately afterwards
            print('Converting ROOT data into DataFrame -- this only needs to be done once!')
            reduceData()
            return
        if opt == '-seed':
            seedNum = int(argv[i+1])
            print('Setting seed to: {}'.format(seedNum))
        if opt == '-debug':
            bDebug = True
            print('Debug mode ON!')
        if opt == '-ds':
            dsNum = str(argv[i+1])
            if dsNum not in ['5b', '6a']:
                print('Error: Dataset {} is not valid! Can only use 5b or 6a!'.format(dsNum))
                return
            print('Setting Resolution to Dataset',dsNum)
        if opt == '-drawPDF':
            bDrawPDF = True
            print('Drawing PDFs only!')
        if opt == '-drawSpectra':
            bFinalSpectra = True
            print('Drawing Final Spectra! MCMC chain must be sampled!')
        if opt == '-diagnostic':
            bDiagnostic = True
            print('Generating Model Diagnostics, MCMC chain needs to be sampled!')
        if opt == '-eff':
            effMode = str(argv[i+1])
            if effMode not in ['Lo90', 'Hi90', 'Sideband', 'IS']:
                print('Error: Dataset {} is not valid! Can only use Lo90, Hi90, Sideband, IS!'.format(effMode))
                return
            print('Setting Resolution to Dataset',dsNum)
        if opt == '-noaxion':
            bFull = False
            print('Building Model without Axion PDF!')
        if opt == '-sample':
            bSample = True
            print('MCMC Sampling ON')
        if opt == '-backend':
            bBackend = True
            backendDir = str(argv[i+1])
            print('Using Backend Directory: {}'.format(backendDir))
        if opt == '-ppc':
            bPPC = True
            print('Generating Posterior Predictive Checks, MCMC chain needs to be sampled!')
        if opt == '-setbinsize':
            axionBinSize = int(argv[i+1])
            print('Setting axion bin size to be {} minutes -- Warning anything different from the default has NOT been vetted'.format(axionBinSize))

    # Build Axion PDF -- Use low edges of PDF to generate bins
    AxionArr, nTimeBins, timeBinLowEdge, removeMask, eeMask, effNorm = convertaxionPDF(startTime, endTime, debug=bDebug, axionBinSize=axionBinSize, dsNum=dsNum, effMode = effMode)
    pdfArrDict['Axion'] = AxionArr

    print('Generated Axion PDF!')
    print('Axion Integral: {} -- Shape: {}'.format(AxionArr.sum(), AxionArr.shape))

    # Load Background data and bin exactly as PDF
    df = pd.read_hdf('{}/Bkg_Spectrum.h5'.format(inDir))
    # Sort by UnixTime -- I forget why I did this...
    df.sort_values(by=['UnixTime'], inplace=True)
    # Select Enriched detectors only, grab numpy array of the values
    dataArr = df[['Energy', 'UnixTime']].loc[(df['isEnr']==1) & (df['UnixTime'] >= startTime) & (df['UnixTime'] <= endTime)].values
    # Bin the data into a 2D histogram -- here we use the exact time bin edges of the axion PDF to make sure everything is binned properly
    pdfArrDict['Data'], xedges, yedges = np.histogram2d(dataArr[:,0], dataArr[:,1], bins=[energyBins, timeBinLowEdge])

    # Build Tritium PDF -- need interpolation because it's 0.2 keV bins!
    # NOTE: We may want to use 0.2 keV bins, the difference may not be noticeable and the sampling will be significantly faster
    # NOTE: I can probably just generate this into 0.1 keV bins at some point but I've been lazy...
    dfTrit = pd.read_hdf('{}/TritSpec.h5'.format(inDir))
    tritEnergy = dfTrit['Energy'].values
    tritSpec = dfTrit['Tritium'].values
    # The tritium calculation has some slightly negative values around 0 (fluctuations), get rid of them!
    tritSpec = np.array([100*x if x >= 0. else 0. for x in tritSpec])
    tritInterp = interp1d(tritEnergy, tritSpec, fill_value='extrapolate')
    tritList = [tritInterp(x) for x in energyBins[:-1]]
    pdfArrDict['Tritium'] = np.array(tritList*nTimeBins).reshape(nTimeBins, len(tritList)).T*eeMask
    pdfArrDict['Bkg'] = np.ones(pdfArrDict['Axion'].shape)*eeMask

    # Generate the Gaussian PDFs
    meansDict = {'Fe55':6.54, 'Zn65':8.98, 'Ge68':10.37}
    gausDict = generateGaus(energyBins[:-1], meansDict, dsNum)
    for name, Arr in gausDict.items():
        pdfArrDict.setdefault(name, np.array(Arr.tolist()*nTimeBins).reshape(nTimeBins, len(Arr)).T)
        pdfArrDict[name] *= eeMask

    # drawPDFs is before applying removeMask so we can make the comparison of removing vs not removing, etc
    if bDrawPDF:
        drawPDFs(pdfArrDict, removeMask, timeBinLowEdge)
        return

    # Remove columns with zeros -- the zeros are from the times when there is no data-taking
    # This will considerably reduce the DOF and speed up the sampling
    pdfArrDict['Axion'] = np.delete(pdfArrDict['Axion'], removeMask, axis=1)
    pdfArrDict['Data'] = np.delete(pdfArrDict['Data'], removeMask, axis=1)
    pdfArrDict['Tritium'] = np.delete(pdfArrDict['Tritium'], removeMask, axis=1)
    pdfArrDict['Fe55'] = np.delete(pdfArrDict['Fe55'], removeMask, axis=1)
    pdfArrDict['Zn65'] = np.delete(pdfArrDict['Zn65'], removeMask, axis=1)
    pdfArrDict['Ge68'] = np.delete(pdfArrDict['Ge68'], removeMask, axis=1)
    pdfArrDict['Bkg'] = np.delete(pdfArrDict['Bkg'], removeMask, axis=1)
    pdfArrDict['Time'] = np.delete(timeBinLowEdge, removeMask)
    effNorm = np.delete(effNorm, removeMask, axis=1)

    print('Generated all PDFs')
    if bDebug:
        print("Data Shape", pdfArrDict['Data'].shape)
        print("Axion Shape", pdfArrDict['Axion'].shape)

    # Flatten 2D arrays into 1D for ease of sampling
    pdfFlatDict['Data'] = pdfArrDict['Data'].flatten()
    pdfFlatDict['Tritium'] = pdfArrDict['Tritium'].flatten()
    pdfFlatDict['Axion'] = pdfArrDict['Axion'].flatten()
    pdfFlatDict['Bkg'] = pdfArrDict['Bkg'].flatten()
    pdfFlatDict['Fe55'] = pdfArrDict['Fe55'].flatten()
    pdfFlatDict['Zn65'] = pdfArrDict['Zn65'].flatten()
    pdfFlatDict['Ge68'] = pdfArrDict['Ge68'].flatten()

    if bDebug:
        for key in pdfFlatDict:
            print('Flattend {} shape: '.format(key), pdfFlatDict[key].shape)
            print('Flattened all PDFs')

    model = constructModel(pdfFlatDict, energyBins, bFull = bFull)
    print('Built model(s)')

    if not bBackend:
        backendDir = '{}/AveragedAxion2_Enr_LowEff_DS5breso_{}_{}'.format(inDir, int(energyThresh), int(energyThreshMax))

    # Sample Here - I didn't include an option to change the number of samples...
    # At LANL, it takes approximately 3 days to sample 16000 samples in one chain
    if bSample:
        with model:
            db = pm.backends.Text(backendDir) # This backend will save the chain into a csv file
            trace = pm.sample(draws=15000, chains=1, n_init=1000, chain_idx=seedNum, seed=seedNum, tune=1000, progressbar=True, trace=db)

    # Perform Diagnostics on the model -- this can only be done AFTER MCMC sampling!
    if bDiagnostic:
        modelDiagnostics(model, pdfArrDict=pdfArrDict, backendDir=backendDir, unNormAxion=AxionArr, effNorm=effNorm)

    if bPPC:
        posteriorChecks(model, pdfArrDict=pdfArrDict, backendDir=backendDir)

    if bFinalSpectra:
        trace = pm.backends.text.load(backendDir, model)
        drawFinalSpectra(trace=trace[nBurn:], pdfDict=pdfArrDict)


def constructModel(pdfDict, energyBins, bFull=True):
    """
        Construct pymc3 model, possibility to create a null hypothesis if bFull is False (no axion signal)
    """
    with pm.Model() as model:
        # Perhaps should switch to a more informative prior in the future
        # These are referred to as free_RVs
        if bFull:
            Axion = pm.HalfFlat("Axion")
        Tritium = pm.HalfFlat("Tritium")
        Bkg = pm.HalfFlat("Bkg")
        Fe55 = pm.HalfFlat("Fe55")
        Zn65 = pm.HalfFlat("Zn65")
        Ge68 = pm.HalfFlat("Ge68")

        # Generate array of deterministic variables (per bin)
        if bFull:
            det = Tritium*pdfDict['Tritium'] + Bkg*pdfDict['Bkg'] + Axion*pdfDict['Axion'] + Fe55*pdfDict['Fe55'] + Zn65*pdfDict['Zn65'] + Ge68*pdfDict['Ge68']
        else:
            det = Tritium*pdfDict['Tritium'] + Bkg*pdfDict['Bkg'] + Fe55*pdfDict['Fe55'] + Zn65*pdfDict['Zn65'] + Ge68*pdfDict['Ge68']

        # Define the likelihood, this is the observed_RVs
        L = pm.Poisson("L", mu=det, observed=pdfDict['Data'])
    return model


def convertaxionPDF(startTime, endTime, axionBinSize = 5, debug = False, dsNum = '5b', effMode = None):
    """
        Converts 2D histogram PDF for axion into numpy array

        axionBinSize determines the minutes per bin

        startTime/endTime are determined by the unix start/stop time of the first/last runs

        Wenqin's Axion PDF is generated for Lambda = 1.
        The units are in counts/kg/day/keV

        Efficiency Correction:
        - Axion PDF * Efficiency * Exposure: counts/kg/day/keV * (kg-day-keV) = Counts
        - For the fit, rescale the eff*exp*axion pdf by maximum value in the PDF = Counts/Norm
        - After the fit, integrate scaled PDF

        returns: NormAxionArr, nTimeBins, timeBinLowEdge, removeMask, eeMask, effNorm

        nTimeBins and timeBinLowEdge: used to bin the data and every other PDF EXACTLY like the axion PDF
        removeMask: 2D (time vs Energy) array of columns where the exposure is 0 (detectors are off during these time bins)
        eeMask: 2D (time vs Energy) array of the exposure*efficiency -- every PDF needs to be multiplied by this!
        effNorm: 2D (time vs Energy) array of the normalization factor (efficiency*exposure/exposure only)
    """

    import ROOT
    inDir = os.environ['LATDIR']+'/data/MCMC'
    fAxion = ROOT.TFile('{}/Axion_averaged_MJD_DS{}reso_Elow0_Ehi20_minutes1_2017_2018.root'.format(inDir, dsNum))

    hday = fAxion.Get('hday')
    hday.RebinX(axionBinSize) # 5 minute bins

    # Find bins to cut off PDF at
    bRepeat, endBin, endBin2 = False, 0, 0
    nHistTime = hday.GetNbinsX()
    endTimeHist = hday.GetXaxis().GetBinLowEdge(nHistTime)
    startBin = hday.GetXaxis().FindBin(startTime) # No check is done on the start bin, assume within 2017-2018

    # If there is enough time bins in the original axion histogram
    if endTimeHist > endTime:
        endBin = hday.GetXaxis().FindBin(endTime)
    # If the end time extends beyond 2018, need to extend the histogram
    else:
        bRepeat = True
        endBin = nHistTime
        endBin2 = hday.GetXaxis().FindBin(endTime - yearInSec)

    AxionList = []

    # Because this is finding the bin edges, need to go beyond the final bin
    # to get the upper end of the last time bin
    timeBinLowEdge = [hday.GetXaxis().GetBinLowEdge(tBin)
                        for tBin in range(hday.GetNbinsX())
                        if tBin >= startBin and tBin <= endBin+1]
    # Here is if we extend the PDF beyond 1 year
    if debug and bRepeat:
        print('Axion PDF Gap Extension (should be binsize):', timeBinLowEdge[-1], hday.GetXaxis().GetBinLowEdge(0)+yearInSec)
    if bRepeat:
        timeBinLowEdge.extend([hday.GetXaxis().GetBinLowEdge(tBin)+yearInSec
                            for tBin in range(hday.GetNbinsX())
                            if tBin >= 0 and tBin <= endBin2+1]) # nBins+1 overlaps with the 0 and 1st bin

    # Total bins is 1 less than the array length of timeBinLowEdge
    nTimeBins = len(timeBinLowEdge)-1

    # Loop through energy bins, grabbing the bin content of each bin and filling into matrix
    # For the loop, start with 1 and go through NbinsY+1 -- this allows for a range up to 20 keV in the fit
    for energyBin in range(1, hday.GetNbinsY()+1):
        # Apply energy cuts on the low edge of the bin -- rounding here is necessary because of ROOT's stupid precision
        if round(hday.GetYaxis().GetBinLowEdge(energyBin),1) < energyThresh: continue
        if round(hday.GetYaxis().GetBinLowEdge(energyBin),1) > energyThreshMax-0.1: continue
        # print("Axion Energy Bin: ", round(hday.GetYaxis().GetBinLowEdge(energyBin),1))
        # Here we want to scale the PDF properly to maintain the units
        # The histogram has a differential rate of c/keV/kg/day
        # We need to divide by 90 because it is averaged across 90 angles (Wenqin's fault)
        # Then multiply by the energy bin size (to get to keV) and divide by 5minutes/1 day for time
        AxionList.append([hday.GetBinContent(timeBin, energyBin) #* binSize * 5./(60*24)/90.
                            for timeBin in range(1, hday.GetNbinsX())
                            if timeBin >= startBin and timeBin <= endBin])
        # Extend Axion histogram
        if bRepeat:
            AxionList[len(AxionList)-1].extend([hday.GetBinContent(timeBin, energyBin) #*binSize * 5./(60*24)/90.
                                for timeBin in range(1, hday.GetNbinsX())
                                if timeBin <= endBin2+1])

    fAxion.Close()
    AxionArr = np.array(AxionList)
    if debug:
        print('Axion Shape', AxionArr.shape, nTimeBins)

    # Efficiency * Exposure mask
    eeMask, expMask, effMask = generateEfficiencyMask(np.array(timeBinLowEdge), energyBins, AxionArr.shape, dsList=dsList, debug=debug, effMode=effMode)

    # Multiply the Axion PDF by the Efficiency*Exposure mask
    NormAxionArr = AxionArr*eeMask
    effNorm = eeMask/expMask

    # Find all columns of full zeros (PDFs don't exist)
    removeMask = np.where(np.any(NormAxionArr, axis=0)==False)[0]

    # Get stuff into counts
    binScaling = (binSize*5./(60*24)/90)

    print('Bin Scaling: {}'.format(binScaling))
    print('Efficiency integral: {}, shape: {}'.format(effMask.sum(axis=0), effMask.sum(axis=0).shape ))
    print('Axion Integral: {} (before exposure) -- {} (after efficiency/before exposure) -- {} (before efficiency) --- {} (after efficiency)'.format(AxionArr.sum()*binScaling, (AxionArr*effMask).sum()*binScaling, (AxionArr*expMask).sum()*binScaling, NormAxionArr.sum()*binScaling))
    print('Axion Integral Ratio: {} -- max: {} -- min {}'.format(((AxionArr*effMask)/AxionArr),np.amax((AxionArr*effMask)/AxionArr), np.amin((AxionArr*effMask)/AxionArr)))
    print('Test Integral: {}'.format((AxionArr*effMask).sum()/AxionArr.sum()))

    return NormAxionArr, nTimeBins, np.array(timeBinLowEdge), removeMask, eeMask, effNorm


def drawFinalSpectra(pdfDict, trace):
    """
        Draws final spectra using mean of posterior
        Also draws time spectra as well
    """
    from matplotlib.ticker import FormatStrFormatter

    rebinSize = 0
    # Here we're rebinning the time spectrum, this is purely for visualization
    # Choose first valid rebin value that will evenly divide the number of time bins
    for rb in range(15, 2001):
        if pdfDict['Data'].shape[1]%rb == 0:
            rebinSize = rb
            break
    # This is done automatically! Don't mess with it
    print('Data Shape Before Rebinning:', pdfDict['Data'].shape[0], pdfDict['Data'].shape[1])
    print('Rebinning time axis by {}, the bin size will be {} minutes'.format(rebinSize, 5*rebinSize))

    rebinShape = (pdfDict['Data'].shape[0], pdfDict['Data'].shape[1]/rebinSize)
    # Convert unix time to date with seconds
    timeBins = np.array(pdfDict['Time'][:-1], dtype='datetime64[s]')
    # Rebin timing bins (here we only select using stepsize = rebinSize rather than summing the times)
    timeBins = timeBins[::rebinSize]

    figSpec, axSpec = plt.subplots(figsize=(10,7))
    figTime, axTime = plt.subplots(figsize=(10,7))

    DataSpec = pdfDict['Data'].sum(axis=1)
    DataSpecErr = np.sqrt(DataSpec)
    # DataSpec2 = rebin(pdfDict['Data'].sum(axis=1), size=2, bPrint=True)
    # DataSpecErr2 = np.sqrt(DataSpec2)
    TimeSpec = rebin(pdfDict['Data'].sum(axis=0), size=rebinSize)
    TimeSpecErr = np.sqrt(TimeSpec)
    # print('Data', TimeSpec, np.amax(TimeSpec), np.amin(TimeSpec))
    axSpec.errorbar(energyBins[:-1], DataSpec, yerr=DataSpecErr, color='black', fmt='o', label='Data') # 0.1 keV
    # axSpec.errorbar(energyBins[:-1:2], DataSpec2, yerr=DataSpecErr2, color='black', fmt='o', label='Data') # 0.2 keV
    axTime.errorbar(timeBins, TimeSpec, yerr=TimeSpecErr, color='black', fmt='o', label='Data')

    print('rebinShape size:', rebinShape)
    print('TimeSpec size:', TimeSpec.shape)
    print('timeBins size:', timeBins.shape)
    # print('Data:', DataSpec.shape, DataSpec2.shape, len(energyBins[:-1:2]))
    # print(DataSpec2, DataSpecErr2)
    print(DataSpec.sum(), TimeSpec.sum())
    # return

    # Create a dummy for the total spectrum
    totalSpec = np.zeros(DataSpec.shape) # Dummy for energy spectrum
    # totalSpec = np.zeros(DataSpec2.shape) # Dummy for energy spectrum
    totalTime = np.zeros(TimeSpec.shape) # Dummy for time spectrum
    totalModel2D = np.zeros((int(rebinShape[0]), int(rebinShape[1]))) # Dummy for energy vs time spectrum

    # Take all the pdfs and flatten time and energy axes
    for key in pdfDict.keys():
        if key == 'Data' or key == 'Time' or key == 'Energy':
            continue
        tempSpec = trace[key].mean()*pdfDict[key].sum(axis=1)
        # tempSpec = trace[key].mean()*rebin(pdfDict[key].sum(axis=1), size=2)
        tempTime = trace[key].mean()*rebin(pdfDict[key].sum(axis=0), size=rebinSize)

        axSpec.plot(energyBins[:-1], tempSpec, label=key, linestyle='--')
        # axSpec.plot(energyBins[:-1:2], tempSpec, label=key, linestyle='--')
        axTime.plot(timeBins, tempTime, label=key, linestyle='--')
        totalSpec += tempSpec
        totalTime += tempTime

    axSpec.plot(energyBins[:-1], totalSpec, c='r', lw='3', label='Total Model')
    # axSpec.plot(energyBins[:-1:2], totalSpec, c='r', lw='3', label='Total Model')
    axSpec.set_title('DS5b, 5c, 6a Enriched')
    axSpec.set_ylabel('Counts/0.1 keV')
    axSpec.set_xlabel('Energy (keV)')
    axSpec.legend()
    plt.tight_layout()
    # figSpec.savefig(os.environ['LATDIR'] + '/plots/AAA_5_20_Spectrum.png')

    # Labels for plots only -- only label 5 bins
    timeLabels = np.linspace(timeBins[0], timeBins[-1], 5, dtype='datetime64[s]')
    timeLabels = np.array(timeLabels, dtype='datetime64[D]')
    axTime.plot(timeBins, totalTime, c='r', lw='3', label='Total Model')
    axTime.set_title('DS5b, 5c, 6a Enriched')
    axTime.set_ylabel('Counts/({} min)'.format(5*rebinSize))
    # axTime.set_xlabel('Time')
    axTime.locator_params(axis='x', nbins=len(timeLabels)) # Draw fewer ticks
    axTime.set_xticklabels(timeLabels, rotation=30)
    axTime.legend()
    plt.tight_layout()
    # figTime.savefig(os.environ['LATDIR'] + '/plots/AAA_5_20_TimeSpec.png')

    # Loop again -- Make sure I don't mess up the other plots
    # The counts vs time plot with binning is very finnicky when rebinning other stuff
    for key in pdfDict.keys():
        if key == 'Data' or key == 'Time':
            continue
        # Make 2D histograms
        temp2D = np.array(trace[key].mean()*pdfDict[key], copy=True)
        totalModel2D += bin_ndarray(temp2D, rebinShape)

    # Make 2D plot of Energy vs Time
    fig2D, (ax2D1, ax2D2) = plt.subplots(ncols=2, figsize=(15,7))
    Data2D = np.array(pdfDict['Data'], copy=True)
    Data2D =  bin_ndarray(Data2D, rebinShape)
    Difference = Data2D - totalModel2D
    print('Total Integral of Data:', Data2D.sum())
    print('Total Integral of Model:', totalModel2D.sum())
    print('Total Integral of Difference:', Difference.sum())

    sns.heatmap(Data2D,  cmap="YlGnBu", ax=ax2D1)
    ax2D1.invert_yaxis()
    ax2D1.set_title('DS5b, 5c, 6a Enriched (Data)')
    ax2D1.set_ylabel('Energy (keV)')
    ax2D1.locator_params(axis='x', nbins=len(timeLabels)) # Draw fewer ticks
    ax2D1.locator_params(axis='y', nbins=len(energyBins[::20])) # Draw fewer ticks
    ax2D1.set_xticklabels(timeLabels, rotation=30)
    ax2D1.set_yticklabels([round(x,1) for x in energyBins[:-1:20]], rotation=0)
    sns.heatmap(totalModel2D, cmap="YlGnBu", ax=ax2D2)
    ax2D2.invert_yaxis()
    ax2D2.set_title('Total Model')
    ax2D2.locator_params(axis='x', nbins=len(timeLabels)) # Draw fewer ticks
    ax2D2.set_xticklabels(timeLabels, rotation=30)
    ax2D2.set_yticklabels([])
    plt.tight_layout()

    # figSpec.savefig('{}/AxionFitResult_Spectrum_Open.png'.format(os.environ['LATDIR']+'/data/MCMC'))
    # figTime.savefig('{}/AxionFitResult_Time_Open.png'.format(os.environ['LATDIR']+'/data/MCMC'))
    # fig2D.savefig('{}/AxionFitResult_2D_Open.png'.format(os.environ['LATDIR']+'/data/MCMC'))
    plt.show()


def drawPDFs(pdfArrDict, removeMask, timeBinLowEdge):
    """
        Helper function to draw PDFs
        pdfArrDict contains all of the PDFs as well as the Data
    """
    timeLabels1 = np.linspace(timeBinLowEdge[0], timeBinLowEdge[-1], 5, dtype='datetime64[s]')
    timeLabels1 = np.array(timeLabels1, dtype='datetime64[D]')

    fig1, ax1 = plt.subplots(ncols=3, nrows=2, figsize=(15,8))
    sns.heatmap(pdfArrDict['Data'], ax=ax1[0,0])
    ax1[0,0].set_title('DS5 Data')
    ax1[0,0].set_ylabel('Energy')
    ax1[0,0].set_xlabel('Time')
    ax1[0,0].invert_yaxis()
    ax1[0,0].locator_params(axis='y', nbins=len(energyBins[::50]))
    ax1[0,0].set_yticklabels([round(x,1) for x in energyBins[::50]], rotation=0)
    ax1[0,0].set_xticklabels([])

    sns.heatmap(pdfArrDict['Tritium'], ax=ax1[0,1])
    ax1[0,1].set_title('Tritium PDF')
    ax1[0,1].set_xlabel('Time')
    ax1[0,1].invert_yaxis()
    ax1[0,1].set_yticklabels([])
    ax1[0,1].set_xticklabels([])

    sns.heatmap(pdfArrDict['Axion'], ax=ax1[0,2])
    ax1[0,2].set_title('Axion PDF')
    ax1[0,2].set_xlabel('Time')
    ax1[0,2].invert_yaxis()
    ax1[0,2].set_yticklabels([])
    ax1[0,2].set_xticklabels([])

    sns.heatmap(pdfArrDict['Fe55'], ax=ax1[1,0])
    ax1[1,0].set_title('Fe55 PDF')
    ax1[1,0].set_ylabel('Energy')
    ax1[1,0].set_xlabel('Time')
    ax1[1,0].invert_yaxis()
    ax1[1,0].locator_params(axis='y', nbins=len(energyBins[::50]))
    ax1[1,0].locator_params(axis='x', nbins=5)
    ax1[1,0].set_yticklabels([round(x,1) for x in energyBins[::50]], rotation=0)
    ax1[1,0].set_xticklabels(timeLabels, rotation=30)

    sns.heatmap(pdfArrDict['Zn65'], ax=ax1[1,1])
    ax1[1,1].set_title('Zn65 PDF')
    ax1[1,1].set_xlabel('Time')
    ax1[1,1].invert_yaxis()
    ax1[1,1].locator_params(axis='x', nbins=5)
    ax1[1,1].set_yticklabels([])
    ax1[1,1].set_xticklabels(timeLabels, rotation=30)

    sns.heatmap(pdfArrDict['Ge68'], ax=ax1[1,2])
    ax1[1,2].set_title('Ge68 PDF')
    ax1[1,2].set_xlabel('Time')
    ax1[1,2].invert_yaxis()
    ax1[1,2].locator_params(axis='x', nbins=5)
    ax1[1,2].set_yticklabels([])
    ax1[1,2].set_xticklabels(timeLabels, rotation=30)
    plt.tight_layout()
    # fig1.savefig('{}/BraggFitPDF_New.png'.format(os.environ['LATDIR']+'/plots/Axion'))

    axionRm = np.delete(pdfArrDict['Axion'], removeMask, axis=1)
    timeRm = np.delete(timeBinLowEdge, removeMask)

    timeLabels2 = np.linspace(timeRm[0], timeRm[-1], 5, dtype='datetime64[s]')
    timeLabels2 = np.array(timeLabels2, dtype='datetime64[D]')

    fig2, (ax2, ax3)= plt.subplots(ncols=2, figsize=(12,6))
    sns.heatmap(pdfArrDict['Axion'], cmap="YlGnBu", ax=ax2)
    ax2.set_title('DS5B AAA PDF (Before Deadtime Removal)')
    # ax2.set_xlabel('Time')
    ax2.set_ylabel('Energy (keV)')
    ax2.invert_yaxis()
    ax2.locator_params(axis='y', nbins=len(energyBins[::25]))
    ax2.locator_params(axis='x', nbins=5)
    ax2.set_yticklabels([round(x,1) for x in energyBins[::25]], rotation=0)
    ax2.set_xticklabels(timeLabels1, rotation=30)

    sns.heatmap(axionRm, cmap="YlGnBu", ax=ax3)
    ax3.set_title('DS5B AAA PDF (After Deadtime Removal)')
    # ax3.set_xlabel('Time')
    ax3.set_ylabel('Energy (keV)')
    ax3.invert_yaxis()
    ax3.locator_params(axis='y', nbins=len(energyBins[::25]))
    ax3.locator_params(axis='x', nbins=5)
    ax3.set_yticklabels([round(x,1) for x in energyBins[::25]], rotation=0)
    ax3.set_xticklabels(timeLabels2, rotation=30)

    plt.tight_layout()
    # fig2.savefig('{}/AxionPDF_Deadtime.png'.format(os.environ['LATDIR']+'/plots/Axion'))
    plt.show()


def generateEfficiencyMask(timeBinLowEdge, energyBins, AxionShape, dsList=None, debug=False, effMode=None):
    """
        Creates mask for efficiency*exposure, will return a number between 0 -- exposure (the weight for each bin)
        for every time and energy bin. Apply this mask to the PDFs to delete columns where the detectors are off

        returns:
            eeMask -- Exposure*Efficnecy mask
            expMask -- Exposure mask only (no efficiency applied)
            effMask -- Efficiency mask only (no exposure applied)
    """

    if not dsList:
        print('Error: DataSet List not specified for efficiency calculation')
        return
    if effMode:
        print('Using {} mode for cut efficiency'.format(effMode))

    import ROOT
    eeMask = np.ones(AxionShape)
    expMask = np.ones(AxionShape) # Exposure only
    effMask = np.ones(AxionShape) # Efficiency only
    if debug:
        print('Start, End Time:', startTime, endTime)
        print('Time Bins:', timeBinLowEdge[0], timeBinLowEdge[-1])
        print('Difference: {} -- {}'.format(startTime - timeBinLowEdge[0], timeBinLowEdge[-1] - endTime))

    timeBinSize = (timeBinLowEdge[1] - timeBinLowEdge[0])
    if debug:
        print('Timing Bin size: ', timeBinSize)

    fEff = ROOT.TFile(os.environ['LATDIR']+'/data/lat-expo-efficiency_Combined.root')
    # Dummy variable for Start time of DS5b to calculate gaps in livetime
    prevEnd = startTime

    for ds in dsList:
        # First build efficiency array -- assumes efficiency for all dataset is constant
        # Also combines efficiency of natural and enriched detectors
        if not effMode:
            hEff = fEff.Get('hDS{}_Norm_Enr'.format(ds))
        else:
            hEff = fEff.Get('hDS{}_Norm_Enr_{}'.format(ds, effMode))
        effArr = np.ones(len(energyBins)-1)
        for eIdx, eBin in enumerate(energyBins[:-1]):
            effArr[eIdx] = hEff.GetBinContent(hEff.FindBin(eBin+binSize/2.))

        if debug:
            print('DS{} Efficiency Array:'.format(ds), effArr, effArr.shape)

        # Create a new map for all runs with the exposure of each run
        # Builds up dictionary of run:exposure for every dataset
        if ds in ['5A', '5B', '5C']:
            dsNum = 5
        else:
            dsNum = int(ds)
        runMatrix = {}
        bkgRanges = bkg.getRanges(dsNum)
        chList = det.getGoodChanList(dsNum, detType='Enr')
        # Loop through channels to build up total runLists
        for ch in chList:
            cpd = int(det.getChanCPD(dsNum,ch))
            inFile = os.environ['LATDIR']+'/data/expLists/DS{}_C{}P{}D{}.txt'.format(ds, *str(cpd))
            try:
                with open(inFile) as f:
                    # Loop through file,
                    for line in f:
                        currArr = np.array(line.split(','), dtype=np.float)
                        # print(currArr)
                        runMatrix.setdefault(int(currArr[0]), 0.)
                        runMatrix[int(currArr[0])] += currArr[1]
            except:
                print("File {} doesn't exist! Skipping".format(inFile))

        if debug:
            expTot = 0
            for idx in runMatrix:
                expTot += runMatrix[idx]
            print('DS{} Total Exposure: {}'.format(ds, expTot))

        # Now loop through the full dataset runList to add in unixtimes and eliminate gaps
        inFile2 = os.environ['LATDIR']+'/data/MCMC/DS{}_RunTimeList.txt'.format(ds)
        # gapMatrix = {}
        totalGap = 0
        with open(inFile2) as f2:
            for line in f2:
                cArr = np.array(line.split(','), dtype=np.int)
                run, sTime, eTime = cArr[0], cArr[1], cArr[2]

                # If the timing between runs is 1 second, it's just a rounding uncertainty
                if sTime - prevEnd < 2: prevEnd = sTime

                # Calculate where the time bins are for the particular run
                prevIndex = np.where(timeBinLowEdge <= prevEnd)[0][-1] # Time bin of end of previous run
                startIndex = np.where(timeBinLowEdge <= sTime)[0][-1] # Time bin of beginning of current run
                stopIndex = np.where(timeBinLowEdge >= eTime)[0][0] # Time bin of end of current run

                # If there is a gap between runs for livetime between runs
                # set scaling to 0 during that time period to remove
                if sTime != prevEnd:
                    # Fills in gapMatrix with livetime gaps, run: end, start, gap size (s)
                    # gapMatrix.setdefault(run, [prevEnd, sTime, sTime-prevEnd])
                    # Set the efficiency for those time bins to be zero
                    if prevIndex == startIndex:
                        eeMask[:, prevIndex] = 0
                        expMask[:, prevIndex] = 0
                        effMask[:, prevIndex] = 0
                    else:
                        eeMask[:, prevIndex:startIndex] = 0
                        expMask[:, prevIndex:startIndex] = 0
                        effMask[:, prevIndex:startIndex] = 0
                    # If there is a gap between runs that's less than the 5 minute bin-size
                    # Not dealing with this for now... < 0.5 kg-day exposure difference
                    if sTime - prevEnd < 300 and debug:
                        print('Time gap for run {} less than 5 min: {} ({}-{})'.format(run,sTime-prevEnd, prevIndex, startIndex))
                    totalGap += sTime-prevEnd
                    # Update end period unixtime to the beginning of the run
                    prevEnd = sTime

                # Otherwise, multiply by efficiency * exposure
                # Update end period unixtime if there is no gap
                # Exposure is split up by number of bins it fills here (nTimeBins * nEnergyBins)
                if run in runMatrix and runMatrix[run] != 0:
                    if startIndex == stopIndex:
                        eeMask[:, startIndex] = effArr*runMatrix[run]/(len(energyBins)-1)
                        expMask[:, startIndex] = runMatrix[run]/(len(energyBins)-1)
                        effMask[:, startIndex] = effArr
                    else:
                        # This part is really long and complicated, basically need to copy (repeat) the effArr to match the size we want to fill inside eeMask, the reshape is just to make the array sizes match

                        # If the shape doesn't match, it's because we're at the last time bin, shrink
                        if eeMask[:, startIndex:stopIndex].shape[1] != stopIndex-startIndex:
                            print('Matching off run {} ({}-{})'.format(run,startIndex, stopIndex) )
                            effMat = np.repeat(effArr,stopIndex-startIndex-1).reshape(len(effArr), stopIndex-startIndex-1)
                        else:
                            effMat = np.repeat(effArr,stopIndex-startIndex).reshape(len(effArr), stopIndex-startIndex)


                        eeMask[:, startIndex:stopIndex] = effMat*runMatrix[run]/(stopIndex-startIndex-1)/(len(energyBins)-1)
                        expMask[:, startIndex:stopIndex] = runMatrix[run]/(stopIndex-startIndex-1)/(len(energyBins)-1) # No efficiency
                        effMask[:, startIndex:stopIndex] = effMat # No exposure
                else:
                    eeMask[:, startIndex:stopIndex] = 0
                    expMask[:, startIndex:stopIndex] = 0
                    effMask[:, startIndex:stopIndex] = 0
                    print('Run {} setting to 0'.format(run))
                    print(eeMask[100, startIndex:stopIndex])

                prevEnd = eTime

        if debug:
            print(eeMask)
            print('DS{} -- Exposure: {} -- Eff Corr Exposure: {} (kg-day)'.format(ds, expMask.sum(), eeMask.sum())) # Sum over one axis at a high energy, should get exposure?
            print('Total time gap for DS{} -- {}'.format(ds, totalGap)) # Residual time gap
            # from mpl_toolkits.mplot3d import Axes3D
            # feff, aeff = plt.subplots(figsize=(8,6))
            # sns.heatmap(eeMask, ax=aeff)
            # feff = plt.figure()
            # aeff = feff.gca(projection='3d')
            # xmesh, ymesh = np.meshgrid(energyBins[:-1], timeBinLowEdge[:-1])
            # aeff.plot_surface(xmesh, ymesh, eeMask.T)
            # plt.show()

    return eeMask, expMask, effMask


def getSigma(energy, dsNum='6a'):

    p0,p1,p2 = 0., 0., 0.
    if dsNum=='0':
        p0 = 0.147; p1=0.0173; p2=0.0003
    elif dsNum=='1':
        p0 = 0.136; p1=0.0174; p2=0.00028
    elif dsNum=='3':
        p0 = 0.162; p1=0.0172; p2=0.000297
    elif dsNum=='4':
        p0 = 0.218; p1=0.015; p2=0.00035
    elif dsNum=='5a':
        p0 = 0.2592; p1=0.2057; p2=0.00030863
    elif dsNum=='5b':
        p0 = 0.1815; p1=0.01705; p2=0.00031527
    elif dsNum=='5c':
        p0 = 0.1361; p1=0.0174; p2=0.0002829
    elif dsNum=='6a':
        p0 = 0.1314; p1=0.01728; p2=0.0002742

    # Use DS6a as default
    else:
        p0 = 0.1314; p1=0.01728; p2=0.0002742

    return np.sqrt(p0*p0 + p1*p1*energy + p2*p2*energy*energy)


def generateGaus(energyList, meansDict, dsNum):
    """
        Generates Gaussian PDFs for X-ray peaks based off of energy array
    """
    outDict = {}
    for name, mean in meansDict.items():
        outDict.setdefault(name, norm.pdf(energyList, loc=mean, scale=getSigma(mean, dsNum=dsNum)))
    return outDict


def modelDiagnostics(model, pdfArrDict=None, backendDir=None, effNorm=None):
    """
        Loads saved traces, draws plots, performs diagnostics on models
        Brian's notes on diagnostics:

        Gelman-Rubin (Rhat) -- tests for lack of convergence by comparing the variance between multiple chains to the variance within one chain.
        If the chains are converged, the variance should be the same (ratio should be 1). Needs multiple chains with differing starting points (different seeds!).

        Auto-correlation -- measures how linearly dependent the current value of the chain is to past values. It essentially tells us how much information is avaliable, if there is a lot of correlation, there's less information than sampled from stationary distribution. The value is between -1 and 1, the lag where the autocorr is ~0 means the number of samples where the values aren't autocorrelated. The amount of autocorrelation determines the effective sample size (n_effective), for setting a 95% credible interval. This number needs to be computed for every parameter in the model

        Effective sample size (effective_n, ESS) -- is an estimate on the effective number of independent samples taking into account auto-correlation (ESS = N_samples/(1 + 2*sum_i(autocorr_i))). This number is not equal to the actual number of samples drawn as the samples are often correlated (since this is a Markov Chain!), this number is also different per parameter in the model.

        Monte Carlo Standard Error (MCSE or mc_error) -- MCSE^2 = Var(f)/ESS

        Highest Posterior Density (HPD) -- The HPD is an estimator that calculates the minimum width of the Bayesian credible interval (BCI).
        If the posterior density is multimodal, the HPD does not result in an interval estimate (thankfully it's not for our case)

    """
    import corner

    # Load from files
    print('Loading Backend from {}'.format(backendDir))
    trace = pm.backends.text.load(backendDir, model)
    # Get summary of trace as a dataframe
    traceBurn = trace[nBurn:]
    dfSum = pm.summary(traceBurn, alpha=0.1, include_transformed=True)
    print(dfSum.head(20))
    AxionIntegral = pdfArrDict['Axion'].sum()


    # This gets printed out from "convertaxionPDF", it's the amount I artificially scale the axion PDF integral
    print('Axion Shape', pdfArrDict['Axion'].shape)
    print('Axion PDF Integral:', AxionIntegral)

    # Loop through and calculate 95% CI for all parameters
    print('Calculating 95% Credible Intervals: ')
    print('Par \t  Mean \t hpd_5 \t hpd_95')
    # Here all the PDFs are normalized in their own way (from model construction section)
    # therefore multiplying the PDF integral by the 95% CI should give number of counts
    print('Efficiency Scaling factor (summed over 1 axis):',1./effNorm.sum(axis=1)[0], 1./effNorm.sum(axis=1)[1])

    bkgIntegral = 0
    axionIntegral = 0
    axionMeanIntegral = 0
    for par in dfSum.index:
        if '_log_' in par:
            continue

        # Integrate PDF
        parInt = (pdfArrDict[par]/effNorm).sum()
        parRaw = pdfArrDict[par].sum()
        # parInt = (pdfArrDict[par]).sum()

        # Multiply integral by interval
        print('{}    {:.2f}    {:.2f}    {:.2f}  {:.2f}'.format(par, parRaw*dfSum.loc[par, 'mean'], parRaw*dfSum.loc[par, 'hpd_5'], parRaw*dfSum.loc[par, 'hpd_95'], parRaw))
        if par == 'Axion':
            # print('Axion: Scaled Integral: {} \t Unscaled Integral: {}'.format(parRaw*dfSum.loc[par, 'hpd_95'], parRaw))
            # print('95% CI limit: {}'.format(parRaw*dfSum.loc[par, 'hpd_95']/parRaw))
            axionMeanIntegral = parRaw*dfSum.loc[par, 'mean']
            axionIntegral = parRaw*dfSum.loc[par, 'hpd_95']
        else:
            bkgIntegral += parRaw*dfSum.loc[par, 'mean']

    print('Data Integral: {}'.format(pdfArrDict['Data'].sum()))
    print('Total Background Integral (Mean): {}'.format(bkgIntegral))
    print('Axion Integral (Mean): {}'.format(axionMeanIntegral))
    print('Axion Integral (95% C.I.): {}'.format(axionIntegral))

    cached = [(var, var.logp_elemwise) for var in model.observed_RVs]

    def logp_vals_point(pt):
        logp_vals = []
        for var, logp in cached:
            logp = logp(pt)
            if var.missing_values:
                logp = logp[~var.observations.mask]
            logp_vals.append(logp.ravel())
        return np.concatenate(logp_vals)

    # points = traceBurn.points()
    # pt = []
    # idx = 0
    # logpArr = np.zeros(5000)
    # for pt in points:
    #     if idx >= 5000:
    #         break
    #     atemp = logp_vals_point(pt)
    #     logpArr[idx] = atemp.sum()
    #     idx += 1
        # print(atemp, atemp.shape)

    # ptMean = {par: dfSum.loc[par, 'mean'] for par in dfSum.index}
    # aMean = logp_vals_point(ptMean).sum()

    # xArr = np.arange(0, 5000, 1)
    # fig0, (ax01, ax02) = plt.subplots(ncols=2, figsize=(12,5))
    # sns.distplot(logpArr, bins=50, kde=False, ax=ax01)
    # ax01.axvline(aMean)
    # ax01.set_xlabel('NLL')
    # ax02.plot(xArr, logpArr)
    # ax02.set_xlabel('Sample')
    # ax02.set_ylabel('NLL')

    # Draw a posterior with the X-axis converted into Axion Counts
    fig11, ax11 = plt.subplots(figsize=(10,6))
    axionTraceArr = trace.get_values('Axion', burn=1000)
    hist, edges = np.histogram(axionTraceArr, bins=5000)
    integralArr = np.cumsum(hist)
    idxPass = np.where(integralArr >= 0.9*hist.sum())[0][0]
    idxPass2 = np.where(integralArr >= 0.95*hist.sum())[0][0]
    print(idxPass, edges[idxPass], pdfArrDict['Axion'].sum()*edges[idxPass])
    print(idxPass2, edges[idxPass2], pdfArrDict['Axion'].sum()*edges[idxPass2])
    print('Axion Mean:', pdfArrDict['Axion'].sum()*dfSum.loc['Axion', 'mean'])
    print('Axion MC sterr:', pdfArrDict['Axion'].sum()*dfSum.loc['Axion', 'mc_error'])

    axionTraceArrScaled = pdfArrDict['Axion'].sum()*axionTraceArr
    sns.distplot(axionTraceArrScaled, kde=False, bins=50, ax=ax11)
    ax11.axvline(pdfArrDict['Axion'].sum()*edges[idxPass])
    ax11.set_title('Axion Posterior')
    ax11.set_xlabel('Axion Counts')
    # fig11.savefig(os.environ['LATDIR']+'/plots/AAA_5_20_AxionPosterior.png')

    # Debate: when converting to lambda,
    # axionLambda = 3600*24*dfSum.loc['Axion', 'hpd_95']/AxionNorm
    # print("Axion Lambda: {} --- g_agg (E-8 GeV): {}".format(axionLambda, np.power(axionLambda, 0.25)))

    # Plot Traces
    fig1, ax1 = plt.subplots(nrows=6,ncols=2,figsize=(16, 18))
    pm.traceplot(traceBurn, ax=ax1)
    # fig1.savefig(os.environ['LATDIR']+'/plots/AAA_5_20_Trace.png')

    # Plot forestplots
    fig2, ax2 = plt.subplots(figsize=(10, 7))
    gs0 = pm.forestplot(traceBurn, varnames=['Axion', 'Tritium', 'Bkg', 'Ge68', 'Fe55', 'Zn65'], alpha=0.1)
    gr0_plot = plt.subplot(gs0[1])
    gr0_plot.set_xlim(0.99, 1.01)
    # fig2.savefig(os.environ['LATDIR']+'/plots/AAA_5_20_Forest.png')

    # fig3, ax3 = plt.subplots(figsize=(10, 7))
    # gs = pm.forestplot(traceBurn, varnames=['Axion'], alpha=0.1)
    # gr_plot = plt.subplot(gs[1])
    # gr_plot.set_xlim(0.99, 1.01)
    # fig3.savefig('AxionFit_5_20_Forest_Axion.png')

    # Plot Posteriors
    # fig4, ax4 = plt.subplots(ncols=3, nrows=3, figsize=(12, 12))
    # pm.plot_posterior(traceBurn, alpha_level=0.1, round_to=5)
    # fig4.savefig(os.environ['LATDIR']+'/plots/AAA_5_20_Posterior.png')


    # print(np.array([traceBurn['Axion'], traceBurn['Tritium']]).T.shape)
    # cornerArr = np.array([traceBurn['Axion'], traceBurn['Tritium'], traceBurn['Bkg'], traceBurn['Ge68'], traceBurn['Fe55'], traceBurn['Zn65']]).T
    # figure = corner.corner(cornerArr,  quantiles=[0.05, 0.95], show_titles=True, title_fmt=".5f", labels=['Axion', 'Tritium', 'Bkg', 'Ge68', 'Fe55', 'Zn65'])

    # Format is hpdDict[nChain][Parameter] = [lower, upper]
    # hpdDict = pm.hpd(trace[nBurn:], alpha=0.1)
    # print(hpdDict[1])


    # fig5, ax5 = plt.subplots(nrows=8, ncols=1, figsize=(15, 7))
    # ax5 = pm.autocorrplot(trace[nBurn:], varnames=['Axion'], max_lag=1000)
    # for i in range(len(ax5)):
        # print(ax5[i], type(ax5[i]))
    # print(pm.autocorr(trace[nBurn:]))
    # pm.densityplot(trace[nBurn:])
    plt.show()


def posteriorChecks(model, pdfArrDict, backendDir=None):
    """
        These functions aren't really used and are for more diagnostics
        Note that with such a large model, running some of these for more than a couple hundred samples will crash the code due to memory purposes.

        In the Bayesian perspective we define a model that we assume is true and we fit it to the data, but what if the model is wrong? Posterior Predictive Checks performs another check on the model by requiring reasonable predictive performance using the model. In standard machine learning, holding out a dataset gives some check of predictive performance but it is only one point, not an ensemble. Cross-validation doesn't give predictive performance, it gives a sense of partitioning performance.

        Posterior Predictive Checks:
            - Simulates replicating data under the fitted model and then comparing these to the observed data
            - This checks for systematic discrepancies between real and simulated data
            - Effectively computes the probability of psuedo-data given the data (comes for free with the MCMC): p(D'|D) = Int( P(D'|w) * P(w|D) * dw )

        Widely-applicable Information Criterion (WAIC):
            - Fully Bayesian criterion for estimating out-of-sample expectation, using the computed log pointwise posterior predictive density (LPPD) and correcting for the effective number of parameters to adjust for overfitting.
            - This is primarilly for the comparison between different models

        Leave-one-out Cross-validation (LOO):
            - Estimate of the out-of-sample predictive fit. In cross-validation, the data are repeatedly partitioned into training and holdout sets, iteratively fitting the model with the former and evaluating the fit with the holdout data. PyMC's implementation of LOO is using Pareto-smoothed importance sampling, it provides and estimate of point-wise out-of-sample prediction accuracy

        Note: out-of-sample is data not used for the fit (ie: making a prediction)

    """

    if not backendDir:
        print('No backend directory specified!')
        return

    # Load from files
    print('Loading Backend from {}'.format(backendDir))
    trace = pm.backends.text.load(backendDir, model)
    traceBurn = trace[15950:]
    # return
    model_waic = pm.waic(traceBurn, model, progressbar=True)
    print(model_waic)
    print(model_waic.WAIC)

    # ppc = pm.sample_ppc(trace[nBurn:], samples=1, size=664, model=model)
    # print(np.asarray(ppc['L'].shape), ppc.keys())
    # _, axppc = plt.subplots(figsize=(12, 7))
    # axppc.hist([n.mean() for n in ppc['L']], bins=19, alpha=0.5)
    # axppc.set(title='Posterior predictive for L', xlabel='L(x)', ylabel='Frequency');
    # plt.show()
    # Nominal: WAIC_r(WAIC=16703.9407014325, WAIC_se=512.1816033806425, p_WAIC=5.251383203562369, var_warn=0)
    # Lo90: WAIC_r(WAIC=16704.64631945686, WAIC_se=512.0447693872607, p_WAIC=5.469910859550034, var_warn=0)
    # IS: WAIC_r(WAIC=16703.581148799192, WAIC_se=512.273151103169, p_WAIC=5.04831880451535, var_warn=0)
    # Sideband: WAIC_r(WAIC=16702.32593988035, WAIC_se=512.1507931617733, p_WAIC=4.98098150964086, var_warn=0)
    return


def rebin(inArr, size=2, bPrint=False):
    """
        Rebin an array
    """
    if not float(inArr.shape[0]/size).is_integer():
        print('Error: Input array length is not evenly divisible by rebin size')
        return
    rebinnedArr = np.zeros(inArr[::size].shape)

    for i in range(size):
        rebinnedArr += inArr[i::size]
        if bPrint:
            print(i, inArr[i::size])
    if bPrint:
        print(rebinnedArr)
    return rebinnedArr


def bin_ndarray(ndarray, new_shape, operation='sum'):
    """
    Bins an ndarray in all axes based on the target shape, by summing or
        averaging.

    Number of output dimensions must match number of input dimensions and
        new axes must divide old ones.

    NOTE: This changes the input array!
    """
    operation = operation.lower()
    if not operation in ['sum', 'mean']:
        raise ValueError("Operation not supported.")
    if ndarray.ndim != len(new_shape):
        raise ValueError("Shape mismatch: {} -> {}".format(ndarray.shape,
                                                           new_shape))
    compression_pairs = [(d, c//d) for d,c in zip(new_shape,
                                                  ndarray.shape)]
    flattened = [l for p in compression_pairs for l in p]
    ndarray = ndarray.reshape(flattened)
    for i in range(len(new_shape)):
        op = getattr(ndarray, operation)
        ndarray = op(-1*(i+1))
    return ndarray


def reduceData():
    import ROOT
    inDir, outDir = os.environ['LATDATADIR']+'/bkg/cut/final95', os.environ['LATDIR']+'/data/MCMC'
    skimTree = ROOT.TChain("skimTree")
    skimTree.Add("{}/final95_DS5B.root".format(inDir))
    skimTree.Add("{}/final95_DS5C.root".format(inDir))
    skimTree.Add("{}/final95_DS6.root".format(inDir))
    theCut = "trapENFCal > 2.4 && trapENFCal < 50"
    nPass = skimTree.Draw('trapENFCal:channel:globalTime:isEnr', theCut, 'goff')
    print ("{} events passed all cuts".format(nPass))
    nEnergy = skimTree.GetV1()
    nCh = skimTree.GetV2()
    nTime = skimTree.GetV3()
    nEnr = skimTree.GetV4()
    nChList = list(int(nCh[n]) for n in range(nPass))
    nEnergyList = list(float(nEnergy[n]) for n in range(nPass))
    nTimeList = list(int(nTime[n]) for n in range(nPass))
    nEnrList = list(int(nEnr[n]) for n in range(nPass))

    df = pd.DataFrame({"Energy":nEnergyList, "Channel":nChList, 'UnixTime':nTimeList, 'isEnr':nEnrList})
    print(df.head(10))
    df.to_hdf('{}/Bkg_Spectrum.h5'.format(outDir), 'skimTree', mode='w', format='table')


def help_options():
    """
        Prints the different options
    """
    print("BigBraggBrand Run Options:")
    print("\t -h or --help: Prints this help message")
    print("\t -backend : Changes the backend directory")
    print("\t -debug: Turns on debug mode, more output/plots")
    print("\t -diagnostic: Generates model diagnostics, only run AFTER MCMC")
    print("\t -drawPDF: Only draws the PDFs")
    print("\t -drawSpectra: Draws final spectra, MCMC chain needs to be sampled")
    print("\t -ds: Sets the energy resolution of the Axion and Gaussian PDFs, default is 5b (can also use 6a)")
    print("\t -eff: Sets efficiency curve, default is None (central), other options are Lo90, Hi90, IS, and Sideband")
    print("\t -noaxion: Builds model without Axion signal")
    print("\t -ppc: Generates posterior predictive checks, only run AFTER MCMC")
    print("\t -reduce: Saves ROOT data into a DataFrame (only run this once!)")
    print("\t -sample: Turns on MCMC sampling")
    print("\t -seed #: Changes seed number (corresponds to MCMC chain number)")
    print("\t -setbinsize: Sets axion bin size (Default/RECOMMENDED 5 minutes)")
    return


if __name__ == '__main__':
    main(sys.argv[1:])
