#!/usr/bin/env python
"""
    Big NUTS Bragg Brand fitting

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

"""

import sys, os, imp, time
import numpy as np
import pymc3 as pm
import matplotlib.pyplot as plt
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

# Define some global parameters
inDir = os.environ['LATDIR']+'/data/MCMC'
dsList = ['5B', '5C', '6']
# dsList = ['5B']
# Energy range (also defines fitting range!)
# Wenqin stops at 12 keV since the axions stop at 12 -- endpoint of trit is 18.
# Scanned from 12 to 19.8 -- makes more sense to go to at least the end of the tritium spectrum
# energyThresh, energyThreshMax, binSize = 2.4, 18, 0.1
energyThresh, energyThreshMax, binSize = 4.0, 18.0, 0.1
energyBins = np.linspace(energyThresh,energyThreshMax, round((energyThreshMax-energyThresh)/0.1)+1)

# The start and end time are for DS5b
startTime = 1485575988
# endTime = 1489762153 # DS5b

# The start and end time are for DS5c
# startTime = 1489773072
# endTime = 1494223870 # DS5c

# The start and end time are for DS5a
# startTime = 1494535909
endTime = 1523845062 # DS6a

totExposure = 2621.72447852 # DS5b, DS5c, DS6
# Fake start times
# startTime = 1483228800
# endTime = 1514678400 # End of 2017
# endTime = 1546300800 # End of 2018
yearInSec = 31536000 # Number of seconds in 1 year
nBurn = 1000 # Number of burn-in samples for the MCMC

def main(argv):

    bDebug, bDrawPDF, bDiagnostic = False, False, False
    bSample, bFull = False, True
    axionBinSize = 5 # This is in units of minutes
    pdfArrDict, pdfFlatDict = {}, {}
    # Set initial seed number here
    seedNum = 1

    if len(argv)==0:
        print("BigBraggBrand Run Options:")
        print("\t -reduce: Saves ROOT data into a DataFrame (only run this once!)")
        print("\t -seed #: Changes seed number (corresponds to MCMC chain number)")
        print("\t -debug: Turns on debug mode, more output/plots")
        print("\t -drawPDF: Only draws the PDFs")
        print("\t -diagnostic: Generates model diagnostics, only run AFTER MCMC")
        print("\t -noaxion: Builds model without Axion signal")
        print("\t -sample: Turns on MCMC sampling")
        print("\t -setbinsize: Sets axion bin size (default 5 minutes)")
        return
    for i,opt in enumerate(argv):
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
        if opt == '-drawPDF':
            bDrawPDF = True
            print('Drawing PDFs only!')
        if opt == '-diagnostic':
            bDiagnostic = True
            print('Generating Model Diagnostics, MCMC chain needs to be sampled!')
        if opt == '-noaxion':
            bFull = False
            print('Building Model without Axion PDF!')
        if opt == '-sample':
            bSample = True
            print('MCMC Sampling ON')
        if opt == '-setbinsize':
            axionBinSize = int(argv[i+1])
            print('Setting axion bin size to be {} minutes'.format(axionBinSize))

    # Build Axion PDF -- Use low edges of PDF to generate bins
    AxionArr, nTimeBins, timeBinLowEdge, removeMask, effMask, effNorm = convertaxionPDF(startTime, endTime, debug=bDebug, axionBinSize = axionBinSize)

    # Normalize Axion PDF by maximum (this is arbitrary for simplifying the sampling)
    AxionNorm = (np.amax(AxionArr))
    pdfArrDict['Axion'] = AxionArr/AxionNorm

    print('Generated Axion PDF')
    print('Axion Integral:', AxionArr.sum())
    print('Axion Normalization Constant: {}'.format(AxionNorm))

    # Load Background data and bin exactly as PDF
    df = pd.read_hdf('{}/Bkg_Spectrum.h5'.format(inDir))
    df.sort_values(by=['UnixTime'], inplace=True)
    # Select Enriched detectors only, grab numpy array of the values
    dataArr = df[['Energy', 'UnixTime']].loc[df['isEnr']==1].values

    # Bin the data into a 2D histogram -- here we use the exact time bin edges of the axion PDF to make sure everything is binned properly
    pdfArrDict['Data'], xedges, yedges = np.histogram2d(dataArr[:,0], dataArr[:,1], bins=[energyBins, timeBinLowEdge])

    # Build Tritium PDF -- need interpolation because it's 0.2 keV bins!
    # NOTE: We may want to use 0.2 keV bins, the difference may not be noticeable and the
    # NOTE: I can probably just generate this into 0.1 keV bins at some point but I've been lazy...
    dfTrit = pd.read_hdf('{}/TritSpec.h5'.format(inDir))
    tritEnergy = dfTrit['Energy'].values
    tritSpec = dfTrit['Tritium'].values
    # The tritium calculation has some slightly negative values around 0 (fluctuations), get rid of them!
    tritSpec = np.array([100*x if x >= 0. else 0. for x in tritSpec])
    tritInterp = interp1d(tritEnergy, tritSpec, fill_value='extrapolate')
    tritList = [tritInterp(x) for x in energyBins[:-1]]
    pdfArrDict['Tritium'] = np.array(tritList*nTimeBins).reshape(nTimeBins, len(tritList)).T*effMask
    pdfArrDict['Bkg'] = np.ones(pdfArrDict['Axion'].shape)*effMask

    # Generate the Gaussian PDFs
    # Fe55 is slightly shifted from 6.54 -- Perhaps should introduce as a constraint of some sort
    meansDict = {'Fe55':6.50, 'Zn65':8.98, 'Ge68':10.37}
    gausDict = generateGaus(energyBins[:-1], meansDict)
    for name, Arr in gausDict.items():
        pdfArrDict.setdefault(name, np.array(Arr.tolist()*nTimeBins).reshape(nTimeBins, len(Arr)).T)
        pdfArrDict[name] *= effMask

    # Create a matrix for energy bin center values -- this was used for efficiency sampling within the model
    # pdfArrDict.setdefault('Energy', np.array((energyBins[:-1]+binSize/2).tolist()*nTimeBins).reshape(nTimeBins, len(energyBins[:-1]) ).T)

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
    # pdfArrDict['Energy'] = np.delete(pdfArrDict['Energy'], removeMask, axis=1)
    effNorm = np.delete(effNorm, removeMask, axis=1)

    # Save un-normalized axion matrix separately for diagnostics
    AxionArr = np.delete(AxionArr, removeMask, axis=1)

    print('Generated all PDFs')
    if bDebug:
        print("Data Shape", pdfArrDict['Data'].shape)
        print("Axion Shape", pdfArrDict['Axion'].shape)

    # Draw PDFs
    if bDrawPDF:
        drawPDFs(pdfArrDict)
        return

    # Flatten 2D arrays into 1D
    pdfFlatDict['Data'] = pdfArrDict['Data'].flatten()
    pdfFlatDict['Tritium'] = pdfArrDict['Tritium'].flatten()
    pdfFlatDict['Axion'] = pdfArrDict['Axion'].flatten()
    pdfFlatDict['Bkg'] = pdfArrDict['Bkg'].flatten()
    pdfFlatDict['Fe55'] = pdfArrDict['Fe55'].flatten()
    pdfFlatDict['Zn65'] = pdfArrDict['Zn65'].flatten()
    pdfFlatDict['Ge68'] = pdfArrDict['Ge68'].flatten()
    # pdfFlatDict['Energy'] = pdfArrDict['Energy'].flatten()

    # 1D plot only
    # pdfFlatDict['Data'] = pdfArrDict['Data'].sum(axis=1)
    # pdfFlatDict['Tritium'] = pdfArrDict['Tritium'][:,0]
    # pdfFlatDict['Axion'] = pdfArrDict['Axion'][:,0]
    # pdfFlatDict['Bkg'] = pdfArrDict['Bkg'][:,0]
    # pdfFlatDict['Fe55'] = pdfArrDict['Fe55'][:,0]
    # pdfFlatDict['Zn65'] = pdfArrDict['Zn65'][:,0]
    # pdfFlatDict['Ge68'] = pdfArrDict['Ge68'][:,0]
    # pdfFlatDict['Energy'] = pdfArrDict['Energy'][:,0]

    if bDebug:
        for key in pdfFlatDict:
            print('Flattend {} shape: '.format(key), pdfFlatDict[key].shape)
            print('Flattened all PDFs')

    model = constructModel(pdfFlatDict, energyBins, bFull = bFull)
    print('Built model(s)')

    # print(trace['Axion'])
    # dfTrace = pm.trace_to_dataframe(trace[nBurn:])
    # print(dfTrace.head())
    backendDir = '{}/AveragedAxion_Enr_{}_{}'.format(inDir, int(energyThresh), int(energyThreshMax))
    # trace = pm.backends.text.load(backendDir, model)
    # drawFinalSpectra(trace=trace[nBurn:], pdfDict=pdfArrDict)
    # Sample Here
    if bSample:
        with model:
            db = pm.backends.Text(backendDir)
            # trace = pm.sample(draws=10000, chains=1, n_init=1500, chain_idx=seedNum, seed=seedNum, tune=1500, progressbar=True, trace=db)
            trace = pm.sample(draws=5000, chains=1, n_init=500, chain_idx=seedNum, seed=seedNum, tune=500, progressbar=True, trace=db)

    # Perform Diagnostics on the model -- this can only be done AFTER MCMC sampling!
    if bDiagnostic:
        modelDiagnostics(model, pdfArrDict=pdfArrDict, backendDir=backendDir, unNormAxion=AxionArr, effNorm=effNorm)

    # dfSum = pm.summary(trace[nBurn:], alpha=0.1)
    # print(dfSum.head())
    # pm.traceplot(trace)
    # plt.show()


class ErrorFnc(pm.Continuous):
    """
        Custom Error Function distribution for pymc3
        This is for efficiency/uncertainty sampling
        Not currently used
    """
    def __init__(self, mu, sd, *args, **kwargs):
        self._mu = tt.as_tensor_variable(mu)
        self._sd = tt.as_tensor_variable(sd)
        # Inherit all the arguments from pm.Continuous
        super(ErrorFnc, self).__init__(*args, **kwargs)

    def logp(self, value):
        """ Log-likelihood of Erf efficiency """
        mu = self._mu
        sigma = self._sd
        return tt.log( 0.5 + (1.+ tt.erf(value-mu)/(sigma*tt.sqrt(2))))


def constructModel(pdfDict, energyBins, bFull=True):
    """
        Construct pymc3 model, possibility to create a null hypothesis if bFull is False (no axion signal)
    """
    with pm.Model() as model:
        # Perhaps should switch to a more informative prior
        if bFull:
            Axion = pm.HalfFlat("Axion")
        Tritium = pm.HalfFlat("Tritium")
        Bkg = pm.HalfFlat("Bkg")
        Fe55 = pm.HalfFlat("Fe55")
        Zn65 = pm.HalfFlat("Zn65")
        Ge68 = pm.HalfFlat("Ge68")

        # Create the detector efficiency deterministic
        # NOTE: Problem with Error function because of the derivative of erf in theano
        # Mu = pm.Normal('Mu', mu = 0.36, sd = 0.5)
        # Sig = pm.Normal('Sig', mu = 1.26, sd = 0.5)
        # eff = 0.5*(1+tt.erf((pdfDict['Energy'] - Mu)/(tt.sqrt(2)*Sig)))

        # Logistic Efficiency
        # Mu = pm.Normal('Mu', mu = 0.591, sd = 0.088)
        # Sig = pm.Normal('Sig', mu = 5.536, sd = 0.103)
        # eff = 1./(1.+tt.exp(-(pdfDict['Energy']-Mu)/Sig))
        # Weibull efficiency
        # amp = pm.Normal('amp', mu = 0.9, sd = 0.01)
        # c = pm.Normal('c', mu = 1.8, sd = 0.2)
        # loc = pm.Normal('loc', mu = -10.8, sd = 1.)
        # scale = pm.Normal('scale', mu = 14.4, sd = 2.)
        # eff = amp*(1.-tt.exp(-tt.pow((pdfDict['Energy']-loc)/scale, c)))

        # Generate array of deterministic variables (per bin)
        if bFull:
            det = Tritium*pdfDict['Tritium'] + Bkg*pdfDict['Bkg'] + Axion*pdfDict['Axion'] + Fe55*pdfDict['Fe55'] + Zn65*pdfDict['Zn65'] + Ge68*pdfDict['Ge68']
        else:
            det = Tritium*pdfDict['Tritium'] + Bkg*pdfDict['Bkg'] + Fe55*pdfDict['Fe55'] + Zn65*pdfDict['Zn65'] + Ge68*pdfDict['Ge68']

        L = pm.Poisson("L", mu=det, observed=pdfDict['Data'])
    return model


def convertaxionPDF(startTime, endTime, axionBinSize = 5, debug = False):
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

    """

    import ROOT
    inDir = os.environ['LATDIR']+'/data/MCMC'
    fAxion = ROOT.TFile('{}/Axion_averaged_MJD_reso_Elow0_Ehi20_minutes1_2017_2018.root'.format(inDir))

    hday = fAxion.Get('hday')
    hday.RebinX(axionBinSize) # 5 minute bins
    # hday.Scale(1./5.) # Scale by 1./5.

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

    if debug and bRepeat:
        print('Axion PDF Gap Extension (should be binsize):', timeBinLowEdge[-1], hday.GetXaxis().GetBinLowEdge(0)+yearInSec)
    if bRepeat:
        timeBinLowEdge.extend([hday.GetXaxis().GetBinLowEdge(tBin)+yearInSec
                            for tBin in range(hday.GetNbinsX())
                            if tBin >= 0 and tBin <= endBin2+1]) # nBins+1 overlaps with the 0 and 1st bin

    # Total bins is 1 less than the array length of timeBinLowEdge
    nTimeBins = len(timeBinLowEdge)-1

    for energyBin in range(hday.GetNbinsY()):
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
    effMask, expMask = generateEfficiencyMask(np.array(timeBinLowEdge), energyBins, AxionArr.shape, dsList = dsList, debug = debug)

    NormAxionArr = AxionArr*effMask
    effNorm = effMask/expMask

    # Find all columns of full zeros (PDFs don't exist)
    removeMask = np.where(np.any(NormAxionArr, axis=0)==False)[0]

    binScaling = (binSize*5./(60*24)/90)

    # print('Axion Eff Norm: ', effNorm, effNorm.shape)
    print('Axion Integral: {} (before exposure) -- {} (before efficiency) --- {} (after efficiency)'.format(AxionArr.sum()*binScaling, (AxionArr*expMask).sum()*binScaling, NormAxionArr.sum()*binScaling))

    return NormAxionArr, nTimeBins, np.array(timeBinLowEdge), removeMask, effMask, effNorm


def drawFinalSpectra(pdfDict, trace):
    """
        Draws final spectra using mean of posterior
        Also draws time spectra as well
    """
    from matplotlib.ticker import FormatStrFormatter

    rebinSize = 0
    # Choose first valid rebin value that will evenly divide the number of time bins
    for rb in range(12, 51):
        if pdfDict['Data'].shape[1]%rb == 0:
            rebinSize = rb
            break
    print('Before Rebinning:', pdfDict['Data'].shape[0], pdfDict['Data'].shape[1])
    print('Rebinning by {}, the bin size will be {} minutes'.format(rebinSize, 5*rebinSize))

    rebinShape = (pdfDict['Data'].shape[0], pdfDict['Data'].shape[1]/rebinSize)
    # Convert unix time to date with seconds
    timeBins = np.array(pdfDict['Time'][:-1], dtype='datetime64[s]')
    timeBinsD = np.array(timeBins, dtype='datetime64[D]')
    # Rebin timing bins (here we only select using stepsize = rebinSize rather than summing the times)
    timeBins = timeBins[::rebinSize]

    figSpec, axSpec = plt.subplots(figsize=(10,7))
    figTime, axTime = plt.subplots(figsize=(10,7))
    fig2D, (ax2D1, ax2D2) = plt.subplots(ncols=2, figsize=(15,7))

    DataSpec = pdfDict['Data'].sum(axis=1)
    DataSpecErr = np.sqrt(DataSpec)
    TimeSpec = rebin(pdfDict['Data'].sum(axis=0), size=rebinSize)
    TimeSpecErr = np.sqrt(TimeSpec)

    axSpec.errorbar(energyBins[:-1], DataSpec, yerr=DataSpecErr, color='black', fmt='o', label='Data')
    axTime.errorbar(timeBins, TimeSpec, yerr=TimeSpecErr, color='black', fmt='o', label='Data')

    # Create a dummy for the total spectrum
    totalSpec = np.zeros(DataSpec.shape) # Dummy for energy spectrum
    totalTime = np.zeros(TimeSpec.shape) # Dummy for time spectrum
    totalModel2D = np.zeros(rebinShape) # Dummy for energy vs time spectrum
    print('rebinShape size:', rebinShape)
    print('TimeSpec size:', TimeSpec.shape)
    print('timeBins size:', timeBins.shape)

    # Take all the pdfs and flatten time and energy axes
    for key in pdfDict.keys():
        if key == 'Data' or key == 'Time' or key == 'Energy':
            continue
        tempSpec = trace[key].mean()*pdfDict[key].sum(axis=1)
        tempTime = trace[key].mean()*pdfDict[key].sum(axis=0)

        print('{} Integral: {:.3f}'.format(key, tempSpec.sum()))

        if key == 'Axion':
            axionRebin = rebin(tempTime, size=rebinSize)
            axTime.plot(timeBins, axionRebin, label=key, linestyle='--')
            totalTime += axionRebin
        else:
            axTime.plot(timeBins, rebinSize*rebin(tempTime, size=rebinSize), label=key, linestyle='--')
            totalTime += rebinSize*rebin(tempTime, size=rebinSize)

        axSpec.plot(energyBins[:-1], tempSpec, label=key, linestyle='--')
        totalSpec += tempSpec

    axSpec.plot(energyBins[:-1], totalSpec, label='Total Model')
    axSpec.set_title('DS5b, 5c, 6a Enriched + Natural')
    axSpec.set_ylabel('Counts/0.1 keV')
    axSpec.set_xlabel('Energy (keV)')
    axSpec.legend()
    plt.tight_layout()

    axTime.plot(timeBins, totalTime, label='Total Model')
    axTime.set_title('DS5b, 5c, 6a Enriched + Natural')
    axTime.set_ylabel('Counts/(140 min)')
    axTime.set_xlabel('Time')
    axTime.locator_params(axis='x', nbins=len(timeBins)/250) # Draw fewer ticks
    axTime.legend()
    plt.tight_layout()

    # Loop again -- Make sure I don't mess up the other plots
    # The counts vs time plot with binning is very finnicky when rebinning other stuff
    for key in pdfDict.keys():
        if key == 'Data' or key == 'Time':
            continue
        # Make 2D histograms
        temp2D = np.array(trace[key].mean()*pdfDict[key], copy=True)
        totalModel2D += bin_ndarray(temp2D, rebinShape)

    # Make 2D plot
    Data2D = np.array(pdfDict['Data'], copy=True)
    Data2D =  bin_ndarray(Data2D, rebinShape)
    Difference = Data2D - totalModel2D
    timeBinsD = timeBinsD[::rebinSize]
    print('Total Integral of Data:', Data2D.sum())
    print('Total Integral of Model:', totalModel2D.sum())
    print('Total Integral of Difference:', Difference.sum())

    sns.heatmap(Data2D, ax=ax2D1)
    ax2D1.invert_yaxis()
    ax2D1.set_title('DS5b, 5c, 6a Enriched + Natural (Data)')
    ax2D1.set_ylabel('Energy (keV)')
    ax2D1.locator_params(axis='x', nbins=len(timeBinsD[::250])) # Draw fewer ticks
    ax2D1.locator_params(axis='y', nbins=len(energyBins[::20])) # Draw fewer ticks
    ax2D1.set_xticklabels(timeBinsD[::250], rotation=30)
    ax2D1.set_yticklabels([round(x,1) for x in energyBins[:-1:20]], rotation=0)
    sns.heatmap(totalModel2D, ax=ax2D2)
    ax2D2.invert_yaxis()
    ax2D2.set_title('Total Model')
    ax2D2.locator_params(axis='x', nbins=len(timeBinsD[::250])) # Draw fewer ticks
    ax2D2.set_xticklabels(timeBinsD[::250], rotation=30)
    ax2D2.set_yticklabels([])
    plt.tight_layout()

    # figSpec.savefig('{}/AxionFitResult_Spectrum_Open.png'.format(os.environ['LATDIR']+'/data/MCMC'))
    # figTime.savefig('{}/AxionFitResult_Time_Open.png'.format(os.environ['LATDIR']+'/data/MCMC'))
    # fig2D.savefig('{}/AxionFitResult_2D_Open.png'.format(os.environ['LATDIR']+'/data/MCMC'))
    plt.show()


def drawPDFs(pdfArrDict):
    """
        Helper function to draw PDFs
        pdfArrDict contains all of the PDFs as well as the Data
    """
    timeLabels = np.linspace(startTime, endTime, 5, dtype='datetime64[s]')
    timeLabels = np.array(timeLabels, dtype='datetime64[D]')
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
    ax1[1,0].set_xticklabels(timeLabels, rotation=50)

    sns.heatmap(pdfArrDict['Zn65'], ax=ax1[1,1])
    ax1[1,1].set_title('Zn65 PDF')
    ax1[1,1].set_xlabel('Time')
    ax1[1,1].invert_yaxis()
    ax1[1,1].locator_params(axis='x', nbins=5)
    ax1[1,1].set_yticklabels([])
    ax1[1,1].set_xticklabels(timeLabels, rotation=50)

    sns.heatmap(pdfArrDict['Ge68'], ax=ax1[1,2])
    ax1[1,2].set_title('Ge68 PDF')
    ax1[1,2].set_xlabel('Time')
    ax1[1,2].invert_yaxis()
    ax1[1,2].locator_params(axis='x', nbins=5)
    ax1[1,2].set_yticklabels([])
    ax1[1,2].set_xticklabels(timeLabels, rotation=50)

    plt.tight_layout()

    fig1.savefig('{}/BraggFitPDF_New.png'.format(os.environ['LATDIR']+'/plots/Axion'))

    fig2, ax2 = plt.subplots(figsize=(10,7))
    sns.heatmap(pdfArrDict['Axion'], ax=ax2)
    ax2.set_title('Axion PDF')
    ax2.set_xlabel('Time')
    ax2.set_ylabel('Energy (keV)')
    ax2.invert_yaxis()
    ax2.locator_params(axis='y', nbins=len(energyBins[::50]))
    ax2.locator_params(axis='x', nbins=5)
    ax2.set_yticklabels([round(x,1) for x in energyBins[::50]], rotation=0)
    ax2.set_xticklabels(timeLabels, rotation=30)
    plt.tight_layout()
    fig2.savefig('{}/BraggAxion_New.png'.format(os.environ['LATDIR']+'/plots/Axion'))
    plt.show()


def generateEfficiencyMask(timeBinLowEdge, energyBins, AxionShape, dsList = None, debug = False):
    """
        Creates mask for efficiency*exposure, will return a number between 0 -- exposure (the weight for each bin)
        for every time and energy bin. Apply this mask to the PDFs to delete columns where the detectors are off
    """
    import ROOT
    effMask = np.ones(AxionShape)
    expMask = np.ones(AxionShape) # Exposure only
    if debug:
        print('Start, End Time:', startTime, endTime)
        print('Time Bins:', timeBinLowEdge[0], timeBinLowEdge[-1])
        print('Difference: {} -- {}'.format(startTime - timeBinLowEdge[0], timeBinLowEdge[-1] - endTime))

    timeBinSize = (timeBinLowEdge[1] - timeBinLowEdge[0])
    if debug:
        print('Timing Bin size: ', timeBinSize)

    fEff = ROOT.TFile(os.environ['LATDIR']+'/data/lat-expo-efficiency_final95_Full.root')
    # Dummy variable for Start time of DS5b to calculate gaps in livetime
    prevEnd = startTime

    for ds in dsList:
        # First build efficiency array -- assumes efficiency for all dataset is constant
        # Also combines efficiency of natural and enriched detectors
        hEff = fEff.Get('hDS{}_Norm_Tot'.format(ds))
        effArr = np.ones(len(energyBins)-1)
        for eIdx, eBin in enumerate(energyBins[:-1]):
            effArr[eIdx] = hEff.GetBinContent(hEff.FindBin(eBin+binSize/2.))

        if debug:
            print('DS{} Efficiency:'.format(ds), effArr, effArr.shape)

        # Create a new map for all runs with the exposure of each run
        if ds in ['5A', '5B', '5C']:
            dsNum = 5
        else:
            dsNum = int(ds)

        # Builds up dictionary of run:exposure for every dataset
        runMatrix = {}
        bkgRanges = bkg.getRanges(dsNum)
        chList = det.getGoodChanList(dsNum)

        # Loop through channels to build up total runLists
        for ch in chList:
            cpd = int(det.getChanCPD(dsNum,ch))
            inFile = os.environ['LATDIR']+'/data/expLists/DS{}_C{}P{}D{}.txt'.format(ds, *str(cpd))
            try:
                with open(inFile) as f:
                    # Loop through file,
                    for line in f:
                        currArr = np.array(line.split(','), dtype=np.float)
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
                        effMask[:, prevIndex] = 0
                        expMask[:, prevIndex] = 0
                    else:
                        effMask[:, prevIndex:startIndex] = 0
                        expMask[:, prevIndex:startIndex] = 0
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
                        effMask[:, startIndex] = effArr*runMatrix[run]/(len(energyBins)-1)
                        expMask[:, startIndex] = runMatrix[run]/(len(energyBins)-1)
                    else:
                        # This part is really long and complicated, basically need to copy (repeat) the effArr to match the size we want to fill inside effMask, the reshape is just to make the array sizes match

                        # If the shape doesn't match, it's because we're at the last time bin, shrink
                        if effMask[:, startIndex:stopIndex].shape[1] != stopIndex-startIndex:
                            print('Matching off run {} ({}-{})'.format(run,startIndex, stopIndex) )
                            effMat = np.repeat(effArr,stopIndex-startIndex-1).reshape(len(effArr), stopIndex-startIndex-1)
                        else:
                            effMat = np.repeat(effArr,stopIndex-startIndex).reshape(len(effArr), stopIndex-startIndex)

                        effMask[:, startIndex:stopIndex] = effMat*runMatrix[run]/(stopIndex-startIndex-1)/(len(energyBins)-1)
                        expMask[:, startIndex:stopIndex] = runMatrix[run]/(stopIndex-startIndex-1)/(len(energyBins)-1) # No efficiency
                else:
                    effMask[:, startIndex:stopIndex] = 0
                    expMask[:, startIndex:stopIndex] = 0
                    print('Run {} setting to 0'.format(run))
                    print(effMask[100, startIndex:stopIndex])

                prevEnd = eTime

        if debug:
            print(effMask)
            print('DS{} -- Exposure: {} -- Eff Corr Exposure: {} (kg-day)'.format(ds, expMask.sum(), effMask.sum())) # Sum over one axis at a high energy, should get exposure?
            print('Total time gap for DS{} -- {}'.format(ds, totalGap))
            # from mpl_toolkits.mplot3d import Axes3D
            # feff, aeff = plt.subplots(figsize=(8,6))
            # sns.heatmap(effMask, ax=aeff)
            # feff = plt.figure()
            # aeff = feff.gca(projection='3d')
            # xmesh, ymesh = np.meshgrid(energyBins[:-1], timeBinLowEdge[:-1])
            # aeff.plot_surface(xmesh, ymesh, effMask.T)
            # plt.show()

    return effMask, expMask


def getSigma(energy, dsNum=1):
    p0,p1,p2 = 0., 0., 0.
    if dsNum==0:
        p0 = 0.147; p1=0.0173; p2=0.0003
    elif dsNum==1:
        p0 = 0.136; p1=0.0174; p2=0.00028
    elif dsNum==3:
        p0 = 0.162; p1=0.0172; p2=0.000297
    elif dsNum==4:
        p0 = 0.218; p1=0.015; p2=0.00035
    elif dsNum==5:
        p0 = 0.2121; p1=0.01838; p2=0.00031137
    else:
        p0 = 0.2121; p1=0.01838; p2=0.00031137
    return np.sqrt(p0*p0 + p1*p1*energy + p2*p2*energy*energy)


def generateGaus(energyList, meansDict):
    """
        Generates Gaussian PDFs for X-ray peaks based off of energy array
    """
    outDict = {}
    for name, mean in meansDict.items():
        outDict.setdefault(name, norm.pdf(energyList, loc=mean, scale=getSigma(mean, 5)))
    return outDict


def modelDiagnostics(model, pdfArrDict=None, backendDir=None, unNormAxion=None, effNorm=None):
    """
        Loads saved traces, draws plots, performs diagnostics on models
        Brian's notes on diagnostics:

        Gelman-Rubin (Rhat) -- tests for lack of convergence by comparing the variance between multiple chains to the variance within one chain.
        If the chains are converged, the variance should be the same (ratio should be 1). Needs multiple chains with differing starting points.

        Auto-correlation -- measures how linearly dependent the current value of the chain is to past values
        It essentially tells us how much information is avaliable, if there is a lot of correlation,
        there's less information than sampled from stationary distribution
        The value is between -1 and 1, the lag where the autocorr is ~0 means the number of samples where the values aren't autocorrelated
        The amount of autocorrelation determines the effective sample size (n_effective), for setting a 95% credible interval

        Effective sample size (effective_n) -- is an estimate on the effective sample size

        Highest Posterior Density (HPD) -- The HPD is an estimator that calculates the minimum width of the Bayesian credible interval (BCI).
        If the posterior density is multimodal, the HPD does not result in an interval estimate (thankfully it's not for our case)

    """
    import corner

    # Load from files
    print('Loading Backend from',backendDir)
    trace = pm.backends.text.load(backendDir, model)
    # Get summary of trace as a dataframe
    dfSum = pm.summary(trace[nBurn:], alpha=0.1)
    print(dfSum.head(10))
    AxionIntegral = pdfArrDict['Axion'].sum()

    # This gets printed out from "convertaxionPDF", it's the amount I artificially scale the axion PDF integral
    AxionNorm = (np.amax(unNormAxion))
    print('Axion PDF Integral:', AxionIntegral)
    print('Normalization: {}'.format(AxionNorm))

    # Loop through and calculate 95% CI for all parameters
    print('Calculating 95% Credible Intervals: ')
    print('Par \t  Mean \t hpd_5 \t hpd_95')
    # Here all the PDFs are normalized in their own way (from model construction section)
    # therefore multiplying the PDF integral by the 95% CI should give number of counts
    for par in dfSum.index:
        # Integrate PDF
        parInt = (pdfArrDict[par]/effNorm).sum()
        # parInt = (pdfArrDict[par]).sum()

        # Multiply integral by interval
        print('{}    {:.2f}    {:.2f}    {:.2f}'.format(par, parInt*dfSum.loc[par, 'mean'], parInt*dfSum.loc[par, 'hpd_5'], parInt*dfSum.loc[par, 'hpd_95']))
        if par == 'Axion':
            print('Axion: Scaled Integral: {} \t Unscaled Integral: {}'.format(parInt*dfSum.loc[par, 'hpd_95'], parInt))
            print('95% CI limit: {}'.format(parInt*dfSum.loc[par, 'hpd_95']/parInt))

    # Debate: when converting to lambda,
    # axionLambda = 3600*24*dfSum.loc['Axion', 'hpd_95']/AxionNorm
    # print("Axion Lambda: {} --- g_agg (E-8 GeV): {}".format(axionLambda, np.power(axionLambda, 0.25)))

    pm.traceplot(trace[nBurn:])
    pm.plot_posterior(trace[nBurn:], alpha_level=0.1, round_to=6)
    # pm.forestplot(trace[nBurn:], varnames=['Axion', 'Tritium', 'Bkg', 'Ge68', 'Fe55', 'Zn65'], alpha=0.1)
    # traceBurn = trace[nBurn:]
    # print(np.array([traceBurn['Axion'], traceBurn['Tritium']]).T.shape)
    # cornerArr = np.array([traceBurn['Axion'], traceBurn['Tritium'], traceBurn['Bkg'], traceBurn['Ge68'], traceBurn['Fe55'], traceBurn['Zn65']]).T
    # figure = corner.corner(cornerArr,  quantiles=[0.05, 0.95], show_titles=True, title_fmt=".5f", labels=['Axion', 'Tritium', 'Bkg', 'Ge68', 'Fe55', 'Zn65'])

    # Format is hpdDict[nChain][Parameter] = [lower, upper]
    # hpdDict = pm.hpd(trace[nBurn:], alpha=0.1)
    # print(hpdDict[1])

    # pm.autocorrplot(trace[nBurn:], max_lag=100)
    # print(pm.autocorr(trace[nBurn:]))
    # pm.densityplot(trace[nBurn:])
    plt.show()


def posteriorChecks(model, trace):
    """
        Performs various posterior checks
        Posterior Predictive Checks:
            - Simulates replicating data under the fitted model and then comparing these to the observed data
            - This checks for systematic discrepancies between real and simulated data

        Widely-applicable Information Criterion (WAIC):
            - Fully Bayesian criterion for estimating out-of-sample expectation, using the computed log pointwise posterior predictive density (LPPD) and correcting for the effective number of parameters to adjust for overfitting.
            - This is primarilly for the comparison between different models

        Leave-one-out Cross-validation (LOO):
            - Estimate of the out-of-sample predictive fit. In cross-validation, the data are repeatedly partitioned into training and holdout sets, iteratively fitting the model with the former and evaluating the fit with the holdout data. PyMC's implementation of LOO is using Pareto-smoothed importance sampling, it provides and estimate of point-wise out-of-sample prediction accuracy

        Note: out-of-sample is data not used for the fit (ie: making a prediction)

    """
    # Posterior Predictive Checks -- Generates 500 toy samples of size 100
    # This is essentially toy MC
    model_ppc = pm.sample_ppc(trace, samples=500, model=model)

    # Widely-applicable Information Criterion (WAIC)
    model_waic = pm.waic(trace[len(trace)-100:], model, progressbar=True)

    # Leave-one-out Cross-validation (LOO)
    model_loo = pm.loo(trace[len(trace)-100:], model, progressbar=True)

    # ppc = pm.sample_ppc(trace[nBurn:], samples=500, model=model)
    # print(np.asarray(ppc['L'].shape), ppc.keys())
    # _, axppc = plt.subplots(figsize=(12, 6))
    # axppc.hist([n.mean() for n in ppc['L']], bins=19, alpha=0.5)
    # axppc.set(title='Posterior predictive for L', xlabel='L(x)', ylabel='Frequency');
    # plt.show()

    # df_comp_WAIC = pm.compare(models = [model, modelBasic], traces = [trace[nBurn:], traceBasic[nBurn]])
    # df_comp_WAIC.head()
    # pm.compareplot(df_comp_WAIC)
    # df_comp_LOO = pm.compare(models = [model, modelBasic], traces = [trace[nBurn:], traceBasic[nBurn]], ic='LOO')
    # df_comp_LOO.head()
    # pm.compareplot(df_comp_LOO)

    # LOO results
    # LOO  pLOO   dLOO weight      SE   dSE warning
    # 0  61479.8  5.68      0   0.94  798.63     0       1
    # 1    61502  4.95  22.12   0.06  798.47  10.2       1

    return model_ppc, model_waic, model_loo


def rebin(inArr, size=2):
    """
        Rebin an array
    """
    if not float(inArr.shape[0]/size).is_integer():
        print('Error: Input array length is not evenly divisible by rebin size')
        return
    rebinnedArr = np.empty(inArr[::size].shape)
    for i in range(size):
        rebinnedArr += inArr[i::size]
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


if __name__ == '__main__':
    main(sys.argv[1:])
