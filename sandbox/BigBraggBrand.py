#!/usr/bin/env python
"""
    Big NUTS Bragg Brand fitting
"""

import os
import numpy as np
import pymc3 as pm
import matplotlib.pyplot as plt
import seaborn as sns
import theano.tensor as tt
from scipy.interpolate import interp1d
from scipy.stats import norm
import pandas as pd
sns.set(style='darkgrid')

# Define some global parameters
inDir = os.environ['LATDIR']+'/data/MCMC'
seedNum, dsNum = 1, 5
# Energy range (also fitting range!)
# Wenqin stops at 12 keV since the axions stop at 12 -- endpoint of trit is 18.
# Scanned from 12 to 19.8
energyThresh, energyThreshMax, binSize = 2.4, 18, 0.1
energyBins = np.linspace(energyThresh,energyThreshMax, round((energyThreshMax-energyThresh)/0.1)+1)
# This start and end time are for DS5b
startTime = 1485575988
endTime = 1489762153
# Fake start times
# startTime = 1483228800
# endTime = 1514678400 # End of 2017
# endTime = 1546300800 # End of 2018
yearInSec = 31536000 # Number of seconds in 1 year
nBurn = 1500

# Save data into dataframe -- if this hasn't been done before
# reduceData(dsNum)

def main():

    pdfArrDict, pdfFlatDict = {}, {}

    # Build Axion PDF -- Use low edges of PDF to generate bins
    AxionArr, nTimeBins, timeBinLowEdge, removeMask = convertaxionPDF(startTime, endTime)

    AxionNorm = (np.amax(AxionArr)/2.)
    # Normalize Axion
    print('Axion Normalization Constant:', AxionNorm)
    pdfArrDict['Axion'] = AxionArr/AxionNorm

    # Load Background data and bin exactly as PDF
    df = pd.read_hdf('{}/DS{}_Spectrum.h5'.format(inDir, dsNum))
    df.sort_values(by=['UnixTime'], inplace=True)
    dataArr = df.loc[:, ['Energy','UnixTime']].values

    #NOTE: Extend data by mirroring it one year later -- this is convoluted but whatever
    # dataArr = np.append(dataArr, np.array([dataArr[:,0], np.add(dataArr[:,1], yearInSec)]).T, axis=0)
    pdfArrDict['Data'], xedges, yedges = np.histogram2d(dataArr[:,0], dataArr[:,1], bins=[energyBins, timeBinLowEdge])

    # Build Tritium PDF -- need interpolation because it's 0.2 keV bins!
    # NOTE: We may want to use 0.2 keV bins, the difference may not be noticeable and the
    dfTrit = pd.read_hdf('{}/TritSpec.h5'.format(inDir))
    tritEnergy = dfTrit['Energy'].values
    tritSpec = dfTrit['Tritium'].values
    tritSpec = np.array([100*x if x >= 0. else 0. for x in tritSpec]) # Get rid of negative values
    tritInterp = interp1d(tritEnergy, tritSpec, fill_value='extrapolate')
    tritList = [tritInterp(x) for x in energyBins[:-1]]
    pdfArrDict['Tritium'] = np.array(tritList*nTimeBins).reshape(nTimeBins, len(tritList)).T

    meansDict = {'Fe55':6.54, 'Zn65':8.98, 'Ge68':10.37}
    gausDict = generateGaus(energyBins[:-1], meansDict)
    for name, Arr in gausDict.items():
        pdfArrDict.setdefault(name, np.array(Arr.tolist()*nTimeBins).reshape(nTimeBins, len(Arr)).T)

    # Create a matrix for energy bin center values
    pdfArrDict.setdefault('Energy', np.array((energyBins[:-1]+binSize/2).tolist()*nTimeBins).reshape(nTimeBins, len(energyBins[:-1]) ).T)


    # Remove columns with zeros
    pdfArrDict['Axion'] = np.delete(pdfArrDict['Axion'], removeMask, axis=1)
    pdfArrDict['Data'] = np.delete(pdfArrDict['Data'], removeMask, axis=1)
    pdfArrDict['Tritium'] = np.delete(pdfArrDict['Tritium'], removeMask, axis=1)
    pdfArrDict['Fe55'] = np.delete(pdfArrDict['Fe55'], removeMask, axis=1)
    pdfArrDict['Zn65'] = np.delete(pdfArrDict['Zn65'], removeMask, axis=1)
    pdfArrDict['Ge68'] = np.delete(pdfArrDict['Ge68'], removeMask, axis=1)
    pdfArrDict['Bkg'] = np.ones(pdfArrDict['Data'].shape)
    pdfArrDict['Time'] = np.delete(timeBinLowEdge, removeMask)
    pdfArrDict['Energy'] = np.delete(pdfArrDict['Energy'], removeMask, axis=1)
    AxionArr = np.delete(AxionArr, removeMask, axis=1)
    print('Generated all PDFs')
    # Draw PDFs
    # drawPDFs(pdfArrDict)
    print("Data Shape", pdfArrDict['Data'].shape)
    print("Axion Shape", pdfArrDict['Axion'].shape)
    # return

    # Flatten 2D arrays into 1D
    pdfFlatDict['Data'] = pdfArrDict['Data'].flatten()
    pdfFlatDict['Tritium'] = pdfArrDict['Tritium'].flatten()
    pdfFlatDict['Axion'] = pdfArrDict['Axion'].flatten()
    pdfFlatDict['Bkg'] = pdfArrDict['Bkg'].flatten()
    pdfFlatDict['Fe55'] = pdfArrDict['Fe55'].flatten()
    pdfFlatDict['Zn65'] = pdfArrDict['Zn65'].flatten()
    pdfFlatDict['Ge68'] = pdfArrDict['Ge68'].flatten()
    pdfFlatDict['Energy'] = pdfArrDict['Energy'].flatten()

    # 1D plot only
    # pdfFlatDict['Data'] = pdfArrDict['Data'].sum(axis=1)
    # pdfFlatDict['Tritium'] = pdfArrDict['Tritium'][:,0]
    # pdfFlatDict['Axion'] = pdfArrDict['Axion'][:,0]
    # pdfFlatDict['Bkg'] = pdfArrDict['Bkg'][:,0]
    # pdfFlatDict['Fe55'] = pdfArrDict['Fe55'][:,0]
    # pdfFlatDict['Zn65'] = pdfArrDict['Zn65'][:,0]
    # pdfFlatDict['Ge68'] = pdfArrDict['Ge68'][:,0]
    # pdfFlatDict['Energy'] = pdfArrDict['Energy'][:,0]

    print('Flattened Data Size:', pdfFlatDict['Data'].shape)
    print('Energy Size', pdfFlatDict['Energy'].shape)
    print('Flattened all PDFs')

    model = constructModel(pdfFlatDict, energyBins)
    # modelBasic = constructBasicModel(pdfFlatDict)
    print('Built model(s)')

    # return

    # print(trace['Axion'])
    # dfTrace = pm.trace_to_dataframe(trace[nBurn:])
    # print(dfTrace.head())
    # drawFinalSpectra(trace=trace[nBurn:], pdfDict=pdfArrDict)
    # modelDiagnostics(model, pdfArrDict=pdfArrDict, backendDir='{}/AveragedAxion_{:d}keV'.format(inDir,energyThreshMax))
    # modelDiagnostics(model, pdfArrDict=pdfArrDict, backendDir='{}/AveragedAxion_WithEff'.format(inDir), unNormAxion=AxionArr)

    # Sample Here
    with model:
        # trace = pm.sample(draws=5000, chains=1, n_init=500, chain_idx=seedNum, seed=seedNum, tune=500, progressbar=True)
        db = pm.backends.Text('{}/AveragedAxion_WithEff'.format(inDir))
        trace = pm.sample(draws=10000, chains=1, n_init=1500, chain_idx=seedNum, seed=seedNum, tune=1500, progressbar=True, trace=db)
    # with modelBasic:
        # dbBasic = pm.backends.Text('{}/AveragedNoAxion'.format(inDir))
        # traceBasic = pm.sample(draws=10000, chains=1, n_init=1500, chain_idx=seedNum, seed=seedNum, tune=1500, progressbar=True, trace=dbBasic)

    # pm.traceplot(trace)
    # plt.show()


class ErrorFnc(pm.Continuous):
    """
        Custom Error Function distribution for pymc3
    """
    def __init__(self, mu, sd, *args, **kwargs):
        self._mu = tt.as_tensor_variable(mu)
        self._sd = tt.as_tensor_variable(sd)
        super(ErrorFnc, self).__init__(*args, **kwargs)

    def logp(self, value):
        """ Log-likelihood of Erf efficiency """
        mu = self._mu
        sigma = self._sd
        return tt.log( 0.5 + (1.+ tt.erf(value-mu)/(sigma*tt.sqrt(2))))



def constructModel(pdfDict, energyBins):
    """
        Construct pymc3 model
    """
    with pm.Model() as model:
        # We using HalfFlat cuz we dgaf -- probably we should switch to a more informative prior
        Axion = pm.HalfFlat("Axion")
        Tritium = pm.HalfFlat("Tritium")
        Bkg = pm.HalfFlat("Bkg")
        Fe55 = pm.HalfFlat("Fe55")
        Zn65 = pm.HalfFlat("Zn65")
        Ge68 = pm.HalfFlat("Ge68")

        # Create the detector efficiency deterministic
        # NOTE: Problem with
        # Mu = pm.Normal('Mu', mu = 0.36, sd = 0.5)
        # Sig = pm.Normal('Sig', mu = 1.26, sd = 0.5)
        # eff = 0.5*(1+tt.erf((pdfDict['Energy'] - Mu)/(tt.sqrt(2)*Sig)))

        # Reparameterize efficiency as logistic function
        Mu = pm.Normal('Mu', mu = 0.591, sd = 0.088)
        Sig = pm.Normal('Sig', mu = 5.536, sd = 0.103)
        eff = 1./(1.+tt.exp(-(pdfDict['Energy']-Mu)/Sig))
        # Generate array of deterministic variables (per bin)
        det = (Tritium*pdfDict['Tritium'] + Bkg*pdfDict['Bkg'] + Axion*pdfDict['Axion'] + Fe55*pdfDict['Fe55'] + Zn65*pdfDict['Zn65'] + Ge68*pdfDict['Ge68'])*eff
        # det = Tritium*pdfDict['Tritium'] + Bkg*pdfDict['Bkg'] + Axion*pdfDict['Axion'] + Fe55*pdfDict['Fe55'] + Zn65*pdfDict['Zn65'] + Ge68*pdfDict['Ge68']

        L = pm.Poisson("L", mu=det, observed=pdfDict['Data'])
    return model


def constructBasicModel(pdfDict):
    """
        Construct model without axions! (Null hypothesis)
    """
    with pm.Model() as model:
        # We using HalfFlat cuz we dgaf -- probably we should switch to a more informative prior
        Tritium = pm.HalfFlat("Tritium")
        Bkg = pm.HalfFlat("Bkg")
        Fe55 = pm.HalfFlat("Fe55")
        Zn65 = pm.HalfFlat("Zn65")
        Ge68 = pm.HalfFlat("Ge68")

        # Generate array of deterministic variables (per bin)
        det = Tritium*pdfDict['Tritium'] + Bkg*pdfDict['Bkg'] + Fe55*pdfDict['Fe55'] + Zn65*pdfDict['Zn65'] + Ge68*pdfDict['Ge68']

        L = pm.Poisson("L", mu=det, observed=pdfDict['Data'])
    return model


def convertaxionPDF(startTime, endTime):
    """
        Converts 2D histogram PDF for axion into numpy array

        Wenqin's Axion PDF is generated for Lambda = 1.
        To convert a result to a lambda:
            1) Take the integral of the pdf and scale by 3600*24 (per day)
            2) Re-scale with any artificial normalizations
            3) Lambda = Limit/Re-scaled Integral

    """
    import ROOT
    inDir = os.environ['LATDIR']+'/data/MCMC'
    fAxion  = ROOT.TFile('{}/Axion_averaged_MJD_reso_Elow0_Ehi20_minutes1_2017_2018.root'.format(inDir))

    hday = fAxion.Get('hday')
    hday.RebinX(5) # 5 minute bins
    # Find bins to cut off PDF at
    startBin = hday.GetXaxis().FindBin(startTime)
    endBin = hday.GetXaxis().FindBin(endTime)
    AxionList = []
    # Because this is finding the bin edges, need to go beyond the final bin
    timeBinLowEdge = [hday.GetXaxis().GetBinLowEdge(time)
                        for time in range(hday.GetNbinsX())
                        if time >= startBin and time <= endBin+1]

    #NOTE: Double the time bins for tests!
    # timeBinLowEdge.extend([hday.GetXaxis().GetBinLowEdge(time)+yearInSec
    #                     for time in range(hday.GetNbinsX())
    #                     if time >= startBin and time <= endBin])

    # Total bins is 1 less
    nTimeBins = len(timeBinLowEdge)-1
    for energyBin in range(hday.GetNbinsY()):
        # Apply energy cuts on the low edge of the bin -- rounding here is necessary because of ROOT's stupid precision
        if round(hday.GetYaxis().GetBinLowEdge(energyBin),1) < energyThresh: continue
        if round(hday.GetYaxis().GetBinLowEdge(energyBin),1) > energyThreshMax-0.1: continue
        # print("Axion Energy Bin: ", round(hday.GetYaxis().GetBinLowEdge(energyBin),1))
        AxionList.append([hday.GetBinContent(time, energyBin)
                            for time in range(hday.GetNbinsX())
                            if time >= startBin and time <= endBin])

    fAxion.Close()
    AxionArr = np.array(AxionList)
    # AxionArr = np.append(AxionArr, AxionArr, axis=1) #NOTE: Debug -- extend pdf for tests
    timeMask = generateLivetimeMask(np.array(timeBinLowEdge), AxionArr.shape)

    # Roughly normalize so the axion scale is the same as the other PDFs
    # NormAxionArr = AxionArr/(np.amax(AxionArr)/2.)*timeMask
    # print('Axion Normalization Constant:', (np.amax(AxionArr)/2.))
    # NormAxionArr = AxionArr/(np.amax(AxionArr)/2.) #NOTE: Debug, remove time mask for tests

    NormAxionArr = AxionArr*timeMask

    # Find all columns of full zeros (PDFs don't exist)
    removeMask = np.where(np.any(NormAxionArr, axis=0)==False)[0]

    return NormAxionArr, nTimeBins, np.array(timeBinLowEdge), removeMask


def drawFinalSpectra(pdfDict, trace):
    """
        Draws final spectra using mean of posterior
        Also draws time spectra as well
    """
    from matplotlib.ticker import FormatStrFormatter

    # Rebin time axis by 28 bins (140 minute bins)
    rebinSize = 28
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
    totalSpec = np.zeros(DataSpec.shape)
    totalTime = np.zeros(TimeSpec.shape)
    totalModel2D = np.zeros(rebinShape)

    # Take all the pdfs and flatten time axis
    for key in pdfDict.keys():
        if key == 'Data' or key == 'Time':
            continue
        tempTime = np.zeros(TimeSpec.shape)
        tempSpec = trace[key].mean()*pdfDict[key].sum(axis=1)
        tempTime = trace[key].mean()*pdfDict[key].sum(axis=0)

        print('{} Integral: {:.3f}'.format(key, tempSpec.sum()))

        if key == 'Axion':
            axionRebin = rebin(tempTime, size=rebinSize)
            axTime.plot(timeBins, axionRebin, label=key, linestyle='--')
            totalTime += axionRebin
        else:
            axTime.plot(timeBins, rebinSize*np.unique(tempTime)*np.ones(TimeSpec.shape), label=key, linestyle='--')
            totalTime += rebinSize*np.unique(tempTime)*np.ones(TimeSpec.shape)

        axSpec.plot(energyBins[:-1], tempSpec, label=key, linestyle='--')
        totalSpec += tempSpec

    axSpec.plot(energyBins[:-1], totalSpec, label='Total Model')
    axSpec.set_title('DS5b Enriched + Natural')
    axSpec.set_ylabel('Counts/0.1 keV')
    axSpec.set_xlabel('Energy (keV)')
    axSpec.legend()
    plt.tight_layout()

    axTime.plot(timeBins, totalTime, label='Total Model')
    axTime.set_title('DS5b Enriched + Natural')
    axTime.set_ylabel('Counts/(140 min)')
    axTime.set_xlabel('Time')
    axTime.locator_params(axis='x', nbins=len(timeBins)/50) # Draw fewer ticks
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
    ax2D1.set_title('DS5b Enriched + Natural (Data)')
    ax2D1.set_ylabel('Energy (keV)')
    ax2D1.locator_params(axis='x', nbins=len(timeBinsD[::50])) # Draw fewer ticks
    ax2D1.locator_params(axis='y', nbins=len(energyBins[::20])) # Draw fewer ticks
    ax2D1.set_xticklabels(timeBinsD[::50], rotation=30)
    ax2D1.set_yticklabels([round(x,1) for x in energyBins[:-1:20]], rotation=0)

    sns.heatmap(totalModel2D, ax=ax2D2)
    ax2D2.invert_yaxis()
    ax2D2.set_title('Total Model')
    ax2D2.locator_params(axis='x', nbins=len(timeBinsD[::50])) # Draw fewer ticks
    ax2D2.set_xticklabels(timeBinsD[::50], rotation=30)
    ax2D2.set_yticklabels([])
    plt.tight_layout()
    # figSpec.savefig('{}/AxionFitResult_Spectrum_2.png'.format(os.environ['LATDIR']+'/data/MCMC'))
    # figTime.savefig('{}/AxionFitResult_Time_2.png'.format(os.environ['LATDIR']+'/data/MCMC'))
    # fig2D.savefig('{}/AxionFitResult_2D_2.png'.format(os.environ['LATDIR']+'/data/MCMC'))
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
    plt.show()
    # fig1.savefig('{}/BraggFitPDF.png'.format(os.environ['LATDIR']+'/data/MCMC'))


def generateLivetimeMask(timeBinLowEdge, AxionShape):
    """
        Creates mask for exposure, will return a number between 0 and 1 (weight for each bin)

    """
    inFile = os.environ['LATDIR']+'/data/MCMC/DS5b_RunTimes.txt'
    livetimeMask = np.ones(AxionShape)
    gapMatrix = {}
    prevEnd = 1485575988 # Start time of DS5b
    totalGap = 0
    try:
        f = open(inFile)
    except:
        print("File doesn't exist!")
    for line in f:
        currentArr = np.array(line.split('\t'), dtype=np.int)
        if currentArr[1] != prevEnd:
            # Fills in gapMatrix with livetime gaps, run: end, start, gap size (s)
            gapMatrix.setdefault(currentArr[0], [prevEnd, currentArr[1], currentArr[1]-prevEnd])
            startIndex = np.where(timeBinLowEdge <= prevEnd)[0][-1]
            stopIndex = np.where(timeBinLowEdge >= currentArr[1])[0][0]
            livetimeMask[:, startIndex:stopIndex] = 0
            totalGap += currentArr[1]-prevEnd
            prevEnd = currentArr[2]
        else:
            prevEnd = currentArr[2]

    return livetimeMask


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


def modelDiagnostics(model, pdfArrDict=None, backendDir=None, unNormAxion=None):
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

    Mu = dfSum.loc['Mu', 'mean']
    Sig = dfSum.loc['Sig', 'mean']
    eff = 1./(1.+np.exp(-(pdfArrDict['Energy']-Mu)/Sig))
    AxionIntegral = pdfArrDict['Axion'].sum()

    # This gets printed out from "convertaxionPDF", it's the amount I artificially scale the axion PDF integral
    # AxionNorm = 11878.384765625
    AxionNorm = (np.amax(unNormAxion)/2.)
    AxionNormEff = (np.amax(unNormAxion*eff)/2.)
    # AxionIntegralEff = (unNormAxion*eff).sum()
    print('Axion PDF Integral:', AxionIntegral)
    # print('Efficiency Corrected Axion PDF Integral:', AxionIntegralEff)
    print('Normalizations: {}  -- {}'.format(AxionNorm, AxionNormEff))

    # Loop through and calculate 90% BCI for all parameters
    print('Calculating 90% intervals: ')
    print('Par \t Mean \t hpd_5 \t hpd_95')
    for par in dfSum.index:
        if par == 'Mu': continue
        if par == 'Sig': continue
        parInt = pdfArrDict[par].sum()
        print('{}, {}, {}, {}'.format(par, parInt*dfSum.loc[par, 'mean'], parInt*dfSum.loc[par, 'hpd_5'], parInt*dfSum.loc[par, 'hpd_95']))


    # Debate: when converting to lambda,
    axionLambda = 3600*24*dfSum.loc['Axion', 'hpd_95']/AxionNorm
    axionLambdaEff = 3600*24*dfSum.loc['Axion', 'hpd_95']/AxionNormEff
    print("Axion Lambda: {} --- g_agg (E-8 GeV): {}".format(axionLambda, np.power(axionLambda, 0.25)))
    print("Axion Lambda (Eff Corrected): {} --- g_agg (E-8 GeV): {}".format(axionLambdaEff, np.power(axionLambdaEff, 0.25)))

    # pm.traceplot(trace[nBurn:])
    pm.plot_posterior(trace[nBurn:], alpha_level=0.1, round_to=6)
    # pm.forestplot(trace[nBurn:], alpha=0.1)
    # traceBurn = trace[nBurn:]
    # print(np.array([traceBurn['Axion'], traceBurn['Tritium']]).T.shape)
    # cornerArr = np.array([traceBurn['Axion'], traceBurn['Tritium'], traceBurn['Bkg'], traceBurn['Ge68'], traceBurn['Fe55'], traceBurn['Zn65'], traceBurn['Mu'], traceBurn['Sig']]).T
    # figure = corner.corner(cornerArr,  quantiles=[0.05, 0.95], show_titles=True, title_fmt=".5f", labels=['Axion', 'Tritium', 'Bkg', 'Ge68', 'Fe55', 'Zn65', 'Mu', 'Sig'])

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


def reduceData(dsNum):
    import ROOT
    inDir, outDir = '/Users/brianzhu/project/cuts/corrfs_rn2', os.environ['LATDIR']+'/data/MCMC'
    skimTree = ROOT.TChain("skimTree")
    skimTree.Add("{}/corrfs_rn-DS{}-*.root".format(inDir, dsNum))
    theCut = "trapENFCal > 0.8 && trapENFCal < 50"
    nPass = skimTree.Draw('trapENFCal:channel:globalTime', theCut, 'goff')
    print ("{} events passed all cuts".format(nPass))
    nEnergy = skimTree.GetV1()
    nCh = skimTree.GetV2()
    nTime = skimTree.GetV3()
    nChList = list(int(nCh[n]) for n in range(nPass))
    nEnergyList = list(float(nEnergy[n]) for n in range(nPass))
    nTimeList = list(int(nTime[n]) for n in range(nPass))

    df = pd.DataFrame({"Energy":nEnergyList, "Channel":nChList, 'UnixTime':nTimeList})
    print(df.head(10))
    df.to_hdf('{}/DS{}_Spectrum.h5'.format(outDir, dsNum), 'skimTree', mode='w', format='table')


if __name__ == '__main__':
    main()
