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

def main():
    dsNum = 5
    pdfArrDict, pdfFlatDict = {}, {}

    # ReduceData(dsNum)
    # Load Background data and bin exactly as PDF
    startTime = 1485575988
    endTime = 1489762153

    # Build Axion PDF -- Use low edges of PDF to generate bins
    pdfArrDict['Axion'], nTimeBins, timeBinLowEdge, removeMask = convertaxionPDF(startTime, endTime)
    # print("Using {} time bins".format(nTimeBins))

    inDir = os.environ['LATDIR']+'/data/MCMC'
    df = pd.read_hdf('{}/DS{}_Spectrum.h5'.format(inDir, dsNum))
    df.sort_values(by=['UnixTime'], inplace=True)
    dataArr = df.loc[:, ['Energy','UnixTime']].values
    energyBins = np.arange(0,20.1,0.1)
    pdfArrDict['Data'], xedges, yedges = np.histogram2d(dataArr[:,0], dataArr[:,1], bins=[energyBins, timeBinLowEdge])

    # Build Tritium PDF -- need interpolation because it's 0.2 keV bins!
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

    # Remove columns with zeros
    pdfArrDict['Axion'] = np.delete(pdfArrDict['Axion'], removeMask, axis=1)
    pdfArrDict['Data'] = np.delete(pdfArrDict['Data'], removeMask, axis=1)
    pdfArrDict['Tritium'] = np.delete(pdfArrDict['Tritium'], removeMask, axis=1)
    pdfArrDict['Fe55'] = np.delete(pdfArrDict['Fe55'], removeMask, axis=1)
    pdfArrDict['Zn65'] = np.delete(pdfArrDict['Zn65'], removeMask, axis=1)
    pdfArrDict['Ge68'] = np.delete(pdfArrDict['Ge68'], removeMask, axis=1)
    pdfArrDict['Bkg'] = np.ones(pdfArrDict['Data'].shape)
    print('Generated all PDFs')

    # Draw PDFs
    # drawPDFs(pdfArrDict)
    dftrace = pd.read_hdf('{}/AxionFitTrace_DS{}.h5'.format(inDir, dsNum))
    drawFinalSpectra(pdfArrDict, dftrace)
    return

    # Flatten 2D arrays into 1D
    pdfFlatDict['Data'] = pdfArrDict['Data'].flatten()
    pdfFlatDict['Tritium'] = pdfArrDict['Tritium'].flatten()
    pdfFlatDict['Axion'] = pdfArrDict['Axion'].flatten()
    pdfFlatDict['Bkg'] = pdfArrDict['Bkg'].flatten()
    pdfFlatDict['Fe55'] = pdfArrDict['Fe55'].flatten()
    pdfFlatDict['Zn65'] = pdfArrDict['Zn65'].flatten()
    pdfFlatDict['Ge68'] = pdfArrDict['Ge68'].flatten()
    print('Flattened all PDFs')

    model = constructModel(pdfFlatDict)
    print('Built model')

    with model:
        trace = pm.sample(draws=75000, chains=1, n_init=1500, progressbar=True)
    dfTrace = pm.trace_to_dataframe(trace, include_transformed=True)
    dfTrace.to_hdf('{}/AxionFitTrace_DS{}_2.h5'.format(inDir, dsNum), "skimTree", mode='w', format='table')
    print(pm.summary(trace))
    fig, axs = plt.subplots(6,2, figsize=(12,12))
    pm.traceplot(trace, ax=axs)
    fig.savefig('{}/AxionFit_TracePlot.png'.format(inDir))

    # drawFinalSpectra(pdfArrDict)
    # plt.show()


def constructModel(pdfDict):
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

        # Generate array of deterministic variables (per bin)
        det = Tritium*pdfDict['Tritium'] + Bkg*pdfDict['Bkg'] + Axion*pdfDict['Axion'] + Fe55*pdfDict['Fe55'] + Zn65*pdfDict['Zn65'] + Ge68*pdfDict['Ge68']

        L = pm.Poisson("L", mu=det, observed=pdfDict['Data'])
    return model


def convertaxionPDF(startTime, endTime):
    """
        Converts 2D histogram PDF for axion into numpy array
    """
    import ROOT
    inDir = os.environ['LATDIR']+'/data/MCMC'
    fAxion  = ROOT.TFile('{}/Axion_plot_v5_INT_phi_det9.00_Elow0.0_Ehi20.0_Esig4.00_NE200_1Minutes_for_SunTreeUTC_Full_hist_minutes1.root'.format(inDir))

    hday = fAxion.Get('hday')
    hday.RebinX(5) # 5 minute bins
    # Find bins to cut off PDF at
    startBin = hday.GetXaxis().FindBin(startTime)
    endBin = hday.GetXaxis().FindBin(endTime)
    nTimeBins = endBin - startBin + 1
    AxionList = []
    timeBinLowEdge = [hday.GetXaxis().GetBinLowEdge(time)
                        for time in range(hday.GetNbinsX())
                        if time >= startBin and time <= endBin+1]
    for energy in range(hday.GetNbinsY()):
        AxionList.append([hday.GetBinContent(time, energy)
                            for time in range(hday.GetNbinsX())
                            if time >= startBin and time <= endBin])
    fAxion.Close()

    AxionArr = np.array(AxionList)

    timeMask = GenerateLivetimeMask(np.array(timeBinLowEdge), AxionArr.shape)
    # AxionArr = AxionArr*timeMask
    # Roughly normalize so the axion scale is the same as the other PDFs
    NormAxionArr = AxionArr/(np.amax(AxionArr)/2.)*timeMask

    removeMask = np.where(np.any(NormAxionArr, axis=0)==False)[0]
    return NormAxionArr, nTimeBins, timeBinLowEdge, removeMask


def generateGaus(energyList, meansDict):
    """
        Generates Gaussian PDFs for X-ray peaks based off of energy array
    """
    outDict = {}
    for name, mean in meansDict.items():
        outDict.setdefault(name, norm.pdf(energyList, loc=mean, scale=GetSigma(mean, 5)))
    return outDict


def ReduceData(dsNum):
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


def drawPDFs(pdfArrDict):
    """
        Helper function to draw PDFs
        pdfArrDict contains
    """
    energyLabels = np.linspace(0, 20, 5)
    startTime = 1485575988
    endTime = 1489762153
    timeLabels = np.linspace(startTime, endTime, 5, dtype='datetime64[s]')
    timeLabels = np.array(timeLabels, dtype='datetime64[D]')
    fig1, ax1 = plt.subplots(ncols=3, nrows=2, figsize=(15,8))
    sns.heatmap(pdfArrDict['Data'], ax=ax1[0,0])
    ax1[0,0].set_title('DS5 Data')
    ax1[0,0].set_ylabel('Energy')
    ax1[0,0].set_xlabel('Time')
    ax1[0,0].invert_yaxis()
    ax1[0,0].locator_params(nbins=5)
    ax1[0,0].set_yticklabels(energyLabels, rotation=0)
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
    ax1[1,0].locator_params(nbins=5)
    ax1[1,0].set_yticklabels(energyLabels, rotation=0)
    ax1[1,0].set_xticklabels(timeLabels, rotation=50)

    sns.heatmap(pdfArrDict['Zn65'], ax=ax1[1,1])
    ax1[1,1].set_title('Zn65 PDF')
    ax1[1,1].set_xlabel('Time')
    ax1[1,1].invert_yaxis()
    ax1[1,1].locator_params(nbins=5)
    ax1[1,1].set_yticklabels([])
    ax1[1,1].set_xticklabels(timeLabels, rotation=50)

    sns.heatmap(pdfArrDict['Ge68'], ax=ax1[1,2])
    ax1[1,2].set_title('Ge68 PDF')
    ax1[1,2].set_xlabel('Time')
    ax1[1,2].invert_yaxis()
    ax1[1,2].locator_params(nbins=5)
    ax1[1,2].set_yticklabels([])
    ax1[1,2].set_xticklabels(timeLabels, rotation=50)

    plt.tight_layout()
    plt.show()
    # fig1.savefig('{}/BraggFitPDF.png'.format(os.environ['LATDIR']+'/data/MCMC'))


def drawFinalSpectra(pdfDict, trace):
    """
        Draws final spectra using mode of posterior
    """
    pdfSpecDict = {}
    figSpec, axSpec = plt.subplots(figsize=(10,7))

    DataSpec = pdfDict['Data'].sum(axis=1)
    axSpec.plot(DataSpec, label='Data', color='black')

    # Create a dummy for the total spectrum
    totalSpec = np.zeros(DataSpec.shape)

    # Take all the pdfs and flatten time axis
    for key in pdfDict.keys():
        if key == 'Data':
            continue
        tempSpec = trace[key].mean()*pdfDict[key].sum(axis=1)
        pdfSpecDict.setdefault(key, tempSpec)
        axSpec.plot(tempSpec, label=key, alpha=0.5, linestyle='--')
        totalSpec += tempSpec

    axSpec.plot(totalSpec, label='Total Model')
    axSpec.legend()
    plt.tight_layout()
    plt.show()


def GenerateLivetimeMask(timeBinLowEdge, AxionShape):
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

def GetSigma(energy, dsNum=1):
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


if __name__ == '__main__':
    main()
