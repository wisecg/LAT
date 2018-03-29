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
    startTime = 148557598
    endTime = 1489762153

    inDir = os.environ['LATDIR']+'/data/MCMC'
    df = pd.read_hdf('{}/DS{}_Spectrum.h5'.format(inDir, dsNum))
    df.sort_values(by=['UnixTime'], inplace=True)
    dataArr = df.loc[:, ['Energy','UnixTime']].values
    # timeBins = np.linspace(startTime, endTime, (endTime-startTime)/300)
    # nTimeBins = 105120
    nTimeBins = 7008
    timeBins = np.linspace(startTime, endTime, nTimeBins)
    energyBins = np.arange(0,20.1,0.1)
    pdfArrDict['Data'], xedges, yedges = np.histogram2d(dataArr[:,0], dataArr[:,1], bins=[energyBins, nTimeBins])


    # Build Tritium PDF -- need interpolation because it's 0.2 keV bins!
    dfTrit = pd.read_hdf('{}/TritSpec.h5'.format(inDir))
    tritEnergy = dfTrit['Energy'].values
    tritSpec = dfTrit['Tritium'].values
    tritSpec = np.array([100*x if x >= 0. else 0. for x in tritSpec])
    tritInterp = interp1d(tritEnergy, tritSpec, fill_value='extrapolate')
    tritList = [tritInterp(x) for x in energyBins[:-1]]
    pdfArrDict['Tritium'] = np.array(tritList*nTimeBins).reshape(nTimeBins, len(tritList)).T

    # Build Axion PDF
    pdfArrDict['Axion'] = convertaxionPDF()

    print('Generated all PDFs')

    meansDict = {'Fe55':6.54, 'Zn65':8.98, 'Ge68':10.37}
    gausDict = generateGaus(energyBins[:-1], meansDict)
    for name, Arr in gausDict.items():
        pdfArrDict.setdefault(name, np.array(Arr.tolist()*nTimeBins).reshape(nTimeBins, len(Arr)).T)

    # Draw PDFs
    # drawPDFs(pdfArrDict)

    # Flatten 2D arrays into 1D
    pdfFlatDict['Data'] = pdfArrDict['Data'].flatten()
    pdfFlatDict['Tritium'] = pdfArrDict['Tritium'].flatten()
    pdfFlatDict['Axion'] = pdfArrDict['Axion'].flatten()
    pdfFlatDict['Bkg'] = np.ones(pdfFlatDict['Data'].shape)
    pdfFlatDict['Fe55'] = pdfArrDict['Fe55'].flatten()
    pdfFlatDict['Zn65'] = pdfArrDict['Zn65'].flatten()
    pdfFlatDict['Ge68'] = pdfArrDict['Ge68'].flatten()
    print('Flattened all PDFs')

    model = constructModel(pdfFlatDict)
    print('Built model')

    with model:
        trace = pm.sample(draws=1000, chains=1, n_init=100)
    dfTrace = pm.trace_to_dataframe(trace, include_transformed=True)
    dfTrace.to_hdf('{}/AxionFitTrace_DS{}.h5'.format(inDir, dsNum), "skimTree", mode='w', format='table')
    pm.traceplot(trace)
    print(pm.summary(trace))


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


def convertaxionPDF():
    """
        Converts 2D histogram PDF for axion into numpy array
    """
    import ROOT
    inDir = os.environ['LATDIR']+'/data/MCMC'
    fAxion  = ROOT.TFile('{}/Axion_plot_v5_INT_phi_det9.00_Elow0.0_Ehi20.0_Esig4.00_NE200_Time525600_for_SunTreeUTC_Neg_hist_minutes1.root'.format(inDir))

    hday = fAxion.Get('hday')
    hday.RebinX(15) # Every 15 minutes
    AxionList = []
    for energy in range(hday.GetNbinsY()):
        AxionList.append([hday.GetBinContent(time, energy) for time in range(hday.GetNbinsX())])

    fAxion.Close()
    return np.array(AxionList)


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

    fig1, ax1 = plt.subplots(ncols=3, nrows=2, figsize=(15,8))
    sns.heatmap(pdfArrDict['Data'], ax=ax1[0,0])
    ax1[0,0].set_title('DS5 Data')
    ax1[0,0].set_ylabel('Energy')
    ax1[0,0].set_xlabel('Time')
    ax1[0,0].invert_yaxis()
    ax1[0,0].set_yticklabels([])
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
    ax1[1,0].set_yticklabels([])
    ax1[1,0].set_xticklabels([])

    sns.heatmap(pdfArrDict['Zn65'], ax=ax1[1,1])
    ax1[1,1].set_title('Zn65 PDF')
    ax1[1,1].set_xlabel('Time')
    ax1[1,1].invert_yaxis()
    ax1[1,1].set_yticklabels([])
    ax1[1,1].set_xticklabels([])

    sns.heatmap(pdfArrDict['Ge68'], ax=ax1[1,2])
    ax1[1,2].set_title('Ge68 PDF')
    ax1[1,2].set_xlabel('Time')
    ax1[1,2].invert_yaxis()
    ax1[1,2].set_yticklabels([])
    ax1[1,2].set_xticklabels([])

    plt.tight_layout()
    plt.show()
    # fig1.savefig('{}/BraggFitPDF.png'.format(inDir))


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