"""
  This script converts ROOT files to Pandas HDF5 format for usage in BigBraggBrand fitter
  Essentially takes the same functions from BigBraggBrand and removes all ROOT dependencies
"""

import sys, os, imp, time, ROOT
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
import seaborn as sns
from scipy.interpolate import interp1d
from scipy.stats import norm
import pandas as pd
sns.set(style='darkgrid')

wl = imp.load_source('waveLibs', '{}/waveLibs.py'.format(os.environ['LATDIR']))
dsi = imp.load_source('dsi', '{}/dsi.py'.format(os.environ['LATDIR']))

bkg = dsi.BkgInfo()
cal = dsi.CalInfo()
det = dsi.DetInfo()

# Define some global parameters
inDir = os.environ['LATDIR']+'/data/MCMC'
dsList = ['5B', '5C', '6']
effList = ['Nominal', 'Lo90', 'IS', 'Sideband']

# The start and end time are for DS5b
startTime = 1485575988
endTime = 1523845062 # DS6a end day (April 16th 2018)
yearInSec = 31536000 # Number of seconds in 1 year
dayInSec = 86400

# Energy range (also defines fitting range!)
# The tritium endpoint is 18 keV, the axion PDF stops at 12 keV.
# We want to extend the PDF to slightly beyond the tritium endpoint in order to capture the flat bkg
energyThresh, energyThreshMax, binSize = 5.0, 20.0, 0.1
energyBins = np.linspace(energyThresh,energyThreshMax, round((energyThreshMax-energyThresh)/binSize)+1)
axionBinSize = 5
bDebug = False
effMode = None

def main():
    # reduceData()
    # generateEfficiency()

    dfList = {}
    # for dsNum in ['5b', '6a']:
    for dsNum in ['5b']:
        AxionArr, timeBinCenter, energyBinCenter = convertaxionPDF(startTime, endTime, bDebug=bDebug, axionBinSize=axionBinSize, dsNum=dsNum, effMode = effMode)

        # The reason we need to take the transpose is because we want the energy bins to be columns and the time bins to be the rows, this is because there is a maximum header size of 64 kB for HDF5 files
        # When we load the files, we will need to take the transpose back
        AxionArr = AxionArr.T

        dfList.setdefault(dsNum, pd.DataFrame(data=AxionArr, columns=energyBinCenter, index=timeBinCenter))
        timeArr = (dfList[dsNum].index.values - axionBinSize/2.*60.).tolist()
        timeArr.append(timeArr[-1]+axionBinSize*60.)
        print(dsNum, dfList[dsNum].shape)
        # print(dfList[dsNum].head(20))
        print('Output',timeArr[0], timeArr[-1])

    # dfList['5b'].to_hdf(inDir + '/AxionAveragedPDFs.h5', key='AAA_DS5b', format='table', complevel=9, mode='w')
    # dfList['6a'].to_hdf(inDir + '/AxionAveragedPDFs.h5', key='AAA_DS6a', format='table', complevel=9, mode='a')



def convertaxionPDF(startTime, endTime, axionBinSize = 5, bDebug = False, dsNum = '5b', effMode = None):
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
    inDir = os.environ['LATDIR']+'/data/MCMC'
    fAxion = ROOT.TFile('{}/Axion_averaged_MJD_DS{}reso_Elow0_Ehi20_minutes1_2017_2018.root'.format(inDir, dsNum))

    hday = fAxion.Get('hday')
    hday.RebinX(axionBinSize) # 5 minute bins -- the PDF is generated in 1 minute bins

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
    if bDebug and bRepeat:
        print('Axion PDF Gap Extension (should be binsize):', timeBinLowEdge[-1], hday.GetXaxis().GetBinLowEdge(0)+yearInSec)
    if bRepeat:
        timeBinLowEdge.extend([hday.GetXaxis().GetBinLowEdge(tBin)+yearInSec
                            for tBin in range(hday.GetNbinsX())
                            if tBin >= 0 and tBin <= endBin2+1]) # nBins+1 overlaps with the 0 and 1st bin

    # Center Bin is just the low edge + bin width/2 -- get rid of the last bin
    timeBinCenter = np.array(timeBinLowEdge[:-1])+axionBinSize/2.*60
    energyBinCenter = energyBins[:-1] + binSize/2

    # Total bins is 1 less than the array length of timeBinLowEdge
    nTimeBins = len(timeBinLowEdge)-1

    # Loop through energy bins, grabbing the bin content of each bin and filling into matrix
    # For the loop, start with 1 and go through NbinsY+1 -- this allows for a range up to 20 keV in the fit
    for energyBin in range(1, hday.GetNbinsY()+1):
        # Apply energy cuts on the low edge of the bin -- rounding here is necessary because of ROOT's stupid precision
        if round(hday.GetYaxis().GetBinLowEdge(energyBin),1) < energyThresh: continue
        if round(hday.GetYaxis().GetBinLowEdge(energyBin),1) > energyThreshMax-0.1: continue
        # print("Axion Energy Bin: ", round(hday.GetYaxis().GetBinLowEdge(energyBin),1))
        AxionList.append([hday.GetBinContent(timeBin, energyBin)
                            for timeBin in range(1, hday.GetNbinsX())
                            if timeBin >= startBin and timeBin <= endBin])
        # Extend Axion histogram
        if bRepeat:
            AxionList[len(AxionList)-1].extend([hday.GetBinContent(timeBin, energyBin)
                                for timeBin in range(1, hday.GetNbinsX())
                                if timeBin <= endBin2+1])

    fAxion.Close()
    AxionArr = np.array(AxionList)
    print('Axion Shape', AxionArr.shape, nTimeBins)
    print(timeBinLowEdge[0], timeBinCenter[0])
    print(timeBinLowEdge[-2], timeBinCenter[-1])
    return AxionArr, timeBinCenter, energyBinCenter


def generateEfficiency():
    fEff = ROOT.TFile(os.environ['LATDIR']+'/data/lat-expo-efficiency_Combined.root')
    effDict = {}
    effDict['Energy'] = energyBins[:-1] + binSize/2
    for ds in dsList:
        for effMode in effList:
            if effMode == 'Nominal':
                hEff = fEff.Get('hDS{}_Norm_Enr'.format(ds))
            else:
                hEff = fEff.Get('hDS{}_Norm_Enr_{}'.format(ds, effMode))
            effArr = np.ones(len(energyBins)-1)
            for eIdx, eBin in enumerate(energyBins[:-1]):
                effArr[eIdx] = hEff.GetBinContent(hEff.FindBin(eBin+binSize/2.))
            effDict.setdefault('DS'+ds + '_'+effMode, effArr)

    df = pd.DataFrame.from_dict(effDict)
    df = df.set_index('Energy')
    print(df.head(10))
    df.to_hdf('{}/CombinedEff.h5'.format(os.environ['LATDIR']+'/data/MCMC'), key='Eff', mode='w', format='table', complevel=9)



def reduceData():
    """
        Saves background and tritium spectrum into pandas dataframes
        This function only needs to be run once!
    """
    inDir, outDir = os.environ['LATDATADIR']+'/bkg/cut/final95', os.environ['LATDIR']+'/data/MCMC'
    inDirTrit = '/projecta/projectdirs/majorana/users/bxyzhu/Axion' # Change this to your own directory!
    bkgFile, tritFile = 'Bkg_Spectrum.h5', 'TritSpec.h5'

    # Check if files already exist
    if os.path.isfile('{}/{}'.format(outDir, bkgFile)):
        print('Error: {} already exists!'.format(bkgFile))
        return
    if os.path.isfile('{}/{}'.format(outDir, tritFile)):
        print('Error: {} already exists!'.format(tritFile))
        return

    import ROOT
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
    df.to_hdf('{}/{}'.format(outDir, bkgFile), 'skimTree', mode='w', format='table')

    f1 = ROOT.TFile('{}/TritSpec.root'.format(inDir))
    tritHist = f1.Get('tritHist')
    tritHist.Scale(1./tritHist.Integral("w"))
    energy = [tritHist.GetBinCenter(xbin) for xbin in range(1, tritHist.GetNbinsX())]
    val = [tritHist.GetBinContent(xbin) for xbin in range(1, tritHist.GetNbinsX())]
    dfTrit = pd.DataFrame({"Energy":energy, "Tritium":val})
    print(dfTrit.head(10))
    dfTrit.to_hdf('{}/{}'.format(outDir, tritFile), 'TritSpec', mode='w', format='table')


if __name__ == '__main__':
    main()
