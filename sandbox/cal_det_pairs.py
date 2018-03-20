# -*- coding: utf-8 -*-
# ¯\_(ツ)_/¯ ¯\_(ツ)_/¯ ¯\_(ツ)_/¯ ¯\_(ツ)_/¯ ¯\_(ツ)_/¯ ¯\_(ツ)_/¯ ¯\_(ツ)_/¯ ¯\_(ツ)_/¯ ¯\_(ツ)_/¯
#!/usr/bin/env python
import sys, os, imp, pathlib, itertools
import numpy as np
import pandas as pd
from statsmodels.stats import proportion
import matplotlib.pyplot as plt
from matplotlib import gridspec
import seaborn as sns
sns.set(style='darkgrid')

# load LAT libraries
# ds = imp.load_source('DataSetInfo',os.environ['LATDIR']+'/DataSetInfo.py')
import DataSetInfo as ds
wl = imp.load_source('waveLibs',os.environ['LATDIR']+'/waveLibs.py')
calInfo = ds.CalInfo()

# load threshold data
import tinydb as db
dsNum, bkgIdx = 5, 83
calDB = db.TinyDB('calDB.json')
pars = db.Query()
thD = ds.getDBRecord("thresh_ds%d_bkgidx%d" % (dsNum, bkgIdx), False, calDB, pars)

# load fitSlo vals for cal run range closest to the run range [22513, 22566]
dsNum, modNum, calIdx = 5, 1, 11  # calIdx 11: [[22568,22635],22568,22841],

def main():

    inDir, outDir = os.environ['LATDIR'], os.environ['LATDIR']+'/plots/CalPairs'
    energyCut = 10

    # getSpecPandas()
    # getSimPandas()
    # loadSpec()
    simMatrix, simdetLabel = loadMatrix('Sim')
    calMatrix, detLabel = loadMatrix('Cal')

    # Draw heatmap
    # fig1, ax1 = plt.subplots(figsize=(10,8))
    # sns.heatmap(calMatrix, cmap="YlGnBu", xticklabels=detLabel, yticklabels=detLabel, ax=ax1)
    # ax1.set_title('DS5 Sim Calibration -- M2 Sum 238 Peak -- Counts <{} keV'.format(energyCut))
    # ax1.set_xlabel('Low Energy Hit')
    # ax1.set_ylabel('High Energy Hit')
    # plt.xticks(rotation=90)
    # plt.yticks(rotation=0)
    # plt.tight_layout()
    # fig1.savefig('{}/CrazyMatrix_Cal_C1.png'.format(outDir))

    # Bar plot showing fractions for each detector
    totCounts = np.sum(calMatrix)
    chCounts = np.sum(calMatrix, axis=1)
    chFrac = chCounts/float(totCounts)
    errLo, errHi = proportion.proportion_confint(chCounts, np.full(chCounts.shape[0], totCounts), alpha=0.1, method='beta')

    simtotCounts = np.sum(simMatrix)
    simchCounts = np.sum(simMatrix, axis=1)
    simchFrac = simchCounts/float(simtotCounts)
    simerrLo, simerrHi = proportion.proportion_confint(simchCounts, np.full(simchCounts.shape[0], simtotCounts), alpha=0.1, method='beta')

    xList = np.linspace(1, len(detLabel), len(detLabel))
    bar_width = 0.35
    # fig2, ax2 = plt.subplots(figsize=(12,7))
    fig2 = plt.figure(figsize=(15,12))
    gs = gridspec.GridSpec(2, 1, height_ratios=[3,1])
    ax2 = plt.subplot(gs[0])
    ax3 = plt.subplot(gs[1])
    # Add 0.5 to each width to shift to the center of the grid
    ax2.bar(xList, chFrac, bar_width, yerr=[chFrac - errLo, errHi - chFrac], label='Calibration Data')
    ax2.bar(xList+bar_width, simchFrac, bar_width, yerr=[simchFrac - simerrLo, simerrHi - simchFrac], label='Calibration Simulation')
    # ax2.set_xticks(np.linspace(1, len(detLabel), len(detLabel)))
    ax2.set_xticks(xList + bar_width/2)
    ax2.set_xticklabels(detLabel, rotation=50)
    ax2.set_title('DS5 Calibration -- M2 Sum 238 Peak -- Fraction of Events <{} keV'.format(energyCut))
    ax2.set_ylabel("Fraction of Events")
    ax2.legend()
    sns.regplot(xList + bar_width/2, simchFrac/chFrac, ax=ax3)
    ax3.set_ylabel('Ration (Sim/Cal)')
    ax3.set_xticks(xList + bar_width/2)
    ax3.set_xticklabels(detLabel, rotation=50)
    plt.tight_layout()
    fig2.savefig('{}/CrazyBar_C1_2.png'.format(outDir))

def convertDF(series):
    """
        Converts a list stored in DataFrame column back into a np.array
        The reason we need to do this is because series.values returns a np.array (of a list)
        So the formatting is wacky
    """
    try:
        return np.asarray(series.values.tolist())
    except:
        print('Input is not a pandas.core.series.Series object! Returning')
        return 0

def simtoCPD(simID):
    """
        Short helper function that converts fWaveformID to CPD
        fWaveformID format: XCCPPDD
    """
    D = simID%10
    P = (simID%1000 - D)/100
    C = (simID%100000 - P*100 - D)/10000
    return int(C*100+P*10+D)

def simtoCh(simID):
    """
        Short helper function that converts fWaveformID to Channel
    """
    CPD, ch = simtoCPD(simID), 0
    for key, val in ds.CPD[dsNum].items():
        if val == CPD and key%2 == 0:
            ch = key
    return int(ch)

def getSpecPandas():
    """
        Get sum and hit spectra w/ threshold cut, without ch. 598 (it's noisy.)
        Here we use "mHT" and "sumET" exclusively.
        Use !EventDC1Bits and trapENFCal > 0.7, all hits.  (Skim file was generated w/ 'dontSkipAnything')
        Save detailed info for events with hits < 2630 keV.
    """
    from ROOT import TFile, TTree
    # inDir, outDir = '/mnt/mjdDisk1/Majorana/data/sandbox/special/lat/cal','/mnt/mjdDisk1/Majorana/users/psz/LAT/data'
    inDir = '/projecta/projectdirs/majorana/users/wisecg/special/lat'
    inPath = pathlib.Path(inDir)
    calInfo = ds.CalInfo()
    runList = calInfo.GetSpecialRuns("longCal",dsNum) # 54 total
    fileList = []
    for run in runList:
        files = inPath.glob('latSkimDS{}_run{}_*.root'.format(dsNum,run))
        fileList.extend([str(f) for f in files]) # Because ROOT doesn't play nice with python paths
    nFiles = len(fileList)
    chList = ds.GetGoodChanListNew(dsNum = dsNum)

    dataList = []
    for iFile, f in enumerate(fileList):
        print("{}/{} {}".format(iFile,nFiles,f))
        tf = TFile(f)
        lTree = tf.Get("skimTree")
        lTree.GetEntry(0)

        for iEnt in range(lTree.GetEntries()):
            lTree.GetEntry(iEnt)
            # Cut pulsers
            if lTree.EventDC1Bits != 0: continue

            idxList = [i for i in range(lTree.channel.size())
                if lTree.channel.at(i) in chList
                and lTree.channel.at(i) != 598
                and lTree.channel.at(i) in thD
                and lTree.trapENFCal.at(i) > thD[lTree.channel.at(i)][0] + 3*thD[lTree.channel.at(i)][1]
                and 0.7 < lTree.trapENFCal.at(i) < 9999
                ]
            if len(idxList) is not 2: continue
            hitE = [lTree.trapENFCal.at(i) for i in idxList]
            mHT, sumET = len(hitE), sum(hitE)

            # Store everything as its own column
            # This is because I'm lazy AF and I didn't want to work with arrays in my dataframe
            # Also, can store this as 'table' here
            if sumET < 2630:
                dataMap = {}
                dataMap['run'] = int(lTree.run)
                dataMap['channel1'] = int(lTree.channel.at(idxList[0]))
                dataMap['channel2'] = int(lTree.channel.at(idxList[1]))
                dataMap['trapENFCal1'] = hitE[0]
                dataMap['trapENFCal2'] = hitE[1]
                dataMap['mHT'] = mHT
                dataMap['sumET'] = sumET
                dataMap['iEvent'] = int(lTree.iEvent)
                dataMap['fitSlo1'] = lTree.fitSlo.at(idxList[0])
                dataMap['fitSlo2'] = lTree.fitSlo.at(idxList[1])
                dataMap['riseNoise1'] = lTree.riseNoise.at(idxList[0])
                dataMap['riseNoise2'] = lTree.riseNoise.at(idxList[1])
                dataMap['CPD1'] = int(100*lTree.C.at(idxList[0])+10*lTree.P.at(idxList[0])+lTree.D.at(idxList[0]))
                dataMap['CPD2'] = int(100*lTree.C.at(idxList[1])+10*lTree.P.at(idxList[1])+lTree.D.at(idxList[1]))
                dataMap['UnixTime'] = lTree.globalTime.GetSec()
                dataList.append(dataMap)

    df = pd.DataFrame.from_dict(dataList)
    print(len(df))
    print(df.head(10))

    # If we want to set columns manually
    # dfCols = ['run', 'channel1', 'channel2', 'trapENFCal1', 'trapENFCal2', 'mHT', 'sumET', 'iEvent', 'fitSlo1', 'fitSlo2', 'riseNoise1', 'riseNoise2', 'CPD1', 'CPD2', 'UnixTime']

    # Write Dataframe to file -- suppress performance warnings
    # There's not enough events to warrant chunk-writing
    import warnings
    warnings.filterwarnings(action="ignore", module="pandas", message="^\nyour performance")
    df.to_hdf('DS{}_Cal_HitData.h5'.format(dsNum), key='skimTree', format = 'table', mode = 'w', complevel=9)

def getSimPandas():
    """
        Get sum and hit events from simulation data
        Here we use the same thresholds and channel selection as the data to be consistent:
        Use threshold cut and without ch 598
        Save detailed info for events with hits < 2630 keV.
    """
    from ROOT import TFile, TTree
    inDir = '/projecta/projectdirs/majorana/sim/MJDG41003GAT/MJDemonstrator/linesource/M1CalSource/A224_Z88'
    inPath = pathlib.Path(inDir)
    files = inPath.glob('processed_MJDemonstrator_linesource_A224_Z88_from_A224_Z88_to_A208_Z81_in_M1CalSource_500000_*.root')
    fileList = []

    fileList.extend([str(f) for f in files]) # Because ROOT doesn't play nice with python paths
    # auxfileList.extend([str(f) for f in auxfiles])
    nFiles = len(fileList)
    print('Total Files: ', nFiles)
    chList = ds.GetGoodChanListNew(dsNum = dsNum)
    cpdList = [ds.CPD[dsNum][ch] for ch in chList]
    dataList = []
    for iFile, f in enumerate(fileList):
        print("{}/{} {}".format(iFile,nFiles,f))
        tf = TFile(f)
        lTree = tf.Get("simTree")
        lTree.GetEntry(0)

        for iEnt in range(lTree.GetEntries()):
            lTree.GetEntry(iEnt)
            ae = lTree.fAnalysisEvent
            sumET = ae.GetTotalEnergy()*1000
            mHT = ae.GetNElements()

            # Initial cut to speed things up
            if mHT != 2: continue

            # Apply threshold and channel selection cuts
            idxList = [i for i in range(mHT)
                if simtoCPD(ae.GetElement(i).GetWaveformID()) in cpdList
                and simtoCPD(ae.GetElement(i).GetWaveformID()) != ds.CPD[5][598]
                and simtoCh(ae.GetElement(i).GetWaveformID()) in thD
                and ae.GetElement(i).GetEnergy()*1000 > thD[simtoCh(ae.GetElement(i).GetWaveformID())][0] + 3*thD[simtoCh(ae.GetElement(i).GetWaveformID())][1]
                and 0.7 < ae.GetElement(i).GetEnergy()*1000 < 9999
                ]
            # If two hits don't remain, skip
            if len(idxList) is not 2: continue
            hitE = [ae.GetElement(i).GetEnergy()*1000 for i in idxList]
            chans = [simtoCh(ae.GetElement(i).GetWaveformID()) for i in idxList]
            cpds = [simtoCPD(ae.GetElement(i).GetWaveformID()) for i in idxList]

            if sumET < 2630:
                dataMap = {}
                dataMap['channel1'] = chans[0]
                dataMap['channel2'] = chans[1]
                dataMap['trapENFCal1'] = hitE[0]
                dataMap['trapENFCal2'] = hitE[1]
                dataMap['mHT'] = mHT
                dataMap['sumET'] = sumET
                dataMap['CPD1'] = cpds[0]
                dataMap['CPD2'] = cpds[1]
                dataList.append(dataMap)

    df = pd.DataFrame.from_dict(dataList)
    print(len(df))
    print(df.head(10))

    # Write Dataframe to file -- suppress performance warnings
    # There's not enough events to warrant chunk-writing
    import warnings
    warnings.filterwarnings(action="ignore", module="pandas", message="^\nyour performance")
    df.to_hdf('DS{}_Sim_HitData.h5'.format(dsNum), key='skimTree', format = 'table', mode = 'w', complevel=9)

def loadSpec():
    """
        Draws various spectra plots and also performs linear regression on scatter plots
    """

    inDir = os.environ['LATDIR']
    df = pd.read_hdf('{}/DS5_LongCal_HitData.h5'.format(inDir))
    df['EMirror'] = 238.63 - df['trapENFCal2']
    dfCut = df.loc[(df['sumET'] > 237) & (df['sumET'] < 240) & (df['CPD1'] < 200) & (df['CPD2'] < 200)]
    # g1 = sns.FacetGrid(df.loc[(df['sumET'] > 237) & (df['sumET'] < 240) & (df['CPD1'] < 200) & (df['CPD2'] < 200)], col='CPD1', hue='CPD2', col_wrap=5, size=10, aspect=1.2).set(xlim=(0, 50), ylim=(0,50))
    # g1 = g1.map(plt.scatter, 'trapENFCal1', 'EMirror').add_legend()
    # g2 = sns.lmplot(x='trapENFCal1', y='EMirror', data=dfCut, col_wrap=5, size=10, col='CPD1').set(xlim=(0, 10), ylim=(0,10))
    # g1 = g1.map(plt.hist, 'trapENFCal1', bins=np.linspace(0, 50, 50)).add_legend()

    # g2.savefig('CrazyPlot_2_M1_Zoomed.png')
    from scipy import stats
    g3 = sns.JointGrid(x="trapENFCal1", y="EMirror", data=dfCut.loc[dfCut['CPD1']==132], xlim=(0,10), ylim=(0,10))
    g3 = g3.plot_joint(sns.regplot)
    g3 = g3.plot_marginals(sns.distplot, kde=False, bins=np.linspace(0,10,10))
    g3 = g3.annotate(stats.pearsonr)
    slope, intercept, r_value, p_value, std_err = stats.linregress(dfCut['trapENFCal1'].loc[dfCut['CPD1']==123].values, dfCut['EMirror'].loc[dfCut['CPD1']==123].values)

    print(slope, intercept, r_value, p_value, std_err)
    # linregress = lambda x,y: stats.linregress(x,y)
    # g3 = g3.annotate(linregress, template="{stat}: {val:.2f}")
    # g3.set('C1P3D2')
    # g3.savefig('CrazyPlot_3_C1P2D3.png')
    # plt.show()

def loadMatrix(dType = 'Cal'):
    """
        Loops through event list and fills a matrix of hit pairs
    """
    inDir, outDir = os.environ['LATDIR'], os.environ['LATDIR']+'/plots/CalPairs'
    df = pd.read_hdf('{}/DS5_{}_HitData.h5'.format(inDir, dType))
    df['EMirror'] = 238.63 - df['trapENFCal2']
    # Select 238 sum peak, select only module 1 detectors
    # This is here right now because I'm a dumbass and I forgot to scale the sum energy
    lowCut, highCut = 0,0
    if dType == 'Sim':
        lowCut, highCut = 0.237, 0.240
    elif dType == 'Cal':
        lowCut, highCut = 237, 240
    dfCut = df.loc[(df['sumET'] > lowCut) & (df['sumET'] < highCut) & (df['CPD1'] < 200) & (df['CPD2'] < 200)]
    # Select module 2 detectors
    # dfCut = df.loc[(df['sumET'] > 237) & (df['sumET'] < 240) & (df['CPD1'] > 200) & (df['CPD2'] > 200)]

    energyCut = 10.

    # Look at detector lists
    detList1 = np.unique(dfCut['CPD1'])
    detList2 = np.unique(dfCut['CPD2'])

    # Combine detector lists with union
    detListFull = np.union1d(detList1, detList2)
    print(detList1)
    print(detList2)
    print(detListFull)
    detLabel = ['C{}P{}D{}'.format(*str(i)) for i in detListFull]

    # Build Matrix of CPD1 vs CPD2
    countMatrix = np.zeros((len(detListFull), len(detListFull)))

    # Loop through event list:
    for index, row in df.iterrows():
        # Get respective indices in matrix
        idx1 = np.where(detListFull == row['CPD1'])
        idx2 = np.where(detListFull == row['CPD2'])

        # If energy condition is met, increment matrix element (count)
        # Row = low, Column = high
        if row['trapENFCal1'] < energyCut:
            countMatrix[idx1, idx2] += 1
        if row['trapENFCal2'] < energyCut:
            countMatrix[idx2, idx1] += 1

    return countMatrix, detLabel

if __name__=="__main__":
    main()
