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
ds = imp.load_source('DataSetInfo',os.environ['LATDIR']+'/sandbox/DataSetInfo.py')
wl = imp.load_source('waveLibs',os.environ['LATDIR']+'/waveLibs.py')
dsi = imp.load_source('dsi',os.environ['LATDIR']+'/dsi.py')
calInfo = ds.CalInfo()
detInfo = dsi.DetInfo()

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

    # Create files with hits
    # getSpecPandas()
    # getSimPandas()
    # loadSpec()
    loadFraction()
    return

    # interceptList, distList, absdistList = loadScatter()
    # print('Mean intercept: ', sum(interceptList)/len(interceptList))
    # print('Mean Distance: ', sum(distList)/len(distList))
    # print('Absolute Distance: ', sum(absdistList)/len(absdistList))

    simMatrix, simdetLabel = loadMatrix('Sim')
    # calMatrix, detLabel = loadMatrix('Cal')

    # Get Total Counts
    # totCounts = np.sum(calMatrix)
    totCountsSim = np.sum(simMatrix)
    # print('Total Counts (Calibration): ', totCounts)
    print('Total Counts (Simulation): ', totCountsSim)

    # Draw heatmap
    # fig1, ax1 = plt.subplots(figsize=(10,8))
    # sns.heatmap(np.round(calMatrix/totCounts*100., 1), cmap="YlGnBu", xticklabels=detLabel, yticklabels=detLabel, annot=True, ax=ax1)
    # ax1.set_title('DS5 Sim Calibration -- M2 Sum 238 Peak -- Percentage of Counts <{} keV'.format(energyCut))
    # ax1.set_xlabel('Low Energy Hit')
    # ax1.set_ylabel('High Energy Hit')
    # plt.xticks(rotation=90)
    # plt.yticks(rotation=0)
    # plt.tight_layout()
    # fig1.savefig('{}/CrazyMatrix_Cal_C1_Percent.png'.format(outDir))
    #
    fig2, ax2 = plt.subplots(figsize=(10,8))
    sns.heatmap(np.round(simMatrix/totCountsSim*100.,1), cmap="YlGnBu", xticklabels=simdetLabel, yticklabels=simdetLabel, annot=True, ax=ax2)
    ax2.set_title('DS5 Sim Calibration -- M2 Sum 238 Peak -- Percentage of Counts <{} keV'.format(energyCut))
    ax2.set_xlabel('Low Energy Hit')
    ax2.set_ylabel('High Energy Hit')
    plt.xticks(rotation=90)
    plt.yticks(rotation=0)
    plt.tight_layout()
    # fig2.savefig('{}/CrazyMatrix_Sim_C1_Percent.png'.format(outDir))


    enrList = detInfo.getGoodChanList(5, mod=1, detType='Enr')
    natList = detInfo.getGoodChanList(5, mod=1, detType='Nat')
    enrCPD = ['C{}P{}D{}'.format(*str(detInfo.getChanCPD(5, i))) for i in enrList]
    natCPD = ['C{}P{}D{}'.format(*str(detInfo.getChanCPD(5, i))) for i in natList]

    simtotCounts = np.sum(simMatrix)
    simchCounts = np.sum(simMatrix, axis=0)
    simchFrac = simchCounts/float(simtotCounts)
    simchFracEnr = np.asarray([simchFrac[i] if simdetLabel[i] in enrCPD else 0. for i in range(len(simchFrac))])
    simchFracNat = np.asarray([simchFrac[i] if simdetLabel[i] in natCPD else 0. for i in range(len(simchFrac))])

    print(simchFracEnr)
    print(simchFracNat)

    xList = np.linspace(1, len(simdetLabel), len(simdetLabel))
    bar_width = 0.35
    fig3, ax3 = plt.subplots(figsize=(12,7))
    ax3.bar(xList, simchFracEnr, bar_width, label='Enriched')
    ax3.bar(xList, simchFracNat, bar_width, label='Natural')
    ax3.set_xticks(xList + bar_width/2)
    ax3.set_xticklabels(simdetLabel, rotation=50)
    ax3.set_title('DS5 Calibration Simulation -- M2 Sum 238 Peak -- Fraction of Events <{} keV'.format(energyCut))
    ax3.set_ylabel("Fraction of Events")
    ax3.legend()
    plt.tight_layout()
    fig3.savefig('{}/CrazyBar_C1_NatEnr.png'.format(outDir))
    plt.show()
    return

    # Bar plot showing fractions for each detector
    chCounts = np.sum(calMatrix, axis=0)
    chFrac = chCounts/float(totCounts)
    errLo, errHi = proportion.proportion_confint(chCounts, np.full(chCounts.shape[0], totCounts), alpha=0.1, method='beta')

    simtotCounts = np.sum(simMatrix)
    simchCounts = np.sum(simMatrix, axis=0)
    simchFrac = simchCounts/float(simtotCounts)
    simerrLo, simerrHi = proportion.proportion_confint(simchCounts, np.full(simchCounts.shape[0], simtotCounts), alpha=0.1, method='beta')

    xList = np.linspace(1, len(detLabel), len(detLabel))
    bar_width = 0.35
    fig2, ax2 = plt.subplots(figsize=(12,7))
    # fig2 = plt.figure(figsize=(15,12))
    # gs = gridspec.GridSpec(2, 1, height_ratios=[3,1])
    # ax2 = plt.subplot(gs[0])
    # ax3 = plt.subplot(gs[1])
    # Add 0.5 to each width to shift to the center of the grid
    ax2.bar(xList, chFrac, bar_width, yerr=[chFrac - errLo, errHi - chFrac], label='Calibration Data')
    ax2.bar(xList+bar_width, simchFrac, bar_width, yerr=[simchFrac - simerrLo, simerrHi - simchFrac], label='Calibration Simulation')
    # ax2.set_xticks(np.linspace(1, len(detLabel), len(detLabel)))
    ax2.set_xticks(xList + bar_width/2)
    ax2.set_xticklabels(detLabel, rotation=50)
    ax2.set_title('DS5 Calibration -- M2 Sum 238 Peak -- Fraction of Events <{} keV'.format(energyCut))
    ax2.set_ylabel("Fraction of Events")
    ax2.legend()
    # sns.regplot(xList + bar_width/2, simchFrac/chFrac, ax=ax3)
    # ax3.set_ylabel('Ration (Sim/Cal)')
    # ax3.set_xticks(xList + bar_width/2)
    # ax3.set_xticklabels(detLabel, rotation=50)
    plt.tight_layout()
    # fig2.savefig('{}/CrazyBar_C1_3.png'.format(outDir))



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
    # inDir = '/projecta/projectdirs/majorana/sim/MJDG41003GAT/MJDemonstrator/linesource/M1CalSource/A224_Z88'
    inDir = '/mnt/mjdDisk1/Majorana/users/psz/MAGE/Processed/5M'
    inPath = pathlib.Path(inDir)
    files = inPath.glob('processed_M1CalSource_A224_Z88_p*.root')
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
            activeNess = [ae.GetElement(i).GetActiveness() for i in idxList]

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
                dataMap['fActiveness1'] = activeNess[0]
                dataMap['fActiveness2'] = activeNess[1]
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
    from scipy import stats
    inDir = os.environ['LATDIR']
    # df = pd.read_hdf('{}/DS5_Cal_HitData.h5'.format(inDir))
    df = pd.read_hdf('{}/DS5_Sim_HitData_G41004.h5'.format(inDir))
    # df['EMirror'] = 238.63 - df['trapENFCal2']
    dfCut = df.loc[(df['sumET'] > 237) & (df['sumET'] < 240) & (df['CPD1'] < 200) & (df['CPD2'] < 200)]

    g1 = g1.map(plt.scatter, 'trapENFCal1', 'EMirror').add_legend()
    g2 = sns.lmplot(x='trapENFCal1', y='EMirror', data=dfCut, col_wrap=5, size=10, col='CPD1').set(xlim=(0, 10), ylim=(0,10))
    g1 = g1.map(plt.hist, 'trapENFCal1', bins=np.linspace(0, 50, 50)).add_legend()

    slope, intercept, r_value, p_value, std_err = stats.linregress(dfCut['trapENFCal1'].loc[dfCut['CPD1']==123].values, dfCut['EMirror'].loc[dfCut['CPD1']==123].values)
    print(slope, intercept, r_value, p_value, std_err)

    # g3 = sns.JointGrid(x="trapENFCal1", y="EMirror", data=dfCut.loc[dfCut['CPD1']==132], xlim=(0,10), ylim=(0,10))
    # g3 = g3.plot_joint(sns.regplot)
    # g3.set_axis_labels('Low Energy Hit (keV)', '238.63 - Pair Hit (keV)')
    # plt.annotate('Slope: {:.2f} \nY-intercept: {:.2f}'.format(slope, intercept), xy=(1, 8))
    # g3 = g3.plot_marginals(sns.distplot, kde=False, bins=np.linspace(0,10,10))
    # g3 = g3.annotate(stats.pearsonr)
    # linregress = lambda x,y: stats.linregress(x,y)
    # g3 = g3.annotate(linregress, template="{stat}: {val:.2f}")
    # g3.set('C1P3D2')
    # g3.savefig('CrazyPlot_3_C1P2D3.png')

    plt.tight_layout()
    plt.show()


def loadScatter():
    """
        Draws various spectra plots and also performs linear regression on scatter plots

    """
    from scipy import stats
    inDir = os.environ['LATDIR']
    df = pd.read_hdf('{}/DS5_Cal_HitData.h5'.format(inDir))
    df['EMirror'] = 238.63 - df['trapENFCal2']
    # df['EMirror'] = 2614.533 - df['trapENFCal2']
    # Make Energy cut
    dfCut = df.loc[(df['sumET'] > 237) & (df['sumET'] < 240) & (df['CPD1'] < 200) & (df['CPD2'] < 200)]
    # dfCut = df.loc[(df['sumET'] > 2610) & (df['sumET'] < 2620) & (df['CPD1'] < 200) & (df['CPD2'] < 200)]

    detList1 = np.unique(dfCut['CPD1'])
    fig1, ax1 = plt.subplots(figsize=(10,7))

    interceptList = []
    distList = []
    absdistList = []
    for CPD in detList1:
        dfCh = dfCut.loc[dfCut['CPD1']==CPD]

        x = dfCh['trapENFCal1'].values
        y = dfCh['EMirror'].values

        slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
        mx = x.mean()
        sx2 = ((x-mx)**2).sum()
        sd_intercept = std_err * np.sqrt(1./len(x) + mx*mx/sx2)
        sd_slope = std_err * np.sqrt(1./sx2)

        print('{} - Slope: {:.3f} +/- {:.3f} Y-int: {:.3f} +/- {:.3f}, std_err: {:.3f}'.format(str(CPD), slope, sd_slope, intercept, sd_intercept, std_err))

        ax1.cla()
        sns.regplot(x='trapENFCal1', y='EMirror', data=dfCh, ax=ax1)
        ax1.set_xlabel('Low Energy Hit (keV)')
        ax1.set_ylabel('238.63 - Pair Hit (keV)')
        # ax1.set_ylabel('2614.53 - Pair Hit (keV)')
        ax1.set_xlim(0, 10)
        ax1.set_ylim(0, 10)
        ax1.set_title('C{}P{}D{} Energy Calibration Check'.format(*str(CPD)))
        ax1.annotate('Slope: {:.2f} \nY-intercept: {:.2f}'.format(slope, intercept), xy=(1, 8))
        plt.tight_layout()
        # fig1.savefig('{}/plots/CalPairs/CalCheck/2615keV_CalCheck_C{}P{}D{}.png'.format(inDir, *str(CPD)))
        interceptList.append(intercept)

        # Find distance between (x,y) and x = y
        xCut = dfCh['trapENFCal1'].loc[dfCh['trapENFCal1'] < 10].values
        yCut = dfCh['EMirror'].loc[dfCh['trapENFCal1'] < 10].values
        absdistance = abs(xCut-yCut)
        distance = xCut-yCut
        if len(distance) > 0:
            distList.append(distance.sum()/len(distance))
            absdistList.append(absdistance.sum()/len(absdistance))
            print('{} - Distance {:.3f} - Absolute Distance {:.3f}'.format(str(CPD), distance.sum()/len(distance), absdistance.sum()/len(absdistance)))
    return interceptList, distList, absdistList


def loadMatrix(dType = 'Sim', chType=None):
    """
        Loops through event list and fills a matrix of hit pairs
    """

    inDir, outDir = os.environ['LATDIR'], os.environ['LATDIR']+'/plots/CalPairs'
    # df = pd.read_hdf('{}/DS5_{}_HitData.h5'.format(inDir, dType))
    df = pd.read_hdf('{}/DS5_{}_HitData_G41004.h5'.format(inDir, dType))
    df['EMirror'] = 238.63 - df['trapENFCal2']
    # Select 238 sum peak, select only module 1 detectors
    # This is here right now because I'm a dumbass and I forgot to scale the sum energy
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


def loadFraction():
    """
        Loops through channels and draws fraction of surface/total events per channel/detector
    """

    inDir = os.environ['LATDIR']
    # df = pd.read_hdf('{}/DS5_Cal_HitData.h5'.format(inDir))
    df1 = pd.read_hdf('{}/DS5_Sim_HitData.h5'.format(inDir))
    df2 = pd.read_hdf('{}/DS5_Sim_HitData_G41004.h5'.format(inDir))
    df = pd.concat([df1, df2])
    # df['EMirror'] = 238.63 - df['trapENFCal2']
    dfCut = df.loc[(df['sumET'] > 237) & (df['sumET'] < 240) & (df['CPD1'] < 200) & (df['CPD2'] < 200)]
    dfSurf1 = dfCut.loc[(dfCut['fActiveness1'] < 1)]
    dfSurf2 = dfCut.loc[(dfCut['fActiveness2'] < 1)]

    totArr = np.concatenate((dfCut['trapENFCal1'].values, dfCut['trapENFCal2'].values))
    surfArr = np.concatenate((dfSurf1['trapENFCal1'].values, dfSurf2['trapENFCal2'].values))

    totHistVals, totHistBins = np.histogram(totArr, bins=50, range=(0, 50))
    surfHistVals, surfHistBins = np.histogram(surfArr, bins=50, range=(0, 50))


    fig1, (ax1, ax2) = plt.subplots(nrows=2, figsize=(10,8))
    ax1.step(totHistBins[:-1]+0.5, totHistVals, label='All Events')
    ax1.step(totHistBins[:-1]+0.5, surfHistVals, label='Transition Layer')
    ax2.step(totHistBins[:-1]+0.5, np.divide(surfHistVals, totHistVals, dtype=float), label='Percentage')
    ax1.set_xlim(0, 50)
    ax2.set_xlim(0, 50)
    ax1.set_title('Module 1 - Fraction of Transition Layer Events')
    ax1.set_ylabel('Counts/keV')
    ax1.set_xlabel('Energy (keV)')
    ax2.set_ylabel('Fraction')
    ax2.set_xlabel('Energy (keV)')
    ax1.legend()
    ax2.legend()

    print(np.sum(totHistVals))
    print(np.sum(surfHistVals))

    # fig1.savefig(inDir + '/plots/DeadLayer/G41003/TransitionFrac_Module2.png')
    # plt.show()

    # Now go through detector lists
    detList1 = np.unique(dfCut['CPD1'])
    detList2 = np.unique(dfCut['CPD2'])
    # Combine detector lists with union
    detListFull = np.union1d(detList1, detList2)

    for det in detListFull:
        ax1.cla()
        ax2.cla()

        dfCh1 = dfCut.loc[dfCut['CPD1'] == det]
        dfCh2 = dfCut.loc[dfCut['CPD2'] == det]
        dfSurfCh1 = dfCh1.loc[dfCh1['fActiveness1'] < 1]
        dfSurfCh2 = dfCh2.loc[dfCh2['fActiveness2'] < 1]

        totArrCh = np.concatenate((dfCh1['trapENFCal1'].values, dfCh2['trapENFCal2'].values))
        surfArrCh = np.concatenate((dfSurfCh1['trapENFCal1'].values, dfSurfCh2['trapENFCal2'].values))

        totHistValsCh, _ = np.histogram(totArrCh, bins=50, range=(0, 50))
        surfHistValsCh, _ = np.histogram(surfArrCh, bins=50, range=(0, 50))
        fracCh = np.divide(surfHistValsCh, totHistValsCh, dtype=float)
        ax1.step(totHistBins[:-1]+0.5, totHistValsCh, label='All Events')
        ax1.step(totHistBins[:-1]+0.5, surfHistValsCh, label='Transition Layer')
        ax2.step(totHistBins[:-1]+0.5, fracCh, label='Percentage')
        ax1.set_xlim(0, 50)
        ax2.set_xlim(0, 50)
        ax1.set_title('C{}P{}D{} - Fraction of Transition Layer Events'.format(*str(det)))
        ax1.set_ylabel('Counts/keV')
        ax1.set_xlabel('Energy (keV)')
        ax2.set_ylabel('Fraction')
        ax2.set_xlabel('Energy (keV)')
        ax1.legend()
        ax2.legend()
        # fig1.savefig(inDir + '/plots/DeadLayer/G41003/TransitionFrac_C{}P{}D{}.png'.format(*str(det)))


if __name__=="__main__":
    main()
