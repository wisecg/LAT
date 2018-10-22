#!/usr/bin/env python
import sys, os, imp, pathlib, itertools
import numpy as np
import pandas as pd
from statsmodels.stats import proportion
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from matplotlib import gridspec
import seaborn as sns
# sns.set(style='darkgrid')
# plt.style.use('../pltReports.mplstyle')

# load LAT libraries
ds = imp.load_source('DataSetInfo',os.environ['LATDIR']+'/sandbox/DataSetInfo.py')
wl = imp.load_source('waveLibs',os.environ['LATDIR']+'/waveLibs.py')
dsi = imp.load_source('dsi',os.environ['LATDIR']+'/dsi.py')
calInfo = ds.CalInfo()
detInfo = dsi.DetInfo()

# load threshold data
import tinydb as db
dsNum, bkgIdx = 5, 83
dbFile = '%s/calDB-v2.json' % (dsi.latSWDir)
# calDB = db.TinyDB('calDB.json')
calDB = db.TinyDB(dbFile)
pars = db.Query()
thD = ds.getDBRecord("thresh_ds%d_bkgidx%d" % (dsNum, bkgIdx), False, calDB, pars)

# load fitSlo vals for cal run range closest to the run range [22513, 22566]
dsNum, modNum, calIdx = 5, 1, 11  # calIdx 11: [[22568,22635],22568,22841],

def main():

    inDir, outDir = os.environ['LATDIR'], os.environ['LATDIR']+'/plots/CalPairs'
    energyCut = 10

    # Create files with hits
    # getSpecPandas()
    # plotM1Spectra()

    # Run these two functions to combine all Enr/Nat detector efficiencies
    #and to calculate uncertainties
    combinedEff(bSave = False, bWriteDB = False)
    # plotMCEffUnc(bSave = False)

    # plotDetSpectra()
    # getSimPandas()
    # loadSpec()
    # loadFraction()
    # plotFraction()
    # plotComptonEdge()
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


def getSpecPandas(dsNum):
    """
        Get sum and hit spectra w/ threshold cut
        Here we use "mHT" and "sumET" exclusively.
        Use !EventDC1Bits and trapENFCal > 0.7, all hits.  (Skim file was generated w/ 'dontSkipAnything')
        Save detailed info for events with hits < 2630 keV.
    """
    from ROOT import TFile, TTree
    inDir, outDir = '/mnt/mjdDisk1/Majorana/data/sandbox/special/lat/cal','/mnt/mjdDisk1/Majorana/users/psz/LAT/data'
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
                # and lTree.channel.at(i) != 598 # This was during DS5, I don't think we need it anymore
                and lTree.channel.at(i) in thD
                # and lTree.trapENFCal.at(i) > thD[lTree.channel.at(i)][0] + 3*thD[lTree.channel.at(i)][1] # Use DB value
                and lTree.trapENFCal.at(i) > lTree.threshKeV.at(i) + 3*lTree.threshSigma.at(i)  # Use threshold from current run
                and 0.7 < lTree.trapENFCal.at(i) < 3000
                ]
            if len(idxList) is not 2: continue
            hitE = [lTree.trapENFCal.at(i) for i in idxList]
            mHT, sumET = len(hitE), sum(hitE)

            # Store everything as its own column
            # This is because I'm lazy AF and I didn't want to work with arrays in my dataframe
            # Also, can store this as 'table' here
            if sumET < 2630:
                dataMap = {}
                dataMap['dsNum'] = int(dsNum)
                dataMap['key'] = key
                dataMap['run'] = int(lTree.run)
                # dataMap['iEvent'] = int(lTree.iEvent)
                dataMap['CPD1'] = int(100*lTree.C.at(idxList[0])+10*lTree.P.at(idxList[0])+lTree.D.at(idxList[0]))
                dataMap['CPD2'] = int(100*lTree.C.at(idxList[1])+10*lTree.P.at(idxList[1])+lTree.D.at(idxList[1]))
                dataMap['channel1'] = int(lTree.channel.at(idxList[0]))
                dataMap['channel2'] = int(lTree.channel.at(idxList[1]))
                dataMap['trapENFCal1'] = hitE[0]
                dataMap['trapENFCal2'] = hitE[1]
                dataMap['sumET'] = sumET
                dataMap['fitSlo1'] = lTree.fitSlo.at(idxList[0])
                dataMap['fitSlo2'] = lTree.fitSlo.at(idxList[1])
                dataMap['riseNoise1'] = lTree.riseNoise.at(idxList[0])
                dataMap['riseNoise2'] = lTree.riseNoise.at(idxList[1])
                dataMap['ToE1'] = lTree.kvorrT.at(idxList[0])/hitE[0]
                dataMap['ToE2'] = lTree.kvorrT.at(idxList[1])/hitE[1]
                dataList.append(dataMap)
        tf.Close()

    df = pd.DataFrame.from_dict(dataList)
    print(len(df))
    print(df.head(10))

    # If we want to set columns manually
    # dfCols = ['run', 'channel1', 'channel2', 'trapENFCal1', 'trapENFCal2', 'mHT', 'sumET', 'iEvent', 'fitSlo1', 'fitSlo2', 'riseNoise1', 'riseNoise2', 'CPD1', 'CPD2', 'UnixTime']

    # Write Dataframe to file -- suppress performance warnings
    # There's not enough events to warrant chunk-writing
    import warnings
    warnings.filterwarnings(action="ignore", module="pandas", message="^\nyour performance")

    chunksize = 100000
    start = 0
    end = chunksize-1
    i = 0
    # for i in len(df):
    if df.shape[0] > end:
        while end < df.shape[0]:
            chunk = df.iloc[start:end]
            try:
                print('Writing to: ', '{}/DS{}_Cal_HitData_{}.h5'.format(outDir,dsNum,i))
                chunk.to_hdf('{}/DS{}_Cal_HitData_{}.h5'.format(outDir,dsNum,i), key='skimTree',
                format='table', mode='w', complevel=9)
            except (Exception) as e:
                print (e)
                print (chunk)
                print (chunk.info())
            start += chunksize
            end += chunksize
            i += 1
    else:
        df.to_hdf('DS{}_Cal_HitData_0.h5'.format(dsNum), key='skimTree', format = 'table', mode = 'w', complevel=9)


def getSimPandas(mod=1):
    """
        Get sum and hit events from simulation data
        Here we use the same thresholds and channel selection as the data to be consistent:
        Use threshold cut and without ch 598
        Save detailed info for events with hits < 2630 keV.
    """
    from ROOT import TFile, TTree
    # inDir = '/projecta/projectdirs/majorana/sim/MJDG41003GAT/MJDemonstrator/linesource/M1CalSource/A224_Z88'
    inDir, outDir = '/mnt/mjdDisk1/Majorana/users/psz/MAGE/Processed/5M', '/mnt/mjdDisk1/Majorana/users/psz/LAT/data'
    inPath = pathlib.Path(inDir)
    files = inPath.glob('DS5_processed_M{}CalSource_A224_Z88_p*.root'.format(mod))
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
            # Update: 8/9/2018 -- Since Micah already added threshold values, don't use threshold here
            idxList = [i for i in range(mHT)
                if simtoCPD(ae.GetElement(i).GetWaveformID()) in cpdList
                # and simtoCPD(ae.GetElement(i).GetWaveformID()) != ds.CPD[5][598]
                and simtoCh(ae.GetElement(i).GetWaveformID()) in thD
                # and ae.GetElement(i).GetEnergy()*1000 > thD[simtoCh(ae.GetElement(i).GetWaveformID())][0] + 3*thD[simtoCh(ae.GetElement(i).GetWaveformID())][1]
                and 0.7 < ae.GetElement(i).GetEnergy()*1000 < 9999 # Basic energy cut
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
    import warnings
    warnings.filterwarnings(action="ignore", module="pandas", message="^\nyour performance")

    chunksize = 10000000
    start = 0
    end = chunksize-1
    i = 0
    # for i in len(df):
    if df.shape[0] > end:
        while end < df.shape[0]:
            chunk = df.iloc[start:end]
            try:
                print('Writing to: ', '{}/DS5_M1_SimHitData_{}.h5'.format(outDir, i))
                chunk.to_hdf('{}/DS5_M1_SimHitData_{}.h5'.format(outDir, i), key='skimTree',
                format='table', mode='w', complevel=9)
            except (Exception) as e:
                print (e)
                print (chunk)
                print (chunk.info())
            start += chunksize
            end += chunksize
            i += 1
    else:
        df.to_hdf('DS5_M1_SimHitData_0.h5', key='skimTree', format = 'table', mode = 'w', complevel=9)

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
    makePlots = False

    inDir = os.environ['LATDIR'] + "/data"
    # df = pd.read_hdf('{}/DS5_Cal_HitData.h5'.format(inDir))
    df1 = pd.read_hdf('{}/DS5_Sim_HitData.h5'.format(inDir))
    df2 = pd.read_hdf('{}/DS5_Sim_HitData_G41004.h5'.format(inDir))
    df = pd.concat([df1, df2])
    # df = df1
    # df['EMirror'] = 238.63 - df['trapENFCal2']
    dfCut = df.loc[(df['sumET'] > 237) & (df['sumET'] < 240) & (df['CPD1'] < 200) & (df['CPD2'] < 200)]
    dfSurf1 = dfCut.loc[(dfCut['fActiveness1'] < 1)]
    dfSurf2 = dfCut.loc[(dfCut['fActiveness2'] < 1)]

    totArr = np.concatenate((dfCut['trapENFCal1'].values, dfCut['trapENFCal2'].values))
    surfArr = np.concatenate((dfSurf1['trapENFCal1'].values, dfSurf2['trapENFCal2'].values))

    # xLo, xHi, xpb = 0, 50, 1 # acceptance study
    xLo, xHi, xpb = 0, 250, 1 # compton edge study
    # xLo, xHi, xpb = 0, 250, 2 # compton edge study, 2 kev bins

    xTot, hTot = wl.GetHisto(totArr, xLo, xHi, xpb, shift=False)
    xSurf, hSurf = wl.GetHisto(surfArr, xLo, xHi, xpb, shift=False)
    hFrac = np.divide(hSurf, hTot, dtype=float)

    if makePlots:
        fig1, (ax1, ax2) = plt.subplots(nrows=2, figsize=(8,6))

        ax1.plot(xTot, hTot, ls='steps', c='k', label="All Events")
        ax1.plot(xSurf, hSurf, ls='steps', c='r', label="Transition Events")
        ax1.axvline(1, c='b', lw=1, label="1.0 keV")
        ax2.plot(xTot, hFrac, ls='steps', c='r')
        ax2.axvline(1, c='b', lw=1, label="1.0 keV")

        ax1.set_xlim(0, 50)
        ax2.set_xlim(0, 50)
        # ax1.set_title('Module 1 - Fraction of Transition Layer Events')
        ax1.set_ylabel('Counts / %.1f keV' % xpb, ha='right', y=1)
        # ax1.set_xlabel('Energy (keV)', ha='right', x=1)
        ax2.set_ylabel('Slow Fraction', ha='right', y=1)
        ax2.set_xlabel('Energy (keV)', ha='right', x=1)
        ax1.legend()
        ax2.legend()
        plt.tight_layout()
        # plt.show()
        plt.savefig('../plots/sim-sloFrac-allDets.pdf')
        # fig1.savefig(inDir + '/plots/DeadLayer/G41003/TransitionFrac_Module2.png')


    # === Now go through detector lists ===
    detList1 = np.unique(dfCut['CPD1'])
    detList2 = np.unique(dfCut['CPD2'])
    # Combine detector lists with union
    detListFull = np.union1d(detList1, detList2)

    # save the histos into this dict
    simHists = {}

    plt.close()
    fig1, (ax1, ax2) = plt.subplots(nrows=2, figsize=(8,6))

    for det in detListFull:
        print("Scanning det", det)
        ax1.cla()
        ax2.cla()

        dfCh1 = dfCut.loc[dfCut['CPD1'] == det]
        dfCh2 = dfCut.loc[dfCut['CPD2'] == det]
        dfSurfCh1 = dfCh1.loc[dfCh1['fActiveness1'] < 1]
        dfSurfCh2 = dfCh2.loc[dfCh2['fActiveness2'] < 1]
        totArrDet = np.concatenate((dfCh1['trapENFCal1'].values, dfCh2['trapENFCal2'].values))
        surfArrDet = np.concatenate((dfSurfCh1['trapENFCal1'].values, dfSurfCh2['trapENFCal2'].values))

        # xTot, hTot = wl.GetHisto(totArrDet, xLo, xHi, xpb, shift=False)
        # xSurf, hSurf = wl.GetHisto(surfArrDet, xLo, xHi, xpb, shift=False)
        # hFrac = np.divide(hSurf, hTot, dtype=float)

        # simHists[str(det)] = [hTot, hSurf, xTot] # save output

        if makePlots:

            ax1.plot(xTot, hTot, ls='steps', c='k', label="All Events, C{}P{}D{}".format(*str(det)))
            ax1.plot(xSurf, hSurf, ls='steps', c='r', label="Transition Events")
            ax1.axvline(1, c='b', lw=1, label="1.0 keV")
            ax2.plot(xTot, hFrac, ls='steps', c='r')
            ax2.axvline(1, c='b', lw=1, label="1.0 keV")

            ax1.set_xlim(0, 50)
            ax2.set_xlim(0, 50)
            ax1.set_ylabel('Counts / %.1f keV' % xpb, ha='right', y=1)
            ax2.set_ylabel('Slow Fraction', ha='right', y=1)
            ax2.set_xlabel('Energy (keV)', ha='right', x=1)
            ax1.legend()
            ax2.legend()
            plt.tight_layout()
            plt.show()
            # fig1.savefig(inDir + '/plots/DeadLayer/G41003/TransitionFrac_C{}P{}D{}.png'.format(*str(det)))
            return

    # save hists for each detector, and the overall total, into one npz file
    # np.savez("../data/efficiency-corr2.npz", hTot, hSurf, xTot, simHists)
    np.savez("../data/efficiency-corr250-2kev.npz", hTot, hSurf, xTot, simHists)


def plotFraction():

    f = np.load('../data/efficiency-corr2.npz')
    hTot, hSurf, xTot, simHists = f['arr_0'], f['arr_1'], f['arr_2'], f['arr_3'].item()

    xLo, xHi, xpb = 0, 50, 1
    hFrac = np.divide(hSurf, hTot, dtype=float)

    fig1, (ax1, ax2) = plt.subplots(nrows=2, figsize=(8,6))

    ax1.plot(xTot, hTot, ls='steps', c='k', label="All Events")
    ax1.plot(xTot, hSurf, ls='steps', c='r', label="Transition Events")
    ax1.axvline(1, c='b', lw=1, label="1.0 keV")
    ax2.plot(xTot, hFrac, ls='steps', c='r')
    ax2.axvline(1, c='b', lw=1, label="1.0 keV")
    ax1.set_xlim(0, 50)
    ax2.set_xlim(0, 50)
    # ax1.set_title('Module 1 - Fraction of Transition Layer Events')
    ax1.set_ylabel('Counts / %.1f keV' % xpb, ha='right', y=1)
    # ax1.set_xlabel('Energy (keV)', ha='right', x=1)
    ax2.set_ylabel('Slow Fraction', ha='right', y=1)
    ax2.set_xlabel('Energy (keV)', ha='right', x=1)
    ax1.legend()
    ax2.legend()
    plt.tight_layout()
    # plt.show()
    plt.savefig('../plots/sim-sloFrac-allDets.pdf')

    for det in simHists:
        if det!='163': continue
        print(det)
        hTot, hSurf, xTot = simHists[det]
        hFrac = np.divide(hSurf, hTot, dtype=float)

        plt.close()
        fig1, (ax1, ax2) = plt.subplots(nrows=2, figsize=(8,6))
        ax1.plot(xTot, hTot, ls='steps', c='k', label="All Events, C{}P{}D{}".format(*str(det)))
        ax1.plot(xTot, hSurf, ls='steps', c='r', label="Transition Events")
        ax1.axvline(1, c='b', lw=1, label="1.0 keV")
        ax2.plot(xTot, hFrac, ls='steps', c='r')
        ax2.axvline(1, c='b', lw=1, label="1.0 keV")
        ax1.set_xlim(0, 50)
        ax2.set_xlim(0, 50)
        ax1.set_ylabel('Counts / %.1f keV' % xpb, ha='right', y=1)
        ax2.set_ylabel('Slow Fraction', ha='right', y=1)
        ax2.set_xlabel('Energy (keV)', ha='right', x=1)
        ax1.legend()
        ax2.legend()
        plt.tight_layout()
        # plt.show()
        plt.savefig("../plots/sim-sloFrac-%s.pdf" % det)

        # return


def plotComptonEdge():

    f = np.load('../data/efficiency-corr250.npz')
    hTot, hSurf, xTot, simHists = f['arr_0'], f['arr_1'], f['arr_2'], f['arr_3'].item()

    xLo, xHi, xpb = 0, 250, 1
    hFrac = np.divide(hSurf, hTot, dtype=float)

    # fig1, (ax1, ax2) = plt.subplots(nrows=2, figsize=(8,6))
    fig1 = plt.figure(figsize=(8,6))
    ax1 = plt.subplot2grid((3,1), (0,0), rowspan=2)
    ax2 = plt.subplot2grid((3,1), (2,0), sharex=ax1)

    ax1.plot(xTot, hTot, ls='steps', c='k', label="All Events")
    ax1.plot(xTot, hSurf, ls='steps', c='r', label="Transition Events")
    ax1.axvline(1, c='b', lw=1, label="1.0 keV")
    ax1.axvline(238, c='g', lw=1, label='238 keV')
    ax2.plot(xTot, hFrac, ls='steps', c='r')
    ax2.axvline(1, c='b', lw=1, label="1.0 keV")
    ax1.axvline(238, c='g', lw=1, label='238 keV')
    ax1.set_xlim(0, 250)
    ax1.set_ylim(ymax=np.amax(hTot)*1.4)
    ax2.set_xlim(0, 250)

    # ax1.set_title('Module 1 - Fraction of Transition Layer Events')
    ax1.set_ylabel('Counts / %.1f keV' % xpb, ha='right', y=1)
    # ax1.set_xlabel('Energy (keV)', ha='right', x=1)
    ax2.set_ylabel('Slow Fraction', ha='right', y=1)
    ax2.set_xlabel('Energy (keV)', ha='right', x=1)
    ax1.legend(loc=1, fontsize=12)
    ax2.legend(loc=1, fontsize=12)
    plt.tight_layout()
    plt.show()
    # plt.savefig('../plots/sim-sloFrac-allDets.pdf')


def plotDetSpectra():
    """
        Draws the detector-by-detector m2s238 spectrum from simulations
    """

    # Want to draw the spectrum of simulated events per detector and compare it with the calibration data
    df = pd.concat([pd.read_hdf('{}/data/DS5_M1_SimHitData_{}.h5'.format(os.environ['LATDIR'], i)) for i in range(2)])
    df = (df.dropna()
            .loc[(df['sumET'] > 237.28) & (df['sumET'] < 239.46)])

    dfg = df.groupby(['CPD1'])

    fig1, ax1 = plt.subplots(figsize=(10,7))
    bins = np.linspace(0, 120, 61)

    for name, g in dfg:
        # highE = g['trapENFCal2'].values
        lowE = g['trapENFCal1'].values

        # print(highE, highE.shape)
        # print(lowE, lowE.shape)

        ax1.cla()
        sns.distplot(lowE, bins=bins, kde=False, label='C{}P{}D{}'.format(*str(name)), ax=ax1)
        fig1.savefig('{}/plots/CalSim/C{}P{}D{}_Spectra.png'.format(os.environ['LATDIR'], *str(name)))


def plotM1Spectra():
    """
        Plot the mHL spectra data with and without cuts
    """
    sns.set(style='ticks')
    df = pd.read_hdf(os.environ['LATDIR']+'/data/CalPairHit_WithDbCut_mH1.h5')

    # fig1, ax1 = plt.subplots(figsize=(10,7))
    # aCuts = df.loc[(df['Pass1'] == True)&(df['cpd1'] > 200), 'trapENFCal1'].values
    # bCuts = df.loc[df['cpd1'] > 200,'trapENFCal1'].values
    # ax1.set_title('Module 2, mHT == 1 Spectrum ({:.2f}% passing cuts)'.format(float(len(aCuts)/len(bCuts))*100))
    # ax1.set_xlabel('Energy (keV)')
    # ax1.set_ylabel('Counts')
    # ax1.legend()
    # plt.tight_layout()

    dfg = df.groupby(['cpd1'])
    cpdList = df['cpd1'].unique()
    fig2, ax2 = plt.subplots(ncols=7, nrows=5, figsize=(17,15))
    ax2 = ax2.flatten()
    idx = 0
    for name, g in dfg:
        # Apply avse right now
        bCuts = g.loc[(df['avse1'] > -1.) ,'trapENFCal1'].values
        aCuts = g.loc[(df['Pass1'] == True) & (df['avse1'] > -1.), 'trapENFCal1'].values
        sns.distplot(aCuts, kde=False, color='b', label='After Cut', ax=ax2[idx])
        sns.distplot(bCuts, kde=False, color='r', label='Before Cut', ax=ax2[idx])
        sns.despine()
        isEnr = detInfo.isEnr(name)
        isEnrStr = ''
        if isEnr:
            isEnrStr = 'Enr'
        else:
            isEnrStr = 'Nat'

        ax2[idx].set(title='C{}P{}D{} ({enr}) ({eff:.1f}%)'.format(*str(name),
                    enr=isEnrStr, eff=100.*len(aCuts)/len(bCuts)))
        idx += 1

    # Turn off the last axes without figures
    for i in range(idx,len(ax2)):
        ax2[i].axis('off')

    plt.tight_layout()
    fig2.savefig('Detector_mH1Eff_avse.png')
    # plt.show()


def combinedEff(bSave = False, bWriteDB = False):
    """
        Utilizes the data produced by lat2 and combines all detectors to calculate one efficiency

        Currently this data only runs out to 50 keV, we can change it to go out farther in energy
    """
    from statsmodels.stats import proportion
    sns.set(style='ticks')

    f = np.load(os.environ['LATDIR']+'/data/lat2-eff-data-95.npz')
    effData = f['arr_0'].item()

    detList = effData.keys()
    xCenterVals = effData['112'][0]
    xVals = effData['112'][7][1:]

    # These are detectors skipped because of low statistics
    #for now I am going to keep all detectors regardless of statistics
    skipList = ['111','211','214','221','261','274']

    # Split between Natural and Enriched efficiencies
    hPassNat, hFullNat = np.zeros(len(xVals)), np.zeros(len(xVals))
    hPassEnr, hFullEnr = np.zeros(len(xVals)), np.zeros(len(xVals))

    for det in detList:
        # Skip the bad stuff
        if det in skipList:
            # print('Skipping', det)
            continue
        isEnr = detInfo.isEnr(det)
        # Pass is the 4th array, total is the 6th
        # print(det, len(effData[det][7]), len(effData[det][4]), len(effData[det][6]))
        # print(effData[det][4], effData[det][6])
        if isEnr:
            hPassEnr += effData[det][4][1:]
            hFullEnr += effData[det][6][1:]
        else:
            hPassNat += effData[det][4][1:]
            hFullNat += effData[det][6][1:]


    hEffEnr = np.nan_to_num(hPassEnr/hFullEnr)
    hEffNat = np.nan_to_num(hPassNat/hFullNat)

    print(len(xVals), xVals)
    print(hPassEnr)
    print(hFullEnr)
    print(hEffEnr)

    # Calculate CI of points based off of stats only
    nat_ci_low, nat_ci_upp = proportion.proportion_confint(hPassNat, hFullNat, alpha=0.1, method='beta')
    enr_ci_low, enr_ci_upp = proportion.proportion_confint(hPassEnr, hFullEnr, alpha=0.1, method='beta')

    fitBnd = ((1,-20,0,0.5),(np.inf,np.inf,np.inf, 0.99))
    initialGuess = [1, -1.6, 2.75, 0.95]
    poptnat, pcovnat = curve_fit(wl.weibull, xVals, hEffNat, p0=initialGuess, bounds=fitBnd)
    poptenr, pcovenr = curve_fit(wl.weibull, xVals, hEffEnr, p0=initialGuess, bounds=fitBnd)

    print(poptenr, pcovenr)
    print(poptnat, pcovnat)

    effNat = wl.weibull(xVals, *poptnat)
    effEnr = wl.weibull(xVals, *poptenr)

    # Generate and save a large dataset into a dataframe
    dfEBFList, dfNBFList = [], []
    dfEList, dfNList = [], []
    xValsFine = np.linspace(0, 200, 2001)

    effBFNat = wl.weibull(xValsFine, *poptnat)
    effBFEnr = wl.weibull(xValsFine, *poptenr)

    dfEBFList.append({xValsFine[i]:effBFEnr[i] for i in range(len(xValsFine))})
    dfNBFList.append({xValsFine[i]:effBFNat[i] for i in range(len(xValsFine))})

    zVal = 1.645
    # Simple Uncertainties
    Sigma_Enr = np.sqrt(np.diagonal(pcovenr))
    Sigma_Nat = np.sqrt(np.diagonal(pcovnat))

    # These are some preliminary uncertainty calculations, the first one is based off of using the diagonal as the uncertainty (ignoring all corrlations)
    # The second one is using a calculation from https://stats.stackexchange.com/questions/135749/confidence-intervals-of-fitted-weibull-survival-function
    # We decided to go with the Toy MC method as it's the most accurate (and doesn't require me to do more math)

    effEnrUpper = wl.weibull(xVals, *(np.array(poptenr) + zVal*Sigma_Enr))
    effEnrLower = wl.weibull(xVals, *(np.array(poptenr) - zVal*Sigma_Enr))
    effNatUpper = wl.weibull(xVals, *(np.array(poptnat) + zVal*Sigma_Nat))
    effNatLower = wl.weibull(xVals, *(np.array(poptnat) - zVal*Sigma_Nat))
    eff2EnrUpper = np.exp( np.log(effEnr)*np.exp(zVal/np.sqrt(np.array(hFullEnr))) )
    eff2EnrLower = np.exp( np.log(effEnr)*np.exp(-1.*zVal/np.sqrt(np.array(hFullEnr))) )
    eff2NatUpper = np.exp( np.log(effNat)*np.exp(zVal/np.sqrt(np.array(hFullNat))) )
    eff2NatLower = np.exp( np.log(effNat)*np.exp(-1.*zVal/np.sqrt(np.array(hFullNat))) )

    # Save stuff into calDB
    dbKey = "fitSlo_Combined_m2s238_eff95"
    cEnr, locEnr, scEnr, ampEnr = poptenr
    cEEnr, locEEnr, scEEnr, ampEEnr = Sigma_Enr
    cNat, locNat, scNat, ampNat = poptnat
    cENat, locENat, scENat, ampENat = Sigma_Nat

    dbVals = {}
    # 0 = Enr, 1 = Nat (similar to gain)
    dbVals[0] = [cEnr, locEnr, scEnr, ampEnr, cEEnr, locEEnr, scEEnr, ampEEnr]
    dbVals[1] = [cNat, locNat, scNat, ampNat, cENat, locENat, scENat, ampENat]

    if bWriteDB:
        print('Writing dbVals:', dbVals)
        dsi.setDBRecord({"key":dbKey, "vals":dbVals}, forceUpdate=True, calDB=calDB, pars=pars)
        print("DB filled:",dbKey)

        fsN = dsi.getDBRecord(dbKey, False, calDB, pars)
        print(fsN)
    return

    fig1, (ax11, ax12) = plt.subplots(ncols=2, figsize=(15,7))
    ax11.set_title('Total Enriched Efficiency')
    ax11.set_xlabel('Energy (keV)')
    ax11.set_ylabel('Efficiency')
    ax11.errorbar(xVals, hEffEnr, yerr=[hEffEnr - enr_ci_low, enr_ci_upp - hEffEnr], color='k', linewidth=0.8, fmt='o', capsize=2)

    ax12.set_title('Total Natural Efficiency')
    ax12.set_xlabel('Energy (keV)')
    ax12.set_ylabel('Efficiency')
    ax12.errorbar(xVals, hEffNat, yerr=[hEffNat - nat_ci_low, nat_ci_upp - hEffNat], color='k', linewidth=0.8, fmt='o', capsize=2)
    sns.despine(ax=ax11)
    sns.despine(ax=ax12)

    # Now Draw all the individual detector efficiencies on the plots in a lighter tone
    # Grab fitSlo efficiency values from DB
    fsD = dsi.getDBRecord("fitSlo_cpd_eff95", False, calDB, pars)

    for det in detList:
        if det in skipList:
            continue
        isEnr = detInfo.isEnr(det)
        cpd = int(det)
        if cpd in fsD:
            c, loc, scale, amp = fsD[cpd][3], fsD[cpd][4], fsD[cpd][5], fsD[cpd][2]
            effCh = wl.weibull(xVals, c, loc, scale, amp)
        else:
            print('Detector {} is not in the DB'.format(det))
            continue
        if isEnr:
            ax11.plot(xVals, effCh, alpha=0.5, label='C{}P{}D{}'.format(*det))
        else:
            ax12.plot(xVals, effCh, alpha=0.5, label='C{}P{}D{}'.format(*det))

    ax11.plot(xVals, effEnr, lw=2, color='r', label='Total Efficiency')
    # ax11.fill_between(xVals, eff2EnrLower, eff2EnrUpper, color='r', alpha=0.2)
    ax12.plot(xVals, effNat, lw=2, color='r', label='Total Efficiency')
    # ax12.fill_between(xVals, eff2NatLower, eff2NatUpper, color='r', alpha=0.2)
    # ax11.set_ylim(0.75, 1.0)
    # ax12.set_ylim(0.75, 1.0)
    ax11.legend(loc='lower right', ncol=3)
    ax12.legend(loc='lower right', ncol=3)
    plt.tight_layout()
    # fig1.savefig('CombinedEfficiency_skipbad_floatX.png')

    fig21, (ax21, ax22) = plt.subplots(ncols=2, figsize=(15,7))
    ax21.set_title('Total Enriched Efficiency')
    ax21.set_xlabel('Energy (keV)')
    ax21.set_ylabel('Efficiency')
    ax21.errorbar(xVals, hEffEnr, yerr=[hEffEnr - enr_ci_low, enr_ci_upp - hEffEnr], color='k', linewidth=0.8, fmt='o', capsize=2)

    ax22.set_title('Total Natural Efficiency')
    ax22.set_xlabel('Energy (keV)')
    ax22.set_ylabel('Efficiency')
    ax22.errorbar(xVals, hEffNat, yerr=[hEffNat - nat_ci_low, nat_ci_upp - hEffNat], color='k', linewidth=0.8, fmt='o', capsize=2)


    # Now generate 50000 toyMC and re-fit the toyMC spectra
    nMC = 50000
    cNatList, locNatList, scaleNatList, ampNatList = [], [], [], []
    cEnrList, locEnrList, scaleEnrList, ampEnrList = [], [], [], []

    for idx in range(nMC):
        if idx%1000 == 0:
            print('Current Loop: {} of {}'.format(idx, nMC))
        ePass = np.random.poisson(hPassEnr)
        eFull = np.random.poisson(hFullEnr)
        nPass = np.random.poisson(hPassNat)
        nFull = np.random.poisson(hFullNat)

        eEff = np.nan_to_num(ePass/eFull)
        nEff = np.nan_to_num(nPass/nFull)

        popte, _ = curve_fit(wl.weibull, xVals, eEff, p0=initialGuess, bounds=fitBnd)
        poptn, _ = curve_fit(wl.weibull, xVals, nEff, p0=initialGuess, bounds=fitBnd)

        effE = wl.weibull(xVals, *popte)
        effN = wl.weibull(xVals, *poptn)

        effEFine = wl.weibull(xValsFine, *popte)
        effNFine = wl.weibull(xValsFine, *poptn)

        dETemp = {xValsFine[i]:effEFine[i] for i in range(len(xValsFine))}
        dNTemp = {xValsFine[i]:effNFine[i] for i in range(len(xValsFine))}

        dfEList.append(dETemp)
        dfNList.append(dNTemp)
        cNatList.append(poptn[0])
        locNatList.append(poptn[1])
        scaleNatList.append(poptn[2])
        ampNatList.append(poptn[3])

        cEnrList.append(popte[0])
        locEnrList.append(popte[1])
        scaleEnrList.append(popte[2])
        ampEnrList.append(popte[3])

        ax21.plot(xVals, effE, color='b', alpha=0.005)
        ax22.plot(xVals, effN, color='b', alpha=0.005)

    if bSave:
        dfBFE = pd.DataFrame.from_dict(dfEBFList)
        dfBFN = pd.DataFrame.from_dict(dfNBFList)

        dfE = pd.DataFrame.from_dict(dfEList)
        dfN = pd.DataFrame.from_dict(dfNList)

        print(dfE.head(10))
        dfBFE.to_hdf(os.environ['LATDIR']+'/data/ToyMCEff.h5', key='EnrBF', format='table', mode='w', complevel=9)
        dfBFN.to_hdf(os.environ['LATDIR']+'/data/ToyMCEff.h5', key='NatBF', format='table', mode='a', complevel=9)
        dfE.to_hdf(os.environ['LATDIR']+'/data/ToyMCEff.h5', key='Enr', format='table', mode='a', complevel=9)
        dfN.to_hdf(os.environ['LATDIR']+'/data/ToyMCEff.h5', key='Nat', format='table', mode='a', complevel=9)

    # Draw these last so they show up!
    ax21.plot(xVals, effEnr, lw=2, color='r', label='Total Efficiency')
    sns.despine()
    ax22.plot(xVals, effNat, lw=2, color='r', label='Total Efficiency')
    sns.despine()
    ax21.set_ylim(0.75, 1.0)
    ax22.set_ylim(0.75, 1.0)
    plt.tight_layout()

    fig22, ax2 = plt.subplots(ncols=4, nrows=2, figsize=(20,17))
    ax2 = ax2.flatten()
    sns.distplot(cEnrList, kde=False, ax=ax2[0])
    ax2[0].set_title('Enriched c')
    sns.distplot(locEnrList, kde=False, ax=ax2[1])
    ax2[1].set_title('Enriched loc')
    sns.distplot(scaleEnrList, kde=False, ax=ax2[2])
    ax2[2].set_title('Enriched scale')
    sns.distplot(ampEnrList, kde=False, ax=ax2[3])
    ax2[3].set_title('Enriched amp')

    sns.distplot(cNatList, kde=False, ax=ax2[4])
    ax2[4].set_title('Natural c')
    sns.distplot(locNatList, kde=False, ax=ax2[5])
    ax2[5].set_title('Natural loc')
    sns.distplot(scaleNatList, kde=False, ax=ax2[6])
    ax2[6].set_title('Natural scale')
    sns.distplot(ampNatList, kde=False, ax=ax2[7])
    ax2[7].set_title('Natural amp')

    plt.tight_layout()
    # fig21.savefig('TotalEfficiency_ToyMC_2_floatX.png')
    # fig22.savefig('TotalEfficiency_ToyMC_Pars_floatX.png')
    # plt.show()

def plotMCEffUnc(bSave = False):
    """
        Loads the toy efficiency data produced from before, does some sanity checks, and calculates the std, etc
    """
    sns.set(style='ticks')

    dfEnr = pd.read_hdf(os.environ['LATDIR']+'/data/ToyMCEff.h5', 'Enr')
    dfNat = pd.read_hdf(os.environ['LATDIR']+'/data/ToyMCEff.h5', 'Nat')

    # Best Fit curve
    dfEnrBF = pd.read_hdf(os.environ['LATDIR']+'/data/ToyMCEff.h5', 'EnrBF')
    dfNatBF = pd.read_hdf(os.environ['LATDIR']+'/data/ToyMCEff.h5', 'NatBF')

    xVals = np.linspace(0, 200, 2001)
    # Here I take the means and compare to the best fit, for the best fit DF, there's only 1 value so it doesn't matter if I take the mean
    EnrEff = dfEnrBF.mean(axis=0).values
    NatEff = dfNatBF.mean(axis=0).values
    toyEnrEff = dfEnr.mean(axis=0).values
    toyNatEff = dfNat.mean(axis=0).values

    toyEnrStd = dfEnr.std(axis=0).values
    toyNatStd = dfNat.std(axis=0).values
    zVal = 1.645
    EnrEffHi = EnrEff + zVal*toyEnrStd
    EnrEffLo = EnrEff - zVal*toyEnrStd
    NatEffHi = NatEff + zVal*toyNatStd
    NatEffLo = NatEff - zVal*toyNatStd

    # Loop through the first couple of columns using iloc and plot the distributions, see if they look gaussian
    #Spoilers: they look Gaussian so we can use the std calculated from Pandas
    fig0, ax0 = plt.subplots(ncols=8, nrows=5, figsize=(15,12))
    ax0 = ax0.flatten()
    bins = np.linspace(0.25, 1.0, 600)
    for i in range(len(ax0)):
        # Start with bin 15 (1 keV)
        idx = i + 11
        # sns.distplot(dfEnr.iloc[idx,:].values, kde=False, bins=bins, ax=ax0[i])
        print('Best Fit mean: {}'.format(EnrEff[idx]))
        # binRange = EnrEff[idx]
        ax0[i].hist(dfEnr.iloc[:,idx].values, bins=bins)
        ax0[i].set_title('{:.1f} keV'.format(0.1*(idx-1)))
        ax0[i].set_xlim(EnrEff[idx] - 5*toyEnrStd[idx], EnrEff[idx] + 5*toyEnrStd[idx])
        ax0[i].axvline(x=EnrEff[idx], color='r')
        ax0[i].axvline(x=EnrEff[idx]-zVal*toyEnrStd[idx], color='r', linestyle='--', alpha=0.7)
        ax0[i].axvline(x=EnrEff[idx]+zVal*toyEnrStd[idx], color='r', linestyle='--', alpha=0.7)
        sns.despine()

    plt.tight_layout()
    fig0.savefig(os.environ['LATDIR']+'/TotalEfficiency_Unc_Dist.png')
    return

    # print(EnrEff)
    # print(toyEnrStd)

    # Plot them side by side to make sure it makes sense
    # they differ by less than 0.1% at each point
    # uncertainties look good
    fig1, (ax11, ax12) = plt.subplots(ncols=2, figsize=(15,7))
    ax11.plot(xVals, EnrEff, 'r', label='Best Fit Efficiency')
    # ax11.fill_between(xVals, EnrEffLo, EnrEffHi, color='r', alpha=0.3, label=r'$1 \sigma$ Uncertainties')
    ax11.fill_between(xVals, EnrEffLo, EnrEffHi, color='r', alpha=0.3, label='90% C.I.')
    ax11.plot(xVals, toyEnrEff, 'b--', label='Toy MC Mean Efficiency')
    sns.despine()
    ax12.plot(xVals, NatEff, 'r', label='Best Fit Efficiency')
    ax12.plot(xVals, toyNatEff, 'b--', label='Toy MC Mean Efficiency')
    # ax12.fill_between(xVals, NatEffLo, NatEffHi, color='r', alpha=0.3, label=r'$1 \sigma$ Uncertainties')
    ax12.fill_between(xVals, NatEffLo, NatEffHi, color='r', alpha=0.3, label='90% C.I.')
    sns.despine()
    ax11.set(title='Total Enriched Efficiency', ylabel='Efficiency', xlabel='Energy (keV)')
    ax12.set(title='Total Natural Efficiency', ylabel='Efficiency', xlabel='Energy (keV)')
    ax11.set_ylim(0.75, 1.0)
    ax12.set_ylim(0.75, 1.0)
    ax11.set_xlim(0,50.)
    ax12.set_xlim(0,50.)
    ax11.legend(loc = 'lower right')
    ax12.legend(loc = 'lower right')
    plt.tight_layout()
    # fig1.savefig(os.environ['LATDIR']+'/TotalEfficiency_Unc.png')
    # plt.show()

    # import ROOT
    # hEnr = ROOT.TH1D()


if __name__=="__main__":
    main()
