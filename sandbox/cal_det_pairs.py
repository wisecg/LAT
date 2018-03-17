# -*- coding: utf-8 -*-
# ¯\_(ツ)_/¯ ¯\_(ツ)_/¯ ¯\_(ツ)_/¯ ¯\_(ツ)_/¯ ¯\_(ツ)_/¯ ¯\_(ツ)_/¯ ¯\_(ツ)_/¯ ¯\_(ツ)_/¯ ¯\_(ツ)_/¯
#!/usr/bin/env python
import sys, os, imp, pathlib, itertools
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(style='darkgrid', context='talk')

# load LAT libraries
ds = imp.load_source('DataSetInfo',os.environ['LATDIR']+'/DataSetInfo.py')
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
# fsD = ds.getDBRecord("fitSlo_ds%d_idx%d_m%d_Peak" % (dsNum, calIdx, modNum), False, calDB, pars)

# load riseNoise vals
# rnSD = ds.getDBRecord("riseNoise_ds%d_idx%d_m%d_SoftPlus" % (dsNum, calIdx, modNum), False, calDB, pars)
# rnCD = ds.getDBRecord("riseNoise_ds%d_idx%d_m%d_Continuum" % (dsNum, calIdx, modNum), False, calDB, pars)


def main():
    getSpecPandas()
    # loadSpec()

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

def getSpecPandas():
    """ Get sum and hit spectra w/ threshold cut, without ch. 598 (it's noisy.)
        Here we use "mHT" and "sumET" exclusively.
        Use !EventDC1Bits and trapENFCal > 0.7, all hits.  (Skim file was generated w/ 'dontSkipAnything')
        Save detailed info for events with hits < 100 keV.
    """
    from ROOT import TFile, TTree
    inDir, outDir = '/mnt/mjdDisk1/Majorana/data/sandbox/special/lat/cal','/mnt/mjdDisk1/Majorana/users/psz/LAT/data'
    inPath = pathlib.Path(inDir)
    calInfo = ds.CalInfo()
    runList = calInfo.GetSpecialRuns("longCal",dsNum) # 54 total
    # runList = runList[:10]
    fileList = []
    for run in runList:
        files = inPath.glob('latSkimDS{}_run{}_*.root'.format(dsNum,run))
        fileList.extend([str(f) for f in files]) # Because ROOT doesn't play nice with python paths
    nFiles = len(fileList)
    chList = ds.GetGoodChanListNew(dsNum = dsNum)

    eCut = 120
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

            # Probably the "right" way to store it
            # if sumET < 2630:
            #     dataMap = {}
            #     dataMap['run'] = int(lTree.run)
            #     dataMap['channel'] = [int(lTree.channel.at(i)) for i in idxList]
            #     dataMap['mHT'] = mHT
            #     dataMap['sumET'] = sumET
            #     dataMap['iEvent'] = int(lTree.iEvent)
            #     dataMap['CPD'] = [int(100*lTree.C.at(i) + 10*lTree.P.at(i) + lTree.D.at(i)) for i in idxList ]
            #     dataMap['trapENFCal'] = hitE
            #     dataMap['fitSlo'] = [lTree.fitSlo.at(i) for i in idxList]
            #     dataMap['riseNoise'] = [lTree.riseNoise.at(i) for i in idxList]
            #     dataMap['UnixTime'] = lTree.globalTime.GetSec()
            #     dataList.append(dataMap)

            # Store everything as its own column
            # This is because I'm lazy AF and I didn't want to work with arrays in my dataframe
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

    # Write Dataframe to file -- suppress performance warnings
    import warnings
    warnings.filterwarnings(action="ignore", module="pandas", message="^\nyour performance")
    chunksize = 500000
    start, end, iChunk = 0, chunksize-1, 0
    while end < df.shape[0]:
        chunk = df.iloc[start:end]
        try:
            chunk.to_hdf('{}/DS{}_LongCal_HitData_{}.h5'.format(outDir, dsNum, iChunk), key='skimTree', data_columns=['run', 'channel', 'mHT', 'iEvent', 'CPD', 'trapENFCal',
            'fitSlo', 'riseNoise', 'UnixTime'], format='fixed', mode='w', complevel=9)
        except (Exception) as e:
            print e
            print chunk
            print chunk.info()
        start += chunksize
        end += chunksize
        iChunk += 1


def loadSpec():
    inDir = '/Users/brianzhu/macros/code/LAT'
    df = pd.read_hdf('{}/DS5_LongCal_HitData_0.h5'.format(inDir))
    print(df.head(10))
    


if __name__=="__main__":
    main()
