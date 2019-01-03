#!/usr/bin/env python3
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
plt.style.use('../clint.mpl')
import dsi
det = dsi.DetInfo()

def main():
    """ lat-expo output in pandas:
    ['ds', 'run', 'det', 'start', 'stop', 'lt', 'rt', 'expo'] """

    # checkExpo()
    # getExpected()
    getCurves()


def checkExpo():
    """ Check exposures of the input file match the LAT unidoc.
    Remember, C2P5D3 was removed because it sucks.
    """
    dsListList = [ # 50 == 5A, etc.
        [0],[1],[2],[3],[4],[50],[51],[52],[6],
        [0,1,2,3,4,50,51,52,6],
        [1,2,3,4,50,51,52],
        [1,2,3,4,50,51,52,6],
        ]

    df = pd.read_hdf("../data/lat-runTimes.h5")

    detList = sorted([int(cpd) for cpd in det.allDetIDs])
    enrList = [cpd for cpd in detList if det.isEnr(cpd)]
    natList = [cpd for cpd in detList if not det.isEnr(cpd)]

    for dsList in dsListList:

        enrExp = df.loc[ df['det'].isin(enrList) & df['ds'].isin(dsList) ]["expo"].sum()
        natExp = df.loc[ df['det'].isin(natList) & df['ds'].isin(dsList) ]["expo"].sum()

        print('dsList:',dsList)
        print("  enrExp: %.2f kg-d, %.2f kg-y,  natExp: %.2f kg-d, %.2f kg-y" % (enrExp,enrExp/365.25,natExp,natExp/365.25))


def getFitResults():
    """ Results from the cosmogenic spectrum fit """

    cols = ["trit","trit-e","210Pb","210Pb-e","68Ge","68Ge-e","68Ga","68Ga-e",
            "65Zn","65Zn-e","55Fe","55Fe-e","54Mn","54Mn-e","flat","flat-e"]
    dfFit = pd.DataFrame(columns=cols)
    dfFit.loc["M1"] = [1817.44, 97.2, 57.89, 10.1, 6.95, 10.9, 11.01, 12.0, 38.13, 13.5, 47.31, 16.3, 0.00, 19.9, 1709.10, 49.0]
    dfFit.loc["M2"] = [302.47, 40.7, 27.32, 6.3, 59.00, 9.3, 0.0, 3.3, 3.28, 5.5, 0.0, 7.0, 12.34, 7.6, 417.77, 24.0]

    return dfFit


def getAbundance(t0=None):
    """ Results from Brandon's activation calculations """

    cols = ["detID","cpd","arrival","abundance-68Ge","surfExpo"]
    dfAbu = pd.DataFrame(columns=cols)
    dfAbu.loc[0] = ["P42537A", 174, "2013/2/12", 58.5, 20.9]
    dfAbu.loc[1] = ["P42538A", 114, "2013/2/12", 56.3, 20.1]
    dfAbu.loc[2] = ["P42538B", 133, "2013/2/12", 67.9, 24.3]
    dfAbu.loc[3] = ["P42573A", 134, "2013/2/12", 67.4, 24.1]
    dfAbu.loc[4] = ["P42573B", 154, "2013/2/12", 67.4, 24.1]
    dfAbu.loc[5] = ["P42574A", 163, "2013/3/26", 60.9, 21.8]
    dfAbu.loc[6] = ["P42574B", 172, "2013/3/26", 56.3, 20.1]
    dfAbu.loc[7] = ["P42574C", 161, "2013/3/26", 59.0, 21.1]
    dfAbu.loc[8] = ["P42575A", 112, "2013/3/26", 60.9, 21.8]
    dfAbu.loc[9] = ["P42575B", 152, "2013/8/13", 82.9, 29.6]
    dfAbu.loc[10] = ["P42661A", 153, "2013/9/8 ", 61.7, 22.0]
    dfAbu.loc[11] = ["P42661B", 162, "2013/9/8 ", 61.7, 22.0]
    dfAbu.loc[12] = ["P42661C", 113, "2013/9/8 ", 61.7, 22.0]
    dfAbu.loc[13] = ["P42662A", 164, "2013/9/8 ", 62.1, 22.2]
    dfAbu.loc[14] = ["P42662B", 173, "2013/9/8 ", 62.1, 22.2]
    dfAbu.loc[15] = ["P42662C", 124, "2013/9/8 ", 62.1, 22.2]
    dfAbu.loc[16] = ["P42664A", 122, "2013/9/8 ", 64.2, 22.9]
    dfAbu.loc[17] = ["P42665A", 123, "2013/12/1", 59.3, 21.2]
    dfAbu.loc[18] = ["P42698A", 132, "2013/12/1", 63.4, 22.6]
    dfAbu.loc[19] = ["P42698B", 111, "2013/12/1", 63.4, 22.6]
    dfAbu.loc[20] = ["P23517A", 261, "2013/9/8 ", 72.2, 25.8]
    dfAbu.loc[21] = ["P42665B", 252, "2013/12/1", 59.3, 21.2]
    dfAbu.loc[22] = ["P42665C", 264, "2013/12/1", 59.3, 21.2]
    dfAbu.loc[23] = ["P42664B", 212, "2014/1/29", 95.9, 34.3]
    dfAbu.loc[24] = ["P42712A", 254, "2014/1/29", 83.2, 29.7]
    dfAbu.loc[25] = ["P42712B", 272, "2014/1/29", 64.3, 23.0]
    dfAbu.loc[26] = ["P42748A", 214, "2015/4/25", 29.9, 10.7]
    dfAbu.loc[27] = ["P42748B", 213, "2015/4/25", 29.9, 10.7]
    dfAbu.loc[28] = ["P42749A", 231, "2015/4/25", 30.8, 11.0]
    dfAbu.loc[29] = ["P42749B", 232, "2015/4/25", 30.8, 11.0]
    dfAbu.loc[30] = ["P42853A", 233, "2015/6/22", 46.5, 16.6]
    dfAbu.loc[31] = ["P42853B", 253, "2015/6/22", 46.5, 16.6]
    dfAbu.loc[32] = ["P42909A", 273, "2015/6/22", 82.5, 29.5]
    dfAbu.loc[33] = ["P42909B", 262, "2015/6/22", 82.5, 29.5]
    dfAbu.loc[34] = ["P42909C", 263, "2015/6/22", 82.5, 29.5]

    dfAbu["t0"] = pd.to_datetime(dfAbu["arrival"])       # make a column w/ an explicit datetime object
    dfAbu["t0u"] = dfAbu["t0"].astype(np.int64) // 10**9 # make the unix timestamp explicit

    if t0 is not None:
        dfAbu['dt'] = t0 - dfAbu["t0u"]

    return dfAbu


def getExpected():
    """ Get the expected 68Ge number of atoms at the beginning of a dsList """

    dsList = [1,2,3,4,50,51,52,6]
    df = pd.read_hdf("../data/lat-runTimes.h5")
    df = df.loc[ df['ds'].isin(dsList) ]
    t0 = df.loc[ df['ds'] == dsList[0] ]['start'].iloc[0]

    t12 = {"68Ge": 0.7414}

    detIDs = det.allDetIDs
    detAMs = det.allActiveMasses
    detList = sorted([int(cpd) for cpd in det.allDetIDs])
    enrList = [cpd for cpd in detList if det.isEnr(cpd)]
    enrDets = set(df.loc[ df['det'].isin(enrList) & df['ds'].isin(dsList) ]['det'])

    dfAbu = getAbundance(t0)

    NTot = 0
    for cpd in sorted(list(enrDets)):

        aMass = detAMs[detIDs[str(cpd)]]/1000
        lmb = np.log(2) / (t12["68Ge"] * 31557600)
        n0 = dfAbu.loc[ dfAbu['cpd']==cpd, "abundance-68Ge" ].values[0]
        dt = dfAbu.loc[ dfAbu['cpd']==cpd, "dt" ].values[0]

        N = aMass * n0 * np.exp(-lmb * dt)
        NTot += N

        print("cpd %d  aMass %.2f  N(t=ds1) %.2f" % (cpd, aMass, N))

    print("NTot:", NTot)


def Fij(tDF, t12, t0):
    """ unix start time, stop time, live time = tDF """
    lmb = np.log(2) / (t12 * 31557600) # decay constant [1/sec]
    eff = tDF[2] / (tDF[1] - tDF[0])   # deadtime factor
    return eff * (np.exp(-lmb * (tDF[0] - t0)) - np.exp(-lmb * (tDF[1] - t0)))


def getCurves():

    dsList = [1,2,3,4,50,51,52,6] # 50=5A, etc.

    df = pd.read_hdf("../data/lat-runTimes.h5")
    df = df.loc[ df['ds'].isin(dsList) ]
    t0 = df.loc[ df['ds'] == dsList[0] ]['start'].iloc[0]

    detIDs = det.allDetIDs
    detAMs = det.allActiveMasses
    detList = sorted([int(cpd) for cpd in det.allDetIDs])

    enrList = [cpd for cpd in detList if det.isEnr(cpd)]
    enrDets = set(df.loc[ df['det'].isin(enrList) & df['ds'].isin(dsList) ]['det'])
    enrMass = sum([ detAMs[detIDs[str(cpd)]] for cpd in enrDets])
    enrExp = df.loc[ df['det'].isin(enrList) & df['ds'].isin(dsList) ]["expo"].sum()

    t12 = {"68Ge": 0.7414}

    dfFit, dfAbu = getFitResults(), getAbundance(t0)

    # print(dfAbu)
    # exit()

    for iso in t12:

        fSum = 0
        expTot = 0
        for cpd in enrList:

            dfCPD = df.loc[ df['det'] == cpd, ['start','stop','lt']]
            if dfCPD.shape[0] == 0: continue
            aMass = detAMs[detIDs[str(cpd)]]

            # weighting factor 1: det active mass / total active mass
            wt1 = detAMs[detIDs[str(cpd)]] / enrMass

            # weighting factor 2: det exposure / total exposure
            expo = dfCPD['lt'].sum()/86400 * aMass/1000
            wt = expo / enrExp
            expTot += expo

            # calculate Fij terms
            eFac = dfCPD.apply(Fij, axis=1, args=(t12[iso], t0))
            fTot = eFac.sum() * wt
            fSum += fTot

            # sanity check stuff
            ts = dfCPD['start'].iloc[0]
            tf = dfCPD['stop'].iloc[-1]
            f0 = Fij([ts, tf, (tf-ts)], t12[iso], t0) * wt
            lt = dfCPD['lt'].sum()
            rt = (dfCPD['stop'] - dfCPD['start']).sum()
            rtCal = tf - ts

            print("%d  fTot %.3f  f0 %.3f  lt/rt %.3f  rt/rtCal %.2f%%  wt1 %.3f  wt2 %.3f  expo %.3f" % (cpd, fTot, f0, lt/rt, 100 * rt/rtCal, wt1, wt, expo))

        N_c = dfFit[iso].sum() # this is all detectors, fit result
        N_0 = N_c / fSum

        print("expo (loop) %.2f, expo (df) %.2f" % (expTot, enrExp))
        print(iso, "N_c: %.1f, N_0: %.2f" % (N_c, N_0))


if __name__=="__main__":
    main()