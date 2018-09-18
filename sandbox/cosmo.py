#!/usr/bin/env python3
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
plt.style.use('../pltReports.mplstyle')
import dsi
det = dsi.DetInfo()

def main():
    """ lat-expo output in pandas:
    ['ds', 'run', 'det', 'start', 'stop', 'lt', 'rt', 'expo'] """

    # checkExpo()
    # getWeights()
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


def getWeights():

    dsList = [1,2,3,4]
    df = pd.read_hdf("../data/lat-runTimes.h5")
    df = df.loc[ df['ds'].isin(dsList) ]

    detIDs = det.allDetIDs
    detAMs = det.allActiveMasses
    detList = sorted([int(cpd) for cpd in det.allDetIDs])
    enrList = [cpd for cpd in detList if det.isEnr(cpd)]
    natList = [cpd for cpd in detList if not det.isEnr(cpd)]
    enrDets = set(df.loc[ df['det'].isin(enrList) & df['ds'].isin(dsList) ]['det'])
    totMass = sum([ detAMs[detIDs[str(cpd)]] for cpd in enrDets])
    detWts = {cpd: detAMs[detIDs[str(cpd)]]/totMass for cpd in enrDets}

    print(totMass/1000)
    print(detWts)
    print(sum([detWts[cpd] for cpd in detWts]))


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

    cols = ["trit","trit-e","210Pb","210Pb-e","68Ge","68Ge-e","68Ga","68Ga-e",
            "65Zn","65Zn-e","55Fe","55Fe-e","54Mn","54Mn-e","flat","flat-e"]
    dfFit = pd.DataFrame(columns=cols)
    dfFit.loc["M1"] = [1817.44, 97.2, 57.89, 10.1, 6.95, 10.9, 11.01, 12.0, 38.13, 13.5, 47.31, 16.3, 0.00, 19.9, 1709.10, 49.0]
    dfFit.loc["M2"] = [302.47, 40.7, 27.32, 6.3, 59.00, 9.3, 0.0, 3.3, 3.28, 5.5, 0.0, 7.0, 12.34, 7.6, 417.77, 24.0]

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
            f0 = Fij([ts, tf, (tf-ts)], t12[iso], t0)
            lt = dfCPD['lt'].sum()
            rt = (dfCPD['stop'] - dfCPD['start']).sum()
            rtCal = tf - ts

            print("%d  fTot %.3f  f0 %.3f  lt/rt %.3f  rt/rtCal %.2f%%  wt1 %.3f  wt2 %.3f  expo %.3f" % (cpd, fTot, f0, lt/rt, 100 * rt/rtCal, wt1, wt, expo))

        N_c = dfFit[iso].sum() # this is all detectors, fit result
        N_0 = N_c / fSum

        print("expo (loop) %.2f, expo (df) %.2f" % (expTot, enrExp))
        print(iso, "N_c: %.1f, N_0: %.2f" % (N_c, N_0))

        # do we need some kinda normalization term?
        # s/t Sum_i^dets m N_i^0 = 1 or something
        # if i had 13 detectors, each 1 kg, then i would have a factor 13*N_0 in the current version of the sum
        # that doesn't seem like the right way to weight it


if __name__=="__main__":
    main()