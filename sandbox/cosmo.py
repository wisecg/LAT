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

    checkExpo()
    getWeights()
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
    print(totMass/1000)




def Fij(tDF, t12, t0):
    """ unix start time, stop time, live time = tDF """
    lmb = np.log(2) / (t12 * 31557600) # decay constant [1/sec]
    eff = tDF[2] / (tDF[1] - tDF[0])   # deadtime factor
    return eff * (np.exp(-lmb * (tDF[0] - t0)) - np.exp(-lmb * (tDF[1] - t0)))


def getCurves():

    dsList = [1,2,3,4,50,51,52,6] # 50=5A, etc.

    detIDs = det.allDetIDs
    detAMs = det.allActiveMasses
    detList = sorted([int(cpd) for cpd in det.allDetIDs])

    t12 = {"68Ge": 0.7414}

    df = pd.read_hdf("../data/lat-runTimes.h5")
    df = df.loc[ df['ds'].isin(dsList) ]
    t0 = df.loc[ df['ds'] == dsList[0] ]['start'].iloc[0]

    for iso in t12:

        expTot = 0
        for cpd in detList:

            dfCPD = df.loc[ df['det'] == cpd, ['start','stop','lt']]
            if dfCPD.shape[0] == 0: continue

            # calculate Fij terms
            eFac = dfCPD.apply(Fij, axis=1, args=(t12[iso], t0))
            fTot = eFac.sum()

            # sanity check stuff
            ts = dfCPD['start'].iloc[0]
            tf = dfCPD['stop'].iloc[-1]
            f0 = Fij([ts, tf, (tf-ts)], t12[iso], t0)
            lt = dfCPD['lt'].sum()
            rt = (dfCPD['stop'] - dfCPD['start']).sum()
            rtCal = tf - ts

            # check exposure
            aMass = detAMs[detIDs[str(cpd)]]
            xp = lt/86400 * aMass/1000
            if det.isEnr(cpd):
                expTot += xp

            print("%d  fTot %.3f  f0 %.3f  lt/rt %.3f, rt/rtCal %.2f %%" % (cpd, fTot, f0, lt/rt, 100 * rt/rtCal))

        print("enr exp tot:",expTot/365.25)

        # we need some kinda normalization term
        # s/t Sum_i^dets m N_i^0 = 1 or something
        # if i had 13 detectors, each 1 kg, then i would have a factor 13*N_0 in the current version of the sum
        # that doesn't seem like the right way to weight it


if __name__=="__main__":
    main()