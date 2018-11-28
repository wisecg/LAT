#!/usr/bin/env python3
import os
import numpy as np
import tinydb as db
import pandas as pd
from statsmodels.stats import proportion
from scipy.optimize import curve_fit
from pprint import pprint
import matplotlib.pyplot as plt
plt.style.use('pltReports.mplstyle')

import dsi
import waveLibs as wl
calDB = db.TinyDB('{}/calDB-v2.json'.format(dsi.latSWDir)) # match LAT's v2 tag
pars = db.Query()
detInfo = dsi.DetInfo()

# Skip these detectors because of low statistics
skipList = ['111', '211', '214', '221', '261', '274']

# Enriched and Natural detector lists, I pulled this from the output of the combined efficiencies -- could probably just use skipList and then check if a detector is Enr or Nat but whatever
enrDetList = [112, 113, 114, 122, 123, 132, 133, 134, 152, 153, 154, 161, 162, 163, 164, 172, 173, 174, 231, 232, 253, 254, 262, 273]
natDetList = [121, 141, 142, 143, 144, 145, 151, 171, 222, 223, 241, 242, 244, 251]

def main():
    """ need to think about how this code should interact w/ lat-expo
    and the spectrum fitting codes.
    """
    # combineSloEff(makePlots=True, writeDB=False, runMC=False, seedNum=1)
    GPXSloEff(makePlots=True, writeDB=False)


def combineSloEff(makePlots=False, writeDB=False, runMC=False, seedNum=1):
    """
    Inputs: m2s238 data (lat2-eff-data-95.npz, from lat2::setSloCut)
    Outputs:
    - npz file w/ efficiency curves
    - combined weibull parameters for enr/nat and toy MC data
    Then in lat-expo:
    - make the efficiency curve per dataset and
      combine w/ the trigger & riseNoise efficiencies
    """

    # load lat2 data
    f = np.load(os.environ['LATDIR'] + '/data/lat2-eff-data-95.npz')
    effData = f['arr_0'].item()
    detList = effData.keys()
    xVals = effData['112'][7][1:]

    # combine individual detector m2s238 histos into overall enr and nat histos
    hPassNat, hAllNat = np.zeros(len(xVals)), np.zeros(len(xVals))
    hPassEnr, hAllEnr = np.zeros(len(xVals)), np.zeros(len(xVals))
    eff5, eff10, cts5, cts10 = {}, {}, {}, {}

    for det in detList:
        if det in skipList:
            continue

        # NOTE: detectors w/ more counts contribute more
        if detInfo.isEnr(det):
            hPassEnr += effData[det][4][1:]
            hAllEnr += effData[det][6][1:]
        else:
            hPassNat += effData[det][4][1:]
            hAllNat += effData[det][6][1:]

        # save efficiencies at 5 and 10 kev for plot 1
        eff5[det] = effData[det][4][4] / effData[det][6][4]
        eff10[det] = effData[det][4][9] / effData[det][6][9]
        cts5[det] = effData[det][4][4]
        cts10[det] = effData[det][4][9]

    hEffEnr = np.nan_to_num(hPassEnr / hAllEnr)
    hEffNat = np.nan_to_num(hPassNat / hAllNat)

    # calculate CI's for each histogram bin
    enrCILo, enrCIHi = proportion.proportion_confint(hPassEnr, hAllEnr, alpha=0.1, method='beta')
    natCILo, natCIHi = proportion.proportion_confint(hPassNat, hAllNat, alpha=0.1, method='beta')

    # ---- fit overall enr/nat efficiency to a weibull function ----
    fitBnd = ((1, -20, 0, 0.5), (np.inf, np.inf, np.inf, 0.99))
    fitInit = [1, -1.6, 2.75, 0.95]
    poptNat, pcovNat = curve_fit(wl.weibull, xVals, hEffNat, p0=fitInit, bounds=fitBnd)
    poptEnr, pcovEnr = curve_fit(wl.weibull, xVals, hEffEnr, p0=fitInit, bounds=fitBnd)
    effNat = wl.weibull(xVals, *poptNat)
    effEnr = wl.weibull(xVals, *poptEnr)

    # ---- fitSlo efficiency uncertainty, method 1 ----
    # use the diagonal as the uncertainty (ignoring correlations)

    zVal = 1.645
    sigmaEnr = np.sqrt(np.diagonal(pcovEnr))
    sigmaNat = np.sqrt(np.diagonal(pcovNat))

    effEnrHi = wl.weibull(xVals, *(np.array(poptEnr) + zVal * sigmaEnr))
    effEnrLo = wl.weibull(xVals, *(np.array(poptEnr) - zVal * sigmaEnr))
    effNatHi = wl.weibull(xVals, *(np.array(poptNat) + zVal * sigmaNat))
    effNatLo = wl.weibull(xVals, *(np.array(poptNat) - zVal * sigmaNat))

    # ---- fitSlo efficiency uncertainty, method 2 ----
    # https://stats.stackexchange.com/questions/135749/\
    # confidence-intervals-of-fitted-weibull-survival-function

    effEnrHi2 = np.exp(np.log(effEnr) * np.exp(zVal / np.sqrt(hAllEnr)))
    effEnrLo2 = np.exp(np.log(effEnr) * np.exp(-1 * zVal / np.sqrt(hAllEnr)))
    effNatHi2 = np.exp(np.log(effNat) * np.exp(zVal / np.sqrt(hAllNat)))
    effNatLo2 = np.exp(np.log(effNat) * np.exp(-1 * zVal / np.sqrt(hAllNat)))

    # ---- run toy MC to get FINAL fitSlo efficiency uncertainty ----
    if runMC:

        np.random.seed(seedNum)

        xLo, xHi = 0, 200
        xCoarse = np.arange(xLo, xHi, 0.1)
        hEnr, hNat = [], [] # store toymc histo efficiencies

        nMC = 10000
        for i in range(nMC):

            if i % 100 == 0:
                wl.update_progress(float(i)/nMC)

            # vary the spectra randomly (toy MC method) and re-fit each one
            ePass = np.random.poisson(hPassEnr)
            nPass = np.random.poisson(hPassNat)
            eAll = np.random.poisson(hAllEnr)
            nAll = np.random.poisson(hAllNat)

            eEff = np.nan_to_num(ePass / eAll)
            nEff = np.nan_to_num(nPass / nAll)

            poptEnr,_ = curve_fit(wl.weibull, xVals, eEff, p0=fitInit, bounds=fitBnd)
            poptNat,_ = curve_fit(wl.weibull, xVals, nEff, p0=fitInit, bounds=fitBnd)

            effCoarseEnr = wl.weibull(xCoarse, *poptEnr)
            effCoarseNat = wl.weibull(xCoarse, *poptNat)

            hEnr.append(effCoarseEnr)
            hNat.append(effCoarseNat)

            # diagnostic plot -- don't delete
            # hScale = np.amax(hAllEnr)
            # plt.plot(xCoarse, effCoarseEnr, '-r')
            # plt.plot(xVals, hAllEnr / hScale, ls='steps', c='k', label='all m2s238 enr evts')
            # plt.plot(xVals, hPassEnr / hScale, ls='steps', c='b', label='orig passing')
            # plt.plot(xVals, ePass / hScale, ls='steps', c='m', label='toyMC variation')
            # plt.axvline(1, c='g', label="1 keV eff: {:.2f}".format(wl.weibull(1, *poptEnr)))
            # plt.xlabel("Energy (keV)", ha='right', x=1)
            # plt.xlim(0, 60)
            # plt.legend()
            # plt.tight_layout()
            # plt.savefig("./plots/toyMCEff.pdf")
            # exit()

        hEnr, hNat = np.vstack(hEnr), np.vstack(hNat)
        toyEffEnr = hEnr.mean(axis=0)
        toyEffNat = hNat.mean(axis=0)
        toyStdEnr = hEnr.std(axis=0)
        toyStdNat = hNat.std(axis=0)
        np.savez("./data/lat-toymc-eff.npz", toyEffEnr, toyEffNat, toyStdEnr, toyStdNat)

        # save results into calDB
        if writeDB:
            dbKey = "fitSlo_Combined_m2s238_eff95"
            dbVals = {0: [*poptEnr, *sigmaEnr], # 0: enr
                      1: [*poptNat, *sigmaNat]} # 1: nat
            print("Writing DB vals for key:", dbKey)
            # pprint(dbVals)
            dsi.setDBRecord({"key":dbKey, "vals":dbVals},
                            forceUpdate=True, calDB=calDB, pars=pars)
            pprint(dsi.getDBRecord(dbKey, False, calDB, pars))
            print("DB filled.")

    # ---------- make some plots ----------
    if makePlots:

        # 1.
        # individual detector efficiency at 5 & 10 keV, vs number of counts
        fig, (p0, p1) = plt.subplots(1, 2, figsize=(10,5))

        nEnr = len([det for det in eff5 if detInfo.isEnr(det)])
        nNat = len(eff5) - nEnr
        cmapEnr = plt.cm.get_cmap('nipy_spectral', nEnr+1)
        cmapNat = plt.cm.get_cmap('jet', nNat+1)
        iEnr, iNat = 0, 0

        for det in eff5:
            if detInfo.isEnr(det):
                iEnr += 1
                p, idx, cm = p0, iEnr, cmapEnr
            else:
                iNat += 1
                p, idx, cm = p1, iNat, cmapNat

            p.plot([eff5[det], eff10[det]], [cts5[det], cts10[det]],
                   '-', c=cm(idx), lw=1, label="C{}P{}D{}".format(*det))
            p.plot(eff5[det], cts5[det], 'v', ms=5, c=cm(idx))
            p.plot(eff10[det], cts10[det], 'o', ms=5, c=cm(idx))

        p0.plot(np.nan, 'kv', ms = 5, label='5 keV')
        p0.plot(np.nan, 'ko', ms = 5, label='10 keV')
        p1.plot(np.nan, 'kv', ms = 5, label='5 keV')
        p1.plot(np.nan, 'ko', ms = 5, label='10 keV')
        p0.legend(ncol=3, fontsize=8)
        p1.legend(ncol=3, fontsize=8)
        p0.set_xlabel("Enr. Efficiency", ha='right', x=1)
        p1.set_xlabel("Nat. Efficiency", ha='right', x=1)
        p0.set_ylabel("Counts (passing, m2s238)", ha='right', y=1)
        plt.tight_layout()
        plt.savefig("./plots/countsVsEff.pdf")
        plt.close()

        # 2.
        # individual and combined detector efficiencies
        fsD = dsi.getDBRecord("fitSlo_cpd_eff95", False, calDB, pars)

        fig, (p0, p1) = plt.subplots(1, 2, figsize=(10,5))
        iEnr, iNat = 0, 0

        for det in eff5:
            if detInfo.isEnr(det):
                iEnr += 1
                p, idx, cm = p0, iEnr, cmapEnr
            else:
                iNat += 1
                p, idx, cm = p1, iNat, cmapNat

            wbPars = fsD[int(det)]
            c, loc, scale, amp = wbPars[3], wbPars[4], wbPars[5], wbPars[2]
            effDet = wl.weibull(xVals, c, loc, scale, amp)

            p.plot(xVals, effDet, alpha=0.4, c=cm(idx), lw=2,
                   label='C{}P{}D{}'.format(*det))

        p0.plot(xVals, effEnr, lw=4, color='k', label='Enr, Combined')
        p1.plot(xVals, effNat, lw=4, color='k', label='Nat, Combined')
        p0.legend(loc=4, ncol=3, fontsize=10)
        p1.legend(loc=4, ncol=3, fontsize=10)
        p0.set_xlabel("Energy (keV)", ha='right', x=1)
        p1.set_xlabel("Energy (keV)", ha='right', x=1)
        p0.set_ylabel("Efficiency", ha='right', y=1)
        plt.tight_layout()
        plt.savefig("./plots/effCombined.pdf")
        plt.close()

        # 3.
        # uncertainties on the combined efficiencies

        fig, (p0, p1) = plt.subplots(1, 2, figsize=(10,5))

        # enriched
        p0.errorbar(xVals, effEnr, yerr=[hEffEnr - enrCILo, enrCIHi - hEffEnr],
                    color='k', lw=1, fmt='o', capsize=2, ms=3,
                    label="Enriched, Combined")

        # NOTE: method 1 swaps the high/low boundaries at the turning point
        # I.E. DO NOT USE!
        # p0.plot(xVals, effEnrHi, 'r-', lw=1, label="Method 1 (hi)")
        # p0.plot(xVals, effEnrLo, 'g-', lw=1, label="Method 1 (lo)")
        # p0.fill_between(xVals, effEnrLo, effEnrHi, color='b', alpha=0.5, label='Method 1')

        # NOTE: method 2 looks like the efficiency uncertainties are too small
        # I.E. DO NOT USE!
        # p0.plot(xVals, effEnrHi2, 'r-', lw=1, label="Method 2 (hi)")
        # p0.plot(xVals, effEnrLo2, 'g-', lw=1, label="Method 2 (lo)")
        # p0.fill_between(xVals, effEnrLo2, effEnrHi2, color='r', alpha=0.5, label='Method 2')

        # Method 3 - uncertainty from Toy MC results
        f = np.load("./data/lat-toymc-eff.npz")
        toyEffEnr, toyEffNat, toyStdEnr, toyStdNat = [f[k] for k in f]
        xLo, xHi = 0, 200
        xCoarse = np.arange(xLo, xHi, 0.1)

        p0.plot(xCoarse, toyEffEnr, c='r', lw=2, label='Toy MC Efficiency')
        effLo = toyEffEnr - zVal * toyStdEnr
        effHi = toyEffEnr + zVal * toyStdEnr
        p0.fill_between(xCoarse, effLo, effHi, color='g', alpha=0.4, label='Toy MC Uncert.')

        # natural
        p1.errorbar(xVals, effNat, yerr=[hEffNat - natCILo, natCIHi - hEffNat],
                    color='k', lw=1, fmt='o', capsize=2, ms=3,
                    label="Natural, Combined")

        p1.plot(xCoarse, toyEffNat, c='r', lw=2, label="Toy MC Efficiency")
        effLo = toyEffNat - zVal * toyStdNat
        effHi = toyEffNat + zVal * toyStdNat
        p1.fill_between(xCoarse, effLo, effHi, color='b', alpha=0.3, label='Toy MC Uncert.')

        p0.set_xlabel("Energy (keV)", ha='right', x=1)
        p1.set_xlabel("Energy (keV)", ha='right', x=1)
        p0.set_ylabel("Efficiency", ha='right', y=1)
        p0.legend()
        p1.legend()
        p0.set_xlim(0, 20)
        p0.set_ylim(0.4, 1)
        p1.set_xlim(0, 20)
        p1.set_ylim(0.4, 1)
        plt.tight_layout()
        plt.savefig("./plots/combinedEff.pdf")



def GPXSloEff(makePlots=False, writeDB=False):
    """
        This function does 2 things:
        1) It calculates the total mHL == 1 efficiency at the 238 keV peak (not really required for anything but Wenqin wanted it)
        2) It uses the mHL == 1 efficiency at the 238 keV peak and performs the sideband method of calculating the cut efficiency.

        The sideband method was suggested by Wenqin as a cross-check to the cut efficiency. Instead of breaking up the m2s238 event pairs and evaluating the efficiency, we keep the pairs in tact and perform a background subtraction on the energy window slightly below the peak. We then can back out the single detector efficiency at low energy by using the mHL==1 efficiency of the 238 keV peak at high energy.

    Requires:
        CalPairHit_WithDbCut.h5 and CalPairHit_WithDbCut_mH1.h5, both generated in lat2.py
        These files are essentially the m2s238 hit data with a Pass/Fail for fitSlo
        Weibull fit parameters of the Combined efficiency

    Writes:
        fitSlo_Sideband_m2s238_eff95 (containing the sideband Weibull fit parameters of Enr and Nat) key to the DB

    """

    df = pd.read_hdf('{}/data/CalPairHit_WithDbCut.h5'.format(os.environ['LATDIR']))
    windowSize = 0.2
    xVals = [round(windowSize*i, 2) for i in range(int(1/windowSize), int((50.+windowSize)/windowSize))]
    fitBnd = ((1,-20,0,0.5),(np.inf,np.inf,np.inf, 0.99)) # eFitHi=30 and these works!
    initialGuess = [1, -1.6, 2.75, 0.95]
    fListPeakE, fListBkgE = [], []
    cListPeakE, cListBkgE = [], []
    fListPeakN, fListBkgN = [], []
    cListPeakN, cListBkgN = [], []

    for idx, er in enumerate(xVals):
        if idx%50 == 0:
            print('Current Energy: {:.2f} of {:.2f}'.format(er, xVals[-1]))
        dFullPeakE, dCutPeakE, dFullBkgE, dCutBkgE, dFullPeakN, dCutPeakN, dFullBkgN, dCutBkgN = runCutVals(df, er, windowSize=windowSize)
        fListPeakE.append(dFullPeakE)
        cListPeakE.append(dCutPeakE)
        fListBkgE.append(dFullBkgE)
        cListBkgE.append(dCutBkgE)
        fListPeakN.append(dFullPeakN)
        cListPeakN.append(dCutPeakN)
        fListBkgN.append(dFullBkgN)
        cListBkgN.append(dCutBkgN)


    # Grab total fitSlo efficiency from DB
    dbKey = "fitSlo_Combined_m2s238_eff95"
    fsN = dsi.getDBRecord(dbKey, False, calDB, pars)
    enrpars = fsN[0]
    natpars = fsN[1]
    EnrEff = wl.weibull(xVals, *(np.array(enrpars[:4])))
    NatEff = wl.weibull(xVals, *(np.array(natpars[:4])))

    # mHL==1 efficiency correction from high energy
    effScaleEnr, effScaleNat = getM1Efficiency()
    print('Scaling Factors: ', effScaleEnr, effScaleNat)

    effCorrE = (np.array(cListPeakE) - np.array(cListBkgE))/(np.array(fListPeakE) - np.array(fListBkgE))
    effCorrE /= effScaleEnr
    enr_ci_low, enr_ci_upp = proportion.proportion_confint(np.array(cListPeakE) - np.array(cListBkgE), np.array(fListPeakE) - np.array(fListBkgE), alpha=0.1, method='beta')

    effCorrN = (np.array(cListPeakN) - np.array(cListBkgN))/(np.array(fListPeakN) - np.array(fListBkgN))
    effCorrN /= effScaleNat
    nat_ci_low, nat_ci_upp = proportion.proportion_confint(np.array(cListPeakN) - np.array(cListBkgN), np.array(fListPeakN) - np.array(fListBkgN), alpha=0.1, method='beta')

    poptEnr, pcovenr = curve_fit(wl.weibull, xVals, effCorrE, p0=initialGuess, bounds=fitBnd)
    poptNat, pcovnat = curve_fit(wl.weibull, xVals, effCorrN, p0=initialGuess, bounds=fitBnd)

    effEFit = wl.weibull(xVals, *poptEnr)
    effNFit = wl.weibull(xVals, *poptNat)

    # Comparison of the parameters
    print(poptEnr, enrpars[:4])
    print(poptNat, natpars[:4])

    sigmaEnr = np.sqrt(np.diagonal(pcovenr))
    sigmaNat = np.sqrt(np.diagonal(pcovnat))

    if writeDB:
        dbKeyFill = "fitSlo_Sideband_m2s238_eff95"
        dbVals = {0: [*poptEnr, *sigmaEnr], # 0: enr
                1: [*poptNat, *sigmaNat]} # 1: nat
        print('Writing dbVals:', dbVals)
        dsi.setDBRecord({"key":dbKeyFill, "vals":dbVals}, forceUpdate=True, calDB=calDB, pars=pars)
        print("DB filled:",dbKeyFill)

        # Get the record
        fsFill = dsi.getDBRecord(dbKeyFill, False, calDB, pars)
        print(fsFill)

    if makePlots:
        fig1, ax1 = plt.subplots(nrows=2, ncols=2, figsize=(15,10))
        ax1 = ax1.flatten()
        ax1[0].errorbar(xVals, effCorrE, yerr=[effCorrE - enr_ci_low/effScaleEnr, enr_ci_upp/effScaleEnr - effCorrE], color='k', linewidth=0.8, fmt='o', alpha=0.75, capsize=2, label='Sideband Method')
        ax1[0].plot(xVals, effEFit, 'b', lw=3, label='Sideband Fit Efficiency')
        ax1[0].plot(xVals, EnrEff, 'r', lw=3, label='Central Fit Efficiency')
        ax1[0].set_title('Enriched Efficiency')
        ax1[0].set_ylabel('Efficiency')
        ax1[0].legend()
        ax1[2].plot(xVals, EnrEff - effEFit)
        ax1[2].set_xlabel('Energy (keV)')
        ax1[2].set_ylabel('Efficiency Difference (Central - Sideband)')

        ax1[1].errorbar(xVals, effCorrN, yerr=[effCorrN - nat_ci_low/effScaleNat, nat_ci_upp/effScaleNat - effCorrN], color='k', linewidth=0.8, fmt='o', alpha=0.75, capsize=2, label='Sideband Method')
        ax1[1].plot(xVals, effNFit, 'b', lw=3, label='Sideband Fit Efficiency')
        ax1[1].plot(xVals, NatEff, 'r', lw=3, label='Central Fit Efficiency')
        ax1[1].set_title('Natural Efficiency')
        ax1[1].set_ylabel('Efficiency')
        ax1[1].legend()
        ax1[3].plot(xVals, EnrEff - effEFit)
        ax1[3].set_xlabel('Energy (keV)')
        ax1[3].set_ylabel('Efficiency Difference (Central - Sideband)')

        plt.tight_layout()
        fig1.savefig('./plots/GPXEfficiencyComparison.png')


def getM1Efficiency():
    df = pd.read_hdf(os.environ['LATDIR']+'/data/CalPairHit_WithDbCut_mH1.h5')
    dfg = df.groupby(['cpd1'])

    dEnrFull, dEnrCut = 0, 0
    dNatFull, dNatCut = 0, 0

    for name, g in dfg:
        bCuts = g['trapENFCal1'].values
        aCuts = g.loc[(df['Pass1'] == True), 'trapENFCal1'].values
        if name in enrDetList:
            dEnrFull += len(bCuts)
            dEnrCut += len(aCuts)
        else:
            dNatFull += len(bCuts)
            dNatCut += len(aCuts)

    effEnr, effNat = float(dEnrCut/dEnrFull), float(dNatCut/dNatFull)
    return np.sqrt(effEnr), np.sqrt(effNat)


def runCutVals(df, eVal=0., windowSize = 2):
    """
        Why do I loop through all detectors? I don't know but in the beginning I did each detector individually and I'm too lazy to change it
    """

    dfg = df.groupby(['cpd1'])

    eMin = round(eVal - windowSize/2, 2)
    eMax = round(eMin + windowSize, 2)
    dFullPeakE, dFullBkgE = 0, 0
    dCutPeakE, dCutBkgE = 0, 0
    dFullPeakN, dFullBkgN = 0, 0
    dCutPeakN, dCutBkgN = 0, 0

    for name, g in dfg:
        valsFull = g['trapENFCal1'].loc[(g['trapENFCal1']>eMin) & (g['trapENFCal1']<eMax)].values + g['trapENFCal2'].loc[(g['trapENFCal1']>eMin) & (g['trapENFCal1']<eMax)].values

        valsCut = g['trapENFCal1'].loc[(g['Pass1']==True) & (g['Pass2']==True) & (g['trapENFCal1']>eMin) & (g['trapENFCal1']<eMax)].values + g['trapENFCal2'].loc[(g['Pass1']==True) & (g['Pass2']==True) & (g['trapENFCal1']>=eMin) & (g['trapENFCal1']<=eMax)].values
        if name in enrDetList:
            dFullPeakE += len(valsFull[(valsFull > 237.28) & (valsFull < 239.46)])
            dCutPeakE += len(valsCut[(valsCut > 237.28) & (valsCut < 239.46)])
            dFullBkgE += len(valsFull[(valsFull > 235) & (valsFull < 237.18)])
            dCutBkgE += len(valsCut[(valsCut > 235) & (valsCut < 237.18)])
        elif name in natDetList:
            dFullPeakN += len(valsFull[(valsFull > 237.28) & (valsFull < 239.46)])
            dCutPeakN += len(valsCut[(valsCut > 237.28) & (valsCut < 239.46)])
            dFullBkgN += len(valsFull[(valsFull > 235) & (valsFull < 237.18)])
            dCutBkgN += len(valsCut[(valsCut > 235) & (valsCut < 237.18)])

    return dFullPeakE, dCutPeakE, dFullBkgE, dCutBkgE, dFullPeakN, dCutPeakN, dFullBkgN, dCutBkgN


if __name__ == "__main__":
    main()
