#!/usr/bin/env python3
import os
import numpy as np
import tinydb as db
from statsmodels.stats import proportion
from scipy.optimize import curve_fit
from pprint import pprint
import matplotlib.pyplot as plt
plt.style.use('pltReports.mplstyle')
import dsi
detInfo = dsi.DetInfo()
import waveLibs as wl
dbFile = '%s/calDB-v2.json' % (dsi.latSWDir) # v2 matches LAT's v2 tag
calDB = db.TinyDB(dbFile)
pars = db.Query()

# Skip these detectors because of low statistics
skipList = ['111', '211', '214', '221', '261', '274']


def main():
    """
    Inputs: m2s238 data (lat2-eff-data-95.npz, from lat2::setSloCut)
    Outputs:
    - npz file w/ efficiency curves
    - combined weibull parameters for enr/nat and toy MC data
    Then in `lat-expo::getEfficiencyROOT`:
    - make the efficiency curve per dataset and
      combine w/ the trigger & riseNoise efficiencies
    """

    getEff(writeDB=False, runMC=True)
    # plotEff() # TODO: move brian's plots in combinedEff here


def getEff(writeDB=False, runMC=False, seedNum=1):
    """ slim version of brian's combinedEff """

    # load lat2 data
    f = np.load(os.environ['LATDIR'] + '/data/lat2-eff-data-95.npz')
    effData = f['arr_0'].item()
    detList = effData.keys()
    xVals = effData['112'][7][1:]

    # combine individual detector m2s238 histos into overall enr and nat histos
    hPassNat, hAllNat = np.zeros(len(xVals)), np.zeros(len(xVals))
    hPassEnr, hAllEnr = np.zeros(len(xVals)), np.zeros(len(xVals))
    for det in detList:
        if det in skipList:
            continue
        if detInfo.isEnr(det):
            hPassEnr += effData[det][4][1:]
            hAllEnr += effData[det][6][1:]
        else:
            hPassNat += effData[det][4][1:]
            hAllNat += effData[det][6][1:]
    hEffEnr = np.nan_to_num(hPassEnr / hAllEnr)
    hEffNat = np.nan_to_num(hPassNat / hAllNat)

    # calculate CI's for each histogram bin
    natCILo, natCIHi = proportion.proportion_confint(hPassNat, hAllNat, alpha=0.1, method='beta')
    enrCILo, enrCIHi = proportion.proportion_confint(hPassEnr, hAllEnr, alpha=0.1, method='beta')

    # fit overall enr/nat efficiency to a weibull function
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

    # ---- run toy MC to get FINAL efficiency uncertainty ----
    if not runMC:
        return

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
    np.savez("./data/lat-eff.npz", toyEffEnr, toyEffNat, toyStdEnr, toyStdNat)

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


    # diagnostic plot
    plt.plot(xCoarse, toyEffEnr, ls='steps', lw=1, c='b', label='enr eff')
    plt.plot(xCoarse, toyEffEnr + zVal * toyStdEnr, ls='steps', lw=1, c='b')
    plt.plot(xCoarse, toyEffEnr - zVal * toyStdEnr, ls='steps', lw=1, c='b')
    plt.plot(xCoarse, toyEffNat, ls='steps', lw=1, c='r', label='nat eff')
    plt.plot(xCoarse, toyEffNat + zVal * toyStdNat, ls='steps', lw=1, c='r')
    plt.plot(xCoarse, toyEffNat - zVal * toyStdNat, ls='steps', lw=1, c='r')
    plt.xlabel("Energy (keV)", ha='right', x=1)
    plt.ylabel("Efficiency", ha='right', y=1)
    plt.legend()
    plt.xlim(0, 15)
    plt.ylim(0.7, 1)
    plt.tight_layout()
    plt.savefig("./plots/combinedEff.pdf")


if __name__ == "__main__":
    main()
