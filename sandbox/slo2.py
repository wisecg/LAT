#!/usr/bin/env python3
import sys, os, time
import numpy as np
import tinydb as db
from scipy.optimize import curve_fit
from statsmodels.stats import proportion

import matplotlib as mpl
import matplotlib.pyplot as plt
plt.style.use('../pltReports.mplstyle')

import waveLibs as wl
import dsi
bkg = dsi.BkgInfo()
cal = dsi.CalInfo()
det = dsi.DetInfo()

# a very important parameter, use 90 or 95
global pctTot
# pctTot = 90
pctTot = 95 # < -- use this one
print("Using pctTot ==",pctTot)

def main():

    # loadSloData()
    # sp_correction()
    # siggenwf()
    # sim_fast_pulse()
    # save_m2s238_hits()
    acceptance_study()


def loadSloData(key):
    """ Load files generated by scanRuns, return data in a dict.
    To avoid confusion, must specify a key from runsCal.json .
    """
    if key not in cal.GetKeys():
        print("Unknown key!")
        return None
    else:
        print("Loading eff data for key:",key)

    # output dict
    eff = {}
    eff["hitE"] = []  # [hitE1, hitE2 , ...] (remove sub-list of input format)
    eff["chan"] = []  # [chan1, chan2 , ...]
    eff["fSlo"] = []  # [fSlo1, fSlo2, ...]
    eff["rise"] = []  # [rise1, rise2, ...]
    eff["run"]  = []  # [run1, run2, ...]
    eff["cIdx"] = []  # [cIdx1, cIdx2, ...]
    eff["spec"] = []  # array of fitSlo histo dicts (i should have used pandas probably)
    eff["specX"] = [] # x values for "spec" histos (all the same)
    for ci in range(cal.GetIdxs(key)):
        eFile = "%s/eff_%s_c%d.npz" % (dsi.effDir, key, ci)
        if not os.path.isfile(eFile):
            print("File not found:",eFile)
            continue
        f = np.load(eFile)
        evtIdx = f['arr_0']          # m2s238 event [[run,iE,cIdx] , ...]
        evtSumET = f['arr_1']        # m2s238 event [sumET , ...]
        evtHitE = f['arr_2']         # m2s238 event [[hitE1, hitE2] , ...]
        evtChans = f['arr_3']        # m2s238 event [[chan1, chan2] , ...]
        thrCal = f['arr_4'].item()   # {ch : [run,thrM,thrS,thrK] for ch in goodList(ds)}
        thrFinal = f['arr_5'].item() # {ch : [thrAvg, thrDev] for ch in goodList(ds)}
        evtCtr = f['arr_6']          # num m2s238 evts
        totCtr = f['arr_7']          # num total evts
        runTime = f['arr_8']         # cal run time
        fSloSpec = f['arr_9'].item() # fitSlo histos (all hits) {ch:[h10, h200, h238] for ch in chList}
        fSloX = f['arr_10']          # xVals for fitSlo histos
        evtSlo = f['arr_11']         # m2s238 event [[fSlo1, fSlo2], ...]
        evtRise = f['arr_12']        # m2s238 event [[rise1, rise2], ...]

        # remove the hit pair
        for i in range(len(evtHitE)):
            eff["hitE"].extend(evtHitE[i])
            eff["chan"].extend(evtChans[i])
            eff["fSlo"].extend(evtSlo[i])
            eff["rise"].extend(evtRise[i])
            eff["run"].extend([evtIdx[i][0], evtIdx[i][0]])
            eff["cIdx"].extend([evtIdx[i][2], evtIdx[i][2]])
        eff["spec"].append(fSloSpec)
        eff["specX"] = fSloX # this doesn't change

    # convert to numpy arrays and return
    for key in eff:
        if key=="spec": continue
        eff[key] = np.asarray(eff[key])
    return eff


def sp_correction():

    from ROOT import TFile, TH1D

    tf = TFile("../data/lat-expo-efficiency-corr.root")
    surfPerc = tf.Get("surfPerc")

    xSp, ySp, xpb = wl.npTH1D(surfPerc)

    plt.plot(xSp, ySp, ls='steps')

    # plt.show()

    # brian's cal_det_pairs data directory
    "/global/projecta/projectdirs/majorana/users/bxyzhu/LATv2"
    # df1 = pd.read_hdf('{}/DS5_Sim_HitData.h5'.format(inDir)) #228M
    # df2 = pd.read_hdf('{}/DS5_Sim_HitData_G41004.h5'.format(inDir)) # 137M


def siggenwf():

    # import seaborn as sns
    # sns.set(style='darkgrid', context='poster')

    f = np.load("../data/siggenwf.npz")

    ts = f['arr_0']
    wf1 = f['arr_1'] # waveform w/ baseline
    wf2 = f['arr_2'] # waveform w/ just 0s as baseline

    bl = np.average(wf1[:200])

    plt.plot(ts, wf1, lw=4, c='b', label="Siggen WF w/ DS0 Baseline")
    plt.plot(ts, wf2 + bl, lw=2, c='r', label="Siggen WF")
    plt.legend(loc=0)
    plt.xlabel("Time (ns)", ha='right', x=1)
    plt.ylabel("Voltage (ADC)", ha='right', y=1)
    plt.tight_layout()
    # plt.show()
    plt.savefig("../plots/lat-siggenwf.pdf")


def sim_fast_pulse():
    """ adapted from sim_wf_study, thx brian """
    import pandas as pd
    import seaborn as sns

    dfSig = pd.read_hdf('../data/SimPSA_P42574A.h5')
    dfSig['Distance'] = dfSig['R']*dfSig['R'] + dfSig['Z']+dfSig['Z']
    dfSig['fitSloShift'] = dfSig['fitSlo'] - 84.5
    dfSig['trapENFCal'] = dfSig['Amp']*0.3959
    dfCut = dfSig.loc[(dfSig['fitSlo'] > 3) & (dfSig['Amp'] > 2)]
    # dfCut = dfSig.loc[(dfSig['Amp'] < 14) & (dfSig['fitSlo'] > 2)]

    # g1 = sns.lmplot(x='Amp', y='fitSlo', data=dfCut, fit_reg=False, size=7, scatter_kws={'alpha':0.3}, legend_out=False)
    # g1.set(yscale='log')
    # plt.subplots_adjust(top=0.95)
    # g1.fig.suptitle('fitSlo vs Amplitude (Simulated Waveforms)')
    # g1.set_axis_labels('Amplitude (ADC)', 'fitSlo')

    # plt.plot()
    # print(dfCut['fitSloShift'])

    plt.plot(dfCut['trapENFCal'], dfCut['fitSloShift'], '.k', ms=2)
    plt.xlabel("Sim. Energy (~keV)", ha='right', x=1)
    plt.ylabel("fitSlo", ha='right', y=1)

    plt.show()
    # plt.show()


def xGaussPDF(x,k,loc,scale,amp):
    """ amp=overall multip factor, loc=mu, sc=sigma, k=tau """
    from scipy.stats import exponnorm
    return amp * exponnorm.pdf(x,k,loc,scale)


def xGaussCDF(x,k,loc,scale,amp):
    from scipy.stats import exponnorm
    return exponnorm.cdf(x,k,loc,scale)


def save_m2s238_hits():

    dsList = [0,1,2,3,4,5]
    # dsList = [1]
    detList = det.allDets
    detIDs = det.allDetIDs

    # this is one of the main output dicts, see lat2
    shiftVals = {} # store the UNshifted 10-200 fitSlo peak vals

    # save the hits in the "low energy" region
    hitE = {cpd:[] for cpd in detList}
    fSlo = {cpd:[] for cpd in detList} # SHIFTED

    # load the slowness data we got from scanRunsSlow()
    for ds in dsList:
        for key in cal.GetKeys(ds):

            # get channels in this DS and map back to CPD
            chList = det.getGoodChanList(ds)
            mod = -1
            if "m1" in key:
                mod = 1
                chList = [ch for ch in chList if ch < 1000]
            if "m2" in key:
                mod = 2
                chList = [ch for ch in chList if ch > 1000]
            cpdList = [det.getChanCPD(ds,ch) for ch in chList]
            chMap = {det.getChanCPD(ds,ch):ch for ch in chList}

            eff = loadSloData(key)
            nCal = cal.GetNCalIdxs(ds,mod)
            shiftVals[key] = {ci:None for ci in range(nCal)}

            # loop over calIdx's
            for ci in range(nCal):

                # save the fs shift value for each cpd/ch in each calIdx, and fill the hit lists for plotting
                shiftVals[key][ci] = {cpd:None for cpd in cpdList}

                for cpd in cpdList:

                    # load fitSlo hist of all cal hits in this channel 10-200 keV
                    ch = chMap[cpd]
                    h1 = eff["spec"][ci][ch][1]
                    x1 = eff["specX"]
                    # shiftSpec[cpd] = np.add(shiftSpec[cpd], h1) # unused
                    if np.sum(h1)==0:
                        # print("ci %d  cpd %d  no counts" % (ci, cpd))
                        shiftVals[key][ci][cpd] = -1
                        continue

                    # get mode (maximum) of the 10-200 hits and save it
                    fMax = x1[np.argmax(h1)]
                    shiftVals[key][ci][cpd] = fMax

                    # fill the hit lists
                    idx = np.where((eff["cIdx"]==ci) & (eff["chan"]==ch))
                    fSlo[cpd].extend([f - fMax for f in eff["fSlo"][idx]])
                    hitE[cpd].extend([e for e in eff["hitE"][idx]])

    for cpd in hitE:
        hitE[cpd] = np.asarray(hitE[cpd])
        fSlo[cpd] = np.asarray(fSlo[cpd])

    np.savez("../data/slo2-m2s238-hits.npz", hitE, fSlo, shiftVals)


def getTot(yVals, pctLo=1, pctHi=91, cutTot=None):
    """ take an array of fitSlo values and find the 90% value """

    # method 1 -- cool xgauss fit
    # fLo, fHi = np.percentile(yVals,pctLo) * 1.2, np.percentile(yVals,pctHi) * 1.2
    # nb = 30
    # fpb = (fHi-fLo)/nb
    # x, h = wl.GetHisto(yVals, fLo, fHi, fpb, shift=False)
    # fMax = x[np.argmax(h)]
    # tau, mu, sig, amp = 10, fMax, 10, 10000
    # bnd = [[0,-np.inf,0,0],[np.inf, np.inf, np.inf, np.inf]]
    # np.seterr(invalid='ignore', over='ignore')
    # try:
    #     popt,_ = curve_fit(xGaussPDF, x, h, p0=(tau, mu, sig, amp), bounds=bnd) # tau, mu, sig, amp
    # except RuntimeError:
    #     popt = 0, 0, 0, 0
    #     pass
    # np.seterr(invalid='warn', over='warn')
    # xFunc = np.arange(fLo, fHi, 0.01)
    # xgInt = xGaussCDF(xFunc, *popt)
    #
    # xgTot = 0
    # idx = np.where(xgInt > 0.9)
    # if len(idx[0])!=0:
    #     xgTot = xFunc[idx][0]
    #
    # if cutTot is not None:
    #     effTot = xGaussCDF(cutTot, *popt)
    #     return xgTot, effTot
    #
    # return xgTot

    # method 2 -- percentiles
    pctPass = np.percentile(yVals, pctTot)
    fShiftLo, fShiftHi, fpbShift = -100, 1000, 1  # this is better ...
    x, h = wl.GetHisto(yVals, fShiftLo, fShiftHi, fpbShift)
    tmp = np.cumsum(h)/np.sum(h)
    idx = np.where(x >= cutTot)
    effTot = tmp[idx][0]
    if cutTot is not None:
        return pctPass, effTot
    return pctPass


def acceptance_study():
    """ 3 fast pulse acceptance plots: ext pulser, simulated wf, and m2s238
    All for C1P6D3.
    Finally, plot the 3 efficiencies on the same plot.
    """
    xLo, xHi = 0, 50
    yLo, yHi = -100, 220

    xPassLo, xPassHi, xpbPass = xLo, xHi, 1     # "low energy" region
    fShiftLo, fShiftHi, fpbShift = -100, 1000, 1  # this is what LAT2 uses
    # fShiftLo, fShiftHi, fpbShift = -100, 20, 1  # debug

    # save the acc/efficiency and plot them below
    accExtPulser, accSimWF, effM2S238, effM2S238Corr = [], [], [], []

    fig = plt.figure(figsize=(8,8))
    p1 = plt.subplot(311)
    p2 = plt.subplot(312)
    p3 = plt.subplot(313)

    # ====== 1. ext pulser ======

    f = np.load("../data/lat-extPulser.npz")
    extData = f['arr_0'].item()  # {run: [pIdx, runTime, extChan, hitE, fSlo]}
    extInfo = {624:[7236, 74.75, 'r'], 688:[7249, 67.25, 'g'], 674:[7220, 65.75, 'b']}

    xVal, yVal = [], []

    fsAll = []
    for run in extData:
        ch = extData[run][2]
        cpd = det.getChanCPD(0, ch)
        if cpd != '163': continue
        if run in [7233]: continue
        hitE = extData[run][3]
        muE = np.average(hitE)
        if muE < 1.: continue
        if muE > xHi: continue
        fSlo = [fs - extInfo[ch][1] for fs in extData[run][4]] # shifted
        fsAll.extend(fSlo)

    # find the (pctTot)% value the way we do it in lat2
    xSloExt, hSloExt = wl.GetHisto(fsAll, fShiftLo, fShiftHi, fpbShift)
    max, avg, std, pct, wid = wl.getHistInfo(xSloExt,hSloExt)
    if pctTot==90:
        fsShiftCutExt = pct[2]
    elif pctTot==95:
        fsShiftCutExt = pct[3]
    else:
        print("WTF are you doing??")
        exit()

    # plot the slices
    for i, run in enumerate(extData):
        ch = extData[run][2]
        cpd = det.getChanCPD(0, ch)
        if cpd != '163': continue
        if run in [7233]: continue

        hitE = extData[run][3]
        fSlo = extData[run][4]
        muE = np.average(hitE)
        sdE = np.std(hitE)
        if muE < 1.: continue
        if muE > xHi: continue

        xE = muE * np.ones(len(fSlo))
        fSlo = [fs - extInfo[ch][1] for fs in extData[run][4]] # shifted
        p1.plot(xE, fSlo, '.b', ms=2)

        xgTot, accTot = getTot(fSlo, cutTot=fsShiftCutExt)
        p1.plot(muE, xgTot, '.r', ms=10)

        xT, yT = wl.GetHisto(fSlo, fShiftLo, fShiftHi, fpbShift)
        fCtr = xT[np.argmax(yT)]
        p1.plot(muE, fCtr, '.c', ms=10)

        xVal.append(xE)
        yVal.append(fSlo)

    p1.plot(np.nan, np.nan, ".b", ms=10, label="Ext.Pulser, C1P6D3")
    p1.axhline(fsShiftCutExt, c='orange', label="%d%% Overall Cut: %.1f" % (pctTot, fsShiftCutExt))

    # try histogram pass/fail method from LAT2
    hitExt, sloExt = [], []
    for i in range(len(xVal)):
        hitExt.extend(xVal[i])
        sloExt.extend(yVal[i])
    hitPass, hitFail = [], []
    for i in range(len(sloExt)):
        if sloExt[i] <= fsShiftCutExt: hitPass.append(hitExt[i])
        else: hitFail.append(hitExt[i])
    hitPass, hitFail = np.asarray(hitPass), np.asarray(hitFail)
    xELow, hPass = wl.GetHisto(hitPass, xPassLo, xPassHi, xpbPass, shift=False)
    xELow, hFail = wl.GetHisto(hitFail, xPassLo, xPassHi, xpbPass, shift=False)
    hTot = np.add(hPass, hFail)
    for i in range(len(xELow)):
        if hTot[i]==0: continue
        effTot = hPass[i]/hTot[i]
        xTot = xELow[i] - xpbPass/2
        accExtPulser.append([xTot, effTot])
    # plt.close()
    # plt.plot(xELow, hTot, ls='steps', c='k')
    # plt.plot(xELow, hPass, ls='steps', c='b')
    # plt.plot(xELow, hFail, ls='steps', c='r')
    # plt.show()
    # return


    # ====== 2. simulated wfs ======

    import pandas as pd
    dfSig = pd.read_hdf('../data/SimPSA_P42574A.h5')
    # dfSig['Distance'] = dfSig['R'] * dfSig['R'] + dfSig['Z'] + dfSig['Z']
    dfSig['trapENFCal'] = dfSig['Amp'] * 0.3959

    # find center value
    xt, yt = wl.GetHisto(dfSig['fitSlo'], fShiftLo, fShiftHi, fpbShift)
    simShift = xt[np.argmax(yt)] # it's 65.5
    dfSig['fitSloShift'] = dfSig['fitSlo'] - simShift
    # plt.close()
    # plt.plot(xt, yt, ls='steps')
    # plt.axvline(simShift, c='r')
    # xt, yt = wl.GetHisto(dfSig['fitSloShift'], fShiftLo, fShiftHi, fpbShift)
    # plt.plot(xt, yt, ls='steps')
    # plt.show()
    # return

    dfCut = dfSig.loc[(dfSig['fitSlo'] > 3) & (dfSig['Amp'] > 2)]

    p2.plot(dfCut['trapENFCal'], dfCut['fitSloShift'], '.b', ms=2, label="")
    p2.plot(np.nan, np.nan, ".b", ms=10, label="Simulated WFs, C1P6D3")

    # find the (pctTot)% value the way we do it in lat2
    xSloSim, hSloSim = wl.GetHisto(dfCut['fitSloShift'].values, fShiftLo, fShiftHi, fpbShift)
    max, avg, std, pct, wid = wl.getHistInfo(xSloSim,hSloSim)
    if pctTot==90:
        fsShiftCutSim = pct[2]
    elif pctTot==95:
        fsShiftCutSim = pct[3]
    else:
        print("WTF are you doing??")
        exit()


    # plot the (pctTot)% values for each slice and the centroids
    # xVals = sorted(list(set(dfCut['trapENFCal'])))
    xVals = np.unique(dfCut['trapENFCal'])
    for xE in xVals:
        dfSlice = dfCut.loc[dfCut['trapENFCal']==xE]
        yVals = dfSlice['fitSloShift']
        xgTot, accTot = getTot(yVals, pctLo=1, pctHi=98, cutTot=fsShiftCutSim)
        p2.plot(xE, xgTot, '.r', ms=10)

        xT, yT = wl.GetHisto(yVals, fShiftLo, fShiftHi, fpbShift)
        fCtr = xT[np.argmax(yT)]
        p2.plot(xE, fCtr, '.c', ms=10)

    p2.axhline(fsShiftCutSim, c='orange', label="%d%% Overall Cut: %.1f" % (pctTot, fsShiftCutSim))

    # try histogram pass/fail method from LAT2
    hitPass, hitFail = [], []
    hitVals, sloVals = dfCut['trapENFCal'].values, dfCut['fitSloShift'].values
    for i in range(len(sloVals)):
        if sloVals[i] <= fsShiftCutSim: hitPass.append(hitVals[i])
        else: hitFail.append(hitVals[i])
    hitPass, hitFail = np.asarray(hitPass), np.asarray(hitFail)
    xELow, hPass = wl.GetHisto(hitPass, xPassLo, xPassHi, xpbPass, shift=False)
    xELow, hFail = wl.GetHisto(hitFail, xPassLo, xPassHi, xpbPass, shift=False)
    hTot = np.add(hPass, hFail)
    for i in range(len(xELow)):
        if hTot[i]==0: continue
        effTot = hPass[i]/hTot[i]
        xTot = xELow[i] - xpbPass/2
        accSimWF.append([xTot, effTot])
    # plt.close()
    # plt.plot(xELow, hTot, ls='steps', c='k')
    # plt.plot(xELow, hPass, ls='steps', c='b')
    # plt.plot(xELow, hFail, ls='steps', c='r')
    # plt.show()
    # return


    # ====== 3. m2s238 data ======

    cpd = '163'
    f2 = np.load("../data/slo2-m2s238-hits.npz")
    hitE, fSlo, shiftVals = f2['arr_0'].item(), f2['arr_1'].item(), f2['arr_2'].item()
    eSlice = np.arange(xLo, xHi+1, 2) # 2 keV slices, 0 - 250
    xVal, yVal = [], []
    for i in range(len(eSlice)-1):
        eLo = eSlice[i]
        eHi = eSlice[i+1]
        idx = np.where((hitE[cpd] >= eLo) & (hitE[cpd] <= eHi))
        xVal.append((eHi-eLo)/2 + eLo)
        yVal.append(fSlo[cpd][idx])

    # find the (pctTot)% value the way we do it in lat2
    xSloData, hSloData = wl.GetHisto(fSlo[cpd], fShiftLo, fShiftHi, fpbShift)
    max, avg, std, pct, wid = wl.getHistInfo(xSloData,hSloData)
    if pctTot==90:
        fsShiftCutData = pct[2]
    elif pctTot==95:
        fsShiftCutData = pct[3]
    else:
        print("WTF are you doing??")
        exit()

    # plot the slices, (pctTot)% dots, and centroid dots
    for i in range(len(xVal)):
        xV = xVal[i] * np.ones(len(yVal[i]))
        p3.plot(xV, yVal[i], ".b", ms=3.)

        xgTot, effTot = getTot(yVal[i], cutTot=fsShiftCutData)
        p3.plot(xVal[i], xgTot, '.r', ms=10)
        # effM2S238.append([xVal[i],effTot])

        xT, yT = wl.GetHisto(yVal[i], fShiftLo, fShiftHi, fpbShift)
        fCtr = xT[np.argmax(yT)]
        p3.plot(xVal[i], fCtr, '.c', ms=10)

    # try histogram pass/fail method from LAT2
    hitPass, hitFail = [], []
    for i in range(len(fSlo[cpd])):
        if fSlo[cpd][i] <= fsShiftCutData: hitPass.append(hitE[cpd][i])
        else: hitFail.append(hitE[cpd][i])
    hitPass, hitFail = np.asarray(hitPass), np.asarray(hitFail)
    xELow, hPass = wl.GetHisto(hitPass, xPassLo, xPassHi, xpbPass, shift=False)
    xELow, hFail = wl.GetHisto(hitFail, xPassLo, xPassHi, xpbPass, shift=False)
    hTot = np.add(hPass, hFail)
    for i in range(len(xELow)):
        if hTot[i]==0: continue
        effTot = hPass[i]/hTot[i]
        xTot = xELow[i] - xpbPass/2
        effM2S238.append([xTot, effTot])

    p3.plot(np.nan, np.nan, ".b", ms=10, label="m2s238 Events, C1P6D3")
    p3.axhline(fsShiftCutData, c='orange', label="%d%% Overall Cut: %.1f" % (pctTot, fsShiftCutData))

    # special -- load slowness fraction
    fS = np.load('../data/efficiency-corr2.npz')
    simHists = fS['arr_3'].item()
    xLo, xHi, xpb = 0, 50, 1
    hTotSim, hSurfSim, xTotSim = simHists[cpd]
    hFracSim = np.divide(hSurfSim, hTotSim, dtype=float)

    for i in range(len(xELow)):
        if hTot[i]==0: continue
        effTotCorr = (hPass[i]/(1-hFracSim[i]))/hTot[i]
        xTot = xELow[i] - xpbPass/2
        effM2S238Corr.append([xTot, effTotCorr])

    # thesis plot, don't delete
    # plt.close()
    # plt.plot(xELow, hTot, ls='steps', c='k', label='total, m2s238 C1P6D3')
    # plt.plot(xELow, hPass, ls='steps', c='b', label='pass')
    # plt.plot(xELow, hFail, ls='steps', c='r', label='fail')
    # plt.plot(xTotSim, hFracSim * hTot, ls='steps', c='m', label='sim, slow fraction of total')
    # plt.plot(xTotSim, hPass/(1-hFracSim), ls='steps', c='g', label='pass, sim-corrected')
    # plt.legend(loc=1, bbox_to_anchor=(0., 0.5, 1, 0.2))
    # plt.xlabel("Energy (keV)", ha='right', x=1)
    # plt.ylabel("Counts", ha='right', y=1)
    # plt.tight_layout()
    # # plt.show()
    # plt.savefig("../plots/lat-sloFrac-%s.pdf" % cpd)
    # return

    p1.plot(np.nan, np.nan, ".r", label="%d%% of slice" % pctTot)
    p2.plot(np.nan, np.nan, ".r", label="%d%% of slice" % pctTot)
    p3.plot(np.nan, np.nan, ".r", label="%d%% of slice" % pctTot)
    p1.plot(np.nan, np.nan, ".c", label="Max of slice")
    p2.plot(np.nan, np.nan, ".c", label="Max of slice")
    p3.plot(np.nan, np.nan, ".c", label="Max of slice")

    p1.legend(loc=1, fontsize=14)
    p2.legend(loc=1, fontsize=14)
    p3.legend(loc=1, fontsize=14)

    p1.set_xlim(xLo, xHi)
    p1.set_ylim(yLo, yHi)
    p2.set_xlim(xLo, xHi)
    p2.set_ylim(yLo, yHi)
    p3.set_xlim(xLo, xHi)
    p3.set_ylim(yLo, 700)

    p1.set_ylabel("fitSlo (shifted)", ha='right', y=1)
    p3.set_xlabel("Energy (keV)", ha='right', x=1)

    p1.set_yticks([-80, 0, 100, 200])
    p2.set_yticks([-80, 0, 100, 200])
    p3.set_yticks([-80, 0, 200, 400, 600])

    plt.tight_layout()
    fig.subplots_adjust(hspace=0.05)
    plt.setp(p1.get_xticklabels(), visible=False)
    plt.setp(p2.get_xticklabels(), visible=False)

    # plt.show()
    plt.savefig('../plots/lat-acceptance-study.png', dpi=300)
    return


    # ========= acceptance / efficiency plot =========
    # return
    plt.close()
    fig = plt.figure()

    accExtPulser = np.asarray(accExtPulser)
    accSimWF = np.asarray(accSimWF)
    effM2S238 = np.asarray(effM2S238)
    effM2S238Corr = np.asarray(effM2S238Corr)
    # plt.plot(accExtPulser[:,0], accExtPulser[:,1], '-r', ls='steps', ms=10, label="Ext. Pulser Acc.")
    # plt.plot(accSimWF[:,0], accSimWF[:,1], '-g', ls='steps', ms=10, label="Sim WF Acc.")
    plt.plot(effM2S238[:,0], effM2S238[:,1], '-c', ls='steps',ms=10, label="m2s238, C1P6D3")
    plt.plot(effM2S238Corr[:,0], effM2S238Corr[:,1], '-m', ls='steps', ms=10, label="Sim-corrected")

    # now you're back to the figure
    plt.xlim(xLo, xHi)
    plt.xlabel("Energy (keV)", ha='right', x=1)
    plt.ylabel("Acceptance", ha='right', y=1)
    plt.legend(loc=4)
    plt.tight_layout()
    # plt.show()
    plt.savefig("../plots/lat-acceptance-vals.pdf")
    return


    # ===== plot the fitSlo distributions for all events with their (pctTot)% cuts =====
    plt.close()

    plt.semilogy(xSloExt, hSloExt/np.amax(hSloExt), c='r', lw=2, ls='steps', label='Ext. Pulser')
    plt.axvline(fsShiftCutExt, c='m', lw=2, label="Ext %d%%: %.1f" % (pctTot, fsShiftCutExt))

    plt.semilogy(xSloSim, hSloSim/np.amax(hSloSim), c='g', lw=2, ls='steps', alpha=0.8, label='Simulated WFs')
    plt.axvline(fsShiftCutSim, c='k', lw=2, label="Simulated %d%%: %.1f" % (pctTot, fsShiftCutSim))

    plt.semilogy(xSloData, hSloData/np.amax(hSloData), c='b', lw=2, ls='steps', alpha=0.7, label='m2s238 Data')
    plt.axvline(fsShiftCutData, c='c', lw=2, label="m2s238 %d%%: %.1f" % (pctTot, fsShiftCutData))

    plt.legend()

    plt.xlabel("fitSlo (shifted)", ha='right', x=1)
    plt.ylabel("Cts (Normalized)", ha='right', y=1)

    plt.tight_layout()
    # plt.show()
    plt.savefig("../plots/lat-acceptance-slo.pdf")


    # === fit the "new (pctTot)%" efficiency to weibull and recreate the "standard" efficiency plot ===
    plt.close()

    # plt.plot(effM2S238[:,0], effM2S238[:,1], '-c', ls='steps',ms=10, label="m2s238 Efficiency")

    # this is the new, corrected efficiency
    idxP = np.where(hPass > 0)
    sloEff = hPass[idxP] / hTot[idxP]
    xEff = xELow[idxP]
    ci_low, ci_upp = proportion.proportion_confint(hPass[idxP], hTot[idxP], alpha=0.1, method='beta')
    idxE = np.where((xEff < xPassHi) & (xEff > 1))
    xEff, sloEff, ci_low, ci_upp = xEff[idxE], sloEff[idxE], ci_low[idxE], ci_upp[idxE]
    xEff -= xpbPass/2.
    eFitHi = 30
    fitBnd = ((1,-20,0,0.5),(np.inf,np.inf,np.inf, 0.99))
    idxF = np.where(xEff <= eFitHi)
    popt, pcov = curve_fit(wl.weibull, xEff[idxF], sloEff[idxF], bounds=fitBnd)
    perr = np.sqrt(np.diag(pcov))
    c, loc, sc, amp = popt
    cE, locE, scE, ampE = perr
    eff1 = wl.weibull(1.,*popt)
    plt.plot(xEff, sloEff, '.b', ms=10., label='C%sP%sD%s' % (cpd[0],cpd[1],cpd[2]))
    plt.errorbar(xEff, sloEff, yerr=[sloEff - ci_low, ci_upp - sloEff], color='k', linewidth=0.8, fmt='none')
    xFunc = np.arange(xPassLo, xPassHi, 0.1)
    plt.plot(xFunc, wl.weibull(xFunc, *popt), 'm-', label=r'Weibull CDF')

    # NOTE: this is wrong, just plotting it to see the diff between old and new
    # plot the old "official" m2s238 efficiency on top -- this is wrong b/c of the wrong fitSlo 90% limit
    # this is a condensed version of lat-plots.py::fitSlo_det_efficiencies
    # cpd = '163'
    # f = np.load('../data/lat2-eff-data.npz')
    # effData = f['arr_0'].item()
    # xPassLo, xPassHi, xpbPass = 0, 50, 1
    # xPassLo, xPassHi = xLo, xHi
    # xEff, sloEff, ci_low, ci_upp = effData[cpd][0], effData[cpd][1], effData[cpd][2], effData[cpd][3]
    # hPass, hFail, hTot, xELow = effData[cpd][4], effData[cpd][5], effData[cpd][6], effData[cpd][7]
    # idxP = np.where(hPass > 0)
    # sloEff = hPass[idxP] / hTot[idxP]
    # xEff = xELow[idxP]
    # ci_low, ci_upp = proportion.proportion_confint(hPass[idxP], hTot[idxP], alpha=0.1, method='beta')
    # idxE = np.where((xEff < xPassHi) & (xEff > 1))
    # xEff, sloEff, ci_low, ci_upp = xEff[idxE], sloEff[idxE], ci_low[idxE], ci_upp[idxE]
    # xEff -= xpbPass/2.
    # b3 = ((0,-15,0,0),(np.inf,np.inf,np.inf,1.)) # this one is working best
    # popt, pcov = curve_fit(wl.weibull, xEff, sloEff, bounds=b3)
    # perr = np.sqrt(np.diag(pcov))
    # c, loc, sc, amp = popt
    # cE, locE, scE, ampE = perr
    # eff1 = wl.weibull(1.,*popt)
    # plt.plot(xEff, sloEff, '.b', ms=10., label='C%sP%sD%s' % (cpd[0],cpd[1],cpd[2]))
    # plt.errorbar(xEff, sloEff, yerr=[sloEff - ci_low, ci_upp - sloEff], color='k', linewidth=0.8, fmt='none')
    # xFunc = np.arange(xPassLo, xPassHi, 0.1)
    # plt.plot(xFunc, wl.weibull(xFunc, *popt), 'r-', label=r'Weibull CDF')


    plt.axvline(1.,color='b', lw=1., label='1.0 keV efficiency: %.2f' % wl.weibull(1.,*popt))
    plt.xlim(xLo, xHi)
    plt.xlabel("Energy (keV)", ha='right', x=1)
    plt.ylabel("Efficiency", ha='right', y=1)
    plt.legend(loc=4)
    plt.tight_layout()
    # plt.show()
    plt.savefig("../plots/lat-efficiency-weibull-%s.pdf" % cpd)


# def sim_compton_edge():



if __name__=="__main__":
    main()