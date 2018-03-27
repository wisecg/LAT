#!/usr/bin/env python
import sys, os, imp, pywt
import scipy.optimize as op
import scipy.special as sp
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import gridspec
import seaborn as sns
sns.set(style='darkgrid')

"""
    Simulated WF study
    -- Generate simulated WFs from SigGen from within the bulk and near the surface of the detector
    -- Add simulated WFs (varying amplitudes) to force trigger baselines
    -- Process with simplified version of LAT (only WF fitting and riseNoise for now)
    -- Study distributions
"""


# load LAT libraries
import DataSetInfo as ds
wl = imp.load_source('waveLibs',os.environ['LATDIR']+'/waveLibs.py')

fig = plt.figure(figsize=(12,7), facecolor='w')
p0 = plt.subplot2grid((6,10), (0,0), colspan=10, rowspan=3) # waveform fit
p1 = plt.subplot2grid((6,10), (3,0), colspan=10, rowspan=1) # residual
p2 = plt.subplot2grid((6,10), (4,0), colspan=2, rowspan=2) # traces
p3 = plt.subplot2grid((6,10), (4,2), colspan=2, rowspan=2)
p4 = plt.subplot2grid((6,10), (4,4), colspan=2, rowspan=2)
p5 = plt.subplot2grid((6,10), (4,6), colspan=2, rowspan=2)
p6 = plt.subplot2grid((6,10), (4,8), colspan=2, rowspan=2)


def main():
    # procSim(5000)
    compareSim()


def compareSim():
    plt.close(fig)

    inDir = os.environ['LATDIR']+'/data'
    # dfSig2_low = pd.read_hdf('{}/SimPSA_sig2_log.h5'.format(inDir))
    dfSig2_low = pd.read_hdf('{}/SimPSA_sig2_low.h5'.format(inDir))
    dfSig2 = pd.read_hdf('{}/SimPSA_sig2.h5'.format(inDir))
    dfSig2['Distance'] = dfSig2['R']*dfSig2['R'] + dfSig2['Z']+dfSig2['Z']
    dfSig2_low['Distance'] = dfSig2_low['R']*dfSig2_low['R'] + dfSig2_low['Z']+dfSig2_low['Z']
    dfSig2_low = dfSig2_low.append(dfSig2, ignore_index=True)
    # dfCut = dfSig2_low.loc[dfSig2_low['Amp'] < 8]
    dfCut = dfSig2_low.loc[(dfSig2_low['Amp'] < 14) & (dfSig2_low['fitSlo'] > 2)]

    # g1 = sns.FacetGrid(dfCut, hue='Distance', size=7)
    # g1 = g1.map(plt.scatter, 'Amp', 'fitSlo').add_legend()
    # g1 = sns.FacetGrid(dfSig2.loc[dfSig2['Amp'] == 6], hue='Dist', size=7)
    # g1 = g1.map(plt.hist, 'fitSlo', bins=np.linspace(0,1000,100), alpha=0.5).add_legend()
    # g1 = g1.map(plt.hist, 'fitSlo', bins=np.linspace(0,2500,250), alpha=0.5).add_legend()
    # g1.set(yscale='log')

    g1 = sns.lmplot(x='Amp', y='fitSlo', data=dfCut, hue='Distance', fit_reg=False, x_jitter=0.3, size=7, scatter_kws={'alpha':0.3}, legend_out=False)
    g1.set(yscale='log')
    plt.subplots_adjust(top=0.95)
    g1.fig.suptitle('fitSlo vs Amplitude (Simulated Waveforms)')
    g1.set_axis_labels('Amplitude (ADC)', 'fitSlo')


    # g2 = sns.FacetGrid(dfCut, hue='Amp', size=7)
    # g2 = g2.map(plt.hist, 'fitSlo', bins=np.linspace(0,2500,250), alpha=0.50).add_legend()
    # g2.set(yscale='log')

    # plt.tight_layout()
    plt.show()
    # g1.savefig('/Users/brianzhu/macros/code/LAT/plots/Systematics/fitSlo/Sim_fitSlo_vs_Amplitude_Scatter.png')
    # g2.savefig('/Users/brianzhu/macros/code/LAT/plots/Systematics/fitSlo/Sim_fitSlo_AmpComparison_log.png')

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


def procSim(nMax = 5000):
    """
        Get sum and hit spectra w/ threshold cut, without ch. 598 (it's noisy.)
        Here we use "mHT" and "sumET" exclusively.
        Use !EventDC1Bits and trapENFCal > 0.7, all hits.  (Skim file was generated w/ 'dontSkipAnything')
        Save detailed info for events with hits < 2630 keV.
    """
    from ROOT import TFile, TTree, MGTWaveform

    # Use the force acquisition baseline + SigGen waveforms
    inDir = os.environ['LATDIR']+'/data'
    f1 = TFile('{}/waveSkimDS0_run4201.root'.format(inDir))
    skimTree = f1.Get('skimTree')

    dfSigGen = pd.read_hdf('{}/SigGen_WFs_low.h5'.format(inDir))

    # Truncates to 2017 samples, cuz we in 2017
    truncLo, truncHi = 8,6
    dataList = []
    # for idx in range(1):
    print('Total Entries: ', skimTree.GetEntries())
    # for idx in range(skimTree.GetEntries()):
    # fig1, ax1 = plt.subplots(figsize=(10,7))
    for idx in range(min(nMax,skimTree.GetEntries())):
        # if idx != 417: continue
        if idx%100 == 0:
            print("Processed Entry ", idx)
        skimTree.GetEntry(idx)
        if skimTree.mH > 1: continue # Only care about mH1 for noise
        wf = MGTWaveform()
        wf = skimTree.MGTWaveforms.at(0)
        signal = wl.processWaveform(wf,truncLo,truncHi)
        data = signal.GetWaveRaw()
        data_blSub = signal.GetWaveBLSub()
        dataTS = signal.GetTS()
        dataBL,dataNoise = signal.GetBaseNoise()
        # Generate one fake WF at each amplitude/r/z combination (30 total)
        for index, row in dfSigGen.iterrows():
            dataMap = {}
            sig = addSigToBaseline(data, row['waveform'])
            sigBL = addSigToBaseline(data_blSub, row['waveform'])
            # Dummy variables
            # fitAmp, TSMax = row['amp'], 1200*10
            dataMap = procEvent(sig, dataTS, sigBL, dataBL, '{}_A{}R{}Z{}'.format(idx, row['amp'], row['r'], row['z']))
            dataMap['Amp'] = row['amp']
            dataMap['R'] = row['r']
            dataMap['Z'] = row['z']
            # dataMap['waveform'] = sigBL.tolist()
            dataList.append(dataMap)

    df = pd.DataFrame.from_dict(dataList)
    print(df.head())
    df.to_hdf('{}/SimPSA_sig2_log.h5'.format(inDir), 'skimTree')


def addSigToBaseline(data_blSub, sigWF):
    padSig = np.pad(sigWF, (len(data_blSub)-len(sigWF),0), 'constant', constant_values=(0,0))
    return data_blSub + padSig


def procEvent(data, dataTS, data_blSub, dataBL, saveStr=''):
    """
        Basically a shortened version of LAT with only riseNoise and fitSlo operating on arrays
    """
    dataMap = {}

    # Calculate trapezoid to get trapENMSample and trapENM
    eTrap = wl.trapFilter(data_blSub, 400, 200, 7200.)
    dataENM = np.amax(eTrap)
    dataTSMax = np.argmax(eTrap)*10 + 6000
    # print(dataENM, dataTSMax, np.argmax(eTrap))
    # wavelet packet transform
    wp = pywt.WaveletPacket(data_blSub, 'db2', 'symmetric', maxlevel=4)
    nodes = wp.get_level(4, order='freq')
    wpCoeff = np.array([n.data for n in nodes],'d')
    wpCoeff = abs(wpCoeff)

    # wavelet parameters
    # First get length of wavelet on the time axis, the scale axis will always be the same
    # due to the number of levels in the wavelet
    wpLength = len(wpCoeff[1,:])
    waveS1 = np.sum(wpCoeff[0:1,1:wpLength//4+1]) # python3 : floor division (//) returns an int
    waveS2 = np.sum(wpCoeff[0:1,wpLength//4+1:wpLength//2+1])
    waveS3 = np.sum(wpCoeff[0:1,wpLength//2+1:3*wpLength//4+1])
    waveS4 = np.sum(wpCoeff[0:1,3*wpLength//4+1:-1])
    waveS5 = np.sum(wpCoeff[2:-1,1:-1])
    S6 = np.sum(wpCoeff[2:9,1:wpLength//4+1])
    S7 = np.sum(wpCoeff[2:9,wpLength//4+1:wpLength//2+1])
    S8 = np.sum(wpCoeff[2:9,wpLength//2+1:3*wpLength//4+1])
    S9 = np.sum(wpCoeff[2:9,3*wpLength//4+1:-1])
    S10 = np.sum(wpCoeff[9:,1:wpLength//4+1])
    S11 = np.sum(wpCoeff[9:,wpLength//4+1:wpLength//2+1])
    S12 = np.sum(wpCoeff[9:,wpLength//2+1:3*wpLength//4+1])
    S13 = np.sum(wpCoeff[9:,3*wpLength//4+1:-1])
    sumList = [S6, S7, S8, S9, S10, S11, S12, S13]
    bcMax = np.max(sumList)
    bcMin = 1. if np.min(sumList) < 1 else np.min(sumList)

    new_wp = pywt.WaveletPacket(data=None, wavelet='db2', mode='symmetric')
    new_wp['aaa'] = wp['aaa'].data
    data_wlDenoised = new_wp.reconstruct(update=False)
    diff = len(data_wlDenoised) - len(data_blSub)
    if diff > 0: data_wlDenoised = data_wlDenoised[diff:]

    amp, mu, sig, tau, bl = dataENM, dataTSMax, 600., -72000., dataBL
    # print('Guess: fitMu {}, fitAmp {}, fitSlo {}, fitTau {}, fitBL {}'.format(mu, amp, sig, tau, bl))
    floats = np.asarray([amp, mu, sig, tau, bl])
    temp = xgModelWF(dataTS, floats)
    MakeTracesGlobal()

    denoisedNoise,_,_ = wl.baselineParameters(data_wlDenoised)
    datas = [dataTS, data_wlDenoised + dataBL, denoisedNoise] # fit wavelet-denoised data w/ Bl added back in
    bnd = ((None,None),(None,None),(2.,None),(-72001.,-71999.),(None,None)) # gets caught much less often.

    # result = op.minimize(lnLike, floats, args=datas, method="L-BFGS-B", options={'disp':True},  bounds=bnd)
    result = op.minimize(lnLike, floats, args=datas, method="L-BFGS-B", bounds=bnd)

    amp, mu, sig, tau, bl = result["x"]
    floats = np.asarray([amp, mu, sig, tau, bl])
    fit = xgModelWF(dataTS, floats)
    fitMu, fitAmp, fitSlo, fitTau, fitBL = mu, amp, sig, tau, bl


    if sig <= 2:
        print('Old: fitMu {}, fitAmp {}, fitSlo {}, sig {}, fitTau {}, fitBL {}'.format(fitMu, fitAmp, fitSlo, sig, fitTau, fitBL))

        amp, mu, sig, tau, bl = dataENM, dataTSMax, np.log(600.), -72000., dataBL
        floats2 = np.asarray([amp, mu, np.log(sig), tau, bl])
        bnd2 = ((None,None),(None,None),(None,None),(-72001.,-71999.),(None,None)) # gets caught much less often.
        result2 = op.minimize(lnLike, floats2, args=datas, method='L-BFGS-B', bounds=bnd2)
        amp, mu, sig, tau, bl = result2["x"]
        floats2 = np.asarray([amp, mu, np.exp(sig), tau, bl])
        fit = xgModelWF(dataTS, floats2)
        fitMu, fitAmp, fitSlo, fitTau, fitBL = mu, amp, np.exp(sig), tau, bl

        p0.cla()
        p0.plot(dataTS,data,color='blue',label='data')
        p0.plot(dataTS,temp,color='orange',label='xgauss guess')
        p0.plot(dataTS,fit,color='red',label='xgauss fit')
        p0.plot(dataTS,datas[1], color='green',label='denoised wf')
        p0.legend(loc='best')
        p1.cla()
        p1.plot(dataTS,data-fit,color='blue',label='residual')
        p1.legend(loc='best')
        p2.cla()
        p2.plot(ampTr[1:],label='amp',color='red')
        p2.legend(loc='best')
        p3.cla()
        p3.plot(muTr[1:],label='mu',color='green')
        p3.legend(loc='best')
        p4.cla()
        p4.plot(sigTr[1:],label='sig',color='blue')
        p4.legend(loc='best')
        p5.cla()
        p5.plot(tauTr[1:],label='tau',color='black')
        p5.legend(loc='best')
        p6.cla()
        p6.plot(blTr[1:],label='bl',color='magenta')
        p6.legend(loc='best')
        plt.tight_layout()
        fig.savefig('/Users/brianzhu/macros/code/LAT/plots/Systematics/fitSlo/ADC2/WF_{}.png'.format(saveStr))

        print('Fixed: fitMu {}, fitAmp {}, fitSlo {}, sig {}, fitTau {}, fitBL {}'.format(fitMu, fitAmp, fitSlo, sig, fitTau, fitBL))
    # print('fitMu {}, fitAmp {}, fitSlo {}, fitTau {}, fitBL {}'.format(fitMu, fitAmp, fitSlo, fitTau, fitBL))

    # find the window of rising edge
    fit_blSub = fit - bl
    fitMaxTime = dataTS[np.argmax(fit_blSub)]
    fitStartTime = dataTS[0]
    idx = np.where(fit_blSub < 0.1)
    if len(dataTS[idx] > 0): fitStartTime = dataTS[idx][-1]
    fitRiseTime50 = (fitMaxTime + fitStartTime)/2.

    # bcMin is 32 samples long in the x-direction.
    # if we make the window half as wide, it'll have the same # of coeff's as bcMin.
    # this is still 'cheating' since we're not summing over the same rows.
    numXRows = wpCoeff.shape[1]
    wpCtrRise = int((fitRiseTime50 - dataTS[0]) / (dataTS[-1] - dataTS[0]) * numXRows)
    wpLoRise = wpCtrRise - 8
    if wpLoRise < 0: wpLoRise = 0
    wpHiRise = wpCtrRise + 8
    if wpHiRise > numXRows: wpHiRise = numXRows
    riseNoise = np.sum(wpCoeff[2:-1,wpLoRise:wpHiRise]) / bcMin

    dataMap['riseNoise'] = riseNoise
    dataMap['fitSlo'] = fitSlo

    return dataMap


def evalGaus(x,mu,sig):
    return np.exp(-((x-mu)**2./2./sig**2.))


def evalXGaus(x,mu,sig,tau):
    """ Ported from GAT/BaseClasses/GATPeakShapeUtil.cc
        Negative tau: Regular WF, high tail
        Positive tau: Backwards WF, low tail
    """
    tmp = (x-mu + sig**2./2./tau)/tau
    # tmp = (x-mu + np.log(sig)**2./2./tau)/tau

    # np.exp of this is 1.7964120280206387e+308, the largest python float value: sys.float_info.max
    fLimit = 709.782

    if all(tmp < fLimit):
        return np.exp(tmp)/2./np.fabs(tau) * sp.erfc((tau*(x-mu)/sig + sig)/np.sqrt(2.)/np.fabs(tau))
        # return np.exp(tmp)/2./np.fabs(tau) * sp.erfc((tau*(x-mu)/np.log(sig) + np.log(sig))/np.sqrt(2.)/np.fabs(tau))
    else:
        # print("Exceeded limit ...")
        # Use an approx. derived from the asymptotic expansion for erfc, listed on wikipedia.
        den = 1./(sig + tau*(x-mu)/sig)
        return sig * evalGaus(x,mu,sig) * den * (1.-tau**2. * den**2.)

        # den = 1./(np.log(sig) + tau*(x-mu)/np.log(sig))
        # return np.log(sig) * evalGaus(x,mu,np.log(sig)) * den * (1.-tau**2. * den**2.)

def xgModelWF(dataTS, floats):
    """ Make a model waveform: Take a timestamp vector, generate an
        xGauss model, normalize to 1, then scale its max value to amp.
    """
    amp, mu, sig, tau, bl = floats
    model = evalXGaus(dataTS,mu,sig,tau)

    if np.isnan(model).any():
        # print("d'oh!",model)
        return(np.zeros(len(model)))
    if np.sum(model)==0:
        # print("dooh!",model)
        return(np.zeros(len(model)))

    # pin max value of function to amp
    model = model * 1./np.sum(model)
    xMax = np.argmax(model)
    model = model * (amp / model[xMax])

    # float the baseline
    model = model + bl
    return model


def MakeTracesGlobal():
    """ This is so 'lnLike' can write to the trace arrays. Has to remain in this file to work. """
    tmp1, tmp2, tmp3, tmp4, tmp5 = [], [], [], [], []
    global ampTr, muTr, sigTr, tauTr, blTr
    ampTr, muTr, sigTr, tauTr, blTr = tmp1, tmp2, tmp3, tmp4, tmp5


def lnLike(floats, *datas):
    """ log-likelihood function: L(A,mu,sig,tau)
    To make this work with op.minimize, 'datas' is passed in as a tuple (the asterisk),
    where the original list is the 1st element.
    """
    amp, mu, sig, tau, bl = floats
    tau = -72000. # manually fix tau

    global ampTr, muTr, sigTr, tauTr, blTr
    ampTr.append(amp)
    muTr.append(mu)
    sigTr.append(sig)
    tauTr.append(tau)
    blTr.append(bl)

    dataTS, data, dataNoise = datas[0][0], datas[0][1], datas[0][2]

    model = xgModelWF(dataTS, [amp, mu, sig, tau, bl])
    lnLike = 0.5 * np.sum ( np.power((data-model)/dataNoise, 2) - np.log( 1 / np.power(dataNoise,2) ) )
    return lnLike


if __name__=="__main__":
    main()
