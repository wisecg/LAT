#!/usr/bin/env python
import sys, time, imp, pywt
import numpy as np
import scipy.special as sp
import scipy.optimize as op
from scipy.signal import butter, lfilter, filtfilt, kaiser
from scipy.ndimage.filters import gaussian_filter
import matplotlib.pyplot as plt
from ROOT import TFile,TTree,TEntryList,MGTWaveform,gDirectory
wl = imp.load_source('waveLibs', '../waveLibs.py')
limit = sys.float_info.max # equivalent to std::numeric_limits::max() in C++

def main(argv):
    """ Investigate the time-domain matched filter, w/ and w/o Kaiser window. """

    scanSpeed = 0.2
    opt1, opt2 = "", ""
    intMode, batMode = False, False
    if (len(argv) >= 1): opt1 = argv[0]
    if (len(argv) >= 2): opt2 = argv[1]
    if "-i" in (opt1, opt2):
        intMode = True
        print "Interactive mode selected."
    if "-b" in (opt1, opt2):
        batMode = True
        print "Batch mode selected."

    inputFile = TFile("../waveSkimDS5_run21975.root")
    waveTree = inputFile.Get("skimTree")
    print "Found",waveTree.GetEntries(),"input entries."

    theCut = inputFile.Get("theCut").GetTitle()
    # theCut += " && Entry$ < 100"
    # theCut += " && trapENFCal > 0.8 && trapENFCal < 2"
    theCut += " && trapENFCal < 6"

    waveTree.Draw(">>elist", theCut, "entrylist")
    elist = gDirectory.Get("elist")
    waveTree.SetEntryList(elist)
    nList = elist.GetN()

    print "Using cut:\n",theCut,"\nFound",nList,"entries passing cuts."

    fig = plt.figure(figsize=(13,7), facecolor='w')
    p0 = plt.subplot(111)
    # p1 = plt.subplot(212)
    # p0 = plt.subplot2grid((6,8), (0,0), colspan=8, rowspan=3) # wf plot
    # p1 = plt.subplot2grid((6,8), (3,0), colspan=8, rowspan=3) # mf plot
    if not batMode: plt.show(block=False)

    # Loop over events
    iList = -1
    while(True):
        iList += 1
        if intMode==True and iList !=0:
            value = raw_input()
            if value=='q': break
            if value=='p': iList -= 2  # previous
            if (value.isdigit()): iList = int(value) # go to entry
        if iList >= elist.GetN(): break

        entry = waveTree.GetEntryNumber(iList);
        waveTree.LoadTree(entry)
        waveTree.GetEntry(entry)
        nChans = waveTree.channel.size()
        nWFs = waveTree.MGTWaveforms.size()
        if (nWFs==0):
            print "Error - nWFs:",nWFs,"nChans",nChans
            continue
        numPass = waveTree.Draw("channel",theCut,"GOFF",1,iList)
        chans = waveTree.GetV1()
        chanList = list(set(int(chans[n]) for n in xrange(numPass)))

        # Loop over hits passing cuts
        hitList = (iH for iH in xrange(nChans) if waveTree.channel.at(iH) in chanList)  # a 'generator expression'
        for iH in hitList:
            run = waveTree.run
            chan = waveTree.channel.at(iH)
            dataENF = waveTree.trapENFCal.at(iH)
            wf = waveTree.MGTWaveforms.at(iH)
            signal = wl.processWaveform(wf)
            data = signal.GetWaveBLSub()
            dataTS = signal.GetTS()
            _,dataNoise = signal.GetBaseNoise()
            dataTSMax = waveTree.trapENMSample.at(iH)*10. - 4000
            dataENM = waveTree.trapENM.at(iH)

            # wavelet packet denoised waveform
            wp = pywt.WaveletPacket(data, 'db2', 'symmetric', maxlevel=4)
            nodes = wp.get_level(4, order='freq')
            waveletYTrans = np.array([n.data for n in nodes],'d')
            waveletYTrans = abs(waveletYTrans)

            # reconstruct waveform w/ only lowest frequency.
            new_wp = pywt.WaveletPacket(data=None, wavelet='db2', mode='symmetric')
            new_wp['aaa'] = wp['aaa'].data
            data_wlDenoised = new_wp.reconstruct(update=False)
            # resize in a smart way
            diff = len(data_wlDenoised) - len(data)
            if diff > 0: data_wlDenoised = data_wlDenoised[diff:]

            # get the noise of the denoised wf
            denoisedNoise,_,_ = wl.baselineParameters(data_wlDenoised)


            # ================= run waveform fitter =================

            amp, mu, sig, tau = dataENM, dataTSMax, 600., -72000.
            floats = np.asarray([amp, mu, sig, tau])
            guess = xgModelWF(dataTS,floats)
            MakeTracesGlobal()

            # datas = [dataTS, data, dataNoise] # fit data
            datas = [dataTS, data_wlDenoised, denoisedNoise] # fit wavelet denoised

            start = time.clock()

            # Powell or Nelder-Mead
            # result = op.minimize(lnLike, floats, args=datas, method="Powell")

            # L-BGFS-B
            bnd = ((None,None),(None,None),(None,None),(None,None)) # A,mu,sig,tau
            opts = {'disp': None,   # None, True to print convergence messages
                    'maxls': 100,   # 20, max line search steps
                    'iprint': -1,   # -1
                    'gtol': 1e-08,  # 1e-05
                    'eps': 1e-08,   # 1e-08
                    'maxiter': 15000,   # 15000
                    'ftol': 2.220446049250313e-09,
                    'maxcor': 10,   # 10
                    'maxfun': 15000}    # 15000

            # numerical gradient - seems trustworthy
            result = op.minimize(lnLike, floats, args=datas, method="L-BFGS-B", options=opts, bounds=bnd)#, jac=lnLikeGrad)

            # analytical gradient - BETA, Don't Use!
            # result = op.minimize(lnLike, floats, args=datas, method="L-BFGS-B", options=opts, bounds=None, jac=lnLikeGrad)

            fitTime = time.clock() - start

            if not result["success"]:
                print "fit 'fail': ", result["message"]
                # errorCode[0] = 1

            amp, mu, sig, tau = result["x"]
            floats = [amp, mu, sig, tau]
            fit = xgModelWF(dataTS, floats)

            # log-likelihood of this fit
            fitLL = result["fun"]

            # chi-square of this fit
            # Textbook is (observed - expected)^2 / expected,
            # but we'll follow MGWFCalculateChiSquare.cc and do (observed - expected)^2 / NDF.
            # NOTE: we're doing the chi2 against the DATA, though the FIT is to the DENOISED DATA.
            fitChi2 = np.sum(np.square(data-fit)) / len(data)


            # =======================================================

            # time-domain matched filter

            # make a longer best-fit waveform s/t it can be shifted L/R.
            matchTS = np.append(dataTS, np.arange(dataTS[-1], dataTS[-1] + 20000, 10)) # add 2000 samples
            match = xgModelWF(matchTS, [amp, mu+10000., sig, tau]) # shift mu accordingly
            match = match[::-1]

            # line up the max of the 'match' (flipped wf) with the max of the best-fit wf
            fitMaxTime = dataTS[np.argmax(fit)]
            matchMaxTime = matchTS[np.argmax(match)]
            matchTS = matchTS + (fitMaxTime - matchMaxTime)

            # windowing function
            window = kaiser(2000, beta=8)
            wMax = np.argmax(window)

            # line up the window with the max of the match wf.
            windowTS = matchTS
            window = np.pad(window, ((len(matchTS)-len(window))/2, (len(matchTS)-len(window))/2), 'constant')
            windowTS = windowTS + (fitMaxTime - matchMaxTime)

            # resize match wf and the window function to have same # samples as dataTS.

            if matchTS[0] <= dataTS[0] and matchTS[-1] >= dataTS[-1]:
                idx = np.where((matchTS >= dataTS[0]) & (matchTS <= dataTS[-1]-5))
                match, matchTS, window = match[idx], matchTS[idx], window[idx]

            elif matchTS[-1] < dataTS[-1]: # too early
                idx = np.where(matchTS >= dataTS[0])
                tmpData, tmpWindow = match[idx], window[idx]
                match = np.pad(tmpData, (0,len(data)-len(tmpData)), 'constant')
                window = np.pad(tmpWindow, (0,len(data)-len(tmpWindow)), 'constant')
                matchTS, windowTS = dataTS, dataTS

            elif matchTS[0] > dataTS[0]: # too late
                idx = np.where(matchTS <= dataTS[-1])
                tmpData, tmpWindow = match[idx], window[idx]
                match = np.pad(tmpData, (len(data)-len(tmpData),0), 'constant')
                window = np.pad(tmpWindow, (len(data)-len(tmpWindow),0), 'constant')
                matchTS, windowTS = dataTS, dataTS

            # get lowpass and bandpass data
            B3, A3 = butter(2,1e6/(1e8/2), btype='lowpass')
            data_lPass = lfilter(B3, A3, data)

            B1,A1 = butter(2, [1e5/(1e8/2),1e6/(1e8/2)], btype='bandpass')
            data_bPass = lfilter(B1, A1, data)

            # make sure match is same length as data, then compute convolution and parameters
            matchMax, matchWidth, matchTime = -888, -888, -888
            smoothMF, windMF = np.zeros(len(dataTS)), np.zeros(len(dataTS))
            if len(match)!=len(data):
                print "array mismatch: len match %i  len data %i " % (len(match),len(data))
                # errorCode[1] = 1
            else:
                # compute match filter parameters
                smoothMF = gaussian_filter(match * data_bPass,sigma=5.)
                # smoothMF = gaussian_filter(match * data, sigma=5.)
                matchMax = np.amax(smoothMF)
                matchTime = matchTS[ np.argmax(smoothMF) ]
                idx = np.where(smoothMF > matchMax/2.)
                if len(matchTS[idx]>1): matchWidth = matchTS[idx][-1] - matchTS[idx][0]


            # matchMaxKW, matchWidthKW, matchTimeKW = -888, -888, -888
            # windMF, windData, windMatch = np.zeros(len(dataTS)), np.zeros(len(dataTS)), np.zeros(len(dataTS))
            # if len(window)!=len(data):
            #     print "window mismatch: len window %d  len data %d" % (len(window),len(data))
            #     # errorCode[8]=1
            # else:
            #     windData = window * data_lPass
            #     windMatch = window * match
            #     windMF = gaussian_filter(windData * windMatch, sigma=5.)
            #     matchMaxKW = np.amax(windMF)
            #     matchTimeKW = matchTS[ np.argmax(windMF) ]
            #     idx = np.where(windMF > matchMaxKW/2.)
            #     if len(matchTS[idx]>1): matchWidthKW = matchTS[idx][-1] - matchTS[idx][0]

            print "%d  ch %d  tENF %.2f  match max %.2f  width %.2f  time %.2f " % (iList, chan, dataENF, matchMax, matchWidth, matchTime)

            # =======================================================
            if batMode: continue

            p0.cla()
            p0.plot(dataTS,data,color='blue',label='data',alpha=0.7)
            p0.plot(dataTS,fit,color='red',label='bestfit',linewidth=3)
            p0.plot(dataTS,data_bPass,color='green',label='bpass',linewidth=2)
            p0.axvline(matchTime,color='orange',label='matchTime',linewidth=2)
            p0.plot(matchTS,smoothMF,color='magenta',label='smoothMF',linewidth=3)
            p0.plot(matchTS,match,color='cyan',label='match',linewidth=2)
            p0.set_xlabel('Time (s)')
            p0.set_ylabel('Voltage (arb)')
            p0.legend(loc=4)
            p0.set_title("Run %d  Entry %d  Channel %d  ENFCal %.2f  matchMax %.2f  matchTime %.2f  matchWidth %.2f" % (run,iList,chan,dataENF,matchMax,matchTime,matchWidth))

            # p1.cla()
            # p1.plot(dataTS, data, color='blue',label='data', alpha=0.5)
            # p1.plot(dataTS, window*amp, color='red', label='window*fitAmp')
            # p1.plot(dataTS, windData, color='orange',label='windData')
            # p1.plot(dataTS, windMatch, color='black',label='windMatch')
            # p1.plot(dataTS, windMF, color='green', label='windMF')
            # p1.legend(loc=4)

            plt.tight_layout()
            plt.pause(scanSpeed)


def evalGaus(x,mu,sig):
    return np.exp(-((x-mu)**2./2./sig**2.))

def evalXGaus(x,mu,sig,tau):
    """ Ported from GAT/BaseClasses/GATPeakShapeUtil.cc
        Negative tau: Regular WF, high tail
        Positive tau: Backwards WF, low tail
    """
    tmp = (x-mu + sig**2./2./tau)/tau
    if all(tmp < limit):
        return np.exp(tmp)/2./np.fabs(tau) * sp.erfc((tau*(x-mu)/sig + sig)/np.sqrt(2.)/np.fabs(tau))
    else:
        print "Exceeded limit ..."
        # Here, exp returns NaN (in C++).  So use an approx. derived from the asymptotic expansion for erfc, listed on wikipedia.
        den = 1./(sig + tau*(x-mu)/sig)
        return sig * gaus(x,mu,sig) * den * (1.-tau**2. * den**2.)

def evalXgGrad(x,amp,mu,sig,tau):
    """ Ported from GAT/BaseClasses/GATPeakShapeUtil.cc """

    gaus = evalGaus(x,mu,sig)
    tailH = evalXGaus(x,mu,sig,tau)
    y = (x-mu)/sig
    sigtau = -sig/tau

    gamp = tailH
    gmu = amp/tau * (-tailH + gaus)
    gsig = amp/tau * (sigtau * tailH - (sigtau - y) * gaus)
    gtau = amp/tau * (-(1. + sigtau * y + sigtau**2.) + tailH + sigtau**2.*gaus) * amp

    return np.asarray((gamp,gmu,gsig,gtau))

def xgModelWF(dataTS, floats):
    """ Make a model waveform: Take a timestamp vector, generate an
        xGauss model, normalize to 1, then scale its max value to amp.
    """
    amp, mu, sig, tau = floats
    model = evalXGaus(dataTS,mu,sig,tau)

    # simple scaling by amp
    # model = model * amp

    # pin max value of function to amp
    model = model * 1./np.sum(model)
    xMax = np.argmax(model)
    model = model * (amp / model[xMax])

    return model

def MakeTracesGlobal():
    """ This is so 'lnLike' can write to the trace arrays. Has to remain in this file to work. """
    tmp1 = tmp2 = tmp3 = tmp4 = np.empty([1,])
    global ampTr, muTr, sigTr, tauTr
    ampTr, muTr, sigTr, tauTr = tmp1, tmp2, tmp3, tmp4

def lnLike(floats, *datas):
    """ log-likelihood function: L(A,mu,sig,tau)
    To make this work with op.minimize, 'datas' is passed in as a tuple (the asterisk),
    where the original list is the 1st element.
    """
    # Add to traces.
    global ampTr, muTr, sigTr, tauTr
    amp, mu, sig, tau = floats
    ampTr = np.append(ampTr, amp)
    muTr = np.append(muTr, mu)
    sigTr = np.append(sigTr, sig)
    tauTr = np.append(tauTr, tau)

    dataTS, data, dataNoise = datas[0][0], datas[0][1], datas[0][2]

    model = xgModelWF(dataTS, floats)
    lnLike = 0.5 * np.sum ( np.power((data-model)/dataNoise, 2) - np.log( 1 / np.power(dataNoise,2) ) )
    return lnLike

def lnLikeGrad(floats, *datas):
    """ grad of log-likelihood function: grad(L(A,mu,sig,tau)) """

    amp, mu, sig, tau = floats
    dataTS, data, dataNoise = datas[0][0], datas[0][1], datas[0][2]

    model = xgModelWF(dataTS, floats)

    # FIXME: do you need to add the original amplitude parameter "h" back in, instead of treating amp == h ?
    grads = evalXgGrad(dataTS,amp,mu,sig,tau) # [4,2016]

    g_amp = 1./dataNoise * np.sum( np.multiply(grads[0],(model - data)) )
    g_mu =  1./dataNoise * np.sum( np.multiply(grads[1],(model - data)) )
    g_sig = 1./dataNoise * np.sum( np.multiply(grads[2],(model - data)) )
    g_tau = 1./dataNoise * np.sum( np.multiply(grads[3],(model - data)) )

    return np.asarray((g_amp, g_mu, g_sig, g_tau))

def numGradient(floats, f, datas, epsilon=1e-8):
    """ Adapted from 'approx_fprime' in scipy's optimize.py, starting at line 649: https://github.com/scipy/scipy/blob/5530bf1d76918c2c49994b431e076723bcf2943b/scipy/optimize/optimize.py
    """
    # f0 = f(*((floats,) + datas)) # original
    f0 = f(floats,datas) # modified by clint to accept 'datas' as a list

    grad = np.zeros((len(floats),), float)
    ei = np.zeros((len(floats),), float)

    for k in range(len(floats)):
        ei[k] = 1.0
        d = epsilon * ei
        f1 = f(floats+d,datas)
        # print "F0:",f0,"F1:",f1, "d:",d
        grad[k] = (f1 - f0) / d[k]
        ei[k] = 0.0

    return grad



if __name__ == "__main__":
    main(sys.argv[1:])
