#!/usr/local/bin/python
import sys, time, imp, pywt
import numpy as np
import scipy.special as sp
import scipy.optimize as op
import matplotlib.pyplot as plt
from ROOT import TFile,TTree,TEntryList,MGTWaveform,gDirectory
wl = imp.load_source('waveLibs', '../waveLibs.py')
limit = sys.float_info.max # equivalent to std::numeric_limits::max() in C++

def main(argv):

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
    theCut += " && trapENFCal > 3"

    waveTree.Draw(">>elist", theCut, "entrylist")
    elist = gDirectory.Get("elist")
    waveTree.SetEntryList(elist)
    nList = elist.GetN()

    print "Using cut:\n",theCut,"\nFound",nList,"entries passing cuts."

    # fit with one xgauss
    fig = plt.figure(figsize=(13,7), facecolor='w')
    p0 = plt.subplot2grid((6,8), (0,0), colspan=8, rowspan=3) # wf plot
    p1 = plt.subplot2grid((6,8), (3,0), colspan=8, rowspan=1) # residual
    p2 = plt.subplot2grid((6,8), (4,0), colspan=2, rowspan=2) # traces
    p3 = plt.subplot2grid((6,8), (4,2), colspan=2, rowspan=2)
    p4 = plt.subplot2grid((6,8), (4,4), colspan=2, rowspan=2)
    p5 = plt.subplot2grid((6,8), (4,6), colspan=2, rowspan=2)

    # fit with two xgauss's
    # fig = plt.figure(figsize=(11,7), facecolor='w')
    # p0 = plt.subplot2grid((8,8), (0,0), colspan=8, rowspan=3) # wf plot
    # p1 = plt.subplot2grid((8,8), (3,0), colspan=8, rowspan=1) # residual
    # p2 = plt.subplot2grid((8,8), (4,0), colspan=2, rowspan=2) # traces
    # p3 = plt.subplot2grid((8,8), (4,2), colspan=2, rowspan=2)
    # p4 = plt.subplot2grid((8,8), (4,4), colspan=2, rowspan=2)
    # p5 = plt.subplot2grid((8,8), (4,6), colspan=2, rowspan=2)
    # p6 = plt.subplot2grid((8,8), (6,0), colspan=2, rowspan=2)
    # p7 = plt.subplot2grid((8,8), (6,2), colspan=2, rowspan=2)
    # p8 = plt.subplot2grid((8,8), (6,4), colspan=2, rowspan=2)

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
            trapENFCal = waveTree.trapENFCal.at(iH)
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

            # fit with one xgauss
            amp, mu, sig, tau = dataENM, dataTSMax, 600., -72000.
            floats = np.asarray([amp, mu, sig, tau])
            guess = xgModelWF(dataTS,floats)
            MakeTracesGlobal()

            # fit with two xgauss's
            # amp, mu, sig, tau, amp2, mu2, sig2 = dataENM, 10000., 600., -72000., 1., 10000., 10.
            # floats = np.asarray([amp,mu,sig,tau,amp2,mu2,sig2])
            # guess = xg2ModelWF(dataTS, floats)
            # MakeTracesGlobal2()

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
            # result = op.minimize(lnLike2, floats, args=datas, method="L-BFGS-B", options=opts, bounds=None)#, jac=lnLikeGrad)

            # analytical gradient - BETA, Don't Use!
            # result = op.minimize(lnLike, floats, args=datas, method="L-BFGS-B", options=opts, bounds=None, jac=lnLikeGrad)

            stop = time.clock()

            if not result["success"]:
                print "fit 'fail': ", result["message"]
                # errorCode[0] = 1

            # fit with one xg
            amp, mu, sig, tau = result["x"]
            floats = [amp, mu, sig, tau]
            fit = xgModelWF(dataTS, floats)

            # fit with two xg's
            # amp, mu, sig, tau, amp2, mu2, sig2 = result["x"]
            # floats = [amp, mu, sig, tau, amp2, mu2, sig2]
            # fit = xg2ModelWF(dataTS, floats)

            # fitchi2 and finalL

            # calculate fitChi2.  Textbook is (observed - expected)^2 / expected,
            # but we'll follow MGWFCalculateChiSquare.cc and do (observed - expected)^2 / NDF.
            # NOTE: we're doing the chi2 against the DATA, though the FIT is to the DENOISED DATA.

            fitChi2 = np.sum(np.square(data - fit)) / len(data)

            finalLL = result["fun"]

            title = "%d cal %.1f  enm %.1f  amp %.2f  mu %.2f  sig %.2f  tau %.2f  chi2 %.2f  spd %.3f" % (iList,trapENFCal,dataENM,amp,mu,sig,tau,fitChi2,stop-start)

            # title = "%d cal %.1f  enm %.1f  amp %.2f  mu %.2f  sig %.2f  tau %.2f\namp2 %.2f  mu2 %.2f  sig2 %.2f  chi2 %.2f  spd %.3f" % (iList,trapENFCal,dataENM,amp,mu,sig,tau,amp2,mu2,sig2,fitChi2,stop-start)

            print title


            # =======================================================
            # plots

            if batMode: continue

            p0.cla()
            p0.plot(dataTS,data,color='blue',label='data')
            p0.plot(dataTS,data_wlDenoised,color='orange',label='wlDenoised')
            # p0.plot(dataTS,guess,color='orange',label='xgauss guess')
            p0.plot(dataTS,fit,color='red',label='xgauss fit')
            p0.set_title(title)
            p0.legend(loc=4)
            p1.cla()
            p1.plot(dataTS,data-fit,color='blue',label='residual')
            p1.legend(loc=1)
            p2.cla()
            p2.plot(ampTr[1:],label='amp',color='red')
            p2.legend(loc=1)
            p3.cla()
            p3.plot(muTr[1:],label='mu',color='green')
            p3.legend(loc=1)
            p4.cla()
            p4.plot(sigTr[1:],label='sig',color='blue')
            p4.legend(loc=1)
            p5.cla()
            p5.plot(tauTr[1:],label='tau',color='black')
            p5.legend(loc=1)
            # p6.cla()
            # p6.plot(amp2Tr[1:],label='amp2',color='magenta')
            # p6.legend(loc=1)
            # p7.cla()
            # p7.plot(mu2Tr[1:],label='mu2',color='cyan')
            # p7.legend(loc=1)
            # p8.cla()
            # p8.plot(sig2Tr[1:],label='sig2',color='orange')
            # p8.legend(loc=1)

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

def xg2ModelWF(dataTS, floats):
    """ Crazy idea: try fitting with TWO xGauss's. """

    amp, mu, sig, tau, amp2, mu2, sig2 = floats

    f1 = evalXGaus(dataTS,mu,sig,tau)
    f2 = amp2 * evalXGaus(dataTS,mu,sig2,tau) # let amp2, mu2, sig2 float
    ftot = f1 + f2

    # pin max value of function to amp
    ftot = ftot * 1./np.sum(ftot)
    xMax = np.argmax(ftot)
    ftot = ftot * (amp / ftot[xMax])

    return ftot

def MakeTracesGlobal2():
    """ This is so 'lnLike' can write to the trace arrays. Has to remain in this file to work. """
    tmp1 = tmp2 = tmp3 = tmp4 = tmp5 = tmp6 = tmp7 = np.empty([1,])
    global ampTr, muTr, sigTr, tauTr, amp2Tr, mu2Tr, sig2Tr
    ampTr, muTr, sigTr, tauTr, amp2Tr, mu2Tr, sig2Tr = tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7

def lnLike2(floats, *datas):
    """ log-likelihood function: L(A,mu,sig,tau,A2,sig2)
    """
    global ampTr, muTr, sigTr, tauTr, amp2Tr, mu2Tr, sig2Tr
    amp, mu, sig, tau, amp2, mu2, sig2 = floats
    ampTr = np.append(ampTr, amp)
    muTr = np.append(muTr, mu)
    sigTr = np.append(sigTr, sig)
    tauTr = np.append(tauTr, tau)
    amp2Tr = np.append(amp2Tr, amp2)
    mu2Tr = np.append(mu2Tr, mu2)
    sig2Tr = np.append(sig2Tr, sig2)

    dataTS, data, dataNoise = datas[0][0], datas[0][1], datas[0][2]

    model = xg2ModelWF(dataTS, floats)
    lnLike = 0.5 * np.sum ( np.power((data-model)/dataNoise, 2) - np.log( 1 / np.power(dataNoise,2) ) )
    return lnLike



if __name__ == "__main__":
    main(sys.argv[1:])
