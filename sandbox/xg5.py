#!/usr/local/bin/python
import sys, time, warnings
import numpy as np
import scipy.special as sp
import scipy.optimize as op
import matplotlib.pyplot as plt
limit = sys.float_info.max # equivalent to std::numeric_limits::max() in C++

def main():
    """ What's better than one xGauss?  TWO xGauss's. """

    file5 = np.load("./data/ds5exampleWaveform5.npz")
    dataTS, data, data5_denoised, dataENM, dataNoise = file5['arr_0'], file5['arr_1'], file5['arr_2'], file5['arr_3'], file5['arr_5']

    amp, mu, sig, tau, amp2, sig2 = dataENM, 10000., 600., -72000., 1., 10.
    floats = np.asarray([amp,mu,sig,tau,amp2,sig2])
    datas = [dataTS, data, dataNoise]

    guess = xg2ModelWF(dataTS, floats)

    print "Running minimizer..."
    MakeTracesGlobal()
    start = time.clock()
    # result = op.minimize(lnLike, floats, args=datas, method="Powell")

    # L-BGFS-B

    bnd = ((None,None),(None,None),(None,None),(None,None))
    opts = {'disp': None,   # None, True to print convergence messages
            'maxls': 100,   # 20, max line search steps
            'iprint': -1,   # -1
            'gtol': 1e-08,  # 1e-05
            'eps': 1e-08,   # 1e-08
            'maxiter': 15000,   # 15000
            'ftol': 2.220446049250313e-09,
            'maxcor': 10,   # 10
            'maxfun': 15000}    # 15000

    result = op.minimize(lnLike, floats, args=datas, method="L-BFGS-B", options=opts, bounds=None)#, jac=lnLikeGrad)
    stop = time.clock()
    print "Time:",stop-start,"seconds."

    amp, mu, sig, tau, amp2, sig2 = result["x"]
    floats = [amp, mu, sig, tau, amp2, sig2]
    fit = xg2ModelWF(dataTS, floats)


    fig = plt.figure(figsize=(13,7), facecolor='w')
    p0 = plt.subplot2grid((6,12), (0,0), colspan=12, rowspan=3) # wf plot
    p1 = plt.subplot2grid((6,12), (3,0), colspan=12, rowspan=1) # residual
    p2 = plt.subplot2grid((6,12), (4,0), colspan=2, rowspan=2) # traces
    p3 = plt.subplot2grid((6,12), (4,2), colspan=2, rowspan=2)
    p4 = plt.subplot2grid((6,12), (4,4), colspan=2, rowspan=2)
    p5 = plt.subplot2grid((6,12), (4,6), colspan=2, rowspan=2)
    p6 = plt.subplot2grid((6,12), (4,8), colspan=2, rowspan=2)
    p7 = plt.subplot2grid((6,12), (4,10), colspan=2, rowspan=2)
    p0.plot(dataTS,data,color='blue',label='data')
    p0.plot(dataTS,guess,color='orange',label='guess')
    p0.plot(dataTS,fit,color='red',label='2xgauss fit')
    p0.legend(loc=4)
    p1.plot(dataTS,data-fit,color='blue',label='residual')
    p1.legend(loc=1)
    p2.plot(ampTr[1:],label='amp',color='red')
    p2.legend(loc=1)
    p3.plot(muTr[1:],label='mu',color='green')
    p3.legend(loc=1)
    p4.plot(sigTr[1:],label='sig',color='blue')
    p4.legend(loc=1)
    p5.plot(tauTr[1:],label='tau',color='black')
    p5.legend(loc=1)
    p6.plot(amp2Tr[1:],label='amp2',color='magenta')
    p6.legend(loc=1)
    p7.plot(sig2Tr[1:],label='sig2',color='orange')
    p7.legend(loc=1)

    plt.tight_layout()
    plt.show()



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
        xGauss model, normalize to 1, then scale its max value to xAmp.
    """
    xAmp, xMu, xSig, xTau = floats
    model = evalXGaus(dataTS,xMu,xSig,xTau)

    # simple scaling by xAmp
    # model = model * xAmp

    # pin max value of function to xAmp
    model = model * 1./np.sum(model)
    xMax = np.argmax(model)
    model = model * (xAmp / model[xMax])

    return model

def MakeTracesGlobal():
    """ This is so 'lnLike' can write to the trace arrays. Has to remain in this file to work. """
    tmp1 = tmp2 = tmp3 = tmp4 = tmp5 = tmp6 = np.empty([1,])
    global ampTr, muTr, sigTr, tauTr, amp2Tr, sig2Tr
    ampTr, muTr, sigTr, tauTr, amp2Tr, sig2Tr = tmp1, tmp2, tmp3, tmp4, tmp5, tmp6


def lnLikeGrad(floats, *datas):
    """ grad of log-likelihood function: grad(L(A,mu,sig,tau)) """

    xAmp, xMu, xSig, xTau = floats
    dataTS, data, dataNoise = datas[0][0], datas[0][1], datas[0][2]

    model = xgModelWF(dataTS, floats)

    # FIXME: do you need to add the original amplitude parameter "h" back in, instead of treating xAmp == h ?
    grads = evalXgGrad(dataTS,xAmp,xMu,xSig,xTau) # [4,2016]

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

    amp, mu, sig, tau, amp2, sig2 = floats

    f1 = evalXGaus(dataTS,mu,sig,tau)
    f2 = amp2 * evalXGaus(dataTS,mu,sig2,tau) # only let amp2 and sig2 float
    ftot = f1 + f2

    # pin max value of function to amp
    ftot = ftot * 1./np.sum(ftot)
    xMax = np.argmax(ftot)
    ftot = ftot * (amp / ftot[xMax])

    return ftot


if __name__ == "__main__":
    main()
