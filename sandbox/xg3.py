#!/usr/bin/env python
import sys, time
import numpy as np
import scipy.special as sp
import scipy.optimize as op
import matplotlib.pyplot as plt
limit = sys.float_info.max # equivalent to std::numeric_limits::max() in C++

def main():

    file5 = np.load("./data/ds5exampleWaveform5.npz")
    dataTS, data, data5_denoised, dataENM, dataNoise = file5['arr_0'], file5['arr_1'], file5['arr_2'], file5['arr_3'], file5['arr_5']

    xAmp, xMu, xSig, xTau = dataENM, 10000., 600., -72000.  # xMu should be from GAT data
    floats = np.asarray([xAmp, xMu, xSig, xTau])
    datas = [dataTS, data, dataNoise]

    guess = xgModelWF(dataTS,floats)

    MakeTracesGlobal()

    print "first lnLike:",lnLike(floats, datas)
    print "guess floats:",xAmp,xMu,xSig,xTau

    # print "Gradients of guess:"
    n = numGradient(floats, lnLike, datas, 1e-5) # 1e-8 is the minimum
    c = lnLikeGrad(floats,datas)
    print "Numerical : %.2e  %.2e  %.2e  %.2e" % (n[0],n[1],n[2],n[3])
    print "Analytic  : %.2e  %.2e  %.2e  %.2e" % (c[0],c[1],c[2],c[3])


    print "Running minimizer..."
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

    print "NUM GRADIENT"
    result = op.minimize(lnLike, floats, args=datas, method="L-BFGS-B", options=opts, bounds=bnd)#, jac=lnLikeGrad)

    # print "ANA GRADIENT"
    # result = op.minimize(lnLike, floats, args=datas, method="L-BFGS-B", options=opts, bounds=None, jac=lnLikeGrad)

    xAmp, xMu, xSig, xTau = result["x"]
    floats = [xAmp, xMu, xSig, xTau]
    fit = xgModelWF(dataTS, floats)
    print "fit floats: ",xAmp,xMu,xSig,xTau
    stop = time.clock()
    print "Time:",stop-start,"seconds."


    print "Making plots ..."
    fig = plt.figure(figsize=(10,7), facecolor='w')
    p0 = plt.subplot2grid((6,8), (0,0), colspan=8, rowspan=3)
    p1 = plt.subplot2grid((6,8), (4,0), colspan=2, rowspan=2)
    p2 = plt.subplot2grid((6,8), (4,2), colspan=2, rowspan=2)
    p3 = plt.subplot2grid((6,8), (4,4), colspan=2, rowspan=2)
    p4 = plt.subplot2grid((6,8), (4,6), colspan=2, rowspan=2)
    p5 = plt.subplot2grid((6,8), (3,0), colspan=8, rowspan=1)

    p0.plot(dataTS,data,color='blue',label='data')
    p0.plot(dataTS,guess,color='orange',label='xgauss guess')
    p0.plot(dataTS,fit,color='red',label='xgauss fit')

    p0.set_title("xAmp %.2f  xMu %.2f  xSig %.2f  xTau %.2f" %(xAmp,xMu,xSig,xTau))
    p0.legend(loc=4)

    p1.plot(xAmpTrace[1:],label='xAmp',color='red')
    p1.legend(loc=1)
    p1.set_xlabel('Fit Steps')
    p2.plot(xMuTrace[1:],label='xMu',color='green')
    p2.legend(loc=1)
    p3.plot(xSigTrace[1:],label='xSig',color='blue')
    p3.legend(loc=1)
    p4.plot(xTauTrace[1:],label='xTau',color='black')
    p4.legend(loc=1)
    p5.plot(dataTS,data-fit,color='blue',label='residual')
    p5.legend(loc=1)
    plt.tight_layout()


    # given the traces, show how the gradients were moving around
    gAmpTrace = gMuTrace = gSigTrace = gTauTrace = np.empty([1,])
    for i in range(len(xAmpTrace)):
        floats = [xAmpTrace[i],xMuTrace[i],xSigTrace[i],xTauTrace[i]]
        c = lnLikeGrad(floats,datas)
        gAmpTrace = np.append(gAmpTrace,c[0])
        gMuTrace = np.append(gMuTrace,c[1])
        gSigTrace = np.append(gSigTrace,c[2])
        gTauTrace = np.append(gTauTrace,c[3])


    fig2 = plt.figure(figsize=(10,2), facecolor='w')
    p6 = plt.subplot2grid((1,4),(0,0))
    p7 = plt.subplot2grid((1,4),(0,1))
    p8 = plt.subplot2grid((1,4),(0,2))
    p9 = plt.subplot2grid((1,4),(0,3))
    p6.plot(gAmpTrace[2:],label='gAmp',color='red')
    p6.legend(loc=1)
    p7.plot(gMuTrace[2:],label='gMu',color='green')
    p7.legend(loc=1)
    p8.plot(gSigTrace[2:],label='gSig',color='blue')
    p8.legend(loc=1)
    p9.plot(gTauTrace[2:],label='gTau',color='black')
    p9.legend(loc=1)
    plt.tight_layout()


    plt.show()



def evalGaus(x,mu,sig):
    return np.exp(-((x-mu)**2./2./sig**2.))

def evalXGaus(x,mu,sig,tau):
    """ Ported from GAT/BaseClasses/GATPeakShapeUtil.cc
        Negative tau: "Regular WF, high tail"
        Positive tau: "Backwards WF, low tail"
    """
    tmp = (x-mu + sig**2./2./tau)/tau
    # print "types:",type(tmp),"tmp is",tmp
    if all(tmp < limit):
        return np.exp(tmp)/2./np.fabs(tau) * sp.erfc((tau*(x-mu)/sig + sig)/np.sqrt(2.)/np.fabs(tau))
    else:
        # in this case exp returns NaN (in C++).  So use an approx. derived from the asymptotic expansion for erfc, listed on wikipedia.
        print "Exceeded limit ..."
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
    tmp1 = tmp2 = tmp3 = tmp4 = np.empty([1,])
    global xAmpTrace, xMuTrace, xSigTrace, xTauTrace
    xAmpTrace, xMuTrace, xSigTrace, xTauTrace = tmp1, tmp2, tmp3, tmp4

def lnLike(floats, *datas):
    """ log-likelihood function: L(A,mu,sig,tau)
    To make this work with op.minimize, 'datas' is passed in as a tuple (the asterisk),
    where the original list is the 1st element.
    """
    # Add to traces.
    global xAmpTrace, xMuTrace, xSigTrace, xTauTrace
    xAmp, xMu, xSig, xTau = floats
    xAmpTrace = np.append(xAmpTrace, xAmp)
    xMuTrace = np.append(xMuTrace, xMu)
    xSigTrace = np.append(xSigTrace, xSig)
    xTauTrace = np.append(xTauTrace, xTau)

    dataTS, data, dataNoise = datas[0][0], datas[0][1], datas[0][2]

    model = xgModelWF(dataTS, floats)
    lnLike = 0.5 * np.sum ( np.power((data-model)/dataNoise, 2) - np.log( 1 / np.power(dataNoise,2) ) )
    return lnLike

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


if __name__ == "__main__":
    main()
