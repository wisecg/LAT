#!/usr/local/bin/python
import sys, time
import numpy as np
import scipy.optimize as op
import scipy.special as sp
import matplotlib.pyplot as plt

def main(argv):
    # checkTemplates()
    plotPartials()

def plotPartials():

    file5 = np.load("./data/ds5exampleWaveform5.npz")
    dataTS, data, data5_denoised, dataENM, dataNoise = file5['arr_0'], file5['arr_1'], file5['arr_2'], file5['arr_3'], file5['arr_5']
    xAmp, xMu, xSig, xTau = dataENM, 10000., 600., 72000.  # xMu should be from GAT data

    consts = [dataTS]
    floats = [xAmp, xMu, xSig, xTau]
    datas = (dataTS, xAmp, xMu, xSig, xTau, data, dataNoise)

    MakeTracesGlobal()
    MakeGradsGlobal()

    # this works
    print "numerical gradient:",approx_fprime(floats, xgLnLike, 1e-08, datas)

    # this doesn't ... yet
    print xgGradLnLike(floats, datas)

    # guess = xgModelWF(consts, floats)
    #
    # fig = plt.figure(figsize=(9,7),facecolor='w')
    # p1 = plt.subplot(211)
    # p1.plot(dataTS,data,color='blue')
    # p1.plot(dataTS,guess,color='orange')
    #
    # pA = xgPartial("A",consts,floats)
    # pMu = xgPartial("mu",consts,floats)
    # pSig = xgPartial("sig",consts,floats)
    # pTau = xgPartial("tau",consts,floats)
    #
    # p2 = plt.subplot(212)
    # p2.plot(dataTS,pA,label='pA')
    # p2.plot(dataTS,pMu,label='pMu')
    # # p2.plot(dataTS,pSig,label='pSig')
    # p2.plot(dataTS,pTau,label='pTau')
    # p2.legend(loc=4)
    #
    # plt.tight_layout()
    # plt.show()

def checkTemplates():

    # =================== Load data =======================

    # load a pysiggen "data waveform"
    pySamp, pyR, pyZ, pyAmp, pyST, pySlo = 2016, 0, 15, 10, 1000, 10
    pyFile = np.load("./data/pysig_test.npz")
    pyData, pyDataTS = pyFile['arr_0'], pyFile['arr_1']

    # load some actual data waveforms (5 and 238 keV, denoised and not-denoised)
    file5 = np.load("./data/ds5exampleWaveform5.npz")
    data5TS, data5, data5_denoised = file5['arr_0'], file5['arr_1'], file5['arr_2'],
    data5ENM, data5Max = file5['arr_3'], file5['arr_4']
    data5Noise, data5DenoisedNoise = file5['arr_5'], file5['arr_6']

    file238 = np.load("./data/ds5exampleWaveform238.npz")
    data238TS, data238, data238_denoised = file238['arr_0'], file238['arr_1'], file238['arr_2']
    data238ENM, data238Max = file238['arr_3'], file238['arr_4']
    data238Noise, data238DenoisedNoise = file238['arr_5'], file238['arr_6']

    # pick the data waveform to use
    # pysiggen
    # data, dataTS, dataENM, dataNoise = pyData, pyDataTS, pyAmp, 1.
    # dataTSMax = np.argmax(data) * 10

    # 5 kev
    data, dataTS, dataENM, dataTSMax, dataNoise = data5, data5TS, data5ENM, data5Max, data5Noise

    # 5 keV denoised
    # data, dataTS, dataENM, dataTSMax, dataNoise = data5_denoised, data5TS, data5ENM, data5Max, data5DenoisedNoise

    # 238 kev
    # data, dataTS, dataENM, dataTSMax, dataNoise = data238, data238TS, data238ENM, data238Max, data238Noise

    # 238 kev denoised
    # data, dataTS, dataENM, dataTSMax, dataNoise = data238_denoised, data238TS, data238ENM, data238Max, data238DenoisedNoise


    # ============ run waveform fitter ====================

    # make initial guesses for our xgauss waveform and save it for plotting
    xAmp, xMu, xSig, xTau = dataENM, 10000., 600., 72000.  # xMu should be from GAT data
    consts = [dataTS]
    floats = [xAmp, xMu, xSig, xTau]

    print "guess floats: ",xAmp,xMu,xSig,xTau
    guess = xgModelWF(consts, floats)

    start = time.clock()
    MakeTracesGlobal()
    MakeGradsGlobal()
    datas = (dataTS, xAmp, xMu, xSig, xTau, data, dataNoise)

    # no-gradient method: Powell / Nelder-Mead
    # result = op.minimize(xgLnLike, floats, args=datas, method="Powell")

    # stock 'op.minimize'
    # opts = {'disp' : True,     # Set to True to print convergence messages.
    #         'maxiter' : None,  # Maximum number of iterations to perform.
    #         'gtol': 1e-8,      # Gradient norm must be less than gtol before successful termination.
    #         'norm': np.inf ,   # Order of norm (Inf is max, -Inf is min).
    #         'eps': 1.4901161193847656e-08} # If jac is approximated, use this value for the step size.
    # result = op.minimize(xgLnLike, floats, args=datas, options=opts, method="BFGS", jac=xgGradLnLike)

    # stock 'l-bfgs-b'
    opt = {'disp':True}
    bnd = ((None,None),(None,None),(None,None),(0.01,1000000))
    # bnd = None
    # jac = xgGradLnLike
    grad = None
    result = op.minimize(xgLnLike, floats, args=datas, method="L-BFGS-B", options=opt, bounds=bnd, jac=grad)

    if not result["success"]: print "fit fail, message:",result["message"]
    xAmp, xMu, xSig, xTau = result["x"]

    floats = [xAmp, xMu, xSig, xTau]
    fit = xgModelWF(consts, floats)


    print "fit floats: ",xAmp,xMu,xSig,xTau
    stop = time.clock()
    print "Time:",stop-start,"seconds."

    # calculate fitChi2.  Textbook is (observed - expected)^2 / expected,
    # but we'll follow MGWFCalculateChiSquare.cc and do (observed - expected)^2 / NDF.
    # fitChi2 = np.sum(np.square(data - fit)) / len(data)
    # print "fitChi2:",fitChi2

    # =====================================================

    # do plotting

    fig = plt.figure(figsize=(10,7), facecolor='w')
    p0 = plt.subplot2grid((6,8), (0,0), colspan=8, rowspan=3) # data & fit
    p1 = plt.subplot2grid((6,8), (4,0), colspan=2, rowspan=2) # trace 1
    p2 = plt.subplot2grid((6,8), (4,2), colspan=2, rowspan=2) # trace 2
    p3 = plt.subplot2grid((6,8), (4,4), colspan=2, rowspan=2) # trace 3
    p4 = plt.subplot2grid((6,8), (4,6), colspan=2, rowspan=2) # trace 4
    p5 = plt.subplot2grid((6,8), (3,0), colspan=8, rowspan=1) # residual

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

    # fig2 = plt.figure(figsize=(10,2), facecolor='w')
    # p6 = plt.subplot2grid((1,4),(0,0))
    # p7 = plt.subplot2grid((1,4),(0,1))
    # p8 = plt.subplot2grid((1,4),(0,2))
    # p9 = plt.subplot2grid((1,4),(0,3))
    # p6.plot(gAmpTrace[1:],color='red')
    # p7.plot(gMuTrace[1:],color='green')
    # p8.plot(gSigTrace[1:],color='blue')
    # p9.plot(gTauTrace[1:],color='black')
    # plt.tight_layout()

    plt.show()

# ==============================================================
def xgFunction(x,mu,sig,tau):
    """ Using first dfn given on wikipedia. Hoping python can handle any overflow. """
    return (1./(2.*tau)) * np.exp( (1./(2.*tau)) * (2.*mu + sig**2./tau - 2.*x) ) * sp.erfc((mu + sig**2./tau - x)/(np.sqrt(2.)*sig))

def xgModelWF(consts, floats):
    """ Make a model waveform: Take a timestamp vector, generate an
        xGauss model, normalize to 1, then scale its max value to xAmp.
    """
    xAmp, xMu, xSig, xTau = floats
    modelTS = consts[0]
    model = xgFunction(modelTS,xMu,xSig,xTau)
    model = model * 1./np.sum(model)
    xMax = np.argmax(model)
    model = model * (xAmp / model[xMax])
    return model

def MakeTracesGlobal():
    """ This is so 'xgLnLike' can write to the trace arrays.
        Has to remain in this file to work.
    """
    tmp1 = tmp2 = tmp3 = tmp4 = np.empty([1,])
    global xAmpTrace, xMuTrace, xSigTrace, xTauTrace
    xAmpTrace, xMuTrace, xSigTrace, xTauTrace = tmp1, tmp2, tmp3, tmp4

def xgLnLike(floats, *datas):
    """ Calculate Log-likelihood of an xGauss function, for use w/ scipy minimizers.
        Can't be moved to another lib unless you don't need the trace arrays.
    """
    global xAmpTrace, xMuTrace, xSigTrace, xTauTrace

    # Can fix parameters here by setting them back to their original guess values
    xAmp, xMu, xSig, xTau = floats
    # if True: xTau = datas[4]
    floats = [xAmp, xMu, xSig, xTau]

    # Add to traces
    xAmpTrace = np.append(xAmpTrace, xAmp)
    xMuTrace = np.append(xMuTrace, xMu)
    xSigTrace = np.append(xSigTrace, xSig)
    xTauTrace = np.append(xTauTrace, xTau)

    # print xAmp, xMu, xSig, xTau

    # Generate the xGauss model waveform
    consts = [datas[0]]
    model = xgModelWF(consts, floats)

    # Find LL of data vs. model
    data, dataNoise = datas[5], datas[6]
    lnLike = 0.5 * np.sum ( np.power((data-model)/dataNoise, 2) - np.log( 1 / np.power(dataNoise,2) ) )
    return lnLike

def MakeGradsGlobal():
    """ This is to plot the value of the gradient vs. step number."""
    tmp1 = tmp2 = tmp3 = tmp4 = np.empty([1,])
    global gAmpTrace, gMuTrace, gSigTrace, gTauTrace
    gAmpTrace, gMuTrace, gSigTrace, gTauTrace = tmp1, tmp2, tmp3, tmp4

def xgGradLnLike(floats, *datas):
    """ Gradient of log-likelihood, returning an array of partials for each point. """
    global gAmpTrace, gMuTrace, gSigTrace, gTauTrace

    xAmp, xMu, xSig, xTau = floats
    # if True: xTau = datas[4]
    # if True: xSig = datas[3]
    floats = [xAmp, xMu, xSig, xTau]

    # print datas[0][1]
    # data, dataNoise = datas[5], datas[6]
    data, dataNoise = datas[0][5], datas[0][6]


    # consts = [datas[0]] # worked for calling it from op.minimize
    consts = [datas[0][0]] # worked for calling xgGradLnLike directly
    model = xgModelWF(consts, floats)

    # this is a good idea ?
    ga = (1./dataNoise) * np.sum( np.multiply(xgPartial("A",consts,floats),(model - data) ))
    gmu = (1./dataNoise) * np.sum( np.multiply(xgPartial("mu",consts,floats),(model - data) ))
    gsig = (1./dataNoise) * np.sum( np.multiply(xgPartial("sig",consts,floats),(model - data) ))
    gtau = (1./dataNoise) * np.sum( np.multiply(xgPartial("tau",consts,floats),(model - data) ))

    print 1./dataNoise, np.sum( np.multiply(xgPartial("A",consts,floats),(model - data) )), np.sum(model-data), np.sum(xgPartial("A",consts,floats))

    # np.set_printoptions(threshold=np.inf)
    # print xgPartial("mu",consts,floats)
    print np.sum(xgPartial("sig",consts,floats))



    # this is a bad idea ?
    # ga = np.sum( xgPartial("A",consts,floats) )
    # gmu = np.sum( xgPartial("mu",consts,floats) )
    # gsig = np.sum( xgPartial("sig",consts,floats) )
    # gtau = np.sum( xgPartial("tau",consts,floats) )

    # print "A %.5f %.4e  mu %.2f %.10e  sig %.2f %.10e  tau %.2f %.10e" % (xAmp,ga,xMu,gmu,xSig,gsig,xTau,gtau)

    # add to traces
    gAmpTrace = np.append(gAmpTrace, ga)
    gMuTrace = np.append(gMuTrace, gmu)
    gSigTrace = np.append(gSigTrace, gsig)
    gTauTrace = np.append(gTauTrace, gtau)

    return np.asarray((ga, gmu, gsig, gtau))

def xgPartial(opt,consts,floats):
    """ Get a numpy array of the requested partial of xGauss."""

    modelTS = consts[0]
    xAmp, xMu, xSig, xTau = floats

    if opt!="A":
        partial = xgGradient(opt,modelTS,xMu,xSig,xTau)
    elif opt=="A":
        floats = [1.,xMu,xSig,xTau] # since xAmp is just a multiplier, the partial is just the model with xAmp==1.
        partial = xgModelWF(consts,floats)

    # normalize
    # partial = (partial / (np.sum(partial)))
    # print "opt",opt,":",partial

    return partial

def xgGradient(opt,x,mu,sig,tau):
    """ Calculate gradient at a single point. """

    A = (-1./(2.*sig**2.)) * (mu + sig**2./tau - x)**2.
    B = (-2.*x + 2.*mu + sig**2./tau) / (2.*tau)
    C = (mu + sig**2./tau - x)/(np.sqrt(2.)*sig)

    if opt=="mu":
        dmu = (-1./(np.sqrt(2. * np.pi)*sig*tau)) * (np.exp(A+B))
        dmu += (1./2.*tau**2.) * (np.exp(B) * sp.erfc(C) )
        return dmu

    if opt=="sig":
        dsig = (-1./(np.pi*tau)) * (np.exp(A+B)) * (-1. * (mu + sig**2./tau - x)/(np.sqrt(2.)*sig**2.) + np.sqrt(2.)/tau )
        dsig += (1./(2*tau**3.)) * (np.exp(B) * sig * sp.erfc(C))
        return dsig

    if opt=="tau":
        dtau = (1./(np.sqrt(2.*np.pi)*tau**3.)) * (np.exp(A+B) * sig)
        dtau += (-1./(2.*tau**2.)) * (np.exp(B) * sp.erfc(C))
        dtau += (1./(2.*tau)) * (np.exp(B) * ( -1.*sig**2./(2.*tau**3.) - (2.*mu + sig**2./tau - 2.*x)/(2.*tau**2.)  ) * sp.erfc(C)  )
        return dtau

def approx_fprime(xk, f, epsilon=1e-8, args=()):
    """ Ripped out of scipy's optimize.py, starting at line 649: https://github.com/scipy/scipy/blob/5530bf1d76918c2c49994b431e076723bcf2943b/scipy/optimize/optimize.py """
    f0 = f(*((xk,) + args))
    grad = np.zeros((len(xk),), float)
    ei = np.zeros((len(xk),), float)
    for k in range(len(xk)):
        ei[k] = 1.0
        d = epsilon * ei
        grad[k] = (f(*((xk + d,) + args)) - f0) / d[k]
        ei[k] = 0.0
    return grad


if __name__ == "__main__":
    main(sys.argv[1:])