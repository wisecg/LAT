#!/usr/local/bin/python
import sys, pywt, imp, time
import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate
import scipy.special as sp
from pysiggen import Detector
from ROOT import TFile,TTree,TEntryList,MGTWaveform,gDirectory
wl = imp.load_source('waveLibs', '../waveLibs.py')
import scipy.optimize as op

def main(argv):

    checkTemplates()
    # fitData(argv)


def checkTemplates():

    # generate a pysiggen "data waveform"
    pySamp, pyR, pyZ, pyAmp, pyST, pySlo = 2016, 0, 15, 10, 1000, 10
    # pyData, pyDataTS = MakeSiggenWaveform(pySamp, pyR, pyZ, pyAmp, pyST, pySlo)
    # np.savez("./data/pysig_test.npz",pyData,pyDataTS)
    pyFile = np.load("./data/pysig_test.npz")
    pyData, pyDataTS = pyFile['arr_0'], pyFile['arr_1']


    # load some actual data waveforms (5 and 238 keV, denoised and not-denoised)
    file5 = np.load("./data/ds5exAmpleWaveform5.npz")
    data5TS, data5, data5_denoised = file5['arr_0'], file5['arr_1'], file5['arr_2'],
    data5ENM, data5Max = file5['arr_3'], file5['arr_4']
    data5Noise, data5DenoisedNoise = file5['arr_5'], file5['arr_6']

    file238 = np.load("./data/ds5exAmpleWaveform238.npz")
    data238TS, data238, data238_denoised = file238['arr_0'], file238['arr_1'], file238['arr_2']
    data238ENM, data238Max = file238['arr_3'], file238['arr_4']
    data238Noise, data238DenoisedNoise = file238['arr_5'], file238['arr_6']


    # pick the data waveform to use
    # pysiggen
    # data, dataTS, dataENM, dataNoise = pyData, pyDataTS, pyAmp, 1.
    # dataTSMax = np.argmax(data) * 10

    # 5 kev
    # data, dataTS, dataENM, dataTSMax, dataNoise = data5, data5TS, data5ENM, data5Max, data5Noise

    # 5 keV denoised
    data, dataTS, dataENM, dataTSMax, dataNoise = data5_denoised, data5TS, data5ENM, data5Max, data5DenoisedNoise

    # 238 kev
    # data, dataTS, dataENM, dataTSMax, dataNoise = data238, data238TS, data238ENM, data238Max, data238Noise

    # 238 kev denoised
    # data, dataTS, dataENM, dataTSMax, dataNoise = data5_denoised, data238TS, data238ENM, data238Max, data238DenoisedNoise


    # generate an xGauss "guess" waveform
    xAmp, xMu, xSig, xTau = dataENM, 10000., 600., 72000. # guesses
    consts = [dataTS]
    floats = [xAmp, xMu, xSig, xTau]
    guess = MakeXGModel(consts, floats)


    # run waveform fitter
    start = time.clock()
    MakeXGTracesGlobal()

    # 1. Powell / Nelder-Mead
    datas = [dataTS, dataENM, 10000., 200., 72000., data, dataNoise]
    # result = op.minimize(findLnLikeXG, floats, args=datas, method="Powell") #method="Nelder-Mead"

    # 2. conjugate gradient
    result = op.fmin_cg(findLnLikeXG, floats, fprime=gradLNLike, args=tuple(datas))
    xAmp, xMu, xSig, xTau = result
    print result

    # if not result["success"]:
        # print "fit fail, message:",result["message"]
    # xAmp, xMu, xSig, xTau = result["x"]
    # floats = [xAmp, xMu, xSig, xTau]
    fit = MakeXGModel(consts, floats)
    stop = time.clock()
    print "Time:",stop-start,"seconds."

    # get residual
    resid = data - fit

    # calculate fitChi2.  Textbook is (observed - expected)^2 / expected,
    # but we'll follow MGWFCalculateChiSquare.cc and do (observed - expected)^2 / NDF.
    # TODO: check this.
    fitChi2 = np.sum(np.square(data - fit)) / len(data)
    print "fitChi2:",fitChi2


    # do plotting stuff

    fig = plt.figure(figsize=(10,7), facecolor='w')
    p0 = plt.subplot2grid((6,8), (0,0), colspan=8, rowspan=3) # data & fit
    p1 = plt.subplot2grid((6,8), (4,0), colspan=2, rowspan=2) # trace 1
    p2 = plt.subplot2grid((6,8), (4,2), colspan=2, rowspan=2) # trace 2
    p3 = plt.subplot2grid((6,8), (4,4), colspan=2, rowspan=2) # trace 3
    p4 = plt.subplot2grid((6,8), (4,6), colspan=2, rowspan=2) # trace 4
    p5 = plt.subplot2grid((6,8), (3,0), colspan=8, rowspan=1) # residual

    p0.plot(dataTS,data,color='blue',label='pysig data')
    p0.plot(dataTS,guess,color='orange',label='xgauss model')
    p0.plot(dataTS,fit,color='red',label='xgauss fit')

    p0.set_title("xAmp %.2f  xMu %.2f  xSig %.2f  xTau %.2f" %(xAmp,xMu,xSig,xTau))
    p0.legend(loc=4)

    p1.plot(xAmpTrace[1:],label='xAmp',color='red')
    p1.legend(loc=4)
    p1.set_xlabel('Fit Steps')
    p2.plot(xMuTrace[1:],label='xMu',color='green')
    p2.legend(loc=4)
    p3.plot(xSigTrace[1:],label='xSig',color='blue')
    p3.legend(loc=4)
    p4.plot(xTauTrace[1:],label='xTau',color='black')
    p4.legend(loc=4)
    p5.plot(dataTS,resid,color='blue',label='residual')
    p5.legend(loc=4)
    plt.tight_layout()
    plt.show()
    return


def xGaussTest(x,mu,sig,tau):
    z = 1/np.sqrt(2) * (sig/tau - (x - mu)/sig)
    if z < 0: return 1, z
    elif z >= 0 and z <= 6.71e7: return 2, z
    elif z > 6.71e7: return 3, z
    else: return -1, z


def xGaussFn(mode,x,h,mu,sig,tau):
    if mode==1:
        return ( h*sig*np.sqrt(np.pi/2.)/tau ) * np.exp( (sig/tau)**2./2. - (x-mu)/tau ) * sp.erfc( ( sig/tau - (x-mu)/sig )/np.sqrt(2.) )
    elif mode==2:
        return ( h*sig*np.sqrt(np.pi/2.)/tau ) * np.exp( (-1./2.)*((x-mu)/sig)**2. ) * sp.erfcx( ( sig/tau - (x-mu)/sig )/np.sqrt(2.) )
    elif mode==3:
        return ( h/(1 - (x-mu)*tau/sig**2.) ) * np.exp( (-1./2.)*((x-mu)/sig)**2. )
    else:
        print "unknown mode!"
        return -1




def fitData(argv):
    """ Compare templates and fitters by actually using them to fit data from a calibration file.
        Make a batch mode and plot the chi2 values for both methods.
    """

    # margs
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

    # Set input file and cuts
    inputFile = TFile("../waveSkimDS5_run21975.root")
    waveTree = inputFile.Get("skimTree")
    print "Found",waveTree.GetEntries(),"input entries."
    theCut = inputFile.Get("theCut").GetTitle()
    # theCut += " && Entry$ < 100"
    theCut += " && trapENFCal > 4.5 && trapENFCal < 5.5"
    print "Using cut:\n",theCut,"\n"
    waveTree.Draw(">>elist", theCut, "entrylist")
    elist = gDirectory.Get("elist")
    waveTree.SetEntryList(elist)
    nList = elist.GetN()
    print "Found",nList,"entries passing cuts."

    # Set figure
    fig = plt.figure(figsize=(10,7), facecolor='w')
    p0 = plt.subplot2grid((6,8), (0,0), colspan=8, rowspan=4) # data & fit
    p1 = plt.subplot2grid((6,8), (4,0), colspan=2, rowspan=2) # trace 1
    p2 = plt.subplot2grid((6,8), (4,2), colspan=2, rowspan=2) # trace 2
    p3 = plt.subplot2grid((6,8), (4,4), colspan=2, rowspan=2) # trace 3
    p4 = plt.subplot2grid((6,8), (4,6), colspan=2, rowspan=2) # trace 4

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
            energy = waveTree.trapENFCal.at(iH)
            wf = waveTree.MGTWaveforms.at(iH)
            signal = wl.processWaveform(wf)
            data = signal.GetWaveBLSub()
            dataTS = signal.GetTS()
            _,dataNoise = signal.GetBaseNoise()
            dataTSMax = waveTree.trapENMSample.at(iH)*10. - 4000
            dataENM = waveTree.trapENM.at(iH)

            print "%d / %d  Run %d  nCh %d  chan %d  trapENF %.1f" % (iList,nList,run,nChans,chan,energy)

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

            np.savez("./data/ds5exampleWaveform5.npz",dataTS,data,data_wlDenoised,dataENM,dataTSMax,dataNoise,denoisedNoise)

            # load and scale template waveform
            # temp, tempTS = tOrig, tOrigTS
            # temp = temp * (dataENM / tAmp)         # scale by amplitudes (not energies)
            # tempTSMax = np.argmax(temp) * 10         # convert to ns
            # tempTS = tempTS - (tempTSMax - dataTSMax)  # align @ max of rising edge
            # tempTSMax = tempTS[np.argmax(temp)]      # get the new max TS after the shifting
            # tempAmp = dataENF                      # amplitude of the wf (ADC), NOT in keV.
            #
            # # set window, guess parameters, pack into lists
            # loWin, hiWin = dataTS[0], dataTS[-1]
            # mt, en, slo = dataTSMax, dataENF, 20.  # guesses
            # InterpFn = interpolate.interp1d(tempTS, temp, kind="linear", copy="False", assume_sorted="True")
            # dataList = [data, dataTS, dataENF, dataTSMax, loWin, hiWin, dataNoise]
            # tempList = [temp, tempTS, tempAmp, tempTSMax]
            # datas = [dataList, tempList, InterpFn]
            # floats = [mt,en,slo]
            #
            # # save the initial guess
            # guess, guessTS = wm.MakeModel(dataList, tempList, floats, fn=InterpFn)
            #
            # # run waveform fitter
            # # nelder-mead (~0.5 sec).
            # MakeTracesGlobal()
            # result = op.minimize(findLnLike, floats, args=datas, method="Nelder-Mead")#,options={"maxiter":10})
            # if not result["success"]:
            #     print "fit 'fail': ", result["message"]
            #     errorCode[0] = 1
            # mt, en, slo = result["x"]
            # nelder, nelderTS = wm.MakeModel(dataList, tempList, [mt,en,slo], fn=InterpFn)
            #
            # # save parameters.  take absval for slowness only.
            # slo = abs(slo)
            # fitMatch[iH], fitE[iH], fitSlo[iH] = mt, en, slo


            # fill the figure
            p0.cla()
            p0.plot(dataTS,data,color='black',label='data',alpha=0.4)
            p0.plot(dataTS,data_wlDenoised,color='blue',label='denoised',alpha=0.8)

            # idx = np.where((tempTS >= dataTS[0]-5) & (tempTS <= dataTS[-1]+5))
            # p0.plot(tempTS[idx],temp[idx],color='orange',label='template')
            # p0.plot(nelderTS,nelder,color='red',label='bestfit',linewidth=3)
            # p0.axvline(fitMatch[iH],color='magenta',label='fitMatch',linewidth=4,alpha=0.5)
            # p0.set_title("Run %d  Entry %d  Channel %d  ENFCal %.2f  fitMatch %.1f  fitE %.2f  fitSlo %.1f" % (run,iEvent,chan,dataENF,mt,en,slo))
            p0.legend(loc=4)
            p1.cla()
            # p1.plot(mtTrace[1:],label='fitMatch',color='red')
            # p1.legend(loc=4)
            # p1.set_xlabel('Fit Steps')
            # p2.cla()
            # p2.plot(enTrace[1:],label='fitE',color='green')
            # p2.legend(loc=4)
            # p3.cla()
            # p3.plot(sloTrace[1:],label='fitSlo',color='blue')
            # p3.legend(loc=4)
            plt.pause(scanSpeed)


def MakeSiggenWaveform(samp,r,z,ene,t0,smooth=1,phi=np.pi/8):
    """ Use pysiggen to generate a waveform w/ arb. ADC amplitude. Thanks Ben. """

    wf_length = samp # in tens of ns.  This sets how long a wf you will simulate
    # r ranges from 0 to detector.detector_radius
    # phi ranges from 0 to np.pi/4                (Clint: Limit is the pcrad, set below)
    # z ranges from 0 to detector.detector_length (Clint: Limit is the pclen, set below)
    # energy sets the waveform amplitude
    # t0 sets the start point of the waveform
    # smooth is a gaussian smoothing parameter
    # default: r, phi, z, energy, t0, smooth = (15, np.pi/8, 15, 80, 10,15)

    # Probably don't mess around with anything in this block
    fitSamples = 1000 # ask me if you think you need to change this (you almost definitely don't)
    timeStepSize = 1 # don't change this, you'll break everything
    # Create a detector model (don't change any of this)
    detName = "../data/conf/P42574A_grad%0.2f_pcrad%0.2f_pclen%0.2f.conf" % (0.05,2.5, 1.65)
    detector =  Detector(detName, timeStep=timeStepSize, numSteps=fitSamples*10./timeStepSize, maxWfOutputLength=5000)
    detector.LoadFieldsGrad("../data/fields_impgrad.npz",pcLen=1.6, pcRad=2.5)
    # Sets the impurity gradient.  Don't bother changing this
    detector.SetFieldsGradIdx(0)

    # First 3 params control the preamp shaping.
    # The 4th param is the RC decay, you can change that if you want.
    # params 5 and 6 are set not to do anything (theyre for the 2 rc constant decay model)
    rc_decay = 72.6 #us
    detector.SetTransferFunction(50, -0.814072377576, 0.82162729751, rc_decay, 1, 1)

    wf_notrap = np.copy(detector.MakeSimWaveform(r, phi, z, ene, t0, wf_length, h_smoothing=smooth))
    timesteps = np.arange(0, wf_length) * 10 # makes it so your plot is in ns

    # Add charge trapping
    # trap_constants = [8, 28,288] # these are in microseconds
    # trap_constant_colors = ["red", "blue", "purple"]
    # for (idx, trap_rc) in enumerate(trap_constants):
    #     detector.trapping_rc = trap_rc
    #     wf = np.copy(detector.MakeSimWaveform(r, phi, z, ene, t0, wf_length, h_smoothing=smooth))
    #     ax0.plot(timesteps, wf, color = trap_constant_colors[idx],  label = "%0.1f us trapping" % trap_rc )
    #     ax1.plot(timesteps, wf_notrap - wf, color = trap_constant_colors[idx])
    #     print "amplitude diff: %f" % ( (np.amax(wf_notrap) - np.amax(wf)) /  np.amax(wf_notrap) )

    return wf_notrap, timesteps


def findLnLikeXG(floats, *datas):
    """ Calculate Log-likelihood of an xGauss function for use w/ scipy minimizers. """
    global xAmpTrace, xMuTrace, xSigTrace, xTauTrace

    # Can fix parameters here by setting them back to their original guess values
    xAmp, xMu, xSig, xTau = floats
    # if True: xTau = datas[4]
    floats = [xAmp, xMu, xSig, xTau]

    # print xAmp, xMu, xSig, xTau

    # Add to traces
    xAmpTrace = np.append(xAmpTrace, xAmp)
    xMuTrace = np.append(xMuTrace, xMu)
    xSigTrace = np.append(xSigTrace, xSig)
    xTauTrace = np.append(xTauTrace, xTau)

    # Generate the xGauss model waveform
    consts = [datas[0]]
    model = MakeXGModel(consts, floats)

    # Find LL of data vs. model
    data, dataNoise = datas[5], datas[6]
    lnLike = 0.5 * np.sum ( np.power((data-model)/dataNoise, 2) - np.log( 1 / np.power(dataNoise,2) ) )
    return lnLike


def findLnLike(floats, datas):
    # SciPy Minimizer.
    #     This could be moved to waveModel.py IF you decide you
    #     don't need to look at the trace arrays anymore.
    #     Also, this version tries to align the MAX time
    #     (not the start time) of the waveforms.
    global mtTrace, enTrace, sloTrace

    # Unpack parameters.
    dataList, tempList, InterpFn = datas
    data, dataTS, dataE, dataTSMax, loWin, hiWin, dataNoise = dataList
    temp, tempTS, tempE, tempTSMax = tempList

    # Can fix parameters here by setting them back to their input guess values
    mt, en, slo = 0, 0, 0
    if len(floats)==1: # basinhopping case
        slo, mt, en = floats[0], dataTSMax, dataE
    if len(floats)==3: # minimize case (no fixing)
        mt, en, slo = floats

    # Make a trace
    mtTrace = np.append(mtTrace,mt)
    enTrace = np.append(enTrace,en)
    sloTrace = np.append(sloTrace,slo)
    # print mt, en, slo

    # Build the model to compare with data
    model, modelTS = wm.MakeModel(dataList, tempList, [mt,en,slo], fn=InterpFn)

    # Find LL of data vs. model
    lnLike = 0.5 * np.sum ( np.power((data-model)/dataNoise, 2) - np.log( 1 / np.power(dataNoise,2) ) )
    return lnLike


def MakeTracesGlobal():
    # This is so 'findLnLike' can write to the trace arrays.
    # It has to remain in this file to work.
    tmp1 = tmp2 = tmp3 = np.empty([1,])
    global mtTrace, enTrace, sloTrace
    mtTrace, enTrace, sloTrace = tmp1, tmp2, tmp3


def MakeXGTracesGlobal():
    """ This is so 'findLnLikeXG' can write to the trace arrays.
        Has to remain in this file to work. """
    tmp1 = tmp2 = tmp3 = tmp4 = np.empty([1,])
    global xAmpTrace, xMuTrace, xSigTrace, xTauTrace
    xAmpTrace, xMuTrace, xSigTrace, xTauTrace = tmp1, tmp2, tmp3, tmp4


def MakeXGModel(consts, floats):
    """ Generate a model waveform using the exponentially modified gaussian (xGauss) function. """

    # take a timestamp vector
    modelTS = consts[0]

    # generate an xGauss model and normalize it
    xAmp, xMu, xSig, xTau = floats
    model = np.zeros(len(modelTS))
    for i in range(len(modelTS)):
        mode,_ = xGaussTest(modelTS[i],xMu,xSig,xTau)
        model[i] = xGaussFn(mode,modelTS[i],1.,xMu,xSig,xTau) # set h==1, not needed
    model = model * 1./np.sum(model)

    # scale the model by setting its max ADC amplitude to "xAmp" - xGauss Energy (uncalibrated)
    xMax = np.argmax(model)
    model = model * (xAmp / model[xMax])

    return model


def gradLNLike(floats, *datas):

    xAmp, xMu, xSig, xTau = floats
    # if True: xTau = datas[4]
    floats = [xAmp, xMu, xSig, xTau]
    data, dataNoise = datas[5], datas[6]

    consts = [datas[0]]
    model = MakeXGModel(consts, floats)

    ga = (1./dataNoise) * np.sum( np.multiply(partialXGArr("A",consts,floats),(model - data) ))
    gmu = (1./dataNoise) * np.sum( np.multiply(partialXGArr("mu",consts,floats),(model - data) ))
    gsig = (1./dataNoise) * np.sum( np.multiply(partialXGArr("sig",consts,floats),(model - data) ))
    gtau = (1./dataNoise) * np.sum( np.multiply(partialXGArr("tau",consts,floats),(model - data) ))

    return np.asarray((ga, gmu, gsig, gtau))


def partialXGArr(opt,consts,floats):
    # get a numpy array of the requested partial of xGauss.

    modelTS = consts[0]
    xAmp, xMu, xSig, xTau = floats
    partial = np.zeros(len(modelTS))

    if opt != "A":
        for i in range(len(modelTS)):
            mode,_ = xGaussTest(modelTS[i],xMu,xSig,xTau)
            partial[i] = xGaussGrad(mode,opt,modelTS[i],xMu,xSig,xTau)

    elif opt=="A":
        floats = [1.,xMu,xSig,xTau] # since xAmp is just a multiplier, the partial is just the model with xAmp==1.
        partial = MakeXGModel(consts,floats)

    return partial

def xGaussGradV1(mode,opt,x,mu,sig,tau):
    """ Computed w/ Mathematica. GEEEZ. """

    gx,gmu,gsig,gtau = 0.,0.,0.,0.
    C = (1/np.sqrt(2.))*(sig/tau - (x-mu)/sig)
    D = np.sqrt(np.pi/2.)

    if mode==1:
        A = (-1./2.)(sig/tau - (x-mu)/sig)**2. + sig**2./(2*tau**2.) - (x-mu)/tau
        B = sig**2./(2*tau**2.) - (x-mu)/tau

        gx = (1./tau) * np.exp(A)
        gx += (-1./tau**2.) * (np.exp(B) * D * sig * sp.erfc(C) )

        gmu = (-1./tau) * np.exp(A)
        gmu += (1./tau**2.) * (np.exp(B) * D * sig * sp.erfc(C) )

        gsig = (-1./tau) * (np.exp(A) * sig * ((x-mu)/sig**2. + 1./tau))
        gsig += (1./tau**3.) * (np.exp(B) * D * sig**2. * sp.erfc(C) )
        gsig += (1./tau) * (np.exp(B) * D * sp.erfc(C) )

        gtau = (1./tau**3.) * (np.exp(A) * sig**2.)
        gtau += (-1./tau**2.) * (np.exp(B) * D * sig * sp.erfc(C) )
        gtau += (1./tau) * (np.exp(B) * D * sig * (-1.*sig**2./tau**3. + (x-mu)/tau**2.) * sp.erfc(C) )

    elif mode==2:
        E = (-1.*(x-mu)**2./(2*sig**2.))
        F = E + (1./2.) * (-1.*(x-mu)/sig + sig/tau)**2.

        gx = (1./tau) * (np.exp(E))
        gx += (1./tau) * (np.exp(F) * D * sig * (-1.*(x-mu)/sig**2.) - (1./sig)*(sig/tau - (x-mu)/sig) * sp.erfc(C) )

        gmu = (-1./tau) * (np.exp(E))
        gmu += (1./tau) * (np.exp(F) * D * sig * ((x-mu)/sig**2. + (1/sig)*(sig/tau - (x-mu)/sig)) * sp.erfc(C) )

        gsig = (-1./tau) * (np.exp(E) * sig * ((x-mu)/sig**2. + 1./tau) )
        gsig += (1./tau) * ( np.exp(F) * D * sp.erfc(C) )
        gsig += (1./tau) * ( np.exp(F) * D * sig * ((x-mu)**2./sig**3. + ((x-mu)/sig**2. + 1./tau) * (sig/tau - (x-mu)/sig)) * sp.erfc(C) )

        gtau = (1./tau**3.) * ( np.exp(E) * sig**2. )
        gtau += (-1./tau**3.) * ( np.exp(E) * D * sig**2. * (sig/tau - (x-mu)/sig) * sp.erfc(C) )
        gtau += (-1./tau**2.) * ( np.exp(F) * D * sp.erfc(C) )

    elif mode==3:
        G = (1 - (x-mu)*tau/sig**2.)

        gx = np.exp(E) * tau / sig**2. * G**2.
        gx -= np.exp(E) * (x-mu) / sig**2. * G

        gmu = -1.* np.exp(E) * tau / sig**2. * G**2.
        gmu += np.exp(E) * (x-mu) / sig**2. * G

        gsig = -2. * np.exp(E) * (x-mu) * tau / sig**3. * G**2.
        gsig += np.exp(E) * (x-mu)**2. / sig**3. * G

        gtau = np.exp(E) * (x-mu) / sig**2. * G**2.

    if opt=="mu":  return gmu
    if opt=="sig": return gsig
    if opt=="tau": return gtau

def xGaussGrad(mode,opt,x,mu,sig,tau):
    # A little bit faster version of the original.

    gmu, gsig, gtau = 0.,0.,0.

    # pre-compute some stuff
    A, B, C, D, E, F, G = 0.,0.,0.,0.,0.,0.,0.
    C = (1/np.sqrt(2.))*(sig/tau - (x-mu)/sig)
    D = np.sqrt(np.pi/2.)
    if mode==1:
        A = (-1./2.) * (sig/tau - (x-mu)/sig)**2. + sig**2./(2*tau**2.) - (x-mu)/tau
        B = sig**2./(2*tau**2.) - (x-mu)/tau
    elif mode==2:
        E = (-1.*(x-mu)**2./(2*sig**2.))
        F = E + (1./2.) * (-1.*(x-mu)/sig + sig/tau)**2.
    elif mode==3:
        G = (1 - (x-mu)*tau/sig**2.)

    # return the numeric value of the partial derivative at the given point
    if opt=="mu":
        if mode==1:
            gmu = (-1./tau) * np.exp(A)
            gmu += (1./tau**2.) * (np.exp(B) * D * sig * sp.erfc(C) )
        elif mode==2:
            gmu = (-1./tau) * (np.exp(E))
            gmu += (1./tau) * (np.exp(F) * D * sig * ((x-mu)/sig**2. + (1/sig)*(sig/tau - (x-mu)/sig)) * sp.erfc(C) )
        elif mode==3:
            gmu = -1.* np.exp(E) * tau / sig**2. * G**2.
            gmu += np.exp(E) * (x-mu) / sig**2. * G
        return gmu

    elif opt=="sig":
        if mode==1:
            gsig = (-1./tau) * (np.exp(A) * sig * ((x-mu)/sig**2. + 1./tau))
            gsig += (1./tau**3.) * (np.exp(B) * D * sig**2. * sp.erfc(C) )
            gsig += (1./tau) * (np.exp(B) * D * sp.erfc(C) )
        elif mode==2:
            gsig = (-1./tau) * (np.exp(E) * sig * ((x-mu)/sig**2. + 1./tau) )
            gsig += (1./tau) * ( np.exp(F) * D * sp.erfc(C) )
            gsig += (1./tau) * ( np.exp(F) * D * sig * ((x-mu)**2./sig**3. + ((x-mu)/sig**2. + 1./tau) * (sig/tau - (x-mu)/sig)) * sp.erfc(C) )
        elif mode==3:
            gsig = -2. * np.exp(E) * (x-mu) * tau / sig**3. * G**2.
            gsig += np.exp(E) * (x-mu)**2. / sig**3. * G
        return gsig

    elif opt=="tau":
        if mode==1:
            gtau = (1./tau**3.) * (np.exp(A) * sig**2.)
            gtau += (-1./tau**2.) * (np.exp(B) * D * sig * sp.erfc(C) )
            gtau += (1./tau) * (np.exp(B) * D * sig * (-1.*sig**2./tau**3. + (x-mu)/tau**2.) * sp.erfc(C) )
        elif mode==2:
            gtau = (1./tau**3.) * ( np.exp(E) * sig**2. )
            gtau += (-1./tau**3.) * ( np.exp(E) * D * sig**2. * (sig/tau - (x-mu)/sig) * sp.erfc(C) )
            gtau += (-1./tau**2.) * ( np.exp(F) * D * sp.erfc(C) )
        elif mode==3:
            gtau = np.exp(E) * (x-mu) / sig**2. * G**2.
        return gtau


def MakeModel(dataList, tempList, params, fn=0, opt=""):
    """ Generate a model waveform from a template and some parameters.

        NOTE: dataSync and tempSync are waveform times we are trying to sync up.
        They don't necessarily refer to the start time or max time of a pulse,
        unless that is specified in a function calling this one.
    """
    # print "entering MakeModel ..."

    # Unpack inputs
    data, dataTS, dataE, dataSync = dataList[0], dataList[1], dataList[2], dataList[3]
    temp, tempTS, tempE, tempSync = tempList[0], tempList[1], tempList[2], tempList[3]
    st, en, slo = params[0], params[1], params[2]

    # Float the template's 'sync' time (can think of it as a start or max time)
    deltaSync = st - tempSync

    # Break the shift of template TS's into two parts:  "unit" and "remainder".
    remShift = 0  #  amount of shift less than 10ns
    if deltaSync < 0:
        remShift = deltaSync % -10.
    else:
        remShift = deltaSync % 10
    unitShift = deltaSync - remShift  # amount of shift greater than 10ns

    # Declare initial model
    model = temp
    modelTS = tempTS + deltaSync

    # Return a guess without windowing
    if opt=="full":
        idxFull = np.where((modelTS >= tempTS[0]) & (modelTS <= tempTS[-1]))
        modelTS = modelTS[idxFull]
        model = temp[idxFull] * (en / tempE)
        return model, modelTS

    # Return a guess without windowing, but with scaling and smoothing
    if opt=="nowindow":
        idxFull = np.where((modelTS >= tempTS[0]) & (modelTS <= tempTS[-1]))
        modelTS = modelTS[idxFull]
        model = temp[idxFull] * (en / tempE)
        model = gaussian_filter(model,sigma=float( slo ))
        return model, modelTS

    # Shift template timestamps and interpolate template ADC values
    if fn!=0:
        tempRemTS = tempTS + remShift   # shift all template TS's by the remainder
        idx = np.where((tempRemTS >= tempTS[0]) & (tempRemTS <= tempTS[-1]))
        tempRemTS = tempRemTS[idx]  # only take indexes which fall within the template's original bounds
        tempInterp = fn(tempRemTS)  # get ADC values corresponding to shifting the template by only the remainder
        tempShiftedTS = tempRemTS + unitShift  # shift timestamps by the "unit shift" (units of 10ns)
        model = tempInterp
        modelTS = tempShiftedTS

    # Window the model
    if dataTS[0] < modelTS[0] or dataTS[-1] > modelTS[-1]:
        print "Model floated out the window.  st %d  tempSync %d  dST %d  loData %d  loModel %d  hiData %d  hiModel %d" % (st,tempSync,deltaSync,dataTS[0],modelTS[0],dataTS[-1],modelTS[-1]) # commented this warning out for LAT
        return np.ones(len(data)),dataTS
    idxFirst = (np.abs(modelTS - dataTS[0])).argmin()
    idxLast = idxFirst + len(dataTS)
    modelTS = modelTS[idxFirst:idxLast]
    model = model[idxFirst:idxLast]

    # Return a guess with windowing & interpolation but nothing else
    if opt=="!fancy":
        model = model * (en / tempE)
        return model, modelTS

    # Float the energy
    model = model * (en / tempE)

    # Float the smoothing
    model = gaussian_filter(model,sigma=float( slo ))

    # Let's make sure modelTS and dataTS have same number of entries ALWAYS.
    if len(modelTS)!=len(dataTS):
        print "array TS mismatch: model %d  data %d  m0 %.0f  m-1 %.0f  d0 %.0f  d-1 %.0f  dST %.0f  tST %.0f  st %.0f" % (len(modelTS),len(dataTS),modelTS[0],modelTS[-1],dataTS[0],dataTS[-1],deltaST,tempST,st)
        return np.ones(len(dataTS)),dataTS

    return model, modelTS


if __name__ == "__main__":
    main(sys.argv[1:])