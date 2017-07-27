#!/usr/local/bin/python
import matplotlib.pyplot as plt
import numpy as np
import scipy.special as sp
from pysiggen import Detector

def main():

    # timestamps
    start, stop, binsPerNS = 0, 20000, 10
    tempTS = np.arange(start, stop, binsPerNS)

    # generate an xgauss wf
    h, mu, sig, tau = 10000., 10000., 200., 72000.
    temp = np.zeros(len(tempTS))
    for i in range(len(tempTS)):
        mode,z = xGaussTest(tempTS[i],mu,sig,tau)
        temp[i] = xGauss(mode,tempTS[i],h,mu,sig,tau)
        if i < 50 or i > len(tempTS)-50: print tempTS[i], mode, z, temp[i]


    # generate a pysiggen wf
    tSamp, tR, tZ, tAmp, tST, tSlo = 2000, 0, 15, 100, 1000, 10
    tOrig, tOrigTS = MakeSiggenWaveform(tSamp, tR, tZ, tAmp, tST, tSlo)
    # result, resultTS = np.zeros(10000), np.zeros(10000)
    # result[5000:] = tOrig
    # resultTS = np.arange(0, 10000) * 10 # ts in ns

    # match


    # plots
    fig = plt.figure(figsize=(8,6),facecolor='w')
    plt.plot(tempTS,temp,color='blue')
    plt.plot(tOrigTS,tOrig,color='red')
    plt.show()

def xGaussTest(x,mu,sig,tau):
    z = 1/np.sqrt(2) * (sig/tau - (x - mu)/sig)
    if z < 0: return 1, z
    elif z >= 0 and z <= 6.71e7: return 2, z
    elif z > 6.71e7: return 3, z
    else: return -1, z

def xGauss(mode,x,h,mu,sig,tau):
    if mode==1:
        return ( h*sig*np.sqrt(np.pi/2.)/tau ) * np.exp( (sig/tau)**2./2. - (x-mu)/tau ) * sp.erfc( ( sig/tau - (x-mu)/sig )/np.sqrt(2.) )
    elif mode==2:
        return ( h*sig*np.sqrt(np.pi/2.)/tau ) * np.exp( (-1./2.)*((x-mu)/sig)**2. ) * sp.erfcx( ( sig/tau - (x-mu)/sig )/np.sqrt(2.) )
    elif mode==3:
        return ( h/(1 - (x-mu)*tau/sig**2.) ) * np.exp( (-1./2.)*((x-mu)/sig)**2. )
    else:
        print "unknown mode!"
        return -1

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



if __name__ == "__main__":
    main()