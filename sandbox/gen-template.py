#!/usr/local/bin/python
from pysiggen import Detector
import numpy as np
import matplotlib.pyplot as plt

def main():

    # Generate a siggen WF
    samp, r, z, ene, t0 = 2016, 0, 15, 3938, 1000
    wave, waveTS = MakeSiggenWaveform(samp,r,z,ene,t0)

    # plots
    # fig = plt.figure(figsize=(8,5), facecolor='w')
    # a1 = plt.subplot(111)
    # a1.set_ylabel("ADC",y=0.95, ha='right')
    # a1.set_xlabel("Time (ns)",x=0.95, ha='right')
    # a1.set_title("Samples %d  r %.1f  z %.1f  ADC %d  startTime %.0f" % (samp,r,z,ene,t0*10))
    # a1.plot(waveTS,wave,color='green',label='Generated WF')
    # a1.axvline(x=t0*10,color='red',label='StartTime')
    # a1.legend(loc=4)
    # plt.show()
    # plt.savefig("./plots/gen-template.pdf")

    # save the template
    # np.savez("./data/genTemplateWF.npz",wave,waveTS,1592.0,t0)

    # LAT Template:
    # templateFile = np.load("./data/lat_template.npz")
    # tOrig, tOrigTS = templateFile['arr_0'], templateFile['arr_1']

    tSamp, tR, tZ, tAmp, tST, tSlo = 5000, 0, 15, 100, 2500, 10
    tOrig, tOrigTS = MakeSiggenWaveform(tSamp, tR, tZ, tAmp, tST, tSlo)
    np.savez("./data/lat_template.npz",tOrig,tOrigTS)


    # TestScaling()
    # FitTemplateEnergy()

def TestScaling():
    """ Worried that scaling a high-e waveform down to low
    energy would give a different tail shape than just generating
    a low energy waveform.  Fortunately, the code below shows
    this is not a problem. """

    samp, r, z, ene, t0 = 2016, 0, 15, 3500, 1000
    wave3500, waveTS = MakeSiggenWaveform(samp,r,z,ene,t0)
    wave3500 = wave3500 * 2 / 3500

    ene = 2
    wave2,_ = MakeSiggenWaveform(samp,r,z,ene,t0)

    fig = plt.figure(figsize=(8,5), facecolor='w')
    a1 = plt.subplot(111)
    a1.set_ylabel("ADC",y=0.95, ha='right')
    a1.set_xlabel("Time (ns)",x=0.95, ha='right')
    a1.plot(waveTS,wave3500,color='green',label='Scaled 3500')
    a1.plot(waveTS,wave2,color='red',label='Unscaled 2')
    a1.legend(loc=4)
    plt.show()


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
    detName = "./data/conf/P42574A_grad%0.2f_pcrad%0.2f_pclen%0.2f.conf" % (0.05,2.5, 1.65)
    detector =  Detector(detName, timeStep=timeStepSize, numSteps=fitSamples*10./timeStepSize, maxWfOutputLength=5000)
    detector.LoadFieldsGrad("./data/fields_impgrad.npz",pcLen=1.6, pcRad=2.5)
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


def FitTemplateEnergy():
    """ The generated waveform is in arbitrary units.
        I want it to be in energy units so that fitting the scale
        of the low energy data waveforms makes sense.

    Idea 1: Use MCMC to fit the scale of the template to one of the
            DEP waveforms, and adjust the template by this amount,
            creating a "fake" DEP waveform.
    Idea 2: Use MGDO/GAT etc to calculate "trapENFCal" for this template
            waveform.

    Heey, why am I doing this?  I can't even be sure that I can trust the
    fit-out scaled down energy.  So let's put this on the back burner.
    Just generate an arbitrary waveform and go with it.
    """

if __name__ == "__main__":
    main()
