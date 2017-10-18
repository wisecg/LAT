#!/usr/bin/env python
from pysiggen import Detector
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec

def plot_waveforms():

    wf_length = 200 #in tens of ns.  This sets how long a wf you will simulate

    #Probably don't mess around with anything in this block
    fitSamples = 1000 #ask me if you think you need to change this (you almost definitely don't)
    timeStepSize = 1 #don't change this, you'll break everything
    #Create a detector model (don't change any of this)
    detName = "conf/P42574A_grad%0.2f_pcrad%0.2f_pclen%0.2f.conf" % (0.05,2.5, 1.65)
    detector =  Detector(detName, timeStep=timeStepSize, numSteps=fitSamples*10./timeStepSize, maxWfOutputLength=5000)
    detector.LoadFieldsGrad("fields_impgrad.npz",pcLen=1.6, pcRad=2.5)
    #sets the impurity gradient.  Don't bother changing this
    detector.SetFieldsGradIdx(0)


    #First 3 params control the preamp shaping.
    #The 4th param is the RC decay, you can change that if you want.
    #params 5 and 6 are set not to do anything (theyre for the 2 rc constant decay model)
    rc_decay = 72.6 #us
    detector.SetTransferFunction(50, -0.814072377576, 0.82162729751, rc_decay, 1, 1)

    #r ranges from 0 to detector.detector_radius
    #phi ranges from 0 to np.pi/4
    #z ranges from 0 to detector.detector_length
    #energy sets the waveform amplitude
    #t0 sets the start point of the waveform
    #smooth is a gaussian smoothing parameter

    r, phi,z, energy, t0, smooth = (15, np.pi/8,15, 80, 10,15)

    trap_constants = [8, 28,288] #these are in microseconds
    trap_constant_colors = ["red", "blue", "purple"]

    # fig1 = plt.figure(0, figsize=(20,10))
    fig1 = plt.figure(0, figsize=(9,7))
    gs = gridspec.GridSpec(2, 1, height_ratios=[4, 1])
    ax0 = plt.subplot(gs[0])
    ax1 = plt.subplot(gs[1], sharex=ax0)
    ax1.set_xlabel("Digitizer Time [ns]")
    ax0.set_ylabel("Voltage [A.U.]")
    ax1.set_ylabel("Residual [A.U]")


    timesteps = np.arange(0, wf_length)*10 #makes it so your plot is in ns


    wf_notrap = np.copy(detector.MakeSimWaveform(r, phi, z, energy, t0, wf_length, h_smoothing=smooth))
    ax0.plot(timesteps, wf_notrap,  color="green", label = "No trapping" )

    for (idx, trap_rc) in enumerate(trap_constants):
        detector.trapping_rc = trap_rc
        wf = np.copy(detector.MakeSimWaveform(r, phi, z, energy, t0, wf_length, h_smoothing=smooth))
        ax0.plot(timesteps, wf, color = trap_constant_colors[idx],  label = "%0.1f us trapping" % trap_rc )
        ax1.plot(timesteps, wf_notrap - wf, color = trap_constant_colors[idx])

        print "amplitude diff: %f" % ( (np.amax(wf_notrap) - np.amax(wf)) /  np.amax(wf_notrap) )

    ax0.legend(loc=4)
    plt.show()
    plt.savefig("wf.pdf")


if __name__=="__main__":
    plot_waveforms()
