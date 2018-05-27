#!/usr/bin/env python

import sys, random
import numpy as np
import matplotlib.pyplot as plt
plt.style.use('../pltReports.mplstyle')

from scipy.ndimage.filters import gaussian_filter
import waveLibs as wl

def main(argv):

    intMode, batMode = False, False
    for i,opt in enumerate(argv):
        if opt == "-i":
            intMode = True
            print("Interactive mode selected.")
        if opt == "-b":
            batMode = True
            print("Bach mode selected.")

    print("Generating signal template ...")
    tSamp, tR, tZ, tAmp, tST, tSlo = 5000, 0, 15, 100, 2500, 10
    # tOrig, tOrigTS = wl.MakeSiggenWaveform(tSamp,tR,tZ,tAmp,tST,tSlo)
    templateFile = np.load("../data/lat_template.npz")
    tOrig, tOrigTS = templateFile['arr_0'], templateFile['arr_1']

    fig = plt.figure(figsize=(10,7),facecolor='w')
    p1 = plt.subplot(211)
    p2 = plt.subplot(212)
    if not batMode: plt.show(block=False)

    # we're going to histogram these
    amps1, amps2, amps3 = [], [], []
    starts1, starts2, starts3 = [], [], []

    print("Starting event loop ...")
    nList = 2000 # how many wf's to generate
    # nList = 100 # how many wf's to generate
    iList = -1
    while True:
        iList += 1
        if intMode==True and iList != 0:
            value = input()
            if value=='q': break        # quit
            if value=='p': iList -= 2   # go to previous
            if (value.isdigit()):
                iList = int(value)      # go to entry number
        elif intMode==False and batMode==False:
            plt.pause(0.00001)          # rapid-draw mode
        if iList >= nList: break        # bail out, goose!

        print("%d/%d" % (iList, nList))


        # make a somewhat convincing template

        temp, tempTS = tOrig, tOrigTS

        maxADCRand = random.uniform(0.1,10) # flat is probably ok here
        temp = temp * (maxADCRand / tAmp)

        slowness = random.uniform(0,20) # should this be gaussian instead of flat?
        temp = gaussian_filter(temp,sigma=float( slowness ))

        # 0 is the mean of the normal distribution you are choosing from
        # 1 is the standard deviation of the normal distribution
        # 100 is the number of elements you get in array noise
        noise = np.random.normal(0, 2, len(tOrig))
        temp += noise




        # trap filter in waveLibs
        # trap = np.zeros(1000)
        # trap = np.append(trap,wl.trapezoidalFilter(temp,rampTime=400, flatTime=200, decayTime=0.))
        # trapMax = np.amax(trap)

        ADCThresh = 2.

        # simple trap filter
        trapTS = tempTS
        trap, trigger, triggerTS = trapFilt(temp,trapThresh=ADCThresh)

        shortTrap, shortTrigger, shortTriggerTS = trapFilt(temp,100,150,0,trapThresh=ADCThresh,padAfter=True)
        walkBackTS,_ = walkBack(shortTrap,thresh=ADCThresh)

        if not trigger:
            amps1.append(maxADCRand)
            starts1.append(triggerTS)
            amps2.append(maxADCRand)
            starts2.append(walkBackTS)
            continue
        amps1.append(maxADCRand)
        starts1.append(triggerTS)
        amps2.append(maxADCRand)
        starts2.append(walkBackTS)

        # make a windowed waveform, the same windowing ORCA uses
        lo = int(triggerTS/10 - 1009)
        hi = int(triggerTS/10 + 1009)
        trig = temp[lo:hi]
        trigTS = tempTS[lo:hi]

        if len(trigTS)==0: continue
        trigTS = trigTS - trigTS[0]

        # asymmetric trap filter
        aTrap, _, _ = asymTrapFilt(temp,trapThresh=ADCThresh)

        # apply a short trap to the WINDOWED waveform
        # start from the max and walk BACK to the threshold AFTER calculating the trap
        windTrap, windTrigger, windTriggerTS = trapFilt(trig,100,150,0,trapThresh=ADCThresh,padAfter=True)
        walkBackTS2,_ = walkBack(windTrap,thresh=ADCThresh)
        if not windTrigger:
            amps3.append(maxADCRand)
            starts3.append(walkBackTS2)
            continue
        amps3.append(maxADCRand)
        starts3.append(walkBackTS2)


        # plooooots
        if not batMode:
            p1.cla()
            p1.plot(tempTS,temp,'b')
            p1.plot(trapTS,trap,'g')
            p1.plot(trapTS,shortTrap,color='red')
            p1.plot(trapTS,aTrap,color='orange')
            p1.axvline(triggerTS,color='green')
            p1.axvline(walkBackTS,color='red')
            p1.set_title("Long Waveform, 50us.  TriggerTS %d" % triggerTS)
            p1.set_xlabel("Time [ns]")
            p1.set_ylabel("ADC")

            p2.cla()
            p2.plot(trigTS,trig,'g')
            p2.plot(trigTS,windTrap,color='magenta')
            p2.axvline(walkBackTS2,color='magenta')
            p2.set_title("Windowed Waveform, 20us.")
            p2.set_xlabel("Time [ns]")
            p2.set_ylabel("ADC")


            ax = plt.gca()
            plt.tight_layout()
            plt.pause(0.000001)


    # Make a scatter plot
    if not intMode:
        fig2 = plt.figure()
        # plt.scatter(amps1,starts1,color='green',s=5,label='long trap - walk up')
        # plt.scatter(amps2,starts2,color='red',s=5,label='short trap - walk back')
        # plt.scatter(amps3,starts3,color='magenta',s=5,alpha=0.8,label='short trap, windowed, walk back')
        plt.scatter(amps3,starts3,color='blue',s=5,alpha=0.8) # thesis plot
        plt.xlabel("Energy (ADC)", ha='right', x=1)
        plt.ylabel("Trigger Time (ns)", ha='right', y=1)
        # plt.legend(loc=4)
        fig2.savefig('../plots/trig-walk.pdf')


def trapFilt(data,ramp=400,flat=200,decay=72,trapThresh=2.,padAfter=False):

    # decay is unused
    # no charge trapping
    # no fixed time pickoff
    # no pole zero correction

    trap = np.zeros(len(data))
    trigger, triggerTS = False, 0

    # compute a moving average
    for i in range(len(data)-1000):
        w1 = ramp
        w2 = ramp+flat
        w3 = ramp+flat+ramp

        # r1 = np.sum(data[i:w1+i])/ramp
        # r2 = np.sum(data[w2+i:w3+i])/ramp
        r1 = np.sum(data[i:w1+i])/100.
        r2 = np.sum(data[w2+i:w3+i])/100.

        if not padAfter:
            trap[i+1000] = r2 - r1
            if trap[i+1000] > trapThresh and not trigger:
                triggerTS = (i + 1000) * 10
                trigger = True
        else:
            trap[i] = r2 - r1
            if trap[i] > trapThresh and not trigger:
                triggerTS = i * 10
                trigger = True

    # if trigger: print("triggered!",triggerTS)
    return trap, trigger, triggerTS

def asymTrapFilt(data,ramp=200,flat=100,fall=40,trapThresh=2.,padAfter=False):

    trap = np.zeros(len(data))
    trigger, triggerTS = False, 0

    # compute a moving average
    for i in range(len(data)-1000):
        w1 = ramp
        w2 = ramp+flat
        w3 = ramp+flat+fall

        r1 = np.sum(data[i:w1+i])/(ramp/10)
        r2 = np.sum(data[w2+i:w3+i])/(fall/10)

        if not padAfter:
            trap[i+1000] = r2 - r1
            if trap[i+1000] > trapThresh and not trigger:
                triggerTS = (i + 1000) * 10
                trigger = True
        else:
            trap[i] = r2 - r1
            if trap[i] > trapThresh and not trigger:
                triggerTS = i * 10
                trigger = True

    # if trigger: print("triggered!",triggerTS)
    return trap, trigger, triggerTS


def walkBack(trap,thresh=2.):

    trapMax = np.argmax(trap)

    foundFirst, triggerTS = False, 0
    for i in range(trapMax,0,-1):
        if trap[i] <= thresh:
            foundFirst = True
            triggerTS = (i+1) * 10
            break

    return triggerTS, foundFirst



if __name__ == "__main__":
    main(sys.argv[1:])