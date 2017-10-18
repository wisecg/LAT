#!/usr/bin/env python
import pywt
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import waveLibs as wl
import waveModel as wm

def main():
    """ Trying to prove to myself what the limits on the start time are.
    OFFICIAL START TIME LIMITS:
        Lowest allowed:  startTime > tempST - tempTS[-1] + hiWin
        Highest allowed: startTime < tempST - tempTS[0] + loWin
    """
    # Load a data and a temp waveform
    npzfile = np.load("./data/tailSlopeInputs.npz")
    dl, tl = npzfile['arr_0'], npzfile['arr_1']
    data, dataTS, dataE, dataST = dl[0], dl[1], dl[2], dl[3]
    temp, tempTS, tempE, tempST = tl[0], tl[1], tl[2], tl[3]


    fig = plt.figure(figsize=(8,5), facecolor='w')
    a1 = plt.subplot(111)
    a1.set_ylabel("ADC",y=0.95, ha='right')
    a1.set_xlabel("Time (ns)",x=0.95, ha='right')

    loWin, hiWin = dataST - 1000, dataST + 4000 # ns
    a1.axvline(x=loWin,color='black')
    a1.axvline(x=hiWin,color='black')
    a1.axvline(x=dataST,color='cyan')

    # original start time
    startTime = tempST
    a1.plot(tempTS,temp,color='green')
    a1.axvline(x=startTime,color='green')

    # lowest possible start time value
    startTime = tempST - tempTS[-1] + hiWin
    deltaST = startTime - tempST # this goes negative
    shiftTS = tempTS + deltaST
    a1.plot(shiftTS,temp,color='blue')
    a1.axvline(x=startTime,color='blue')

    # highest possible start time value
    startTime = tempST - tempTS[0] + loWin
    deltaST = startTime - tempST # this goes positive
    shiftTS = tempTS + deltaST
    a1.set_title("loWin %d  hiWin %d  red shiftTS[0] %d  %d  %d" % (loWin,hiWin,shiftTS[0], loWin-shiftTS[0], tempTS[0]) )
    a1.plot(shiftTS,temp,color='red')
    a1.axvline(x=startTime,color='red')




    # a1.legend(loc=4)
    plt.show()
    # plt.savefig("./plots/tailSlope.pdf")

if __name__ == "__main__":
    main()
