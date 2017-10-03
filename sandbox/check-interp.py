#!/usr/bin/env python
import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt

def main():

    # Load a data and a temp waveform
    npzfile = np.load("./data/mcmcInputs.npz")
    dl, tl = npzfile['arr_0'], npzfile['arr_1']
    data, dataTS, dataE, dataRT = dl[0], dl[1], dl[2], dl[3]
    temp, tempTS, tempE, tempRT = tl[0], tl[1], tl[2], tl[3]

    # Create the interp function for the template
    fn = interpolate.interp1d(tempTS, temp, kind="linear", copy="False", assume_sorted="True")

    # Set some guess parameters
    rt, en, slo = dataRT-105, dataE, 5

    # Float the absolute rise time ("rt").
    deltaRT = rt - tempRT       # will usually be negative

    # Break the shift of template TS's into two parts:  "unit" and "remainder".
    remShift = 0  #  amount of shift less than 10ns
    if deltaRT < 0: remShift = deltaRT % -10.
    else: remShift = deltaRT % 10
    unitShift = deltaRT - remShift  # amount of shift greater than 10ns

    # Shift template timestamps and interpolate ADC values
    tempRemTS = tempTS + remShift   # shift all template TS's by the remainder
    idx = np.where((tempRemTS >= tempTS[0]) & (tempRemTS <= tempTS[-1]))
    tempRemTS = tempRemTS[idx]  # only take indexes which fall within the template's original bounds
    tempInterp = fn(tempRemTS)  # get ADC values corresponding to shifting the template by only the remainder
    tempShiftedTS = tempRemTS + unitShift  # shift timestamps by the "unit shift - units of 10ns"

    # NOTE: plot tempRemTS, tempInterp : template shifted only by remainder, w/ correct interpolation
    #       plot tempShiftedTS, tempInterp : template shifted by full amount, w/ interpolation correcting for the remainder

    # Window the model
    model = tempInterp
    modelTS = tempShiftedTS
    if dataTS[0] < modelTS[0] or dataTS[-1] > modelTS[-1]:
        print "Model floated out the window.  loData %d  loModel %d  hiData %d  hiModel %d" % (dataTS[0],modelTS[0],dataTS[-1],modelTS[-1])
        return np.ones(len(data)),dataTS
    idxFirst = (np.abs(modelTS - dataTS[0])).argmin()
    idxLast = idxFirst + len(dataTS)
    modelTS = modelTS[idxFirst:idxLast]
    model = model[idxFirst:idxLast]

    # plotz
    fig = plt.figure(figsize=(7,7), facecolor='w')

    a1 = plt.subplot(211)
    a1.set_xlabel("Time (ns)")
    a1.set_title("RT Template (ns) %.1f  RT Model %.1f  Delta %.1f  Remainder %.1f" % (tempRT,rt,deltaRT,remShift))
    a1.plot(tempTS,temp,color='blue',label='Unshifted Template')
    a1.plot(modelTS,model,color='red',label='Shifted, Windowed')
    a1.legend(loc=4)

    a2 = plt.subplot(212)
    a2.set_xlabel("Time (ns)")
    idx1 = np.where((tempTS >= 10300) & (tempTS <= 10500))
    idx2 = np.where((tempRemTS >= 10300) & (tempRemTS <= 10500))
    a2.plot(tempTS[idx1],temp[idx1],'o',color='blue',markersize=1,label='Unshifted Template')
    a2.plot(tempRemTS[idx2],tempInterp[idx2],'o',color='red',markersize=1,label="Shifted by Remainder")
    a2.axvline(x = 10400, color='green',linewidth=0.5,alpha=0.5)    # draw a marker on an unshifted TS
    a2.legend(loc=4)

    plt.tight_layout()
    plt.subplots_adjust(hspace=0.35)
    plt.show(block=False)
    plt.savefig("./plots/check-interp.pdf")


def version2():
    """ This mostly works but seems really bulky.  Also,
    the np.wheres make it hard to keep the model the same size as the data.
    """

    # Load a data and a temp waveform
    npzfile = np.load("./data/wfArrays.npz")
    wave,waveTS,temp,tempTS,tDiffGuess,scaleGuess = npzfile['a'],npzfile['b'],npzfile['c'],npzfile['d'],npzfile['e'],npzfile['f']

    # Create the interp function for the temp
    InterpFn = interpolate.interp1d(tempTS, temp, kind="linear", copy="False", assume_sorted="True")

    # Start with the total amount we want to shift by (mcmc parameter)
    bigOffset = -2005.678

    # Get the remainder of the shift that's smaller than 10ns
    remShift = bigOffset % 10
    print "remShift: ",remShift

    # Shift all timestamps by the remainder
    # Only take values of temp TS's within the original limits (otherwise interpolation will fail)
    tempShiftedTS = tempTS + remShift
    idxShifted = np.where((tempShiftedTS >= tempTS[0]) & (tempShiftedTS <= tempTS[-1]))
    tempShiftedTS = tempShiftedTS[idxShifted]

    # ** Make a new interpolated set of ADC values
    tempInterp = InterpFn(tempShiftedTS)

    # ** Finally, shift the timestamps by the original offset.
    tempShiftedFullTS = tempShiftedTS + bigOffset


    # plotz
    fig = plt.figure(figsize=(9,7), facecolor='w')

    # Zoom in on the remainder shift
    a1 = plt.subplot(211)
    idx2 = np.where((tempTS >= 10500) & (tempTS <= 11000))
    a1.plot(tempTS[idx2],temp[idx2],'o',markersize=1,color='blue')    # should stay in place

    idx3 = np.where((tempShiftedTS >= 10500) & (tempShiftedTS <= 11000))
    a1.plot(tempShiftedTS[idx3],tempInterp[idx3],'o',markersize=1,color='red') # should be shifted by the remainder

    a1.axvline(x = 10700, color='green',linewidth=0.5,alpha=0.5)  # draw a marker on an unshifted TS

    # Zoom out for the overall shift
    a2 = plt.subplot(212)
    a2.plot(tempTS, temp, color='blue')  # should stay in place

    a2.plot(tempShiftedFullTS, tempInterp, color='red') # should be shifted by the "big offset"

    plt.show()


def version1():
    """This doesn't actually shift the interpolated waveform.
    Kept it in case I get too confused in main(). """
    # Load a data and a temp waveform
    npzfile = np.load("./data/wfArrays.npz")
    wave,waveTS,temp,tempTS,tDiffGuess,scaleGuess = npzfile['a'],npzfile['b'],npzfile['c'],npzfile['d'],npzfile['e'],npzfile['f']

    InterpFn = interpolate.interp1d(tempTS, temp, kind="linear", copy="False", assume_sorted="True")

    rt = 1000
    tempTSRT = tempTS + rt

    idx = np.where((tempTSRT >= tempTS[0]) & (tempTSRT <= tempTS[-1]))
    newTemp = temp[idx]
    newTempTS = tempTSRT[idx]

    interpTemp = InterpFn(newTempTS)

    print "temp  ",newTemp[:10]
    print "old time ",tempTS[idx][:10]
    print " "
    print "adcInterp",interpTemp[:10]
    print "timestamp",newTempTS[:10]

    # manually check interp
    def linInterp(x,x0,x1,y0,y1):
        return y0 + (x-x0)*(y1-y0)/(x1-x0)
    manual = [linInterp(newTempTS[i],tempTS[idx][i],tempTS[idx][i+1],newTemp[i],newTemp[i+1]) for i in range(10)]
    print('\n'.join('{}: {}'.format(*k) for k in enumerate(manual)))

    # idx2 = np.where((newTempTS >= 9000) & (newTempTS <= 13000))
    fig = plt.figure(figsize=(9,5), facecolor='w')
    a1 = plt.subplot(111)
    a1.plot(newTempTS,newTemp,color='red')
    a1.plot(newTempTS,interpTemp,color='blue')
    plt.show()


if __name__ == "__main__":
    main()
