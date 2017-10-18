#!/usr/bin/env python
import sys
import numpy as np
# from pymc import TruncatedNormal, Normal, HalfNormal, deterministic, Uniform
from scipy.ndimage.filters import gaussian_filter
from scipy import interpolate

"""
# This function is perfectly fine, but I can't get pymc to work on PDSF.  Fuck you, PDSF.
def TemplateModel(dataList, dataNoise, tempList):
    # PyMC Minimizer.

    # NOTE: do you still want to refer to it as "startTime", given what you learned below?

    # Unpack inputs
    data, dataE, dataST, loWin, hiWin = dataList[0], dataList[2], dataList[3], dataList[4], dataList[5]
    temp, tempTS, tempST = tempList[0], tempList[1], tempList[3]

    # Find limits on the start time for this template and this data waveform
    minST = tempST - tempTS[-1] + hiWin
    maxST = tempST - tempTS[0] + loWin

    # Set up stochastic variables
    startTime = TruncatedNormal('startTime', mu=dataST, tau=np.power(100.,-2),a=minST,b=maxST)
    energy = TruncatedNormal('energy',mu=dataE, tau=np.power(1.,-2), a=-10, b=3000) # let it go negative to spot bad events?
    slowness = TruncatedNormal('slowness', mu=1., tau=np.power(150.,-2), a=0, b=200)

    # Set up ADC interpolation function
    InterpFn = interpolate.interp1d(tempTS, temp, kind="linear", copy="False", assume_sorted="True")

    # Deterministic variable is the signal model
    @deterministic(plot=False, name="TemplateModel")
    def SignalModel(rt=startTime, en=energy, slo=slowness): #, blo=blOffset): #, sl1=tailSlope1): #, sl2=tailSlope2):
        params = [rt,en,slo] #,blo] #,sl1]
        model,_ = MakeModel(dataList, tempList, params, fn=InterpFn)
        return model

    # Final model for waveform
    baseline_observed = Normal('baseline_observed', mu=SignalModel, tau=np.power(dataNoise,-2), value=data, observed=True)

    # Only return objects pymc can recognize (don't just do 'return locals()')
    # Model definition must be in terms of Stochastics, Deterministics, Potentials and Containers.
    outDict = {'SignalModel':SignalModel, 'startTime':startTime, 'energy':energy, 'slowness':slowness, 'baseline_observed':baseline_observed}
    return outDict
"""

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
