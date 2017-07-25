#!/usr/local/bin/python
import numpy as np
from scipy import fft, ifft, interpolate
from scipy.optimize import minimize
import matplotlib.pyplot as plt
import pymc
import waveLibs as wl
import waveModel as wm

def main():
    """ Get a cool SciPy minimizer working and compare to MCMC. """

    # Load data and template
    npzfile = np.load("./data/optimumInputs.npz")
    rl, tl = npzfile['arr_0'], npzfile['arr_1']
    wave, waveTS, dataE, dataST = rl[0], rl[1], rl[2], rl[3]
    temp, tempTS, tempE, tempST = tl[0], tl[1], tl[2], tl[3]

    # Window the fit around rising edge - start time calculator method
    loWin, hiWin = dataST - 1000, dataST + 4000 # ns
    if loWin < waveTS[0] or hiWin > waveTS[-1]:
        print "Window out of range!  dataST: %.1f  loWin %.1f  hiWin %.1f" % (dataST,loWin,hiWin)
    idx = np.where((waveTS >= loWin) & (waveTS <= hiWin))
    data = wave[idx]
    dataTS = waveTS[idx]

    # Pack into lists
    dataNoise = 2.  # just a guess - 1 sigma baseline adc values
    rawList = [wave, waveTS, dataE, dataST]
    dataList = [data, dataTS, dataE, dataST, loWin, hiWin, dataNoise]
    tempList = [temp, tempTS, tempE, tempST]

    # Recreate the guess and the guess's rising edge
    guessFull, guessFullTS = wm.MakeModel(rawList, tempList, [dataST,dataE,1.], opt="full")
    guess, guessTS = wm.MakeModel(dataList, tempList, [dataST,dataE,1.], opt="!fancy")

    # Make an "almost complete" guess - no MCMC
    # st, en, slo = dataST-100, dataE, 5
    # InterpFn = interpolate.interp1d(tempTS, temp, kind="linear", copy="False", assume_sorted="True")
    # model, modelTS = wm.MakeModel(dataList, tempList, [st,en,slo], fn=InterpFn)

    # Fit with MCMC and get best-fit parameters
    numSteps, burnIn = 3000, 1800  # default: 10000, 5000.  fast: 3000, 1800  long test: 20000,10000
    wfModel = wm.TemplateModel( dataList, dataNoise, tempList )
    M = pymc.MCMC( pymc.Model( wfModel ) )
    M.use_step_method(pymc.Metropolis, M.startTime, proposal_sd=100., proposal_distribution='Normal')
    M.use_step_method(pymc.Metropolis, M.energy, proposal_sd=1., proposal_distribution='Normal')
    M.use_step_method(pymc.Metropolis, M.slowness, proposal_sd=100., proposal_distribution='Normal')
    M.sample(iter=numSteps, verbose=0)
    st = np.median(M.trace("startTime")[:])
    en = np.median(M.trace("energy")[:])
    slo = np.median(M.trace("slowness")[:])
    InterpFn = interpolate.interp1d(tempTS, temp, kind="linear", copy="False", assume_sorted="True")
    model, modelTS = wm.MakeModel(dataList, tempList, [st,en,slo], fn=InterpFn)
    print "MCMC:",st,en,slo

    # Fit with SciPy minimizer
    MakeTracesGlobal() # creates 3 global arrays: startTrace, enTrace, sloTrace
    floats = [dataST, dataE, 1]
    print "Minimizer guesses:",floats
    datas = [dataList,tempList,InterpFn]
    result = minimize(findLnLike, floats, args=datas, method="Nelder-Mead")
    st, en, slo = result["x"]
    print "Minimizer: %.1f  %.1f  %.1f  Success: %s.  %s" % (st, en, slo, result["success"], result["message"])
    minimizer, minimizerTS = wm.MakeModel(dataList, tempList, [st,en,slo], fn=InterpFn)

    # plots
    fig = plt.figure(figsize=(11,7), facecolor='w')
    p1 = plt.subplot2grid((6,7), (0,0), colspan=4, rowspan=2) # original
    p2 = plt.subplot2grid((6,7), (2,0), colspan=4, rowspan=3) # rising edge
    p3 = plt.subplot2grid((6,7), (0,4), colspan=3, rowspan=2 )           # trace 1
    p4 = plt.subplot2grid((6,7), (2,4), colspan=3, rowspan=2, sharex=p3) # trace 2
    p5 = plt.subplot2grid((6,7), (4,4), colspan=3, rowspan=2, sharex=p3) # trace 3

    # p1 = plt.subplot(211)
    p1.set_title("Energy %.1f keV  Start Time %.0f ns" % (dataE, dataST))
    p1.set_ylabel("ADC [A.U.]",y=0.95, ha='right')
    p1.set_xlabel("Time (ns)",x=0.95, ha='right')
    p1.plot(waveTS,wave,color='blue',alpha=0.8,label='Data WF')
    p1.plot(guessFullTS,guessFull,color='orange',alpha=0.8,label='Guess WF')
    p1.axvline(x=dataST,color='green')
    p1.legend(loc=4)

    # p2 = plt.subplot(212)
    p2.plot(dataTS, data, color='blue',label='Data')
    p2.plot(guessTS, guess, color='orange',label='Guess')
    p2.plot(modelTS, model, color='red',linewidth=4,alpha=0.8,label='MCMC')
    p2.plot(minimizerTS, minimizer, color='cyan',linewidth=1,label='Nelder-Mead')
    p2.legend(loc=4)

    p3.cla()
    p3.set_title("startTime %.1f  Energy %.2f  Slow %.1f" % (st,en,slo))
    p3.plot(stTrace[1:])
    p3.set_ylabel('startTime')

    p4.cla()
    p4.plot(enTrace[1:])
    p4.set_ylabel('energy')

    p5.cla()
    p5.plot(sloTrace[1:])
    p5.set_ylabel('slowness')

    plt.tight_layout()
    plt.subplots_adjust(hspace=0.35)
    # plt.show(block=False)
    # plt.show()
    plt.savefig("./plots/minimizer-test.pdf")


def MakeTracesGlobal():
    """ This is so 'findLnLike' can write to the trace arrays """
    tmp1 = tmp2 = tmp3 = np.empty([])
    global stTrace, enTrace, sloTrace
    stTrace, enTrace, sloTrace = tmp1, tmp2, tmp3


def findLnLike(floats, datas):
    global stTrace, enTrace, sloTrace

    # Unpack all parameters
    st, en, slo = floats
    dataList, tempList, InterpFn = datas
    data, dataTS, dataE, dataST, loWin, hiWin, dataNoise = dataList
    temp, tempTS, tempE, tempST = tempList

    # Make a trace
    stTrace = np.append(stTrace,st)
    enTrace = np.append(enTrace,en)
    sloTrace = np.append(sloTrace,slo)

    # Build the model to compare with data
    model, modelTS = wm.MakeModel(dataList, tempList, floats, fn=InterpFn)

    # Return NLL
    lnLike = -0.5 * np.sum ( np.power((data-model)/dataNoise, 2) - np.log( 1 / np.power(dataNoise,2) ) )
    return -1.0 * lnLike

if __name__ == "__main__":
    main()
