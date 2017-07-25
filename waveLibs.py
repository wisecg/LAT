#!/usr/local/bin/python
import sys, pywt, random
import numpy as np
from scipy.fftpack import fft
from scipy.optimize import curve_fit
from scipy.ndimage.filters import gaussian_filter
from scipy import signal as sg
# from pysiggen import Detector
from ROOT import TTree, std, TH1D, TH2D
from ROOT import MGTWaveform
"""
A collection of 'useful' crap.
C. Wiseman, B. Zhu
v1. 20 March 2017
v2. 21 June 2017
"""

def H1D(tree,bins,xlo,xhi,drawStr,cutStr,xTitle="",yTitle=""):
    nameStr = str(random.uniform(1.,2.))
    h1 = TH1D(nameStr,nameStr,bins,xlo,xhi)
    tree.Project(nameStr,drawStr,cutStr)
    if xTitle!="": h1.GetXaxis().SetTitle(xTitle)
    if yTitle!="": h1.GetYaxis().SetTitle(yTitle)
    return h1

def H2D(tree,xbins,xlo,xhi,ybins,ylo,yhi,drawStr,cutStr,xTitle="",yTitle=""):
    nameStr = str(random.uniform(1.,2.))
    h2 = TH2D(nameStr,nameStr,xbins,xlo,xhi,ybins,ylo,yhi)
    tree.Project(nameStr,drawStr,cutStr)
    if xTitle!="": h2.GetXaxis().SetTitle(xTitle)
    if yTitle!="": h2.GetYaxis().SetTitle(yTitle)
    return h2

def Get1DBins(hist,xmin,xmax):
    bx1 = hist.FindBin(xmin)
    bx2 = hist.FindBin(xmax)
    return bx1, bx2


def Get2DBins(hist,xmin,xmax,ymin,ymax):
    bx1 = hist.GetXaxis().FindBin(xmin)
    bx2 = hist.GetXaxis().FindBin(xmax)
    by1 = hist.GetYaxis().FindBin(ymin)
    by2 = hist.GetYaxis().FindBin(ymax)
    return bx1, bx2, by1, by2

def npTH1D(hist,opt=""):
    bins = hist.GetNbinsX()
    xArr = np.zeros(bins)
    yArr = np.zeros(bins)
    for i in range(bins):
        ctr = hist.GetXaxis().GetBinCenter(i)
        if opt=="i": xArr[i] = int(ctr)
        else: xArr[i] = ctr
        yArr[i] = hist.GetBinContent(i)
    return xArr,yArr

def integFunc(arr):
    integ = np.zeros(len(arr))
    sum = 0
    for i in range(0,len(arr)):
        sum+=arr[i]
        integ[i] = sum
    return integ

def GetIntegralPoints(hist):
    x_h0, y_h0 = npTH1D(hist)
    int_h0 = integFunc(y_h0)

    idx99 = np.where(int_h0 > 0.99)
    idx95 = np.where(int_h0 > 0.95)
    idx90 = np.where(int_h0 > 0.90)
    idx01 = np.where(int_h0 > 0.01)
    idx05 = np.where(int_h0 > 0.05)

    val99 = x_h0[idx99][0]
    val95 = x_h0[idx95][0]
    val01 = x_h0[idx01][0]
    val05 = x_h0[idx05][0]
    val90 = x_h0[idx90][0]

    return val99,val95,val01,val05,val90

def SetPars(f1,parList):
    for idx, par in enumerate(parList):
        f1.SetParameter(idx,par)

def GetPars(f1):
    parList = []
    npar = f1.GetNpar()
    for i in range(npar):
        parList.append(f1.GetParameter(i))
    return parList


def gauss_function(x, a, x0, sigma):
    """ Just a simple gaussian. """
    return a * np.exp(-(x - x0)**2 / (2 * sigma**2))


def tailModelExp(t, a, b):
    return a * np.exp(-1.0 * t / b)


def tailModelPol(t,a,b,c,d):
    return a + b * t + c * np.power(t,2) + d * np.power(t,3)


def rootToArray(hist, xLow, xHigh):
    """ Take a ROOT TH1D histogram and a range, get a numpy array.
    Note on plotting:
    Can't use np.histogram, the hist has already been done.
    So make a plot that fakes it with :
        xaxis = np.arange(xlow, xhi, step)
        plt.bar(xaxis,hist,width=step)  # 1: fake it with bar
        plt.plot(xaxis,hist,ls='steps-post')  # 2: fake it with plot
    """
    binLow = hist.FindBin(xLow)
    binHigh = hist.FindBin(xHigh)
    loopArray = range(binLow, binHigh )
    npArray = np.empty_like(loopArray, dtype=np.float)
    for (iArray, iBin) in enumerate(loopArray):
        npArray[iArray] = np.float(hist.GetBinContent(iBin))
    return npArray


def th2Array(hist,xbins,xlo,xhi,ybins,ylo,yhi):
    """ Take a ROOT TH2D histogram and get a numpy array.
        Input the same values as the TH2 constructor.
    """
    bx1,bx2,by1,by2 = Get2DBins(hist,xlo,xhi,ylo,yhi)

    print "xbins %d  bx1 %d  bx2 %d" % (xbins,bx1,bx2)

    arr = np.zeros((xbins,ybins),dtype=float)
    for i in range(bx1-1, bx2-1): # it's zero indexed
        for j in range(by1-1, by2-1):
            arr[i,j] = hist.GetBinContent(i,j)

    return arr




def baselineParameters(signalRaw):
    """ Finds basic parameters of baselines using first 500 samples. """
    rms = 0
    slope = 0
    baselineMean = 0
    baselineAveSq = 0
    for x in xrange(0,500):
        wfX = signalRaw[x]
        baselineMean += wfX
        baselineAveSq += wfX*wfX
    baselineMean /= 500.
    baselineAveSq /= 500.
    rms = np.sqrt( baselineAveSq - baselineMean*baselineMean);
    return rms, slope, baselineMean


def findBaseline(signalRaw):
    """ Find the average starting baseline from an MJD waveform. """
    (hist, bins) = np.histogram(signalRaw[:200], bins=np.arange(-8000, 8000, 1))
    fitfunc = lambda p, x: p[0] * np.exp(-0.5 * ((x - p[1]) / p[2])**2) + p[3]
    errfunc = lambda p, x, y: (y - fitfunc(p, x))
    c = np.zeros(3)
    try:
        c,_ = curve_fit(gauss_function, bins[:-1], hist, p0=[np.amax(hist), bins[np.argmax(hist)], 5])
    except:
        print "baseline fit failed! using simpler method."
        c[2],_,c[1] = baselineParameters(signalRaw)
    mu, std = c[1], c[2]
    return mu, std


def fourierTransform(signalRaw):
    """ Simple FFT for a waveform """
    yf = fft(signalRaw)
    T = 1e-8 # Period
    N = len(yf)
    xf = np.linspace(0.0, 1.0/(2.0*T), N/2)
    yret = 2.0/N*np.abs(yf[:N/2]) # FT spectrum
    ypow = np.multiply(yret, yret) # Power spectrum
    return xf, yret, ypow


def waveletTransform(signalRaw, level=4, wavelet='db2', order='freq'):
    """ Use PyWavelets to do a wavelet transform. """
    wp = pywt.WaveletPacket(signalRaw, wavelet, 'symmetric', maxlevel=level)
    nodes = wp.get_level(level, order=order)
    yWT = np.array([n.data for n in nodes], 'd')
    yWT = abs(yWT)
    return wp.data, yWT


def trapezoidalFilter(signalRaw, rampTime=400, flatTime=200, decayTime=0.):
    """ Apply a trap filter to a waveform. """
    baseline = 0.
    decayConstant = 0.
    norm = rampTime
    if decayTime != 0:
        decayConstant = 1./(np.exp(1./decayTime) - 1)
        norm *= decayConstant

    trapOutput = np.linspace(0, len(signalRaw), num=len(signalRaw), dtype=np.double)
    fVector = np.linspace(0, len(signalRaw), num=len(signalRaw), dtype=np.double)
    fVector[0] = signalRaw[0] - baseline
    trapOutput[0] = (decayConstant+1.)*(signalRaw[0] - baseline)
    scratch = 0.
    for x in xrange(1,len(signalRaw)):
        scratch = signalRaw[x] - (
            signalRaw[x-rampTime] if x >= rampTime else baseline) - (
            signalRaw[x-flatTime-rampTime] if x >= (flatTime+rampTime) else baseline) + (
            signalRaw[x-flatTime-2*rampTime] if x >= (flatTime+2*rampTime) else baseline)
        if decayConstant != 0:
            fVector[x] = fVector[x-1] + scratch
            trapOutput[x] = trapOutput[x-1] + fVector[x] + decayConstant*scratch
        else:
            trapOutput[x] = trapOutput[x-1] + scratch

    # Normalize and resize output
    for x in xrange(2*rampTime+flatTime, len(signalRaw)):
        trapOutput[x-(2*rampTime+flatTime)] = trapOutput[x]/norm
    trapOutput.resize( (len(signalRaw) - (2*rampTime+flatTime), 1))
    return trapOutput


def wfDerivative(signalRaw,sp=10.):
    """ Take a derivative of a waveform numpy array.
    Adapted from $MGDODIR/Transforms/MGWFBySampleDerivative
    y[n] = ( x[n+1] - x[n] )/sp
    where sp is the sampling period of the waveform.
    """
    signalDeriv = np.zeros(len(signalRaw))
    for i in xrange(len(signalRaw)-1):
        signalDeriv[i] = (signalRaw[i+1] - signalRaw[i])/sp
    return signalDeriv


def MGTWFFromNpArray(npArr):
    """ Convert a numpy array back into an MGTWaveform. """
    vec = std.vector("double")()
    for adc in npArr: vec.push_back(adc)
    mgtwf = MGTWaveform()
    mgtwf.SetData(vec)
    return mgtwf

"""
# This works FINE on my machine but not on PDSF.  Damn you, PDSF.
def MakeSiggenWaveform(samp,r,z,ene,t0,smooth=1,phi=np.pi/8):
    # Use pysiggen to generate a waveform w/ arb. ADC amplitude. Thanks Ben.

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
"""


class processWaveform:
    """ Handy class for auto-processing waveforms into various numpy arrays. """
    def __init__(self, wave, removeNSamples=2, N=1, Wn=0.08, lo=50, hi=400, opt=''):

        self.waveMGT = wave                            # input an MGTWaveform object
        self.offset = self.waveMGT.GetTOffset()        # time offset [ns]
        vec = wave.GetVectorData()
        npArr = np.fromiter(vec, dtype=np.double)      # raw numpy array

        hist = wave.GimmeUniqueHist()                  # get timestamp limits and make an array
        self.start = hist.GetXaxis().GetXmin() + 5     # add 5 ns to make it start at 0 (default is -5.)
        self.stop = hist.GetXaxis().GetXmax() + 5.
        self.binsPerNS = (self.stop - self.start) / hist.GetNbinsX()
        ts = np.arange(self.start,self.stop,self.binsPerNS)
        removeSamples = [npArr.size - i for i in xrange(1,removeNSamples+1)]
        self.ts = np.delete(ts,removeSamples)
        self.waveRaw = np.delete(npArr,removeSamples)   # force the size of the arrays to match

        # self.baseAvg, self.noiseAvg = findBaseline(self.waveRaw) # get baseline and baseline RMS
        self.noiseAvg,_,self.baseAvg = baselineParameters(self.waveRaw)

        self.waveBLSub = self.waveRaw
        self.waveBLSub[:] = [x - self.baseAvg for x in self.waveRaw] # subtract the baseline value

        # self.lastZero, self.lastZeroTime, self.loWin, self.hiWin = 0,0,0,0
        # self.waveEdge, self.tsEdge = np.empty([1]), np.empty([1])
        self.waveletYOrig, self.waveletYTrans = np.empty([1]), np.empty([1])

        if (opt == 'full'):

            self.waveletYOrig, self.waveletYTrans = waveletTransform(self.waveBLSub) # wavelet transform

            self.b, self.a = sg.butter(N, Wn)
            self.waveFilt = sg.filtfilt(self.b, self.a, self.waveBLSub)  # do a smoothing

            # self.zeros = np.asarray(np.where(self.waveBLSub < 0.1))
            # zeros = np.asarray(np.where(self.waveFilt < 0.01))   # make a list of negative points ("zeros")
            # if (zeros.size > 0): self.lastZero = zeros[0,-1]     # find the last zero crossing
            # self.lastZeroTime = ts[self.lastZero]
            # self.loWin, self.hiWin = self.lastZero-lo, self.lastZero+hi  # indexes, not times
            # waveEnds = np.concatenate((np.arange(0, self.loWin), np.arange(self.hiWin, self.waveRaw.size)), axis=0)
            # self.waveEdge = np.delete(self.waveBLSub, waveEnds) # rising edge of the waveform
            # self.tsEdge = np.delete(self.ts, waveEnds)          # rising edge timestamps

            # self.waveFTX, self.waveFTY, self.waveFTPow = fourierTransform(self.waveBLSub) # fourier transform

            # self.waveTrap = trapezoidalFilter(self.waveBLSub) # trapezoidal filter
            # self.waveFiltTrap = trapezoidalFilter(self.waveFilt) # trapezoidal filter on the smoothed waveform


    # constants
    def GetOffset(self): return self.offset
    def GetStartTime(self): return self.start
    def GetStopTime(self): return self.stop
    def GetBins(self): return self.binsPerNS
    def GetBaseNoise(self): return self.baseAvg, self.noiseAvg
    def GetWindowIndex(self): return self.loWin, self.hiWin
    # def GetBaselineRMS(self): return self.baselineRMS
    # def GetBaselineSlope(self): return self.baselineSlope

    # arrays
    def GetWaveEdge(self): return self.waveEdge
    def GetTSEdge(self): return self.tsEdge
    def GetTS(self): return self.ts
    def GetWaveRaw(self): return self.waveRaw
    def GetWaveBLSub(self): return self.waveBLSub
    def GetWaveFilt(self): return self.waveFilt
    def GetLastZeroIndex(self): return self.lastZero
    def GetLastZeroTime(self): return self.lastZeroTime
    # def GetFourier(self): return self.waveFTX, self.waveFTY, self.waveFTPow
    def GetWavelet(self): return self.waveletYOrig, self.waveletYTrans
    def GetTrapezoid(self): return self.waveTrap
    # def GetFiltTrapezoid(self): return self.waveFiltTrap


class latBranch:
    # this was a good idea that ROOT will not play nice with.  Damn you, ROOT.

    def __init__(self, name, outTree, val="double"):
        self.vec = std.vector(val)
        # fails -- "TypeError: can not resolve method template call for 'Branch'""
        self.branch = outTree.Branch(name,self.vec)

    def reset(self,nChans,val=-999):
        self.vec.assign(nChans,val)

    def assign(self,idx,val):
        self.vec[idx] = val

    def fill(self):
        self.branch.Fill()

    # # initialize
    # waveS0 = std.vector("double")()
    # bWaveS0 = outTree.Branch("waveS0", waveS0)
    #
    # # reset for every tree entry
    # waveS0.assign(nChans,-999)
    #
    # # assign a value for a particular hit
    # waveS0[iH] = np.sum(waveletYTrans[2:-1,1:-1])
    #
    # # fill the branch at event level
    # bWaveS0.Fill()

