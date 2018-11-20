#!/usr/bin/env python
""" A collection of 'useful' LAT routines.
    C. Wiseman, B. Zhu
"""
import sys, pywt, random, os, glob
import numpy as np
import tinydb as db
from scipy.fftpack import fft
from scipy.optimize import curve_fit
from scipy.ndimage.filters import gaussian_filter
from scipy import signal as sg
import scipy.special as sp

limit = sys.float_info.max # equivalent to std::numeric_limits::max() in C++
homePath = os.path.expanduser('~')
bgDir = homePath + "/project/bg-lat"
calDir = homePath + "/project/cal-lat"


def update_progress(progress, runNumber=None):
    """ adapted from from https://stackoverflow.com/a/15860757 """
    barLength = 20  # Modify this to change the length of the progress bar
    status = ""
    if isinstance(progress, int):
        progress = float(progress)
    if not isinstance(progress, float):
        progress = 0
        status = "error: progress var must be float\r\n"
    if progress < 0:
        progress = 0
        status = "Halt...\r\n"
    if progress >= 1:
        progress = 1
        status = "Done...\r\n"
    block = int(round(barLength * progress))
    if runNumber is None:
        text = "\rPROGRESS : [{}] {:0.3f}% {}".format(
            "#" * block + "-" * (barLength - block), progress * 100, status)
    else:
        text = "\rPROGRESS : [{}] {:0.3f}% {} (Run {})".format(
            "#" * block + "-" * (barLength - block), progress * 100, status,
            runNumber)
    sys.stdout.write(text)
    sys.stdout.flush()


def nPol(x, *par):
    y = 0
    for i, p in enumerate(par):
        y += p*x**i
    return y


def expFunc(x, a, tau, b):
    return a * np.exp(x/tau) + b


def oneExp(x, a, c):
    # in the style of roofit, which only has one parameter per exp
    return a * np.exp(x*c)


def twoExp(x, a, c1, c2):
    # in the style of roofit, which only has one parameter per exp
    return a * (np.exp(x*c1) + np.exp(x*c2))


def gaus(x, mu, sig, amp):
    """ Normalize this to 1 so that amp = nCts. """
    return amp * np.exp(-(x - mu)**2 / (2 * sig**2)) / (sig * np.sqrt(2*np.pi))


def nGaus(x, *par):
    y = np.zeros_like(x)
    for i in range(0, len(par), 3):
        mu = par[i]
        sig = par[i+1]
        amp = par[i+2]
        y = y + gaus(x, mu, sig, amp)
    return y


def sig_ae(E,m):
    """ Axioelectric cross section.
    E, m are in units of keV.  must multiply this result by sig_pe. """
    beta = (1 - m**2./E**2.)**(1./2)
    return (1 - (1./3.)*beta**(2./3.)) * (3. * E**2.) / (16. * np.pi * (1./137.) * 511.**2. * beta)


def logistic(x,mu,sig,amp,skew):
    from scipy.stats import genlogistic
    return amp * genlogistic.cdf(x, skew, mu, sig)


def weibull(x,c,loc,scale,amp):
    from scipy.stats import weibull_min
    return amp * weibull_min.cdf(x,c,loc,scale)


def xgauss(x,k,loc,scale,amp):
    from scipy.stats import exponnorm
    return amp * exponnorm.cdf(x,k,loc,scale)


def linFunc(x, m, b):
    return m*x + b


def pol1(x, a, b, c):
    return a*x**2 + b*x + c


def logisticFunc(x,mu,sig,amp):
    return amp/(1.+np.exp(-(x-mu)/sig))


def erFunc(x,mu,sig,amp):
    from scipy.special import erf
    return amp*0.5*(1 + erf( (x - mu)/(sig*np.sqrt(2) ) ))
    # from scipy.special import expit
    # return expit((x-mu)/sig)


def getHistInfo(x,h):
    """ Computes max, mean, width, percentiles of a numpy
    array based histogram , w/ x values 'x' and counts 'h'. """
    if np.sum(h)==0:
        return 0, 0, 0, [0,0,0,0], 0

    max = x[np.argmax(h)]
    avg = np.average(x, weights=h/np.sum(h))
    std = np.sqrt(np.average((h-max)**2, weights=h)/np.sum(h))
    pct = []
    for p in [5, 10, 90, 95]:
        tmp = np.cumsum(h)/np.sum(h)*100
        idx = np.where(tmp > p)
        pct.append(x[idx][0])
    wid = pct[2]-pct[0]
    return max, avg, std, pct, wid


def niceList(lst, fmt="%.2f", dtype="f"):
    """ trick to make list printing prettier """
    tmp = []
    for l in lst:
        if dtype=="f":
            tmp.append(float(str(fmt % l)))
        elif dtype=="i":
            tmp.append(int(str(fmt % l)))
        else:
            tmp.append(str(fmt % l))
    return tmp


def GetHisto(npArr, xLo, xHi, xpb, nb=None, shift=True, wts=None):
    """ This returns a histogram w/ shifted bins.  For use with: plt.plot(x, y, ls='steps') """
    if nb is None: nb = int((xHi-xLo)/xpb)
    y, x = np.histogram(npArr, bins=nb, range=(xLo, xHi), weights=wts)
    y = np.insert(y, 0, 0, axis=0)
    if shift: x = x-xpb/2.
    return x, y


def H1D(tree,bins,xlo,xhi,drawStr,cutStr,xTitle="",yTitle="",Title=None, Name=None):
    from ROOT import TH1D
    nameStr, titleStr = "", ""
    if Name == None: nameStr = str(random.uniform(1.,2.))
    else: nameStr = str(Name)
    if Title == None: titleStr = str(random.uniform(1.,2.))
    else: titleStr = str(Title)
    h1 = TH1D(nameStr,titleStr,bins,xlo,xhi)
    tree.Project(nameStr,drawStr,cutStr)
    if xTitle!="": h1.GetXaxis().SetTitle(xTitle)
    if yTitle!="": h1.GetYaxis().SetTitle(yTitle)
    return h1


def H2D(tree,xbins,xlo,xhi,ybins,ylo,yhi,drawStr,cutStr,xTitle="",yTitle="",Title=None, Name=None):
    from ROOT import TH2D
    nameStr, titleStr = "", ""
    if Name == None: nameStr = str(random.uniform(1.,2.))
    else: nameStr = str(Name)
    if Title == None: titleStr = str(random.uniform(1.,2.))
    else: titleStr = str(Title)
    h2 = TH2D(nameStr,titleStr,xbins,xlo,xhi,ybins,ylo,yhi)
    tree.Project(nameStr,drawStr,cutStr)
    if xTitle!="": h2.GetXaxis().SetTitle(xTitle)
    if yTitle!="": h2.GetYaxis().SetTitle(yTitle)
    return h2


class processWaveform:
    """ Auto-processes waveforms into python-friendly formats.

    NOTE: DS2 (multisampling) waveform bug (cf. wave-skim.cc):
       -> All samples from regular and aux waveforms work correctly with GetVectorData in PyROOT.
       -> Combining them into an MJTMSWaveform and then casting to MGTWaveform causes the first
          4 samples to be set to zero with GetVectorData in PyROOT.  (They are actually nonzero.)
       -> This problem doesn't exist on the C++ side. There's got to be something wrong with the python wrapper.
       -> The easiest workaround is just to set remLo=4 for MS data.
    """
    def __init__(self, wave, remLo=0, remHi=2):
        from ROOT import MGTWaveform
        # initialize
        # self.waveMGT = wave
        self.offset = wave.GetTOffset()
        self.period = wave.GetSamplingPeriod()
        self.length = wave.GetLength()
        vec = wave.GetVectorData()
        npArr = np.fromiter(vec, dtype=np.double, count=self.length)
        ts = np.arange(self.offset, self.offset + self.length * self.period, self.period) # superfast!

        # resize the waveform
        removeSamples = []
        if remLo > 0: removeSamples += [i for i in range(0,remLo+1)]
        if remHi > 0: removeSamples += [npArr.size - i for i in range(1,remHi+1)]
        self.ts = np.delete(ts,removeSamples)
        self.waveRaw = np.delete(npArr,removeSamples)   # force the size of the arrays to match

        # compute baseline and noise
        self.noiseAvg,_,self.baseAvg = baselineParameters(self.waveRaw)
        self.waveBLSub = np.copy(self.waveRaw) - self.baseAvg

    # constants
    def GetOffset(self): return self.offset
    def GetBaseNoise(self): return self.baseAvg, self.noiseAvg

    # arrays
    def GetTS(self): return self.ts
    def GetWaveRaw(self): return self.waveRaw
    def GetWaveBLSub(self): return self.waveBLSub


def GetVX(tree, bNames, theCut="", showVec=True):
    """ Sick of using GetV1234(), let's try to write a versatile tree parser.
    - This isn't quite as fast as Draw, but it can draw more branches and still do entry lists.
    - If this gets too complicated I'll have to figure out how to use root_numpy.
    - NOTE: If the tree uses vectors, only HITS PASSING CUTS will be returned for each event.
      EXCEPT, if you cut ONLY on a non-vector branch (say mH), then Iteration$ is always 0 (WRONG).
      It's not my fault!  The Draw() that creates the entry list is messed up.
      HACKY FIX: add "channel > 0". It should always be true.  Who knows how tree->Scan does it right.
    """
    from ROOT import TChain, MGTWaveform

    # make sure tree has entries
    if tree.GetEntries()==0:
        print("No entries found.")
        return None

    # make sure branches exist
    missingBranches = False
    branchList = tree.GetListOfBranches()
    for br in bNames:
        if br == "Entry$" or br=="Iteration$": continue
        if br not in branchList:
            print("ERROR, couldn't find branch:",br)
            missingBranches = True
    if missingBranches:
        return None

    # this is what we'll return: {branchName:[numpy array of vals]}
    tVals = {br:[] for br in bNames}

    # create an entry list of hits passing cuts
    evtList = False
    if theCut != "":
        evtList = True
        if showVec:
            theCut += " && channel > 0"
        nPass = tree.Draw("Entry$:Iteration$",theCut,"GOFF")
        evtList = tree.GetV1()
        itrList = tree.GetV2()
        evtList = [int(evtList[idx]) for idx in range(nPass)]
        itrList = [int(itrList[idx]) for idx in range(nPass)]
        # for idx in range(len(evtList)):
            # print(evtList[idx],itrList[idx])

        # match the iteration numbers to the entry
        evtSet = sorted(set(evtList))
        evtDict = {evt:[] for evt in evtSet}
        for idx in range(len(itrList)):
            evtDict[evtList[idx]].append(itrList[idx])
        treeItr = evtSet
    else:
        treeItr = range(tree.GetEntriesFast())
        # treeItr = range(7) # debug, obvs

    # activate only branches we want, and detect vector types
    sizeVec = None
    brTypes = {br:0 for br in bNames}
    tree.SetBranchStatus('*',0)
    for name in bNames:
        if name=="Entry$" or name=="Iteration$": continue
        tree.SetBranchStatus(name,1)
        if "vector" in tree.GetBranch(name).GetClassName():
            brTypes[name] = "vec"
            sizeVec = name

    # loop over the tree
    for iEvt in treeItr:
        tree.GetEvent(iEvt)
        nHit = 1 if sizeVec is None else getattr(tree,sizeVec).size()
        branchItr = evtDict[iEvt] if evtList else range(nHit)
        for name in tVals:
            for iH in branchItr:
                if brTypes[name] == "vec":
                    if name=="MGTWaveforms": # mgtwaveforms are weird, they don't persist in memory
                        wf = getattr(tree,name)[iH]
                        sig = processWaveform(wf)
                        tVals[name].append(sig)
                    else:
                        tVals[name].append(getattr(tree,name)[iH])
                    continue
                elif name == "Entry$":
                    tVals[name].append(iEvt)
                    continue
                elif name == "Iteration$":
                    tVals[name].append(iH)
                    continue
                else:
                    tVals[name].append(getattr(tree,name))

    # convert lists to np arrays (except waveforms) and return
    # tValsNP = {}
    # for br in bNames:
    #     tValsNP[br] = tVals[br] if br=="MGTWaveforms" else np.asarray(tVals[br])
    # return tValsNP

    tVals = {br:np.asarray(tVals[br]) for br in bNames}
    return tVals


def getDetPos(run):
    """ Load position info for all enabled channels in a run. """
    from ROOT import GATDataSet
    gds = GATDataSet(run)
    chMap = gds.GetChannelMap()
    chSet = gds.GetChannelSettings()
    enabledIDs = chSet.GetEnabledIDList()
    enabledIDs = [enabledIDs[idx] for idx in range(enabledIDs.size())]
    detPos = {}
    for ch in enabledIDs:
        hglg = "H" if ch%2==0 else "L"
        pos = "%sD%d-%s" % (chMap.GetString(ch,"kStringName"), chMap.GetInt(ch,"kDetectorPosition"),hglg)
        detPos[ch] = pos
    return detPos


def getChan(crate,card,gain):
    """ In hex, 0x(crate+1)(card)(gain). """
    return(((crate+1)*2 << 8) + (card << 4) + gain)


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


def npTH1D(hist, opt="up"):
    """ By default, xArr is the bin centers.
    If you have e.g. an energy histogram,
    you want to shift up by half a bin width.
    (ROOT does that by default).
    """
    bins = hist.GetNbinsX()
    xArr = np.zeros(bins+1)
    yArr = np.zeros(bins+1)
    for i in range(bins+1):
        ctr = hist.GetXaxis().GetBinCenter(i)
        if opt=="i": xArr[i] = int(ctr)
        else: xArr[i] = ctr
        yArr[i] = hist.GetBinContent(i)
    xpb = xArr[1] - xArr[0]

    if opt=="up":
        xArr = xArr + xpb/2.

    return xArr, yArr, xpb


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


def evalXGaus(x,mu,sig,tau):
    """ Updated with correct limit, 1 May 2018.  (Correct version in LAT already.)
    Ported from GAT/BaseClasses/GATPeakShapeUtil.cc
    Negative tau: Regular WF, high tail
    Positive tau: Backwards WF, low tail
    """
    tmp = (x-mu + sig**2./2./tau)/tau

    # np.exp of this is 1.7964120280206387e+308, the largest python float value: sys.float_info.max
    fLimit = 709.782

    if all(tmp < fLimit):
        return np.exp(tmp)/2./np.fabs(tau) * sp.erfc((tau*(x-mu)/sig + sig)/np.sqrt(2.)/np.fabs(tau))
    else:
        # print("Exceeded limit ...")
        # Use an approx. derived from the asymptotic expansion for erfc, listed on wikipedia.
        den = 1./(sig + tau*(x-mu)/sig)
        return sig * evalGaus(x,mu,sig) * den * (1.-tau**2. * den**2.)


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

    print("xbins %d  bx1 %d  bx2 %d" % (xbins,bx1,bx2))

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
    for x in range(0,500):
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
        print("baseline fit failed! using simpler method.")
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


def trapFilter(signalRaw, rampTime=400, flatTime=200, decayTime=0.):
    """ Apply a trap filter to a waveform. """
    baseline = 0.
    decayConstant = 0.
    norm = rampTime
    if decayTime != 0:
        decayConstant = 1./(np.exp(1./decayTime) - 1)
        norm *= decayConstant

    trapOutput = np.zeros_like(signalRaw)
    fVector = np.zeros_like(signalRaw)
    scratch = np.zeros_like(signalRaw)

    fVector[0] = signalRaw[0] - baseline
    trapOutput[0] = (decayConstant+1.)*(signalRaw[0] - baseline)

    wf_minus_ramp = np.zeros_like(signalRaw)
    wf_minus_ramp[:rampTime] = baseline
    wf_minus_ramp[rampTime:] = signalRaw[:len(signalRaw)-rampTime]

    wf_minus_ft_and_ramp = np.zeros_like(signalRaw)
    wf_minus_ft_and_ramp[:(flatTime+rampTime)] = baseline
    wf_minus_ft_and_ramp[(flatTime+rampTime):] = signalRaw[:len(signalRaw)-flatTime-rampTime]

    wf_minus_ft_and_2ramp = np.zeros_like(signalRaw)
    wf_minus_ft_and_2ramp[:(flatTime+2*rampTime)] = baseline
    wf_minus_ft_and_2ramp[(flatTime+2*rampTime):] = signalRaw[:len(signalRaw)-flatTime-2*rampTime]

    scratch = signalRaw - (wf_minus_ramp + wf_minus_ft_and_ramp + wf_minus_ft_and_2ramp )

    if decayConstant != 0:
        fVector = np.cumsum(fVector + scratch)
        trapOutput = np.cumsum(trapOutput +fVector+ decayConstant*scratch)
    else:
        trapOutput = np.cumsum(trapOutput + scratch)

    # Normalize and resize output
    trapOutput[:len(signalRaw) - (2*rampTime+flatTime)] = trapOutput[2*rampTime+flatTime:]/norm
    trapOutput.resize( (len(signalRaw) - (2*rampTime+flatTime)))
    return trapOutput


def wfDerivative(signalRaw,sp=10.):
    """ Take a derivative of a waveform numpy array.
    Adapted from $MGDODIR/Transforms/MGWFBySampleDerivative
    y[n] = ( x[n+1] - x[n] )/sp
    where sp is the sampling period of the waveform.
    """
    signalDeriv = np.zeros(len(signalRaw))
    for i in range(len(signalRaw)-1):
        signalDeriv[i] = (signalRaw[i+1] - signalRaw[i])/sp
    return signalDeriv


def MGTWFFromNpArray(npArr):
    """ Convert a numpy array back into an MGTWaveform. """
    from ROOT import MGTWaveform, std
    vec = std.vector("double")()
    for adc in npArr: vec.push_back(adc)
    mgtwf = MGTWaveform()
    mgtwf.SetData(vec)
    return mgtwf


def generateSimBaseline():
    """ Generates a fake baseline from a force triggered power spectrum, based off of WC's thesis """
    from ROOT import TFile,MGTWaveformFT

    fNoiseFile = TFile("./data/noise_spectra_DS0.root")
    noiseFT = fNoiseFile.Get("noiseWF_ch624_P42574A");
    tempFT = np.zeros(noiseFT.GetDataLength(), dtype=complex)
    # Skip 0th component (DC component)
    for i in range(1, noiseFT.GetDataLength()):
        sigma = np.sqrt(noiseFT.At(i).real()*noiseFT.At(i).real()+noiseFT.At(i).imag()*noiseFT.At(i).imag())/np.sqrt(2.)
        tempFT[i] = np.complex(sigma*np.random.random_sample(), sigma*np.random.random_sample())

    # This WF is actually 2030 samples because the force trigger WF is longer
    tempWF = np.fft.irfft(tempFT)
    # Cut off first 50 and last 100 samples, symmetrically pad
    # Padding with an asymmetric number of samples removes an artifact in the power spectrum
    tempWF = np.pad(tempWF[50:-100], (50, 100), mode='symmetric')
    # Seems to be off by a scale of around pi
    tempWF = tempWF*np.pi
    # Cut off first 6 and last 6 samples to make the length 2018
    return tempWF[6:-6]


def peakdet(v, delta, x = None):
    """
    Converted from MATLAB script at http://billauer.co.il/peakdet.html
    Returns two arrays

    function [maxtab, mintab]=peakdet(v, delta, x)
    # PEAKDET Detect peaks in a vector
    #        [MAXTAB, MINTAB] = PEAKDET(V, DELTA) finds the local
    #        maxima and minima ("peaks") in the vector V.
    #        MAXTAB and MINTAB consists of two columns. Column 1
    #        contains indices in V, and column 2 the found values.
    #
    #        With [MAXTAB, MINTAB] = PEAKDET(V, DELTA, X) the indices
    #        in MAXTAB and MINTAB are replaced with the corresponding
    #        X-values.
    #
    #        A point is considered a maximum peak if it has the maximal
    #        value, and was preceded (to the left) by a value lower by
    #        DELTA.
    # Eli Billauer, 3.4.05 (Explicitly not copyrighted).
    # This function is released to the public domain; Any use is allowed.
    """
    maxtab, mintab = [], []

    if x is None:
        x = np.arange(len(v))
    v = np.asarray(v)
    if len(v) != len(x):
        sys.exit('Peak Finder Error: Input vectors v and x must have same length')
    if not np.isscalar(delta):
        sys.exit('Peak Finder Error: Input argument delta must be a scalar')
    if delta <= 0:
        sys.exit('Peak Finder Error: Input argument delta must be positive')

    mn, mx, mnpos, mxpos = np.Inf, -np.Inf, np.NaN, np.NaN
    lookformax = True

    for i in np.arange(len(v)):
        this = v[i]
        if this > mx:
            mx = this
            mxpos = x[i]
        if this < mn:
            mn = this
            mnpos = x[i]
        if lookformax:
            if this < mx-delta:
                maxtab.append((mxpos, mx))
                mn = this
                mnpos = x[i]
                lookformax = False
        else:
            if this > mn+delta:
                mintab.append((mnpos, mn))
                mx = this
                mxpos = x[i]
                lookformax = True
    return np.array(maxtab), np.array(mintab)


def GetPeaks(hist, xvals, thresh):
    """ Wrapper for peakdet function. Returns two numpy arrays."""
    scan,_ = peakdet(hist, thresh)
    pk1, ct1 = [], []
    for i in range(len(scan)):
        pk1.append( xvals[ int(scan[i][0]) ] )
        ct1.append( scan[i][1] )
    pk1, ct1 = np.asarray(pk1), np.asarray(ct1)
    return pk1, ct1


def peakModel238(ene,amp,mu,sig,c1):
    """ Gaussian plus a flat background component. """
    return gauss_function(ene,amp,mu,sig) + c1


def peakModel238240(x,a1,c0,mu,sig,c1,tau,c2,a2,b):
    """ Combined 238 + 240 peak model.
        Source: http://nucleardata.nuclear.lu.se/toi/radSearch.asp
    Parameter [guess value] [description]
    Pb212 (238.6322 keV) - Gaussian + xGauss backward step + erfc forward step
        a1  100.   overall amplitude
        c0  1.     gaussian amplitude
        mu  pkFit  238 peak mean in raw 'x' units, input from peakFit
        sig 0.4    width
        c1  10.    xgauss amplitude
        tau 100.   decay constant of xgauss function
        c2  0.1    erfc amplitude
    Ra224 (240.9866 keV) - Gaussian only (it's a much smaller contribution)
        a2  0.1*a1 240 peak amplitude
    Linear BG
        # m  -0.001  slope
        b  5.      offset
    """
    f0 = gauss_function(x,c0,mu,sig)
    f1 = evalXGaus(x,mu,sig,tau)
    f2 = sp.erfc((x-mu)/sig/np.sqrt(2.))
    mu240 = 240.9866 / (238.6322 / mu) # peak/gain
    f3 = gauss_function(x,a2,mu240,sig)
    return a1 * (f0 + c1*f1 + c2*f2) + f3 + b


def peakModel238_2(x,a1,c0,mu,sig,c1,tau,c2,b):
    """ See above. """
    f0 = gauss_function(x,c0,mu,sig)
    f1 = evalXGaus(x,mu,sig,tau)
    f2 = sp.erfc((x-mu)/sig/np.sqrt(2.))
    return a1 * (f0 + c1*f1 + c2*f2) + b


def walkBackT0(trap, timemax=10000., thresh=2., rmin=0, rmax=1000, forward=False):
    """
        Leading Edge start time -- walk back or forward from a maximum to threshold
        Times are returned in ns
    """
    minsample = np.amax([0,rmin])
    maxsample = np.amin([len(trap),rmax])
    trapMax = np.argmax(trap[minsample:maxsample])
    foundFirst, triggerTS = False, 0
    sampleArr = []
    if forward:
        sampleArr = range(trapMax,maxsample)
    else:
        sampleArr = range(trapMax,minsample,-1)
    for idx, i in enumerate(sampleArr):
        # If passed threshold, means we walked past the sample that crossed
        if trap[i] <= thresh:
            foundFirst = True
            # Interpolate between the current and previous sample if the difference isn't zero
            if (trap[sampleArr[idx]]-trap[sampleArr[idx-1]]) != 0:
                triggerTS = ((thresh-trap[sampleArr[idx]]) * (sampleArr[idx]-sampleArr[idx-1])/(trap[sampleArr[idx]]-trap[sampleArr[idx-1]]) + i)*10
            else: triggerTS = (i+1)*10
            break
    # Save-guards if the t0 goes out of range for picking off on the large trapezoid
    if triggerTS >= timemax:
        return timemax, foundFirst
    elif triggerTS <= 0:
        return 0., foundFirst
    else:
        return triggerTS, foundFirst


def constFractiont0(trap, frac=0.1, delay=200, thresh=0., rmin=0, rmax=1000):
    """
        Constant Fraction start time
        1) Invert the signal
        2) Delay and sum the original + inverted
        3) Walk back from maximum to a threshold (usually zero crossing)
        Times are returned in ns
    """
    minsample = np.amax([0,rmin])
    maxsample = np.amin([len(trap)-delay,rmax])
    invertTrap = np.multiply(trap, -1.*frac)
    summedTrap = np.add(invertTrap[delay:], trap[:-delay])
    trapMax = np.argmax(summedTrap[0:1000])
    foundFirst, triggerTS = False, 0
    for i in range(trapMax,0,-1):
        if summedTrap[i] <= thresh:
            foundFirst = True
            if (summedTrap[i+1]-summedTrap[i]) != 0:
                triggerTS = ((thresh-summedTrap[i])*((i+1)-i)/(summedTrap[i+1]-summedTrap[i]) + i)*10
            else: triggerTS = (i+1)*10
            break
    return triggerTS, foundFirst


def asymTrapFilter(data,ramp=200,flat=100,fall=40,padAfter=False):
    """ Computes an asymmetric trapezoidal filter """
    trap = np.zeros(len(data))
    for i in range(len(data)-1000):
        w1 = ramp
        w2 = ramp+flat
        w3 = ramp+flat+fall
        r1 = np.sum(data[i:w1+i])/(ramp)
        r2 = np.sum(data[w2+i:w3+i])/(fall)
        if not padAfter:
            trap[i+1000] = r2 - r1
        else:
            trap[i] = r2 - r1
    return trap


def MakeSiggenWaveform(samp,r,z,ene,t0,smooth=1,phi=np.pi/8):
    # This works FINE on my machine but not on PDSF.  Damn you, PDSF.
    # Use pysiggen to generate a waveform w/ arb. ADC amplitude. Thanks Ben.
    from pysiggen import Detector

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
    detName = "../data/conf/P42574A_grad%0.2f_pcrad%0.2f_pclen%0.2f.conf" % (0.05, 2.5, 1.65)
    detector = Detector(detName, timeStep=timeStepSize, numSteps=fitSamples*10./timeStepSize, maxWfOutputLength=5000)
    # detector.LoadFieldsGrad("../data/fields_impgrad.npz",pcLen=1.6, pcRad=2.5)

    data = np.load("../data/fields_impgrad.npz")
    print(data.keys())
    wpArray  = data['wpArray']
    # # efld_rArray = data['efld_rArray']
    # # efld_zArray = data['efld_zArray']
    # gradList = data['gradList']
    print(wpArray.shape)
    return





    detector.LoadFieldsGrad("../data/fields_impgrad.npz")
    return

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
    #     print("amplitude diff: %f" % ( (np.amax(wf_notrap) - np.amax(wf)) /  np.amax(wf_notrap) ))

    return wf_notrap, timesteps
