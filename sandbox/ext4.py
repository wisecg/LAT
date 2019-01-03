#!/usr/bin/env python3
import numpy as np
from scipy.optimize import curve_fit

import matplotlib.pyplot as plt
plt.style.use('../clint.mpl')
from matplotlib.colors import LogNorm

import dsi
bkg = dsi.BkgInfo()
det = dsi.DetInfo()
cal = dsi.CalInfo()
import waveLibs as wl
import tinydb as db

def main():

    # get_ext_pulser_data()
    # plot_ext_pulser()
    ext_pulser_width()


def get_ext_pulser_data():
    """ Adapted from ext2.py::getEff
    Just gets data and save into a npz file.
    """
    from ROOT import TChain, GATDataSet
    import glob

    # this is the output
    extData = {} # {run: [pIdx, runTime, extChan, hitE, fSlo]}

    for pIdx in [19,20,21]:
    # for pIdx in [19]:

        extPulserInfo = cal.GetSpecialList()["extPulserInfo"]
        attList = extPulserInfo[pIdx][0] # unused
        extChan = extPulserInfo[pIdx][-1]
        syncChan = wl.getChan(0,10,0) # 672

        runList = cal.GetSpecialRuns("extPulser",pIdx)
        for run in runList:

            # elogs: "20 Hz, 150 second runs"
            gds = GATDataSet(run)
            runTime = gds.GetRunTime() # sec
            # pulseRate = 20 # Hz

            fList = glob.glob(dsi.specialDir+"/lat/latSkimDS0_run%d_*.root" % run)
            tt = TChain("skimTree")
            for f in fList: tt.Add(f)

            tCut = "(channel==%d || channel==%d) && mH==2" % (syncChan, extChan) # enforce correct sync
            n = tt.Draw("trapENFCal:channel:fitSlo",tCut,"goff")
            hitE, chan, fSlo = tt.GetV1(), tt.GetV2(), tt.GetV3()
            hitE = np.asarray([hitE[i] for i in range(n) if chan[i]==extChan])
            fSlo = np.asarray([fSlo[i] for i in range(n) if chan[i]==extChan])

            if len(hitE)==0:
                continue

            extData[run] = [pIdx, runTime, extChan, hitE, fSlo]

            tt.Reset()

    # save output
    # np.savez("./data/lat-extPulser.npz",extData)


def plot_ext_pulser():
    """ The ext pulser data doesn't have the same centroid as the physics data
    since it wasn't tuned quite right.  So we find the centroid of a higher-E set for each channel
    and plot the other ones relative to that. """

    ds, cIdx, bIdx, sIdx = 0, 33, 75, 0 # runs 6887-6963, closest to this one
    # calDB = db.TinyDB('./calDB.json')
    # pars = db.Query()
    # fsD = dsi.getDBRecord("fitSlo_%s_idx%d_m2s238" % ("ds0_m1", cIdx), False, calDB, pars)
    # thD = dsi.getDBRecord("thresh_ds%d_bkg%d_sub%d" % (ds, bIdx, sIdx), False, calDB, pars)

    f = np.load("./data/lat-extPulser.npz")
    extData = f['arr_0'].item()  # {run: [pIdx, runTime, extChan, hitE, fSlo]}

    # get centroids
    # for run in extData:
    #     ch = extData[run][2]
    #     muE = np.mean(extData[run][3])
    #     stdE = np.std(extData[run][3])
    #     muFS = np.mean(extData[run][4])
    #     stdFS = np.std(extData[run][4])
    #     print("run %d  ch %d  E %.2f pm %.2f  FS %.2f pm %.2f" % (run, ch, muE, stdE, muFS, stdFS))
    # run 7236  ch 624  E 54.78 pm 0.14  FS 74.80 pm 1.41
    # run 7249  ch 688  E 46.83 pm 0.11  FS 67.14 pm 1.81
    # run 7220  ch 674  E 71.53 pm 1.30  FS 67.56 pm 95.30
    # runMap = {7236:624, 7249:688, 7220:674}

    # cent = {}
    # for run in extData:
    #     if run not in runMap: continue
    #     ch = extData[run][2]
    #     fLo, fHi, fpb = 50, 100, 0.5
    #     x, hist = wl.GetHisto(extData[run][4], fLo, fHi, fpb)
    #     fMax = x[np.argmax(hist)]
    #     cent[ch] = fMax
    #     # plt.plot(x, hist, ls='steps')
    #     # plt.axvline(fMax, c='r')
    #     # plt.show()

    # this is the result of the above block
    centMap = {624: 74.75, 688: 67.25, 674: 65.75}

    fig = plt.figure(figsize=(8,6))
    p1 = plt.subplot(211)
    p2 = plt.subplot(212)

    # hitE, fSlo = [], []
    cols = {624:'r', 688:'g', 674:'b'}
    for run in extData:
        ch = extData[run][2]
        hitE = extData[run][3]
        # fSlo = extData[run][4] # unshifted
        fSlo = [fs - centMap[ch] for fs in extData[run][4]] # shifted
        p1.plot(hitE, fSlo, '.', c=cols[ch], ms=0.5, alpha=0.7)

        fLo, fHi, fpb = -50, 200, 2
        x, hist = wl.GetHisto(fSlo, fLo, fHi, fpb)
        fMax = x[np.argmax(hist)]
        muE = np.mean(hitE)
        p2.plot(muE, fMax, ".", c=cols[ch], ms=10)

    for ch in cols:
        cpd = det.getChanCPD(ds, ch)
        p1.plot(np.nan, np.nan, '.', c=cols[ch], label="C%sP%sD%s" % (cpd[0],cpd[1],cpd[2]))
        p2.plot(np.nan, np.nan, '.', c=cols[ch], label="C%sP%sD%s" % (cpd[0],cpd[1],cpd[2]))

    p1.set_ylim(-100, 150) # shifted
    # p1.set_ylim(0, 200) # unshifted

    p1.set_xlim(0, 20)
    p2.set_xlim(0, 20)
    p2.legend()
    p2.set_xlabel("Energy (keV)", ha='right', x=1)
    p2.set_ylabel("fitSlo Centroid", ha='right', y=1)
    p1.set_ylabel("fitSlo", ha='right', y=1)

    plt.tight_layout()
    # plt.show()
    plt.savefig("./plots/lat-extPulser-centroid.pdf")


def ext_pulser_width():

    f = np.load("../data/lat-extPulser.npz")
    extData = f['arr_0'].item()  # {run: [pIdx, runTime, extChan, hitE, fSlo]}

    extInfo = {624:[7236, 74.75, 'r'], 688:[7249, 67.25, 'g'], 674:[7220, 65.75, 'b']}

    def xGaussPDF(x,k,loc,scale,amp):
        """ amp=overall multip factor, loc=mu, sc=sigma, k=tau """
        from scipy.stats import exponnorm
        return amp * exponnorm.pdf(x,k,loc,scale)

    def xGaussCDF(x,k,loc,scale,amp):
        from scipy.stats import exponnorm
        return exponnorm.cdf(x,k,loc,scale)

    fsMax, fsLo, fsHi, fsE, fsR = [], [], [], [], []

    xgFits = {}

    for i, run in enumerate(extData):

        ch = extData[run][2]
        cpd = det.getChanCPD(0, ch)

        # if ch != 674: continue
        # if i != 14: continue
        if run in [7233]: continue

        hitE = extData[run][3]
        fSlo = extData[run][4]
        idx = np.where(np.isfinite(hitE))
        hitE, fSlo = hitE[idx], fSlo[idx]

        muE = np.average(hitE)
        sdE = np.std(hitE)
        if muE < 1.: continue
        if muE > 100: continue

        runTime = extData[run][1]/1e9 # elogs: 20 Hz, 150 sec runs
        pulseRate = 20 # Hz
        nExp = int(runTime * pulseRate)
        nTot = len(hitE)
        tEff = nTot / nExp

        # fSlo = extData[run][4] # unshifted
        fSlo = [fs - extInfo[ch][1] for fs in extData[run][4]] # shifted

        fLo, fHi = np.percentile(fSlo,1) * 1.2, np.percentile(fSlo,97) * 1.2
        if fHi > 50: fHi = 50

        nb = 100
        fpb = (fHi-fLo)/nb

        x, h = wl.GetHisto(fSlo, fLo, fHi, fpb, shift=False)
        fMax = x[np.argmax(h)]
        pct = []
        for p in [10, 90]:
            tmp = np.cumsum(h)/np.sum(h)*100
            idx = np.where(tmp > p)
            pct.append(x[idx][0])
        wid = pct[1]-pct[0]

        fsMax.append(fMax)
        fsLo.append(pct[0])
        fsHi.append(pct[1])
        fsE.append(muE)
        fsR.append(run)

        tau, mu, sig, amp = 10, fMax, 10, 10000
        np.seterr(invalid='ignore', over='ignore')
        popt,_ = curve_fit(xGaussPDF, x, h, p0=(tau, mu, sig, amp)) # tau, mu, sig, amp
        np.seterr(invalid='warn', over='warn')
        xgFits[run] = [muE, popt]

        # print("%d  %d  %d  lo %.2f  hi %.2f  fpb %.2f  E %-5.1f  rt %.1f  %d/%d (%.3f)  FS  10%% %-5.1f  max %-5.1f  90%% %-5.1f  W %.1f" % (i, run, ch, fLo, fHi, fpb, muE, runTime, nTot, nExp, tEff, pct[0], fMax, pct[1], wid))
        # plt.plot(x, h, ls='steps-mid', c=extInfo[ch][2], label="C%sP%sD%s, E %.1f  Max %.1f  Wid %.1f" % \
        #     (cpd[0],cpd[1],cpd[2],muE, fMax,wid))
        # xFunc = np.arange(fLo, fHi, 0.01)
        # tau, mu, sig, amp = popt
        # xgFit = xGaussPDF(xFunc, *popt)
        # xgMax = xFunc[np.argmax(xgFit)]
        # plt.plot(xFunc, xgFit, '-b', label = "xGaus, tau %.1f  mu %.1f  sig %.1f  amp %.1f" % (tau, mu, sig, amp))
        # # get CDF
        # # xgInt = xGaussCDF(xFunc, *popt)
        # # distMax = np.amax(xgFit)
        # # plt.plot(xFunc, xgInt, '-k', label = "xGauss CDF")
        # # plt.axvline(xgMax, c='m', label='xg max %.1f' % xgMax)
        # # plt.axvline(fMax, c='k')
        # plt.axvline(pct[0], c='b')
        # plt.axvline(pct[1], c='g')
        # plt.xlabel("fitSlo (shifted)", ha='right', x=1)
        # plt.legend(loc=1, fontsize=10)
        # plt.tight_layout()
        # plt.show()
        # return

    # ======== plot the fast pulse envelope ========

    fig = plt.figure(figsize=(8,7))
    p1 = plt.subplot(211)
    p2 = plt.subplot(212)

    # set the starting 90% value
    # i 12 muE 18.511 run 7251
    run90 = 7251
    chan90 = extData[run90][2]
    hitE90 = extData[run90][3]
    fSlo90 = [fs - extInfo[chan90][1] for fs in extData[run90][4]] # shifted
    fLo, fHi, fpb = -20, 20, 0.1
    x, h = wl.GetHisto(fSlo90, fLo, fHi, fpb, shift=False)
    fMax = x[np.argmax(h)]
    tmp = np.cumsum(h)/np.sum(h)*100
    idx = np.where(tmp > 90)
    pct90 = x[idx][0]
    # plt.plot(x, h, ls='steps')
    # plt.axvline(fMax,c='r')
    # plt.axvline(pct90,c='g')
    # plt.show()
    # return
    p1.axhline(pct90, c='orange')

    # plot the acceptance of each CDF
    for i, run in enumerate(xgFits):

        muE = xgFits[run][0]
        if muE > 30: continue

        # print(i, muE, run)

        popt = xgFits[run][1]
        tau, mu, sig, amp = popt

        yFunc = np.arange(-80, 150, 0.1)

        # get the cdf and suppress the errors
        np.seterr(invalid='ignore', over='ignore')
        yCDF = xGaussCDF(yFunc, *popt)
        yCDF[ np.where((np.isnan(yCDF)) | (yCDF < 1e-5)) ] = 0 # fix bad values
        np.seterr(invalid='warn', over='warn')
        # plt.plot(yFunc, yCDF)
        # plt.show()
        # return

        # plot each CDF as a color gradient

        y1 = yFunc[np.where(yCDF > 0.005)][0]
        y2 = mu
        y3 = yFunc[np.where(yCDF > 0.95)][0]
        y90 = yFunc[np.where(yCDF > 0.9)][0]



        nPts = 50
        lineX = muE * np.ones(nPts)
        lineY = np.arange(y1, y3, (y3-y1)/nPts)
        if len(lineY) > nPts: lineY = lineY[:nPts]

        cmap = plt.cm.get_cmap('inferno',nPts)
        for i in range(nPts):
            col = xGaussCDF(lineY[i], *popt)
            p1.plot(lineX[i:i+2], lineY[i:i+2], c=cmap(col))

        p1.plot(muE, y2, ".k")
        p1.plot(muE, y90, '.r', ms=10)

        # calculate the acceptance for the given 90% value
        acc90 = xGaussCDF(pct90, *popt)
        print("%i  pct90 %.2f  muE %.2f  acc90 %.2f" % (i, pct90, muE, acc90))

        p2.plot(muE, acc90, '.k')


    p2.set_xlabel("Energy (keV)", ha='right', x=1)
    p1.set_ylabel("fitSlo (shifted)", ha='right', y=1)
    p2.set_ylabel("acc90", ha='right', y=1)

    plt.tight_layout()
    plt.show()


if __name__=="__main__":
    main()