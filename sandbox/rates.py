# -*- coding: utf-8 -*-
#!/usr/bin/env python3
import os, imp
import numpy as np
import matplotlib.pyplot as plt
plt.style.use('../pltReports.mplstyle')
from matplotlib.colors import LogNorm

dsi = imp.load_source('dsi',os.environ['LATDIR']+'/dsi.py')
wl = imp.load_source('waveLibs',os.environ['LATDIR']+'/waveLibs.py')
bkg = dsi.BkgInfo()
det = dsi.DetInfo()
cal = dsi.CalInfo()

from ROOT import TFile, TChain, TTree

def main():

    dcrRates()
    # lowRates()


def getExpoLTD(dsList, useBlind=False):
    """ get detector exposures from parsing the livetime document tables"""

    with open("{}/data/ds06-exp.txt".format(os.environ['LATDIR'])) as f:
        lines = f.readlines()

    detList = det.allDets
    totExpo = {cpd:0 for cpd in detList}

    for line in lines[2:]:

        tmp = line.rstrip().split()
        if len(tmp)==0: continue
        if tmp[0]=="DS":
            dsNum = int(tmp[1])
            ds = str("5%s" % tmp[2]) if tmp[1]=="5" else int(tmp[1])
            isBlind = True if "Blind" in tmp else False
            # print("New DS:",ds,"Blind?",isBlind)
            continue

        if ds not in dsList: continue
        if isBlind and not useBlind: continue
        if not tmp[0].isdigit(): continue

        chan = int(tmp[0])
        if chan%2==1:
            # print("LG chan %d, cpd %s" % (chan, cpd))
            chan -= 1 # hack to combine LG exposure with the HG exposure -- only in DS1 open, ok cgw 21/8/2018
            continue
        cpd = det.getChanCPD(dsNum, chan)
        expo = float(tmp[5])
        totExpo[cpd] += expo

    enrExp, natExp = 0, 0
    for cpd in totExpo:
        # print(cpd, totExpo[cpd]/365.25)
        if det.isEnr(cpd): enrExp += totExpo[cpd]
        else: natExp += totExpo[cpd]

    # print("Enr exp (kg-y): %.2f  Nat exp (kg-y): %.2f" % (enrExp/365.25, natExp/365.25))
    # print(d, "Enr exp (kg-d): %.2f  Nat exp (kg-d): %.2f" % (enrExp, natExp))

    return totExpo


def dcrRates():

    tt = TChain("skimTree")
    # tt.Add("/Users/brianzhu/project/skim/light/lightDS*")
    tt.Add("/Users/wisecg/project/bkg/highE/*")

    # all detectors in the skims
    n = tt.Draw("(C*100+P*10+D)","","goff")
    detList = tt.GetV1()
    detList = [str(int(detList[i])) for i in range(n)]
    detList = sorted(list(set(detList)))

    # hits w/ positive dcr (alphas), 1--5.5 MeV
    n1 = tt.Draw("trapENFCal:dcr99:(C*100+P*10+D)","trapENFCal>1000 && trapENFCal<5500 && avse>-1 && dcr99 > 0 && isGood","goff")
    hit1, dcr1, det1 = tt.GetV1(), tt.GetV2(), tt.GetV3()
    hit1 = np.asarray([hit1[i] for i in range(n1)])
    dcr1 = np.asarray([dcr1[i] for i in range(n1)])
    det1 = np.asarray([str(int(det1[i])) for i in range(n1)])

    # events 2.7--5.5 MeV passing DCR, these are probably alphas too
    n2 = tt.Draw("trapENFCal:dcr99:(C*100+P*10+D)","trapENFCal>2700 && trapENFCal<5500 && avse>-1 && dcr99 < 0 && isGood","goff")
    hit2, dcr2, det2 = tt.GetV1(), tt.GetV2(), tt.GetV3()
    hit2 = np.asarray([hit2[i] for i in range(n2)])
    dcr2 = np.asarray([dcr2[i] for i in range(n2)])
    det2 = np.asarray([str(int(det2[i])) for i in range(n2)])

    # combine types together
    hit1 = np.append(hit1, hit2)
    dcr1 = np.append(dcr1, dcr2)
    det1 = np.append(det1, det2)

    # debug - 2d dcr vs energy histogram
    # xLo, xHi, xpb = 500, 8000, 5
    # yLo, yHi, ypb = -0.003, 0.005, 0.0001
    # nbx = int((xHi-xLo)/xpb)
    # nby = int((yHi-yLo)/ypb)
    # plt.hist2d(hit1, dcr, bins=[nbx, nby], range=[[xLo, xHi],[yLo,yHi]], norm=LogNorm(), cmap='jet' )
    # plt.show()

    # get detector exposures
    detExpo = getExpoLTD([0,1,2,3,4,"5A","5B","5C",6], useBlind=True)

    # print alpha rates for each detector
    for cpd in detList:

        expo = detExpo[cpd]/365.25
        idx = np.where(det1 == cpd)
        nAlphas = len(idx[0])
        rAlphas = nAlphas/expo if expo > 0 else -1

        if nAlphas > 0 and expo == 0:
            print("Warning, %d alphas w/ 0 exposure" % nAlphas)

        # print("det %s  nA %-6d  exp %-6.2f  r_A %.3f (cts/kg-y)" % (str(cpd), nAlphas, expo, rAlphas))
        print("%s,%d,%.2f,%.3f" % (str(cpd), nAlphas, expo, rAlphas))


def lowRates():
    """ ds1-5c, detector rate at 46 vs rate in 0-5, marking which ones are enriched """

    from ROOT import TFile, TChain, TTree, gROOT
    gROOT.ProcessLine("gErrorIgnoreLevel = 3001;")
    pctTot = 95

    f = np.load("%s/data/lat-expo-efficiency-all-e%d.npz" % (dsi.latSWDir, pctTot))
    xEff0 = f['arr_0']
    detEff, detExpo = f['arr_13'].item(), f['arr_14'].item()
    detList = det.allDets

    dsList = [1,2,3,4,"5A","5B","5C",6]

    eff = {cpd:np.zeros(len(xEff0)) for cpd in detList}
    expo = {cpd:0 for cpd in detList}

    tt = TChain("skimTree")
    for ds in dsList:
        tt.Add("%s/bkg/cut/final%dt/final%dt_DS%s.root" % (dsi.dataDir, pctTot, pctTot, ds))
        for cpd in detList:
            eff[cpd] = np.add(eff[cpd], detEff[ds][cpd])
            expo[cpd] += detExpo[ds][cpd]

    exit()

    r5n, r46n, pbn, pbUn, r5Un = [], [], [], [], []
    r5e, r46e, pbe, pbUe, r5Ue = [], [], [], [], []

    # for cpd in detList[:10]:
    for cpd in detList:
        if expo[cpd]==0: continue
        if cpd in ['254']: continue # not enough counts in peak
        isEnr = det.isEnr(cpd)

        # draw the data
        tCut = "C==%s && P==%s && D==%s" % (cpd[0], cpd[1], cpd[2])
        xLo, xHi, xpb = 1, 50, 0.5

        n = tt.Draw("trapENFCal",tCut,"goff")
        if n == 0: continue
        hitE = tt.GetV1()
        hitE = [hitE[i] for i in range(n)]
        x, hCts = wl.GetHisto(hitE, xLo, xHi, xpb, shift=False)

        thisEff = eff[cpd]
        thisExp = expo[cpd]

        effNorm = (np.amax(thisEff)/365.25) / (thisExp/365.25)
        thisEff = effNorm * np.divide(thisEff, np.amax(thisEff))
        idxE = np.where((xEff0 <= xHi) & (xEff0 >= xLo))
        xEff, thisEff = xEff0[idxE], thisEff[idxE]
        hSpec = np.divide(hCts, thisExp * xpb) # scale by exposure and binning to get cts/(keV kg d)
        hErr = np.asarray([np.sqrt(hBin /(thisExp * xpb)) for hBin in hSpec]) # statistical error in each bin

        # get the efficiency-corrected rates

        # 1--5 keV
        eLo, eHi = 1, 5
        idxE = np.where((xEff >= eLo) & (xEff <= eHi))
        effCorr = 1 - 1 / ((xEff[1] - xEff[0]) * np.sum(thisEff[idxE]))
        idxR = np.where((x >= eLo) & (x <= eHi))
        rate5 = np.sum(xpb * hSpec[idxR])/ (eHi-eLo) / effCorr
        rate5Unc = rate5 * np.sqrt(np.sum(hCts[idxR])) / (np.sum(hCts[idxR])) / effCorr
        if rate5 < 0: continue

        if isEnr:
            r5e.append(rate5)
            r5Ue.append(rate5Unc)
        else:
            r5n.append(rate5)
            r5Un.append(rate5Unc)

        # 46--47 keV
        eLo, eHi = 45.5, 47.5
        idxE = np.where((xEff >= eLo) & (xEff <= eHi))
        effCorr = 1 - 1 / ((xEff[1] - xEff[0]) * np.sum(thisEff[idxE]))
        idxR = np.where((x >= eLo) & (x <= eHi))
        rate46 = np.sum(xpb * hSpec[idxR])/ (eHi-eLo) / effCorr
        rate46Unc = rate46 * np.sqrt(np.sum(hCts[idxR])) / (np.sum(hCts[idxR])) / effCorr

        if isEnr:
            r46e.append(rate46)
        else:
            r46n.append(rate46)

        # try to plot the 46.5 -- maybe do sideband analysis to subtract the bkg, and skip the eff corr (it's 95%)
        tCut = "C==%s && P==%s && D==%s" % (cpd[0], cpd[1], cpd[2])
        xLo, xHi, xpb = 40, 55, 0.2
        n = tt.Draw("trapENFCal",tCut,"goff")
        if n == 0: continue
        hitE = tt.GetV1()
        hitE = [hitE[i] for i in range(n)]
        x, hCts = wl.GetHisto(hitE, xLo, xHi, xpb, shift=False)
        hSpec = np.divide(hCts, thisExp * xpb)

        pLo, pHi = 46, 47

        idxB = np.where((x <= pLo) | (x >= pHi))
        rateBkg = np.sum(xpb * hSpec[idxB]) /  (pLo-xLo + xHi-pHi)
        rateBkgU = rateBkg * np.sqrt(np.sum(hCts[idxB]))/np.sum(hCts[idxB])

        idxP = np.where((x >= pLo) & (x <= pHi))
        ratePk = np.sum(xpb * hSpec[idxP]) / (pHi-pLo)

        ratePkU = ratePk * np.sqrt(np.sum(hCts[idxP]))/np.sum(hCts[idxP])

        pkBkg = ratePk / rateBkg
        pkBkgU = pkBkg * np.sqrt( (ratePkU/ratePk)**2 + (rateBkgU/rateBkg)**2 )

        if isEnr:
            pbe.append(pkBkg)
            pbUe.append(pkBkgU)
        else:
            pbn.append(pkBkg)
            pbUn.append(pkBkgU)

        print("cpd %s  r5: %.3f ± %.3f  r46: %.3f ± %.3f  P %.4f  B %.4f  P/B %.4f" % (cpd, rate5, rate5Unc, rate46, rate46Unc, rateBkg, ratePk, pkBkg))

        # plt.axhline(rateBkg, c='g')
        # plt.axhline(ratePk, c='r')
        # plt.step(x, hSpec)
        # plt.show()
        # exit()



    # plot the rates

    fig, (p1, p2) = plt.subplots(1,2, figsize=(9,5))

    p1.plot(r5e, pbe, ".b", label="Enriched")
    p1.errorbar(r5e, pbe, yerr=pbUe, xerr=r5Ue, color='k', linewidth=0.8, fmt='none')
    p1.set_xlabel(r"$\mathregular{r_{1-5}}$ [cts/kg-d]", ha='right', x=1)
    p1.set_ylabel(r"$\mathregular{r_{46.5}\ /\ r_{bkg}}$", ha='right', y=1)
    p1.legend(loc=1)

    p2.plot(r5n, pbn, ".r", label="Natural")
    p2.errorbar(r5n, pbn, yerr=pbUn, xerr=r5Un, color='k', linewidth=0.8, fmt='none')
    p2.set_xlabel(r"$\mathregular{r_{1-5}}$ [cts/kg-d]", ha='right', x=1)
    p2.set_ylabel(r"$\mathregular{r_{46.5}\ /\ r_{bkg}}$", ha='right', y=1)
    p2.legend(loc=1)

    plt.tight_layout()
    # plt.show()
    plt.savefig("%s/plots/rates46.pdf" % dsi.latSWDir)


if __name__=="__main__":
    main()
