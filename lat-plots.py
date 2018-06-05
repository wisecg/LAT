#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
plt.style.use('./pltReports.mplstyle')

import dsi
bkg = dsi.BkgInfo()
det = dsi.DetInfo()
import waveLibs as wl

def main():

    # spec()
    spec_vs_cpd()


def spec():
    from ROOT import TFile, TChain, TTree

    dsList = [0,1,2,3,4,"5A","5B","5C"]

    tt = TChain("skimTree")
    enrExp, natExp = 0, 0
    for ds in dsList:
        inFile = "%s/bkg/cut/final/final_DS%s.root" % (dsi.dataDir, ds)
        tf = TFile(inFile)
        enrExp += float(tf.Get("enrExp (kg-d)").GetTitle())
        natExp += float(tf.Get("natExp (kg-d)").GetTitle())
        tf.Close()
        tt.Add(inFile)

    print("enrExp %.2f  natExp %.2f, ds" % (enrExp/365.25, natExp/365.25),dsList)

    fig = plt.figure(figsize=(8,7))
    p0 = plt.subplot(211)
    p1 = plt.subplot(212)
    xLo, xHi, xpb = 0, 20, 0.1

    # natural
    tCut = "!isEnr"
    n = tt.Draw("trapENFCal",tCut,"goff")
    hitE = tt.GetV1()
    hitE = [hitE[i] for i in range(n)]
    x, h1 = wl.GetHisto(hitE, xLo, xHi, xpb, shift=False)
    h1 = np.divide(h1, natExp*xpb) # scale by exposure and binning
    p0.plot(x, h1, 'b', ls='steps', lw=2, label="Natural")
    # p0.set_xlabel("Energy (keV)", ha='right', x=1)
    p0.set_ylabel("Counts/keV-kg-d", ha='right', y=1)
    p0.legend(loc=1)

    # enriched
    tCut = "isEnr"
    n = tt.Draw("trapENFCal",tCut,"goff")
    hitE = tt.GetV1()
    hitE = [hitE[i] for i in range(n)]
    x, h1 = wl.GetHisto(hitE, xLo, xHi, xpb, shift=False)
    h1 = np.divide(h1, enrExp*xpb) # scale by exposure and binning
    p1.plot(x, h1, 'b', ls='steps', lw=2, label="Enriched")
    p1.set_xlabel("Energy (keV)", ha='right', x=1)
    # p1.set_ylabel("Counts/keV-kg-d", ha='right', y=1)
    p1.legend(loc=1)

    plt.tight_layout()
    plt.show()


def spec_vs_cpd():
    from ROOT import TFile, TChain, TTree

    dsList = [0,1,2,3,4,"5A","5B","5C"]

    tt = TChain("skimTree")
    enrExp, natExp = 0, 0
    for ds in dsList:
        inFile = "%s/bkg/cut/final/final_DS%s.root" % (dsi.dataDir, ds)
        tf = TFile(inFile)
        enrExp += float(tf.Get("enrExp (kg-d)").GetTitle())
        natExp += float(tf.Get("natExp (kg-d)").GetTitle())
        tf.Close()
        tt.Add(inFile)
    print("enrExp %.2f  natExp %.2f, ds" % (enrExp/365.25, natExp/365.25),dsList)

    # the channel map changes across datasets, so we have to plot by CPD

    cpdList = []
    for ds in dsList:
        dsTmp = int(ds[0]) if isinstance(ds,str) else ds
        chTmp = det.getGoodChanList(dsTmp)
        for ch in chTmp:
            cpdList.append(int(det.getChanCPD(dsTmp,ch)))
    cpdList = sorted(list(set(cpdList)))

    enrList = [cpd for cpd in cpdList if det.allDetIDs[str(cpd)] > 100000]
    cpdEnrMap = {enrList[i]:i for i in range(len(enrList))}
    # enrLabels = enrList

    natList = [cpd for cpd in cpdList if det.allDetIDs[str(cpd)] < 100000]
    cpdNatMap = {natList[i]:i for i in range(len(natList))}
    # natLabels =


    fig = plt.figure(figsize=(8,7))
    p0 = plt.subplot(211)
    p1 = plt.subplot(212)
    xLo, xHi, xpb = 0, 20, 0.2

    # natural
    tCut = "!isEnr"
    yLo, yHi = 0, len(natList)
    nbx, nby = int((xHi-xLo)/xpb), len(natList)

    n = tt.Draw("trapENFCal:C:P:D",tCut,"goff")
    hitE, hitC, hitP, hitD = tt.GetV1(), tt.GetV2(), tt.GetV3(), tt.GetV4()
    hitE = [hitE[i] for i in range(n)]
    hitCPD = [cpdNatMap[ int("%d%d%d" % (hitC[i],hitP[i],hitD[i])) ] for i in range(n)]

    hNat,_,_,im0 = p0.hist2d(hitE, hitCPD, bins=[nbx, nby], range=[[xLo,xHi],[yLo,yHi]], cmap='jet')
    # p0.set_xlabel("Energy (keV)", ha='right', x=1.)
    p0.set_xticks(np.arange(xLo, xHi+1, 1.0))
    p0.set_ylabel("CPD, Natural", ha='right', y=1.)
    p0.set_yticks(np.arange(0, len(natList))+0.5)
    p0.set_yticklabels(natList, fontsize=8)

    # enriched
    tCut = "isEnr"
    yLo, yHi = 0, len(enrList)
    nbx, nby = int((xHi-xLo)/xpb), len(enrList)

    n = tt.Draw("trapENFCal:C:P:D",tCut,"goff")
    hitE, hitC, hitP, hitD = tt.GetV1(), tt.GetV2(), tt.GetV3(), tt.GetV4()
    hitE = [hitE[i] for i in range(n)]
    hitCPD = [ cpdEnrMap[int("%d%d%d" % (hitC[i],hitP[i],hitD[i])) ] for i in range(n)]

    hEnr,_,_,im1 = p1.hist2d(hitE, hitCPD, bins=[nbx, nby], range=[[xLo,xHi],[yLo,yHi]], cmap='jet')
    p1.set_xlabel("Energy (keV)", ha='right', x=1.)
    p1.set_xticks(np.arange(xLo, xHi+1, 1.0))
    p1.set_ylabel("CPD, Enriched", ha='right', y=1.)
    p1.set_yticks(np.arange(0, len(enrList))+0.5)
    p1.set_yticklabels(enrList, fontsize=8)

    cb0 = fig.colorbar(im0, ax=p0)
    cb1 = fig.colorbar(im1, ax=p1)
    # cb0.set_label('cts', ha='right', rotation=270, labelpad=20)


    plt.tight_layout()
    # plt.show()
    plt.savefig("./plots/lat-spec-vs-cpd.pdf")

    # ========= get detector counts =========

    eLo, eHi = 10, 11
    bLo, bHi = int(eLo/xpb), int(eHi/xpb)

    print("Natural Counts, %d-%d keV" % (eLo, eHi))
    for cpd in natList:
        col = cpdNatMap[cpd]
        print(cpd, int(np.sum(hNat[bLo:bHi,col])))

    print("Enriched Counts, %d-%d keV" % (eLo, eHi))
    for cpd in enrList:
        col = cpdEnrMap[cpd]
        print(cpd, int(np.sum(hEnr[bLo:bHi,col])))



if __name__=="__main__":
    main()