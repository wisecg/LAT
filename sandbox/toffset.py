#!/usr/bin/env python3
import numpy as np

import matplotlib.pyplot as plt
plt.style.use('../clint.mpl')

import dsi
det = dsi.DetInfo()

import waveLibs as wl

def main():

    # printOffsets()
    # plotChanOrphans()
    # getOffset()
    plotOffset()


def printOffsets():

    from ROOT import GATDataSet, TChain

    ds, cIdx, run = 1, 41, 13387

    det = dsi.DetInfo()
    # pMons = det.getPMon(ds) # this cal run isn't getting these
    # print(pMons)
    # return
    cpdList = det.dets["M1"]
    cpdPairs = {cpd:0 for cpd in cpdList}
    cpdOrphanHG = {cpd:0 for cpd in cpdList}
    cpdOrphanLG = {cpd:0 for cpd in cpdList}

    skim = TChain("skimTree")
    skim.Add("/global/homes/w/wisecg/project/cal/skim/skimDS1_run13387_low.root")
    # n = skim.Draw("iEvent:tOffset:channel","mH==1 && tOffset > 100","goff")
    n = skim.Draw("iEvent:tOffset:channel","mH==1","goff")
    iEnt, tOff, chan = skim.GetV1(), skim.GetV2(), skim.GetV3()
    iEnt = [int(iEnt[i]) for i in range(n)]
    tOff = [tOff[i] for i in range(n)]
    chan = [chan[i] for i in range(n)]

    gds = GATDataSet(run)
    gat = gds.GetGatifiedChain(False)

    iNWF, hitChans, tOffs, hitEs = [], [], [], []
    for iE in iEnt:
        gat.GetEntry(iE)
        iN = gat.channel.size()
        iNWF.append(iN)

        hChs = [int(gat.channel.at(j)) for j in range(iN)]
        hitChans.append(hChs)
        tOffs.append([int(gat.tOffset.at(j)) for j in range(iN)])

        hits = [gat.trapENFCal.at(j) for j in range(iN)]
        hitEs.append(wl.niceList(hits))

        pairs = []
        for ch in hChs:
            if ch+1 in hChs: pairs.extend([ch, ch+1])
        hgOrphans = [ch for ch in hChs if ch+1 not in hChs and ch not in pairs and ch%2==0]
        lgOrphans = [ch for ch in hChs if ch-1 not in hChs and ch not in pairs and ch%2==1]
        # print(hChs, "pairs:", pairs, "hg orph", hgOrphans, "lg orph", lgOrphans)

        for ch in pairs:
            if ch%2==1: continue
            cpd = det.getChanCPD(ds,ch)
            cpdPairs[cpd] += 1

        for ch in hgOrphans:
            cpd = det.getChanCPD(ds,ch)
            cpdOrphanHG[cpd] += 1

        for ch in lgOrphans:
            cpd = det.getChanCPD(ds,ch-1)
            cpdOrphanLG[cpd] +=1

    # nLim = 200 if n > 200 else n
    # for i in range(nLim):
        # print("%d  %-5d  %-4d  gNWF %d" % (chan[i], iEnt[i], tOff[i], iNWF[i]), hitChans[i], tOffs[i], hitEs[i])

    np.savez("../data/toffset-orphans.npz", cpdPairs, cpdOrphanHG, cpdOrphanLG)


def plotChanOrphans():

    f = np.load("../data/toffset-orphans.npz")
    cpdPairs, cpdOrphanHG, cpdOrphanLG = f['arr_0'].item(), f['arr_1'].item(), f['arr_2'].item()

    detList = det.dets["M1"]
    detMap = {i:detList[i] for i in range(len(detList))}
    detIdx = np.arange(0, len(detList), 1)

    detCtsP, detCtsOHG, detCtsOLG, detCtsTot = [], [], [], []
    for cpd in detList:
        detCtsP.append(cpdPairs[cpd])
        detCtsOHG.append(cpdOrphanHG[cpd])
        detCtsOLG.append(cpdOrphanLG[cpd])
        detCtsTot.append(cpdPairs[cpd]+cpdOrphanHG[cpd]+cpdOrphanLG[cpd])

    # 610/611 is C1P3D2
    # print(det.getChanCPD(1,610))

    # get bad/veto only?
    # maybe color code by good/bad channels?
    # goodCh = det.getGoodChanList(1)
    # goodDet = [det.getChanCPD(1, ch) for ch in goodCh]

    detCtsP, detCtsOHG, detCtsOLG, detCtsTot = np.asarray(detCtsP), np.asarray(detCtsOHG), np.asarray(detCtsOLG), np.asarray(detCtsTot)
    detIdx = np.asarray(detIdx)
    idx = np.where(detCtsTot > 0)
    detIdx2 = np.arange(0, len(detIdx[idx]), 1)

    nTot = [detCtsP[idx][i] +detCtsOHG[idx][i] + detCtsOLG[idx][i] for i in range(len(detIdx2))]
    pctP = np.asarray([100*detCtsP[idx][i]/nTot[i] for i in range(len(detIdx2))])
    pctH = np.asarray([100*detCtsOHG[idx][i]/nTot[i] for i in range(len(detIdx2))])
    pctL = np.asarray([100*detCtsOLG[idx][i]/nTot[i] for i in range(len(detIdx2))])

    width = 0.3       # the width of the bars: can also be len(x) sequence

    # plt.bar(detIdx2, pctL, width, color='red', log=True, label="LG Orphan")
    # plt.bar(detIdx2, pctH, width, color='green', log=True, bottom=pctL, label="HG Orphan")
    # plt.bar(detIdx2, pctP, width, color='blue', log=True, bottom=pctL+pctH, label="Paired")

    plt.bar(detIdx2, detCtsTot[idx], width*2, color='blue', alpha=0.5, log=True, label="Total Counts")
    plt.bar(detIdx2-width/2, detCtsOLG[idx], width, color='red', log=True, label="LG Orphan")
    plt.bar(detIdx2+width/2, detCtsOHG[idx], width, color='green', log=True, label="HG Orphan")
    # plt.bar(detIdx2, detCtsP[idx], 0.8, color='blue', log=True, label="Paired")


    plt.ylim(ymin=0.7)


    xticks = np.arange(0, len(detIdx2))
    plt.xticks(xticks)
    xlabels = [detMap[i] for i in detIdx[idx]]
    plt.gca().set_xticklabels(xlabels, fontsize=12)

    plt.xlabel("CPD", ha='right', x=1)
    plt.ylabel("Counts", ha='right', y=1)

    plt.legend()
    plt.tight_layout()

    # plt.show()
    plt.savefig("../plots/toffset-pairs-run13387.pdf")


def getOffset():
    """ we're looking for the difference in the offsets between HG and LG waveforms. """

    from ROOT import GATDataSet, TChain

    ds, cIdx, run = 1, 41, 13387

    det = dsi.DetInfo()
    pMons = det.getPMon(ds)

    # skim = TChain("skimTree")
    # skim.Add("/global/homes/w/wisecg/project/cal/skim/skimDS1_run13387_low.root")
    gds = GATDataSet(run)
    gat = gds.GetGatifiedChain(False)

    n = gat.Draw("Entry$:tOffset:channel","trapENFCal < 250","goff")
    iEnt, tOff, chan = gat.GetV1(), gat.GetV2(), gat.GetV3()
    iEnt = np.asarray([int(iEnt[i]) for i in range(n)])
    tOff = np.asarray([tOff[i] for i in range(n)])
    chan = np.asarray([chan[i] for i in range(n)])

    chList = sorted(list(set(chan)))
    chList = [ch for ch in chList if ch%2==0]
    chDiff = {ch:[] for ch in chList}
    chNoPair = {ch:[] for ch in chList}

    eList = sorted(list(set(iEnt))) # we do this b/c the list isn't insanely huge
    for iE in eList:
        idx = np.where(iEnt==iE)
        tOffTmp = tOff[idx]
        chTmp = chan[idx]

        for i, ch in enumerate(chTmp):
            if ch%2==1: continue

            iM = np.where(chTmp==ch+1)
            if len(iM[0])==0:
                chNoPair[ch].append(tOffTmp[i])
                continue

            diff = tOffTmp[i] - tOffTmp[iM] # HG - LG
            chDiff[ch].append(tOffTmp[i] - tOffTmp[iM])

    for ch in chList:
        chDiff[ch] = np.asarray(chDiff[ch])

    np.savez("../data/tOff-%d.npz" % run, chDiff, chNoPair)


def plotOffset():

    ds, run = 1, 13387

    det = dsi.DetInfo()

    f = np.load("../data/tOff-%d.npz" % run)
    chDiff = f['arr_0'].item()
    chNoPair = f['arr_1'].item()

    tLo, tHi, tpb = -5000, 500, 10

    xTot, hTot = wl.GetHisto([], tLo, tHi, tpb)

    # remember, diffs are HG - LG

    # ==== 1. plot tOffset_HG - tOffset_LG ====
    cmap = plt.cm.get_cmap('jet',len(chDiff)+1)
    for i, ch in enumerate(chDiff):
        # print(ch, len(chDiff[ch]))
        x, hDiff = wl.GetHisto(chDiff[ch], tLo, tHi, tpb)
        hTot = np.add(hTot, hDiff)
        cpd = det.getChanCPD(ds, ch)
        plt.semilogy(x, hDiff, ls='steps', lw=2, c=cmap(i), alpha=0.5, label="C%sP%sD%s" % (cpd[0],cpd[1],cpd[2]))

    p = 99.9
    tmp = np.cumsum(hTot)/np.sum(hTot)*100
    idx = np.where(tmp > p)
    x99 = xTot[idx][0]
    plt.plot(xTot, hTot, "k", ls='steps', label="Total, run %d" % run)
    plt.axvline(x99, c='r', lw=5, label="99.9%% value: %d" % x99)
    plt.legend(loc=2, ncol=3, fontsize=12)
    plt.xlabel("HG-LG tOffset (10ns)", ha='right', x=1)
    plt.ylabel("Counts", ha='right', y=1)
    plt.tight_layout()
    # plt.show()
    plt.savefig("../plots/tOffset-run%d.pdf" % run)

    # ==== 2. plot tOffset_HG of any HG hits w/o a paired LG hit ====
    plt.close()

    tLo, tHi, tpb = 0, 10000, 50

    xTot, hTot = wl.GetHisto([], tLo, tHi, tpb)

    for i, ch in enumerate(chNoPair):

        cpd = det.getChanCPD(ds, ch)
        if cpd in ['173', '112']: continue
        print(ch, cpd)

        x, hNoPair = wl.GetHisto(chNoPair[ch], tLo, tHi, tpb)
        hTot = np.add(hTot, hNoPair)

        plt.semilogy(x, hNoPair, ls='steps', lw=2, c=cmap(i), alpha=0.7, label="C%sP%sD%s" % (cpd[0],cpd[1],cpd[2]))

    pctTot = 99
    tmp = np.cumsum(hTot)/np.sum(hTot)*100

    idx = np.where(tmp > pctTot)
    xPctVal = xTot[idx][0]
    plt.plot(xTot, hTot, "k", ls='steps', label="Total, run %d" % run)
    plt.axvline(xPctVal, c='r', lw=5, label="%d%% value: %d" % (pctTot, xPctVal))
    plt.legend(loc=1, ncol=3, fontsize=10)
    plt.xlabel("Unpaired tOffset_HG (10ns)", ha='right', x=1)
    plt.ylabel("Counts", ha='right', y=1)
    plt.tight_layout()
    # plt.show()
    plt.savefig("../plots/tOffset-unpaired-run%d.pdf" % run)



if __name__=="__main__":
    main()