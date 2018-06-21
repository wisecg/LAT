#!/usr/bin/env python3
import numpy as np

import matplotlib.pyplot as plt
plt.style.use('../pltReports.mplstyle')

import dsi
import waveLibs as wl

def main():

    # printOffsets()
    # getOffset()
    plotOffset()


def printOffsets():

    from ROOT import GATDataSet, TChain

    ds, cIdx, run = 1, 41, 13387

    det = dsi.DetInfo()
    pMons = det.getPMon(ds)
    print(pMons)
    return

    skim = TChain("skimTree")
    skim.Add("/global/homes/w/wisecg/project/cal/skim/skimDS1_run13387_low.root")
    n = skim.Draw("iEvent:tOffset:channel","mH==1 && tOffset > 100","goff")
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
        hitChans.append([int(gat.channel.at(j)) for j in range(iN)])
        tOffs.append([int(gat.tOffset.at(j)) for j in range(iN)])

        hits = [gat.trapENFCal.at(j) for j in range(iN)]
        hitEs.append(wl.niceList(hits))


    nLim = 200 if n > 200 else n

    for i in range(nLim):
        print("%d  %-5d  %-4d  gNWF %d" % (chan[i], iEnt[i], tOff[i], iNWF[i]), hitChans[i], tOffs[i], hitEs[i])


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

    tLo, tHi, tpb = -2000, 500, 10

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