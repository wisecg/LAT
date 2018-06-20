#!/usr/bin/env python3
import numpy as np
import glob
import dsi
import waveLibs as wl

def main():

    # printOffsets()
    diffOffset()
    loadDiff()


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


def diffOffset():
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

    eList = sorted(list(set(iEnt))) # we do this b/c the list isn't insanely huge
    for iE in eList:
        idx = np.where(iEnt==iE)
        tOffTmp = tOff[idx]
        chTmp = chan[idx]

        for i, ch in enumerate(chTmp):
            if ch%2==1: continue

            iM = np.where(chTmp==ch+1)
            if len(iM[0])==0: continue

            diff = tOffTmp[i] - tOffTmp[iM] # HG - LG
            chDiff[ch].append(tOffTmp[i] - tOffTmp[iM])

    for ch in chList:
        chDiff[ch] = np.asarray(chDiff[ch])

    np.savez("../data/tOff-%d.npz" % run, chDiff)


def loadDiff():

    run = 13387
    f = np.load("../data/tOff-%d.npz" % run)
    chDiff = f['arr_0'].item()

    print(len(chDiff))
    for ch in chDiff:
        print(ch, len(chDiff[ch]))


if __name__=="__main__":
    main()