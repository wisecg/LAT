#!/usr/bin/env python
import sys, time, os, ROOT
import numpy as np
from scipy.stats import ks_2samp
import seaborn as sns
import matplotlib.pyplot as plt
sns.set_style('darkgrid')

def main(argv):

    # dsList, module = [0,1,2,3,5,6], 1
    dsList, module = [4,5,6], 2

    theCut = "isGood && !wfDCBits && !(isLNFill1 && C==1) && !(isLNFill2&&C==2) && isEnr && !muVeto && mHL==1 && C==%d && globalTime > 0"%(module)

    nuCut = theCut + "&&trapENFCalC>1000&&trapENFCalC<1400&&avse>-1&&dcr99<0"
    dcrCut = theCut + "&&trapENFCalC>2350&&trapENFCalC<3350&&avse>-1&&dcr99>=0"
    alphaCut = theCut + "&&trapENFCalC>4000&&trapENFCalC<8000&&avse>-1"

    h2nList = []
    halphaList = []
    pval = []

    sortnutot,sortnublindtot = [], []
    pnutot, pnublindtot = [], []
    sortalphatot, sortalphablindtot = [], []
    palphatot, palphablindtot = [], []
    ksresultTot, ksblindresultTot  = [], []
    kseresultTot = []
    endTime, endTimeblind = [], []
    sortetot,sorteblindtot = [], []
    petot, peblindtot = [], []

    print "Using cut for two nu: ", nuCut
    print "Using cut for alpha: ", alphaCut
    print "Using cut for dcr: ", dcrCut

    for idx,ds in enumerate(dsList):
        print "Scanning dataset", ds
        skimOpen = ROOT.TChain("skimTree")
        skimOpenAlpha = ROOT.TChain("skimTree")
        skimOpenEnergy = ROOT.TChain("skimTree")
        if ds == 5:
            skimOpen.Add("/Users/brianzhu/project/skim/GAT-v02-00-51-g69c5025/skimDS%d_*.root"%(ds))
            skimOpenAlpha.Add("/Users/brianzhu/project/skim/GAT-v02-00-51-g69c5025/skimDS%d_*.root"%(ds))
            skimOpenEnergy.Add("/Users/brianzhu/project/skim/GAT-v02-00-51-g69c5025/skimDS%d_*.root"%(ds))

            skimOpen.Add("/Users/brianzhu/project/skim/GAT-v01-07/skimDS{}_*.root".format(ds))
            skimOpenAlpha.Add("/Users/brianzhu/project/skim/GAT-v01-07/skimDS{}_*.root".format(ds))
            skimOpenEnergy.Add("/Users/brianzhu/project/skim/GAT-v01-07/skimDS{}_*.root".format(ds))

            skimOpen.Add("/Users/brianzhu/project/skim/GAT-v02-00-66-gf078278/skimDS{}_*.root".format(ds))
            skimOpenAlpha.Add("/Users/brianzhu/project/skim/GAT-v02-00-66-gf078278/skimDS{}_*.root".format(ds))
            skimOpenEnergy.Add("/Users/brianzhu/project/skim/GAT-v02-00-66-gf078278/skimDS{}_*.root".format(ds))

        elif ds == 6:
            skimOpen.Add("/Users/brianzhu/project/skim/GAT-v02-00-66-gf078278/skimDS{}_*.root".format(ds))
            skimOpenAlpha.Add("/Users/brianzhu/project/skim/GAT-v02-00-66-gf078278/skimDS%d_*.root"%(ds))
            skimOpenEnergy.Add("/Users/brianzhu/project/skim/GAT-v02-00-66-gf078278/skimDS%d_*.root"%(ds))

            skimOpen.Add("/Users/brianzhu/project/skim/GAT-v02-00-74-g794b9f8/skimDS{}_*.root".format(ds))
            skimOpenAlpha.Add("/Users/brianzhu/project/skim/GAT-v02-00-74-g794b9f8/skimDS%d_*.root"%(ds))
            skimOpenEnergy.Add("/Users/brianzhu/project/skim/GAT-v02-00-74-g794b9f8/skimDS%d_*.root"%(ds))

        else:
            skimOpen.Add("/Users/brianzhu/project/skim/GAT-v01-07/skimDS%d_*.root"%(ds))
            skimOpenAlpha.Add("/Users/brianzhu/project/skim/GAT-v01-07/skimDS%d_*.root"%(ds))
            skimOpenEnergy.Add("/Users/brianzhu/project/skim/GAT-v01-07/skimDS%d_*.root"%(ds))

        n2nu = skimOpen.Draw("globalTime", nuCut, "goff")
        nalpha = skimOpenAlpha.Draw("globalTime", alphaCut, "goff")
        ne = skimOpenEnergy.Draw("globalTime", dcrCut, "goff")

        nu = skimOpen.GetV1()
        alpha = skimOpenAlpha.GetV1()
        e = skimOpenEnergy.GetV1()
        nuList = list(int(nu[n]) for n in xrange(n2nu))
        alphaList = list(int(alpha[n]) for n in xrange(nalpha))
        eList = list(int(e[n]) for n in xrange(ne))
        sortnu = np.sort(nuList)
        sortalpha = np.sort(alphaList)
        sorte = np.sort(eList)

        pnu = 1. * np.arange(len(nuList)) / (len(nuList) - 1)
        palpha = 1. * np.arange(len(alphaList)) / (len(alphaList) - 1)
        pe = 1. * np.arange(len(eList)) / (len(eList) - 1)

        ksresult = ks_2samp(nuList, alphaList)
        ksresultTot.append(ksresult)

        kseresult = ks_2samp(nuList, eList)
        kseresultTot.append(kseresult)

        # Save the ending time of each dataset (maximum time of last event of the 3 distributions)
        endTime.append( np.amax( [sortnu[-1], sorte[-1]] ) )
        sortnutot.extend(sortnu)
        pnutot.extend([x+idx for x in pnu])
        sortalphatot.extend(sortalpha)
        palphatot.extend([x+idx for x in palpha])

        sortetot.extend(sorte)
        petot.extend([x+idx for x in pe])


    # This is CDFs if we combine all of the datasets
    # pnuSumTot = 1. * np.arange(len(sortnutot)) / (len(sortnutot) - 1)
    # palphaSumTot = 1. * np.arange(len(sortalphatot)) / (len(sortalphatot) - 1)
    # ksSumResult = ks_2samp(sortnutot, sortgetot)

    nuDate = np.asarray(sortnutot, dtype='datetime64[s]')
    alphaDate = np.asarray(sortalphatot, dtype='datetime64[s]')
    dcrDate = np.asarray(sortetot, dtype='datetime64[s]')
    endDate = np.asarray(endTime, dtype='datetime64[s]')

    fig1, (a11, a12) = plt.subplots(nrows=2, figsize=(15,10))
    a11.step(nuDate, pnutot, color='green', label=r"$2\nu\beta\beta$")
    a11.step(alphaDate, palphatot, color='blue', label="Alphas (4-8 MeV)")
    a12.step(nuDate, pnutot, color='green', label=r"$2\nu\beta\beta$")
    a12.step(dcrDate, petot, color='blue', label="DCR rejected")

    a11.set_xlabel("Date")
    a12.set_xlabel("Date")
    a11.set_ylabel("CDF")
    for idx, ds in enumerate(dsList[:-1]):
        a11.axvline(endDate[idx], color='red', alpha=0.5, linestyle=':')
        a12.axvline(endDate[idx], color='red', alpha=0.5, linestyle=':')
    a11.set_title("Module {} Eriched, All DS -- Alphas (4 - 8 MeV)".format(module))
    a12.set_title("Module {} Eriched, All DS -- DCR Rejected Alphas (2.35 - 3.35 MeV)".format(module))
    a11.legend(loc=4)
    labelText = ""
    for idx, ksres in enumerate(ksresultTot):
        labelText = labelText + "DS%d -- KS statistic: %.3f   p-value: %.3f \n"%(dsList[idx], ksres[0], ksres[1])
    a11.text(0.05, 0.95, labelText, transform=a11.transAxes, fontsize=12, verticalalignment='top' )

    a12.legend(loc=4)
    labelText = ""
    for idx, ksres in enumerate(kseresultTot):
        labelText = labelText + "DS%d -- KS statistic: %.3f   p-value: %.3f \n"%(dsList[idx], ksres[0], ksres[1])
    a12.text(0.05, 0.95, labelText, transform=a12.transAxes, fontsize=12, verticalalignment='top' )
    plt.tight_layout()

    # fig2, a2 = plt.subplots(figsize=(10,6))
    # a2.step(np.asarray(sortetot), petot, color='green', label="Open Data")
    # a2.step(np.asarray(sorteblindtot), peblindtot, color='blue', label="Blind Data")
    # a2.set_xlabel("Energy (keV)")
    # a2.set_ylabel("CDF")
    # a2.legend(loc=4)
    # labelText = ""
    # for idx, ksres in enumerate(kseresultTot):
    #     labelText = labelText + "DS%d -- KS statistic: %.3f   p-value: %.3f \n"%(dsList[idx], ksres[0], ksres[1])
    #
    # a2.text(0.05, 0.95, labelText, transform=a2.transAxes, fontsize=12, verticalalignment='top' )
    # plt.tight_layout()

    plt.show()


if __name__ == "__main__":
	main(sys.argv[1:])
