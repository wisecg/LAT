#!/usr/common/usg/software/python/2.7.9/bin/python
import sys, time, os, ROOT
import numpy as np
from scipy.stats import ks_2samp
import seaborn as sns
from matplotlib import pyplot as plt
sns.set_style('whitegrid')

def main(argv):

    dsList = [0, 1, 2, 3, 5]
    # dsList = [4, 5]

    theCut = "gain==0 && isGood && !wfDCBits && !(isLNFill1 && C==1) && !(isLNFill2&&C==2) && isNat && !muVeto && mHL==1 && C==1"

    nuCut = theCut + "&&trapENFCal>1000&&trapENFCal<1400&&avse>-1&&dcr99<0"
    alphaCut = theCut + "&&trapENFCal>2000&&trapENFCal<3000&&avse>-1&&dcr99>=0"
    geCut = theCut + "&& trapENFCal > 9.7 && trapENFCal < 11.1 && trapETailMin < 0.5"

    h2nList = []
    halphaList = []
    pval = []

    sortnutot = []
    pnutot = []
    sortalphatot = []
    palphatot = []
    sortgetot = []
    pgetot = []
    ksresultTot = []

    print "Using cut for two nu: ", nuCut
    print "Using cut for alpha: ", alphaCut
    print "Using cut for Ge68: ", geCut

    for idx,ds in enumerate(dsList):
        print "Scanning dataset", ds

        skimTree = ROOT.TChain("skimTree")
        skimTree.Add("/Users/brianzhu/project/skim/GAT-v01-06-134-g3f44fab/skimDS%d_*.root"%(ds))

        skimTree2 = ROOT.TChain("skimTree")
        skimTree2.Add("/Users/brianzhu/project/skim/GAT-v01-06-134-g3f44fab/skimDS%d_*.root"%(ds))

        skimTree3 = ROOT.TChain("skimTree")
        skimTree3.Add("/Users/brianzhu/project/skim/GAT-v01-06-134-g3f44fab/skimDS%d_*.root"%(ds))

        n2nu = skimTree.Draw("globalTime", nuCut, "goff")
        nalpha = skimTree2.Draw("globalTime", alphaCut, "goff")
        nge = skimTree3.Draw("globalTime", geCut, "goff")

        nu = skimTree.GetV1()
        alpha = skimTree2.GetV1()
        ge = skimTree3.GetV1()

        nuList = list(set(int(nu[n]) for n in xrange(n2nu)))
        alphaList = list(set(int(alpha[n]) for n in xrange(nalpha)))
        geList = list(set(int(ge[n]) for n in xrange(nge)))
        sortnu = np.sort(nuList)
        sortalpha = np.sort(alphaList)
        sortge = np.sort(geList)

        pnu = 1. * np.arange(len(nuList)) / (len(nuList) - 1)
        palpha = 1. * np.arange(len(alphaList)) / (len(alphaList) - 1)
        pge = 1. * np.arange(len(geList)) / (len(geList) - 1)

        ksresult = ks_2samp(nuList, alphaList)
        # ksresult = ks_2samp(nuList, geList)
        ksresultTot.append(ksresult)

        print "Statistic:",ksresult[0], " p-value:", ksresult[1]

        sortnutot.extend(sortnu)
        pnutot.extend([x+idx for x in pnu])

        sortalphatot.extend(sortalpha)
        palphatot.extend([x+idx for x in palpha])

        sortgetot.extend(sortge)
        pgetot.extend([x+idx for x in pge])


    fig = plt.figure(figsize=(9,6))
    a1 = plt.subplot(111)
    a1.plot(sortnutot, pnutot, color='green', label="2nbb")
    a1.plot(sortalphatot, palphatot, color='blue', label="DCR rejected")
    a1.plot(sortgetot, pgetot, color='purple', label="Ge68")
    a1.set_xlabel("UnixTime")
    a1.set_ylabel("CDF")
    plt.title("All DS (M1)")
    a1.legend(loc=4)
    labelText = ""
    for idx, ksres in enumerate(ksresultTot):
        labelText = labelText + "DS%d -- KS statistic: %.3f p-value: %.3f \n"%(dsList[idx], ksres[0], ksres[1])
    a1.text(0.05, 0.95, labelText, transform=a1.transAxes, fontsize=12, verticalalignment='top' )
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
	main(sys.argv[1:])
