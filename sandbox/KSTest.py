#!/usr/common/usg/software/python/2.7.9/bin/python
import sys, time, os, ROOT
import numpy as np
from scipy.stats import ks_2samp
import seaborn as sns
from matplotlib import pyplot as plt
sns.set_style('whitegrid')

def main(argv):


    # dsList, module = [0, 1, 2, 3, 5], 1
    dsList, module = [4, 5], 2

    theCut = "gain==0 && isGood && !wfDCBits && !(isLNFill1 && C==1) && !(isLNFill2&&C==2) && isNat && !muVeto && mHL==1 && C==%d"%(module)

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
    endTime = []

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

        # ksresult = ks_2samp(nuList, alphaList)
        ksresult = ks_2samp(nuList, geList)
        ksresultTot.append(ksresult)

        # Save the ending time of each dataset (maximum time of last event of the 3 distributions)
        endTime.append( np.amax( [sortnu[-1], sortalpha[-1], sortge[-1]] ) )

        sortnutot.extend(sortnu)
        pnutot.extend([x+idx for x in pnu])

        sortalphatot.extend(sortalpha)
        palphatot.extend([x+idx for x in palpha])

        sortgetot.extend(sortge)
        pgetot.extend([x+idx for x in pge])

    # This is CDFs if we combine all of the datasets
    pnuSumTot = 1. * np.arange(len(sortnutot)) / (len(sortnutot) - 1)
    palphaSumTot = 1. * np.arange(len(sortalphatot)) / (len(sortalphatot) - 1)
    pgeSumTot = 1. * np.arange(len(sortgetot)) / (len(sortgetot) - 1)

    ksSumResult = ks_2samp(sortnutot, sortgetot)

    nuDate = np.asarray(sortnutot, dtype='datetime64[s]')
    alphaDate = np.asarray(sortalphatot, dtype='datetime64[s]')
    geDate = np.asarray(sortgetot, dtype='datetime64[s]')
    endDate = np.asarray(endTime, dtype='datetime64[s]')

    fig = plt.figure(figsize=(9,6))
    a1 = plt.subplot(111)
    # Using r in front changes to latex
    # a1.plot(nuDate, pnutot, color='green', label=r"$2\nu\beta\beta$")
    # a1.plot(alphaDate, palphatot, color='blue', label="DCR rejected")

    a1.plot(nuDate, pnuSumTot, color='green', label=r"$2\nu\beta\beta$")
    a1.plot(geDate, pgeSumTot, color='purple', label=r"$^{68}$Ge")

    a1.set_xlabel("Date")
    a1.set_ylabel("CDF")
    for idx, ds in enumerate(dsList[:-1]):
        a1.axvline(endDate[idx], color='red', alpha=0.5, linestyle=':')

    plt.title("All DS (M%d Natural)"%(module))
    a1.legend(loc=4)
    # labelText = ""
    # for idx, ksres in enumerate(ksresultTot):
        # labelText = labelText + "DS%d -- KS statistic: %.3f p-value: %.3f \n"%(dsList[idx], ksres[0], ksres[1])
    labelText = "KS statistic: %.3f p-value: %.3f \n"%(ksSumResult[0], ksSumResult[1])
    a1.text(0.05, 0.95, labelText, transform=a1.transAxes, fontsize=12, verticalalignment='top' )
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
	main(sys.argv[1:])
