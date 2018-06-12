#!/usr/bin/env python
"""
    This script compares 46 keV events by detector/total and alpha events (both full peak alphas and rejected by DCR) using open data from DS0-6
"""
import os, imp, ROOT
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.misc import factorial
ds = imp.load_source('dsi',os.environ['LATDIR']+'/dsi.py')
sns.set(style='darkgrid', context='talk')

detInfo = ds.DetInfo()

def main():
    outDir = os.environ['LATDIR'] + '/plots/AlphaRate'
    dsList, module = [6], 1
    # dsList, module = [4,5,6], 2
    # dsList, module = [0, 1, 2, 3, 4,5,6], 2
    # dsList, module = [4], 2

    theCut = "isGood && !wfDCBits && !(isLNFill1 && C==1) && !(isLNFill2&&C==2) && !muVeto && mHL==1 && !(channel==592 && run>=11348 && run<=11488) && !((channel & 0xFFFFFFF0) == 0x2B0 && run >= 4239 && run <= 4436)"
    # pbCut = theCut + "&&trapENFCalC>45.5&&trapENFCalC<47.5"
    # pbCut = theCut + "&& avse>-1 && dcr99<0 && ((trapENFCalC>=1950 && trapENFCalC<=2094) || (trapENFCalC>=2127 && trapENFCalC<=2195) || (trapENFCalC>=2212 && trapENFCalC<=2350))" # Now this stands for ROI
    # pbCut = theCut + "&& avse>-1 && dcr99<0 && ((trapENFCalC>=1950 && trapENFCalC<=2034) || (trapENFCalC>=2044 && trapENFCalC<=2099) || (trapENFCalC>=2109 && trapENFCalC<=2113) || (trapENFCalC>=2123 && trapENFCalC<=2199) || (trapENFCalC>=2209 && trapENFCalC<=2350))" # Now this stands for ROI

    # pbCut = theCut + "&& avse>-1 && dcr99<0 && ((trapENFCalC>=1950 && trapENFCalC<=2099) || (trapENFCalC>=2109 && trapENFCalC<=2113) || (trapENFCalC>=2123 && trapENFCalC<=2199) || (trapENFCalC>=2209 && trapENFCalC<=2350))" # Now this stands for ROI
    pbCut = theCut + "&& avse>-1 && dcr99<0 && trapENFCalC>=1500 && trapENFCalC<=2500" # Now

    dcrCut = theCut + "&&trapENFCalC>2000&&trapENFCalC<4000&&avse>-1&&dcr99>=0"
    alphaCut = theCut + "&&trapENFCalC>4000&&trapENFCalC<8000"

    pbDict = {}
    dcrDict = {}
    alphaDict = {}

    bar_width = 0.20
    fig1, ax1 = plt.subplots(figsize=(15,7))

    pbTot = {}
    dcrTot = {}
    alphaTot = {}

    for idx,ds in enumerate(dsList):
        ax1.cla()

        print ("Scanning dataset {}, module {}".format(ds, module))
        skimOpen = ROOT.TChain("skimTree")
        if ds == 5:
            # DS5a and DS5b
            skimOpen.Add("/Users/brianzhu/project/skim/GAT-v01-07/skimDS{}_*.root".format(ds))
            # DS5b
            # for i in range(80, 113):
                # skimOpen.Add("/Users/brianzhu/project/skim/GAT-v01-07/skimDS{}_{}.root".format(ds,i))
            # DS5c
            skimOpen.Add("/Users/brianzhu/project/skim/GAT-v02-00-51-g69c5025/skimDS{}_*.root".format(ds))
            skimOpen.Add("/Users/brianzhu/project/skim/GAT-v02-01/skimDS{}_*.root".format(ds))
        elif ds == 6:
            skimOpen.Add("/Users/brianzhu/project/skim/GAT-v02-00-66-gf078278/skimDS{}_*.root".format(ds))
            skimOpen.Add("/Users/brianzhu/project/skim/GAT-v02-01/skimDS{}_*.root".format(ds))
        elif ds == 1 or ds == 2:
            skimOpen.Add("/Users/brianzhu/project/skim/GAT-v02-01/skimDS{}_*.root".format(ds))
            skimOpen.Add("/Users/brianzhu/project/skim/GAT-v01-07/skimDS{}_*.root".format(ds))
        else:
            skimOpen.Add("/Users/brianzhu/project/skim/GAT-v01-07/skimDS{}_*.root".format(ds))

        pbDict.setdefault(ds, getTreeList(skimOpen, ['trapENFCalC', 'channel'], pbCut, ds))
        dcrDict.setdefault(ds, getTreeList(skimOpen, ['trapENFCalC', 'channel'], dcrCut, ds))
        alphaDict.setdefault(ds, getTreeList(skimOpen, ['trapENFCalC', 'channel'], alphaCut, ds))

        pbch, pbchCount = np.unique(pbDict[ds]['channel'], return_counts=True)
        dcrch, dcrchCount = np.unique(dcrDict[ds]['channel'], return_counts=True)
        alphach, alphachCount = np.unique(alphaDict[ds]['channel'], return_counts=True)

        # This is for some empty datasets
        pbch = np.asarray(pbch, dtype=np.str)
        dcrch = np.asarray(dcrch, dtype=np.str)
        alphach = np.asarray(alphach, dtype=np.str)

        # Make arrays the same  -- this is probably the most inefficient way at doing this
        pbDiff = np.setdiff1d(dcrch, pbch)
        pbidx = np.searchsorted(pbch, pbDiff)
        pbch = np.insert(pbch, pbidx, pbDiff)
        pbchCount = np.insert(pbchCount, pbidx, np.zeros(len(pbDiff)))

        alphaDiff = np.setdiff1d(dcrch, alphach)
        alphaidx = np.searchsorted(alphach, alphaDiff)
        alphach = np.insert(alphach, alphaidx, alphaDiff)
        alphachCount = np.insert(alphachCount, alphaidx, np.zeros(len(alphaDiff)))

        dcrDiff = np.setdiff1d(pbch, dcrch)
        dcridx = np.searchsorted(dcrch, dcrDiff)
        dcrch = np.insert(dcrch, dcridx, dcrDiff)
        dcrchCount = np.insert(dcrchCount, dcridx, np.zeros(len(dcrDiff)))

        alphaDiff2 = np.setdiff1d(dcrch, alphach)
        alphaidx2 = np.searchsorted(alphach, alphaDiff2)
        alphach = np.insert(alphach, alphaidx2, alphaDiff2)
        alphachCount = np.insert(alphachCount, alphaidx2, np.zeros(len(alphaDiff2)))

        dcrDiff2 = np.setdiff1d(alphach, dcrch)
        dcridx2 = np.searchsorted(dcrch, dcrDiff2)
        dcrch = np.insert(dcrch, dcridx2, dcrDiff2)
        dcrchCount = np.insert(dcrchCount, dcridx2, np.zeros(len(dcrDiff2)))

        pbDiff2 = np.setdiff1d(alphach, pbch)
        pbidx2 = np.searchsorted(pbch, pbDiff2)
        pbch = np.insert(pbch, pbidx2, pbDiff2)
        pbchCount = np.insert(pbchCount, pbidx2, np.zeros(len(pbDiff2)))

        pbLen = len(pbch)
        dcrLen = len(dcrch)
        alphaLen = len(alphach)
        maxLength = max(pbLen, dcrLen, alphaLen)


        xList = np.linspace(1, maxLength, maxLength)
        ax1.bar(xList-bar_width, pbchCount, bar_width, label='Background (350 keV window)')
        ax1.bar(xList, dcrchCount, bar_width, label='DCR Rejected (2 - 4 MeV)')
        ax1.bar(xList+bar_width, alphachCount, bar_width, label='Alpha (4 - 8 MeV)')
        ax1.set_xticks(xList)
        ax1.set_xticklabels(pbch)
        ax1.set_title('DS{} Background Comparison'.format(ds))
        ax1.set_xlabel('Channel')
        ax1.set_ylabel('Counts')
        ax1.legend()
        plt.tight_layout()
        # fig1.savefig(outDir+'/DS{}_BkgComparison.png'.format(ds))

        # print('Channels:',chList)
        print("Background (350 keV window): ", len(pbDict[ds]['trapENFCalC']))
        print("DCR Rejected (2-4 MeV): ", len(dcrDict[ds]['trapENFCalC']))
        print("Alphas (4-8 MeV)", len(alphaDict[ds]['trapENFCalC']))
        for E, Ch in zip(pbDict[ds]['trapENFCalC'], pbDict[ds]['channel']):
            print(Ch, E)

        for ch, count in zip(pbch, pbchCount):
            pbTot.setdefault(ch, []).append(count)
        # plt.show()

    print(pbDict)
    print(pbTot)
    countArr = []
    for ch, countList in pbTot.items():
        print(ch, sum(countList))
        countArr.append(sum(countList))
    print(countArr)

    fig2, ax2 = plt.subplots(figsize=(10,6))
    entries, bin_edges, patches = ax2.hist(countArr,bins=np.arange(0,10,1))
    bin_middles = 0.5*(bin_edges[1:] + bin_edges[:-1])
    ax2.set_title('DS5 and 6 Enriched Poisson Fit (350 keV window)')
    ax2.set_ylabel('Counts/Detector')
    ax2.set_xlabel('Counts')
    # fit with curve_fit
    parameters, cov_matrix = curve_fit(poisson, bin_middles, entries)
    print(parameters)
    print(cov_matrix)
    x_plot = np.linspace(0, 10, 500)
    plt.plot(x_plot, poisson(x_plot, *parameters), 'r-', lw=2)
    plt.tight_layout()
    plt.show()


def getTreeList(tree, parList=None, cut=None, ds=None):
    """
        If event list is too large, might not work. Might work better with a generator
    """
    eventDic = {}
    if parList == None or cut == None:
        print("Parameter or Cut empty")
        return
    nPar = len(parList)
    if nPar > 4:
        print('Too Many Parameters!')
        return
    nEvents = tree.Draw(':'.join(parList), cut, 'goff')
    for idx, par in enumerate(parList):
        if idx == 0:
            parN = tree.GetV1()
        elif idx == 1:
            parN = tree.GetV2()
        elif idx == 2:
            parN = tree.GetV3()
        elif idx == 3:
            parN = tree.GetV4()
        if par == "blah":
            continue
        #     chList = list(parN[n] for n in range(nEvents))
        #     cpdList = ['C{}P{}D{}'.format(*str(detInfo.getChanCPD(ds,ch)))
        #                 if ch%2==0 else 'C{}P{}D{}'.format(*str(detInfo.getChanCPD(ds, ch-1)))
        #                 for ch in chList ]
        #     eventDic.setdefault(par, cpdList)

        else:
            eventDic.setdefault(par, list(parN[n] for n in range(nEvents)))

    return eventDic

# poisson function, parameter lamb is the fit parameter
def poisson(k, lamb, a):
    return a*(lamb**k/factorial(k)) * np.exp(-lamb)


if __name__ == "__main__":
	main()
