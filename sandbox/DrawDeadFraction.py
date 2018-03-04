#!/usr/bin/env python
import os, math, ROOT, glob
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
sns.set(style='darkgrid', context='poster')

if __name__ == '__main__':

    # inDir = "/global/projecta/projectdirs/majorana/users/mbuuck/sim/MJDWithGrahamTDL_byDetESmearing/5.0_TDL/0.75_transition_point/0.50_transition_level/MJDemonstrator/linesource/M1CalSource/A224_Z88"
    # inDir = '/mnt/mjdDisk1/Majorana/users/psz/MAGE/micah'
    inDir = os.environ['SLURM_TMP']

    fileList = sorted(glob.glob("%s/*.root" % inDir))
    print("nFiles:",len(fileList))

    simTree = ROOT.TChain("simTree")
    for f in fileList: simTree.Add(f)

    theCut = 'fNWaveforms==2 && fTotalEnergy/fActiveness > 0.237 && fTotalEnergy/fActiveness < 0.24'
    # theCut = 'fNWaveforms==2'

    bins, bmin, bmax = 60,0,240
    ebins = np.linspace(bmin,bmax,bins)

    nPass1 = simTree.Draw('fEnergy*1000/fActiveness', theCut, 'goff')
    nEnergy1 = simTree.GetV1()
    Total = list(float(nEnergy1[n]) for n in range(nPass1))

    nPass2 = simTree.Draw('fEnergy*1000/fActiveness', theCut + '&& fActiveness == 1', 'goff')
    nEnergy2 = simTree.GetV1()
    Bulk = list(float(nEnergy2[n]) for n in range(nPass2))

    nPass3 = simTree.Draw('fEnergy*1000/fActiveness', theCut + '&& fActiveness < 1', 'goff')
    nEnergy3 = simTree.GetV1()
    Transition = list(float(nEnergy3[n]) for n in range(nPass3))

    fig1, ax1 = plt.subplots(figsize=(10,7))
    ax1.set_yscale('log')
    sns.distplot(Total, bins=ebins, kde=False, label='Total', ax=ax1)
    sns.distplot(Bulk, bins=ebins, kde=False, label='Bulk', ax=ax1)
    sns.distplot(Transition, bins=ebins, kde=False, label='Transition', ax=ax1)
    # sns.distplot(Total, kde=False, label='Total', ax=ax1)
    # sns.distplot(Bulk, kde=False, label='Bulk', ax=ax1)
    # sns.distplot(Transition, kde=False, label='Transition', ax=ax1)
    ax1.set_xlabel('Energy (keV)')
    ax1.legend()

    hTransition, hTransitionEdge = np.histogram(Transition, bins=bins, range=(bmin,bmax))
    hTotal, hTotalEdge = np.histogram(Total, bins=bins, range=(bmin,bmax))
    print (hTransition)
    print (hTotal)
    print (np.divide(hTransition,hTotal, dtype=float) )
    print(type(hTransition), hTransition.shape, ebins.shape)
    fig2, ax2 = plt.subplots(figsize=(10,7))
    # ax2.set_yscale('log')
    sns.regplot(x=ebins, y=100*np.divide(hTransition,hTotal, dtype=float), ax=ax2)
    ax2.set_ylabel('Percentage of Transition Layer Events')
    ax2.set_xlabel('Energy (keV)')
    plt.tight_layout()
    # plt.show()
    fig1.savefig('TransitionFractionHist.png')
    fig2.savefig('TransitionFractionPercentage.png')
