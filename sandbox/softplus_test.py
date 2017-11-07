#!/usr/bin/env python
import ROOT
import numpy as np
from scipy.optimize import curve_fit
from matplotlib import pyplot as plt
import seaborn as sns
sns.set(style='whitegrid', context='talk')

"""
    Test of softplus function to fit riseNoise distribution

"""

def softplus(x, a, b, c, d):
    """
        A = Y-offset, B = Slope , C = X shift for flatness, D = Curvature
    """
    return a + b*np.log(1+np.exp((x - c)/d) )

def calibres(x, a, b, c):
    return np.sqrt(a+b*x+c*x*x)

if __name__ == "__main__":
    skimTree = ROOT.TChain("skimTree")
    # Pick arbitrary range
    # skimTree.Add("~/project/cal-lat/latSkimDS1_run940*_*.root")
    # skimTree.Add("~/project/cal-lat/latSkimDS1_run941*_*.root")
    # skimTree.Add("~/project/cal-lat/latSkimDS1_run945*_*.root")
    # skimTree.Add("~/project/cal-lat/latSkimDS1_run946*_*.root")
    # skimTree.Add("~/project/cal-lat/latSkimDS1_run949*_*.root")
    # skimTree.Add("~/project/cal-lat/latSkimDS1_run950*_*.root")
    skimTree.Add("~/project/cal-lat/latSkimDS5_run2118*_*.root")

    theCut = "isGood && !muVeto && !(C==1&&isLNFill1) && !(C==2&&isLNFill2) && C!=0&&P!=0&&D!=0 && gain == 0"

    # DS1 good channels
    # chList = [578, 580, 582, 592, 598, 600, 608, 610, 626, 632, 640, 648, 664, 672, 690, 692]

    # DS5 good channels
    chList = [584, 592, 598, 608, 610, 614, 624, 626, 628, 632, 640, 648, 658, 660, 662, 672, 678, 680, 688, 690, 694]
    # chList = [626]
    # chList = [672]
    for ch in chList:
        nPass1 = skimTree.Draw("trapENFCalC:riseNoise",  theCut + "&& trapENFCal>5 && trapENFCal<50 && channel==%d"%(ch), "goff")
        nCutArray1, nCutArray2, nCutArray3 = [], [], []
        if nPass1 != 0:
            nEnergy1 = skimTree.GetV1()
            nCut1 = skimTree.GetV2()
            nCutList1 = list(float(nCut1[n]) for n in xrange(nPass1))
            nEnergyList1 = list(float(nEnergy1[n]) for n in xrange(nPass1))
            nCutArray1 = [[x,y] for x,y in zip(nCutList1, nEnergyList1) if x > np.percentile(nCutList1, 5) and x < np.percentile(nCutList1, 85)]

        nPass2 = skimTree.Draw("trapENFCalC:riseNoise",  theCut+"&& trapENFCal>50 && trapENFCal<150 && channel==%d"%(ch), "goff")
        if nPass2 != 0:
            nEnergy2 = skimTree.GetV1()
            nCut2 = skimTree.GetV2()
            nCutList2 = list(float(nCut2[n]) for n in xrange(nPass2))
            nEnergyList2 = list(float(nEnergy2[n]) for n in xrange(nPass2))
            nCutArray2 = [[x,y] for x,y in zip(nCutList2, nEnergyList2) if x > np.percentile(nCutList2, 5) and x < np.percentile(nCutList2, 90)]

        nPass3 = skimTree.Draw("trapENFCalC:riseNoise",  theCut+"&& trapENFCal>150 && trapENFCal<240 && channel==%d"%(ch), "goff")
        if nPass3 != 0:
            nEnergy3 = skimTree.GetV1()
            nCut3 = skimTree.GetV2()
            nCutList3 = list(float(nCut3[n]) for n in xrange(nPass3))
            nEnergyList3 = list(float(nEnergy3[n]) for n in xrange(nPass3))
            nCutArray3 = [[x,y] for x,y in zip(nCutList3, nEnergyList3) if x > np.percentile(nCutList3, 5) and x < np.percentile(nCutList3, 90)]

        nCutArray = np.asarray(nCutArray1 + nCutArray2 + nCutArray3)
        if len(nCutArray) == 0:
            print "No events, skipping"
            continue
        # popt, _ = curve_fit(softplus, nEnergyList, nCutList)
        popt, _ = curve_fit(softplus, nCutArray[:,1], nCutArray[:,0], p0=[1., 0.005, 10, 0.5] , bounds = ((-10,0,-10,0),(10,10,150,100)))
        nTest = np.linspace(0, 250, 1000)
        print popt
        yFit = softplus(nTest, *popt)
        # g2 = sns.JointGrid(x=np.array(nEnergyList), y=np.array(nCutList), size=10, space=0.2)
        g2 = sns.JointGrid(x=nCutArray[:,1], y=nCutArray[:,0], size=10, space=0.2)
        g2.plot_joint(sns.kdeplot)
        plt.plot(nTest, yFit, "-", color='red')
        # plt.plot(nTest, yFit2, "-", color='blue')
        plt.title('Y-Offset: %.2f  Slope: %.3f  X-Shift: %.2f  Curvature: %.2f'%(popt[0], popt[1], popt[2], popt[3]))
        g2.ax_marg_y.set_axis_off()
        _ = g2.ax_marg_x.hist(np.array(nEnergyList1 + nEnergyList2 + nEnergyList3), alpha=0.8, bins=np.linspace(0, 250, 500))
        # plt.show()
        g2.savefig("/Users/brianzhu/macros/code/LAT/plots/SoftPlus/Test4SoftPlus_riseNoise_ch%d.png"%(ch))
