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
    return a + b*np.log(1+np.exp((x - c)/d) )

def calibres(x, a, b, c):
    return np.sqrt(a+b*x+c*x*x)

if __name__ == "__main__":
    skimTree = ROOT.TChain("skimTree")
    # Pick arbitrary range
    skimTree.Add("~/project/cal-lat/latSkimDS1_run941*_*.root")

    # DS1 good channels
    # chList = [578, 580, 582, 592, 598, 600, 608, 610, 626, 632, 640, 648, 664, 672, 690, 692]
    chList = [626]
    for ch in chList:
        nPass = skimTree.Draw("trapENFCalC:riseNoise", "trapENFCal>5 && trapENFCal<200 && (mHL==1 || mHL==2) && isGood && !muVeto && !(C==1&&isLNFill1) && !(C==2&&isLNFill2) && C!=0&&P!=0&&D!=0 && gain == 0 && channel==%d"%(ch), "goff")
        nEnergy = skimTree.GetV1()
        nCut = skimTree.GetV2()
        nCutList = list(float(nCut[n]) for n in xrange(nPass))
        nEnergyList = list(float(nEnergy[n]) for n in xrange(nPass))
        popt, _ = curve_fit(softplus, nEnergyList, nCutList)
        popt2, _ = curve_fit(calibres, nEnergyList, nCutList)

        nTest = np.linspace(0, 250, 2000)
        print popt
        print popt2
        yFit = softplus(nTest, *popt)
        yFit2 = calibres(nTest, *popt2)
        g2 = sns.JointGrid(x=np.array(nEnergyList), y=np.array(nCutList), size=10, space=0.2)
        g2.plot_joint(sns.kdeplot)
        plt.plot(nTest, yFit, "-", color='red')
        plt.plot(nTest, yFit2, "-", color='blue')
        plt.title('a: %.1f  b: %.1f  c: %.1f  d: %.1f'%(popt[0], popt[1], popt[2], popt[3]))
        g2.ax_marg_y.set_axis_off()
        _ = g2.ax_marg_x.hist(np.array(nEnergyList), color='g', alpha=0.8, bins=np.linspace(0, 250, 500))
        plt.show()
        g2.savefig("SoftPlus_riseNoise_ch%d.png"%(ch))
