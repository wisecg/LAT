#!/usr/bin/env python
import ROOT
import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
import waveLibs as wl
import DataSetInfo as ds
import seaborn as sns
sns.set(style='whitegrid', context='talk')

"""
    Test fitting of wfstd
"""

def ralphFit(x, a, b, c, d, e, f):
    return np.sqrt(np.power(a+ b*x+ c*x*x+ d*x*x*x+ e*x*x*x*x,2)+f*f)

def ralphUpper(x, a, b, c, d, e, f, g, h):
    return np.sqrt(np.power(a+ b*x+ c*x*x+ d*x*x*x+ e*x*x*x*x,2)+f) + g + h*x

def ralphLower(x, a, b, c, d, e, f, g, h):
    return np.sqrt(np.power(a+b*x+c*x*x+d*x*x*x+e*x*x*x*x,2)+f) - g - h*x

if __name__ == "__main__":
    cInfo = ds.CalInfo()
    dsNum, subNum, modNum = 1, 1, 1
    pathToInput = '/projecta/projectdirs/majorana/users/wisecg/cal-lat'
    skimTree = ROOT.TChain("skimTree")
    calList = cInfo.GetCalList("ds%d_m%d" % (dsNum, modNum), subNum, runLimit=10)
    for i in calList: skimTree.Add("%s/latSkimDS%d_run%d_*" % (pathToInput, dsNum, i))
    theCut = "isGood && !muVeto && !(C==1&&isLNFill1) && !(C==2&&isLNFill2) && C!=0&&P!=0&&D!=0 && gain == 0"

    chList = ds.GetGoodChanList(dsNum)
    wfD = wl.getDBCalRecord("wfstd_ds%d_idx%d_mod%d" % (dsNum, subNum, modNum))
    tD = wl.getDBCalRecord("thresh_ds%d_bkgidx%d" % (dsNum, subNum))
    ROOT.gStyle.SetOptStat(0)
    for ch in chList:
        nPass = skimTree.Draw("trapENFCal:wfstd",  theCut + "&& trapENFCal<250 && channel==%d"%(ch), "goff")
        print "Channel %d, %d events"%(ch, nPass)
        if nPass != 0:
            nEnergy = skimTree.GetV1()
            nCut = skimTree.GetV2()
            nCutList = list(float(nCut[n]) for n in xrange(nPass))
            nEnergyList = list(float(nEnergy[n]) for n in xrange(nPass))
            nCutArray = np.asarray([[x,y] for x,y in zip(nCutList, nEnergyList) if x > np.percentile(nCutList, 1) and x < np.percentile(nCutList, 99)])
            nCutArray = np.asarray([[x,y] for x,y in zip(nCutList, nEnergyList)])
        if len(nCutArray) == 0:
            print "No events, skipping"
            continue

        popt, _ = curve_fit(ralphFit, nCutArray[:,1], nCutArray[:,0], p0=[1., 0.1, 0.1, 0.1, 0., tD[ch][1]], bounds=((-10,0,0,0,0,0), (10,10,10,10,10,20)))
        # print popt
        nTest = np.linspace(0, 250, 25000)
        yFit = ralphFit(nTest, *popt)

        g2 = sns.JointGrid(x=nCutArray[:,1], y=nCutArray[:,0], size=10, space=0.2)
        g2.plot_joint(sns.kdeplot)
        g2.ax_marg_y.set_axis_off()
        g2.ax_marg_x.set_axis_off()
        plt.plot(nTest, yFit, "-", color='red', label='Fit')
        # plt.plot(nTest, yFit2, "-", color='blue', label='Lower cut')
        plt.title('DS%d calidx%d ch%d'%(dsNum, subNum, ch))
        plt.legend()
        g2.savefig("/projecta/projectdirs/majorana/users/bxyzhu/LATv2/plots/wfstd/wfstd_band_kde_ch%d.png"%(ch))

        c1 = ROOT.TCanvas("c%d"%(ch),"c%d"%(ch),800,600)
        h1 = ROOT.TH2D("h%d"%(ch), "h%d"%(ch), 100,0,10,100,0,10)
        skimTree.Project('h%d'%(ch), 'wfstd:trapENFCal', theCut + "&& channel==%d"%(ch))
        h1.SetTitle('DS%d calidx%d ch %d'%(dsNum, subNum, ch))
        h1.GetXaxis().SetTitle('Energy (keV)')
        h1.GetYaxis().SetTitle('wfstd')
        print "Parameters", popt

        fUpper = ROOT.TF1("fUpper","sqrt(pow([0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x,2)+[5])+[6]+[7]*x",0,300)
        fLower = ROOT.TF1("fLower","sqrt(pow([0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x,2)+[5])-[6]-[7]*x",0,300)
        fUpper.SetParameters(popt[0],popt[1],popt[2],popt[3],popt[4],popt[5], 0.3, 0.003)
        fLower.SetParameters(popt[0],popt[1],popt[2],popt[3],popt[4],popt[5], 0.3, 0.003)

        c1.SetLogz()
        h1.SetMinimum(1.)
        h1.Draw('colz')
        fUpper.Draw('SAME')
        fLower.Draw('SAME')
        c1.SaveAs("/projecta/projectdirs/majorana/users/bxyzhu/LATv2/plots/wfstd/wfstd_ch%d_DS%d_calidx%d.png"%(ch, dsNum, subNum))
