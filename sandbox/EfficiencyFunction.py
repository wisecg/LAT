import ROOT
import numpy as np
from matplotlib import pyplot as plt
import waveLibs as wl
import DataSetInfo as ds
import math

"""
    Example of applying trigger efficiency to threshold
"""

if __name__ == "__main__":
    ROOT.gStyle.SetOptStat(0)
    dsNum = 2
    chList = ds.GetGoodChanList(dsNum)
    if dsNum==5: # remove 692 and 1232 (both beges, so who cares)
        chList = [584, 592, 598, 608, 610, 614, 624, 626, 628, 632, 640, 648, 658, 660, 662, 672, 678, 680, 688, 690, 694, 1106, 1110, 1120, 1124, 1128, 1170, 1172, 1174, 1176, 1204, 1208, 1298, 1302, 1330, 1332]

    nRanges = [0, ds.dsMap[dsNum]]
    if dsNum==5: nRanges[0] = 80 # exclude DS-5A

    threshDict = {}
    dummyExposure = 1.
    for bkgidx in range(nRanges[0], nRanges[1]+1):
        tD = wl.getDBCalRecord("thresh_ds%d_bkgidx%d" % (dsNum, bkgidx))
        # Get Analysis Threshold Dictionary here
        # Get Exposure Dictionary here
        for idx, ch in enumerate(chList):
            if ch in tD.keys():
                mu, sigma = tD[ch][0], tD[ch][1]
                if mu > 5 or mu == 0:
                    print "Threshold for ch%d bkgidx%d zero or too high, skipping"%(ch, bkgidx)
                    continue

                threshFnc = ROOT.TF1("fEff_%d_%d_%d"%(dsNum, ch,bkgidx),"0.5*(1+TMath::Erf((x-[0])/(TMath::Sqrt(2)*[1])))", 0, 2)
                threshFnc.SetParameters(mu, abs(sigma))

                # If cut efficiency functions exist, add as TF1 here

                # Easiest to just add as histogram... is this the right way to go? In the future we can replace with a giganto-function?
                h3 = ROOT.TH1D("hCh%d_Bkgidx%d"%(ch, bkgidx), "", 20000,0,2)
                for i in range(h3.GetNbinsX()+1):
                    if i < math.ceil(mu*10000): continue # Round up to set analysis threshold
                    h3.SetBinContent(i, dummyExposure)

                # Scale here
                h3.Multiply(threshFnc)

                if ch in threshDict.keys():
                    threshDict[ch].Add(h3)
                else:
                    threshDict[ch] = h3.Clone("Efficiency_ch%d"%(ch))


    # Draw here
    histTot = ROOT.TH1D("DS%dTot"%(dsNum), "DS%d Total"%(dsNum), 20000,0,2)
    for ch, hist in threshDict.iteritems():
        c1 = ROOT.TCanvas("c%d"%(ch), "c%d"%(ch), 800,600)
        histTot.Add(hist)
        hist.Draw()
        c1.SaveAs("TestDS2Eff_ch%d.pdf"%(ch))

    c1 = ROOT.TCanvas("cTot", "cTot", 800,600)
    histTot.GetXaxis().SetTitle("Energy (keV)")
    histTot.GetYaxis().SetTitle("Exposure (kg-yr)")
    histTot.Draw()
    c1.SaveAs("TestD2Eff_Tot.pdf")
