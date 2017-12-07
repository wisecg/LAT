import ROOT
import numpy as np
from matplotlib import pyplot as plt
import waveLibs as wl
import DataSetInfo as ds
import seaborn as sns
sns.set(style='whitegrid', context='talk')

"""
    Saves cut acceptance histograms
    This script is confusing as fk sometimes, whatever
"""

def SaveCalHistograms():

    outDir = '/projecta/projectdirs/majorana/users/bxyzhu/LATv2/plots/spectra'
    calDir = '/projecta/projectdirs/majorana/users/wisecg/cal-lat'
    cutDir = '/projecta/projectdirs/majorana/users/bxyzhu/cuts'
    cInfo = ds.CalInfo()
    bins,lower,upper = 250,0,250
    # Basic cut
    mNum = 2
    cuts = "gain==0 && mHL==%d && isGood && !muVeto && !(C==1&&isLNFill1) && !(C==2&&isLNFill2) && C!=0&&P!=0&&D!=0"%(mNum)
    skimTree = ROOT.TChain("skimTree")
    skimCut = ROOT.TChain("skimTree")
    dsList = [0, 1, 2, 3, 4, 5]
    # dsList = [1]
    outFile = ROOT.TFile(outDir + "/CalAcceptance_wfstd_mHL%d.root"%(mNum), "RECREATE")
    outFile.cd()
    # Total histogram (all datasets)

    hCutTotal = ROOT.TH1D("hCutTotal", "", bins,lower,upper)
    hFullTotal = ROOT.TH1D("hFullTotal", "", bins,lower,upper)
    # Total histogram for Dataset
    hDSTotal = []
    hDSFullTotal = []
    dMissingCh = [0, 1, 2, 3, 4, 5]
    dThreshCutCh = [0, 1, 2, 3, 4, 5]
    for iDS, dsNum in enumerate(dsList):
        hDSTotal.append(ROOT.TH1D())
        hDSTotal[iDS] = ROOT.TH1D("hDS%d"%(dsNum), "", bins,lower,upper)
        hDSFullTotal.append(ROOT.TH1D())
        hDSFullTotal[iDS] = ROOT.TH1D("hDS%d_Full"%(dsNum), "", bins,lower,upper)
        nMods = [1]
        if dsNum == 4: nMods = [2]
        if dsNum == 5: nMods = [1, 2]
        for modNum in nMods:
            chList = ds.GetGoodChanList(dsNum)
            if dsNum==5 and modNum==1: # remove 692 and 1232 (both beges, so who cares)
                chList = [584, 592, 598, 608, 610, 614, 624, 626, 628, 632, 640, 648, 658, 660, 662, 672, 678, 680, 688, 690, 694]
            elif dsNum==5 and modNum==2:
                chList = [1106, 1110, 1120, 1124, 1128, 1170, 1172, 1174, 1176, 1204, 1208, 1298, 1302, 1330, 1332]
            # Total channel histograms, split by dataset
            hChTotalList = []
            hChFullTotalList = []
            for idx, ch in enumerate(chList):
                hChTotalList.append(ROOT.TH1D())
                hChTotalList[idx] = ROOT.TH1D("hDS%d_Ch%d"%(dsNum, ch), "", bins,lower,upper)
                hChFullTotalList.append(ROOT.TH1D())
                hChFullTotalList[idx] = ROOT.TH1D("hDS%d_Ch%d_Full"%(dsNum, ch), "", bins,lower,upper)

            nRanges = [0, len(cInfo.master['ds%d_m%d'%(dsNum, modNum)])]
            for calidx in range(nRanges[0], nRanges[1]):
                print "Drawing DS%d calidx%d mod%d"%(dsNum, calidx, modNum)
                hCutList, hFullList = [], []
                skimTree.Reset()
                calList = cInfo.GetCalList("ds%d_m%d" % (dsNum, modNum), calidx, runLimit=10)
                for run in calList:
                    skimTree.Add("%s/latSkimDS%d_run%d_*.root"%(calDir, dsNum, run))

                for idx, ch in enumerate(chList):
                    # Reset Tree every calidx + ch
                    skimCut.Reset()
                    skimCut.Add("%s/calwf/calwfstd-DS%d-%d-ch%d.root"%(cutDir, dsNum, calidx, ch))
                    if skimCut.GetEntries() == 0:
                        hCutList.append(ROOT.TH1D())
                        hFullList.append(ROOT.TH1D())
                        print "Channel %d, idx %d has no entries, skipping"
                        continue

                    hCutList.append(ROOT.TH1D())
                    hFullList.append(ROOT.TH1D())
                    # Add additional cut here for channel
                    hCutList[idx] = wl.H1D(skimCut,bins,lower,upper, "trapENFCal", cuts+"&& channel==%d"%(ch), Title="hDS%d_Ch%d_%d"%(dsNum,ch,calidx), Name="hDS%d_Ch%d_%d"%(dsNum, ch,calidx))
                    hFullList[idx] = wl.H1D(skimTree,bins,lower,upper, "trapENFCal", cuts+"&&channel==%d"%(ch),Title="hFullDS_%d_Ch%d_%d"%(dsNum,ch,calidx), Name="hFullDS%d_Ch%d_%d"%(dsNum,ch,calidx))

                    # Only write channel specific if there are counts
                    if hCutList[idx].Integral() > 0:
                        # Add ch+calidx histograms to ch total histogram if
                        hCutList[idx].Write()
                        hFullList[idx].Write()
                        hChTotalList[idx].Add(hCutList[idx])
                        hChFullTotalList[idx].Add(hFullList[idx])

                        # Add individual ch+calidx histograms to total histogram and total DS histogram
                        hCutTotal.Add(hCutList[idx])
                        hFullTotal.Add(hFullList[idx])
                        hDSTotal[iDS].Add(hCutList[idx])
                        hDSFullTotal[iDS].Add(hFullList[idx])

        # Write Total channel histograms
        for idx, ch in enumerate(chList):
            if hChTotalList[idx].Integral() > 0:
                hChTotalList[idx].Write()
                hChFullTotalList[idx].Write()
            else:
                print "Channel %d has no entries!"%(ch)

        # Write total DS histograms
        hDSTotal[iDS].Write()
        hDSFullTotal[iDS].Write()

    # Write total histogram and close
    hFullTotal.Write()
    hCutTotal.Write()
    outFile.Close()
    return 0

def CutAcceptance():
    """
        This only works for histograms with the naming convention created by SaveCalHistograms
    """
    inDir = '/projecta/projectdirs/majorana/users/bxyzhu/LATv2/plots/spectra'
    dsList = [0, 1, 2, 3, 4, 5]
    # dsList = [1]
    cutName = 'AllCuts'
    mNum = 2
    f1 = ROOT.TFile("%s/CalAcceptance_%s_mHL%d.root"%(inDir, cutName, mNum))
    for dsNum in dsList:
        fig, ax = plt.subplots(figsize=(10,6))
        ax.cla()
        fig2, ax2 = plt.subplots(figsize=(10,6))
        ax2.cla()
        cutDict = {}
        fullDict = {}

        # Loop over all histograms in file
        for key in f1.GetListOfKeys():
            # Check if histogram is in keys
            if "hDS%d"%(dsNum) in key.GetName():
                # Skip empty histograms
                if f1.Get(key.GetName()).Integral() == 0: continue

                # Skip calidx specific histograms
                if 'FullDS' in key.GetName(): continue
                elif "Full" not in key.GetName() and len(key.GetName().split('_')) > 2: continue

                # Save only digits in name -- format is "hDS%d_Ch%d_%d"
                name = ''.join(c for c in key.GetName() if c.isdigit())
                name = name[1:] # Skips the dsNum
                # name = ''.join(c for c in key.GetName().split('_') if c.isdigit())

                # If a number exists, otherwise it will be empty
                if name:
                    nameStr = str(ds.CPD[dsNum][int(name)])
                    # nameList.append("C%sP%sD%s (%s)"%(nameStr[0], nameStr[1], nameStr[2], name)) # Name
                    nameStr = "C%sP%sD%s (%s)"%(nameStr[0], nameStr[1], nameStr[2], name)
                    # Histograms without cuts
                    f1.Get(key.GetName()).Rebin(5)
                    if "Full" in key.GetName():
                        fullDict[nameStr] = []
                        for xbin in range(f1.Get(key.GetName()).GetNbinsX()):
                            fullDict[nameStr].append(f1.Get(key.GetName()).GetBinContent(xbin))

                    # Histograms with cuts
                    else:
                        cutDict[nameStr] = []
                        for xbin in range(f1.Get(key.GetName()).GetNbinsX()):
                            cutDict[nameStr].append(f1.Get(key.GetName()).GetBinContent(xbin))

        # print fullDict
        # print cutDict
        idx = 0
        for nameStr in cutDict.keys():
            print 'Drawing DS%d %s'%(dsNum, nameStr)
            # resize to every 10 bins (keV)
            dCut = np.asarray(cutDict[nameStr], dtype=np.float)
            dFull = np.asarray(fullDict[nameStr], dtype=np.float)
            accept = dCut/dFull
            energy = np.linspace(0, 250, len(accept))
            ax.plot(energy, accept*100, label=nameStr, color = plt.cm.tab20c(idx))
            ax.set_title('DS%s Fraction of Events Retained (%s)'%(dsNum, cutName))
            ax.set_ylabel("Fraction of Events Retained/(5 keV) (%)")
            ax.set_xlabel("Energy (keV)")

            ax2.plot(energy[:int(len(accept)*100./250):], accept[:int(len(accept)*100./250):]*100, label=nameStr, color = plt.cm.tab20c(idx))
            ax2.set_title('DS%s Fraction of Events Retained (%s)'%(dsNum, cutName))
            ax2.set_ylabel("Fraction of Events Retained/(5 keV) (%)")
            ax2.set_xlabel("Energy (keV)")

            idx += 1
        plt.tight_layout()
        # Shrink current axis and put legend on the side
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width*0.80, box.height])
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), ncol=1)
        fig.savefig('%s/CutAcc_DS%d_%s_mHL%d.png'%(inDir, dsNum, cutName, mNum))

        box2 = ax2.get_position()
        ax2.set_position([box2.x0, box2.y0, box2.width*0.80, box2.height])
        ax2.legend(loc='center left', bbox_to_anchor=(1, 0.5), ncol=1)
        fig2.savefig('%s/CutAcc_DS%d_%s_mHL%d_Zoom.png'%(inDir, dsNum, cutName, mNum))

        ROOT.gStyle.SetOptStat(0)
        c1 = ROOT.TCanvas("c%d"%(dsNum),"c%d"%(dsNum),800,600)
        c1.SetLogy()
        f1.Get('hDS%d_Full'%(dsNum)).SetLineColor(ROOT.kBlack)
        f1.Get('hDS%d_Full'%(dsNum)).SetTitle('DS%d Calibration(%s)'%(dsNum, cutName))
        f1.Get('hDS%d_Full'%(dsNum)).GetXaxis().SetTitle('Energy (keV)')
        f1.Get('hDS%d_Full'%(dsNum)).GetYaxis().SetTitle('Counts')
        f1.Get('hDS%d_Full'%(dsNum)).SetMinimum(100.)
        f1.Get('hDS%d'%(dsNum)).SetLineColor(ROOT.kRed)
        f1.Get('hDS%d_Full'%(dsNum)).Draw()
        f1.Get('hDS%d'%(dsNum)).Draw("SAME")
        leg1 = ROOT.TLegend(0.35, 0.7, 0.65, 0.88)
        leg1.Clear()
        leg1.AddEntry('hDS%d_Full'%(dsNum), "DS%d Basic Cuts"%(dsNum), "l")
        leg1.AddEntry('hDS%d'%(dsNum), "DS%d with %s"%(dsNum, cutName), "l")
        leg1.Draw()
        c1.SaveAs('%s/CalReducedSpectra_DS%d_%s_mHL%d.pdf'%(inDir, dsNum,cutName,mNum))

    return 0


def CutCounts():
    dsList = [0, 1, 2, 3, 4, 5]
    f1 = ROOT.TFile("/Users/brianzhu/macros/code/LAT/plots/spectra/PrelimSpectra/BkgHisto_AllCuts.root")
    for dsNum in dsList:
        Sum12 = []
        Sum25 = []
        Sum520 = []
        nameList = []

        # Loop over all histograms in file
        for key in f1.GetListOfKeys():
            # Check if histogram is in keys
            if "hDS%d"%(dsNum) in key.GetName():
                if f1.Get(key.GetName()).Integral() == 0: continue
                # Save only digits in name
                name = ''.join(c for c in key.GetName() if c.isdigit())
                name = name[1:]
                # If a number exists, otherwise it will be empty
                if name:
                    nameStr = str(ds.CPD[dsNum][int(name)])
                    nameList.append("C%sP%sD%s (%s)"%(nameStr[0], nameStr[1], nameStr[2], name)) # Name
                # nameList.append(name) # Number
                    Sum12.append(f1.Get(key.GetName()).Integral(f1.Get(key.GetName()).FindBin(1.), f1.Get(key.GetName()).FindBin(2.), 'W')/1.)
                    Sum25.append(f1.Get(key.GetName()).Integral(f1.Get(key.GetName()).FindBin(2.), f1.Get(key.GetName()).FindBin(5.), 'W')/3.)
                    Sum520.append(f1.Get(key.GetName()).Integral(f1.Get(key.GetName()).FindBin(5.), f1.Get(key.GetName()).FindBin(20.), 'W')/15.)

        xList = np.linspace(1, len(nameList), len(nameList))
        bar_width = 0.25
        fig, ax = plt.subplots(figsize=(10,6))
        ax.cla()
        # Add 0.5 to each width to shift to the center of the grid
        plt.bar(xList+0.5-bar_width, Sum12, bar_width, label='1 - 2 keV')
        plt.bar(xList+0.5, Sum25, bar_width, label='2 - 5 keV')
        plt.bar(xList+0.5+bar_width, Sum520, bar_width, label='5 - 20 keV')
        ax.set_xticks(np.linspace(1, len(nameList), len(nameList)))
        # ax.set_xticklabels(nameList, rotation='vertical')
        ax.set_xticklabels(nameList, rotation=50)
        ax.set_title('DS%d %s'%(dsNum, Type))
        ax.set_ylabel("Rate (counts/keV)")
        # ax.set_xlabel("Energy")
        ax.legend()
        plt.tight_layout()
        # Move up plot so that labels will fit
        plt.subplots_adjust(bottom=0.2)
        fig.savefig('CountsAllCuts_DS%d_%s'%(dsNum, Type))

    return 0

if __name__ == "__main__":
    SaveCalHistograms()
    #CutAcceptance()
