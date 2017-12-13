import ROOT
import numpy as np
from matplotlib import pyplot as plt
import waveLibs as wl
import DataSetInfo as ds
import math
import seaborn as sns
sns.set(style='whitegrid', context='talk')

"""
    Sample for getting analysis thresholds using Ralph's algorithm:
    0) Make histograms (0.1 keV binning) for all calidx/bkgidx combinations with cuts
    1) Start with threshold for each bkgidx
    2) Walk up until first 0 bin in the bkg histogram

"""

def SaveHistogramsIDX():

    outDir = '/projecta/projectdirs/majorana/users/bxyzhu/LATv2/plots/spectra'
    calDir = '/projecta/projectdirs/majorana/users/wisecg/cal-lat'
    bkgDir = '/projecta/projectdirs/majorana/users/wisecg/bg-lat'
    bkgcutDir = '/projecta/projectdirs/majorana/users/wisecg/cuts'
    calcutDir = '/projecta/projectdirs/majorana/users/bxyzhu/cuts'
    cInfo = ds.CalInfo()
    bins,lower,upper = 2500,0,250
    # Basic cut
    mNum = 1
    cuts = "gain==0 && mHL==%d && isGood && !muVeto && !(C==1&&isLNFill1) && !(C==2&&isLNFill2) && C!=0&&P!=0&&D!=0"%(mNum)
    skimTreeCal = ROOT.TChain("skimTree")
    skimTreeBkg = ROOT.TChain("skimTree")
    skimCutCal = ROOT.TChain("skimTree")
    skimCutBkg = ROOT.TChain("skimTree")
    dsList = [0, 1, 2, 3, 4, 5]
    # dsList = [2]
    outFile = ROOT.TFile(outDir + "/AThresh_mHL%d.root"%(mNum), "RECREATE")
    outFile.cd()

    for iDS, dsNum in enumerate(dsList):
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
            nRangesCal = [0, len(cInfo.master['ds%d_m%d'%(dsNum, modNum)])]
            nRangesBkg = [0, ds.dsMap[dsNum]]
            for calidx in range(nRangesCal[0], nRangesCal[1]):
                print "Drawing DS%d calidx%d mod%d"%(dsNum, calidx, modNum)
                hCutList, hFullList = [], []
                skimTreeCal.Reset()
                calList = cInfo.GetCalList("ds%d_m%d" % (dsNum, modNum), calidx, runLimit=10)
                for run in calList:
                    skimTreeCal.Add("%s/latSkimDS%d_run%d_*.root"%(calDir, dsNum, run))

                for idx, ch in enumerate(chList):
                    # Reset Tree every calidx + ch
                    skimCutCal.Reset()
                    skimCutCal.Add("%s/calfs_rn/calfs_rn-DS%d-%d-ch%d.root"%(calcutDir, dsNum, calidx, ch))
                    if skimCutCal.GetEntries() == 0:
                        hCutList.append(ROOT.TH1D())
                        hFullList.append(ROOT.TH1D())
                        print "Channel %d, calidx %d has no entries, skipping"%(ch, calidx)
                        continue

                    hCutList.append(ROOT.TH1D())
                    hFullList.append(ROOT.TH1D())
                    # Add additional cut here for channel
                    hCutList[idx] = wl.H1D(skimCutCal,bins,lower,upper, "trapENFCal", cuts+"&& channel==%d"%(ch), Title="hCalDS%d_Ch%d_CalIdx%d"%(dsNum,ch,calidx), Name="hDS%d_Ch%d_CalIdx%d"%(dsNum, ch,calidx))
                    hFullList[idx] = wl.H1D(skimTreeCal,bins,lower,upper, "trapENFCal", cuts+"&&channel==%d"%(ch),Title="hFullDS_%d_Ch%d_CalIdx%d"%(dsNum,ch,calidx), Name="hCalFullDS%d_Ch%d_CalIdx%d"%(dsNum,ch,calidx))

                    # Only write channel specific if there are counts
                    if hCutList[idx].Integral() > 0:
                        # Add ch+calidx histograms to ch total histogram if
                        hCutList[idx].Write()
                        hFullList[idx].Write()

            for bkgidx in range(nRangesBkg[0], nRangesBkg[1]+1):
                print "Drawing DS%d bkgidx%d mod%d"%(dsNum, bkgidx, modNum)
                hCutList, hFullList = [], []
                skimTreeBkg.Reset()
                skimTreeBkg.Add("%s/latSkimDS%d_%d_*.root"%(bkgDir, dsNum, bkgidx))

                for idx, ch in enumerate(chList):
                    # Reset Tree every bkgidx + ch
                    skimCutBkg.Reset()
                    skimCutBkg.Add("%s/fs_rn/fs_rn-DS%d-%d-ch%d.root"%(bkgcutDir, dsNum, bkgidx, ch))
                    if skimCutBkg.GetEntries() == 0:
                        hCutList.append(ROOT.TH1D())
                        hFullList.append(ROOT.TH1D())
                        print "Channel %d, bkgidx %d has no entries, skipping"%(ch, bkgidx)
                        continue
                    hCutList.append(ROOT.TH1D())
                    hFullList.append(ROOT.TH1D())
                    # Add additional cut here for channel
                    hCutList[idx] = wl.H1D(skimCutBkg,bins,lower,upper, "trapENFCal", cuts+"&& channel==%d"%(ch), Title="hBkgDS%d_Ch%d_BkgIdx%d"%(dsNum,ch,bkgidx), Name="hDS%d_Ch%d_BkgIdx%d"%(dsNum, ch,bkgidx))
                    hFullList[idx] = wl.H1D(skimTreeBkg,bins,lower,upper, "trapENFCal", cuts+"&&channel==%d"%(ch),Title="hFullDS_%d_Ch%d_BkgIdx%d"%(dsNum,ch,bkgidx), Name="hBkgFullDS%d_Ch%d_BkgIdx%d"%(dsNum,ch,bkgidx))

                    # Only write channel specific if there are counts
                    if hCutList[idx].Integral() > 0:
                        # Add ch+bkgidx histograms to ch total histogram if
                        hCutList[idx].Write()
                        hFullList[idx].Write()

    # Write total histogram and close
    outFile.Close()
    return 0

def GetAnaylsisThreshold(dsNum=2):
    inDir = '/Users/brianzhu/macros/code/LAT/plots/spectra'
    f1 = ROOT.TFile("%s/AThresh_DS%d_mHL1.root"%(inDir, dsNum))

    bkgDict = {}
    for key in f1.GetListOfKeys():
        histName = key.GetName()
        if f1.Get(histName).Integral() == 0: continue
        name = ''.join(c for c in histName.split('_')[1] if c.isdigit())
        idx = ''.join(c for c in histName.split('_')[2] if c.isdigit())
        ch = int(name)
        bkgidx = int(idx)
        if 'BkgIdx' in histName and 'BkgFull' not in histName:
            if bkgidx not in bkgDict.keys(): bkgDict[bkgidx] = {}
            if ch not in bkgDict.keys(): bkgDict[bkgidx][ch] = []
            for xbin in range(f1.Get(histName).GetNbinsX()):
                # Skip first 7 bins, -0.05 to 0.65 keV
                if f1.Get(histName).GetBinCenter(xbin) < 0.70: continue
                if f1.Get(histName).GetBinCenter(xbin) > 50.: continue
                bkgDict[bkgidx][ch].append(f1.Get(histName).GetBinContent(xbin))

    for bkgidx in bkgDict.keys():
        tD = wl.getDBCalRecord("thresh_ds%d_bkgidx%d" % (dsNum, bkgidx))
        for ch in bkgDict[bkgidx].keys():
            x = np.array(bkgDict[bkgidx][ch])
            noWall = True
            if np.cumsum(x[:50])[-1] > 5: noWall = False
            print x
            # Method 1: Find first 0 walking up in energy from threshold
            # Works if no noise wall
            thresh1 = tD[ch][0]
            if noWall:
                thresh1 = tD[ch][0]
                thresh2 = 0.1*np.where(x==0)[0][0]+0.7 # Add 0.65 or 0.7 here?
                aThresh = max(thresh1, thresh2)
                print "Ch%d, idx%d -- Thresh: %.1f -- AThresh: %.1f"%(ch, bkgidx, thresh1, thresh2)
            # Method 2: If there's a noise wall, find the maximum index and then start the walk from there
            elif not noWall:
                amax = np.argmax(x[:50])
                thresh3 = 0.1*(np.where(x[amax:]==0)[0][0]+amax) + 0.7
                print "Noise Wall Exists -- Ch%d, idx%d -- AThresh: %.1f"%(ch, bkgidx, thresh3)


if __name__ == "__main__":
    # SaveHistogramsIDX()
    GetAnaylsisThreshold(3)
