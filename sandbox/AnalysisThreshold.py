#!/usr/bin/env python
import os, math, ROOT
import numpy as np
from matplotlib import pyplot as plt
import DataSetInfo as ds
import waveLibs as wl
import seaborn as sns
import Exposure as ex
sns.set(style='darkgrid', context='talk')

"""
    Sample for getting analysis thresholds using Ralph's algorithm:
    0) Make histograms (0.1 keV binning) for all calidx/bkgidx combinations with cuts
    1) Start with threshold for each bkgidx
    2) Walk up until first 0 bin in the bkg histogram

"""

def main():
    # SaveHistogramsIDX()
    # GetAnaylsisThreshold(4, True)
    # GenerateCorrectedSpectra(dsNum = 3, dType = 'isEnr')
    # dsNumList = [0,1,2,3,4,5]
    # dTypeList = ['isEnr', 'isNat']
    dsNumList = [1]
    dTypeList = ['isNat']
    for dsNum in dsNumList:
        for dType in dTypeList:
            GenerateCorrectedSpectra(dsNum, dType)
            # MovingAve(dsNum, dType)



def NoiseFunc(x, a, tau, sigma):
    return 1./(4*tau)*np.exp(-(x*x)/(2*sigma*sigma))


def SaveHistogramsIDX(dType='Bkg', binsize = 0.1, lower = 0, upper = 250):
    """
        Saves background or calibration data before and after cuts into histograms
    """

    outDir = '/projecta/projectdirs/majorana/users/bxyzhu/LATv2/plots/spectra'
    calDir = '/projecta/projectdirs/majorana/users/wisecg/cal-lat'
    bkgDir = '/projecta/projectdirs/majorana/users/wisecg/bg-lat'
    bkgcutDir = '/projecta/projectdirs/majorana/users/wisecg/cuts'
    calcutDir = '/projecta/projectdirs/majorana/users/bxyzhu/cuts'
    cInfo = ds.CalInfo()
    bins = int((upper-lower)/binsize)
    # Basic cut
    mNum = 1
    cuts = "gain==0 && mHL=={} && isGood && !muVeto && !(C==1&&isLNFill1) && !(C==2&&isLNFill2) && C!=0&&P!=0&&D!=0"%(mNum)
    skimTreeCal = ROOT.TChain("skimTree")
    skimTreeBkg = ROOT.TChain("skimTree")
    skimCutCal = ROOT.TChain("skimTree")
    skimCutBkg = ROOT.TChain("skimTree")
    dsList = [0, 1, 2, 3, 4, 5]
    outFile = ROOT.TFile(outDir + "/AThresh_mHL{}.root"%(mNum), "RECREATE")
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
            nRangesCal = [0, len(cInfo.master['ds{}_m{}'.format(dsNum, modNum)])]
            nRangesBkg = [0, ds.dsMap[dsNum]]

            if dType == 'Cal':
                for calidx in range(nRangesCal[0], nRangesCal[1]):
                    print ("Drawing DS{} calidx{} mod{}".format(dsNum, calidx, modNum))
                    hCutList, hFullList = [], []
                    skimTreeCal.Reset()
                    calList = cInfo.GetCalList("ds{}_m{}" .format (dsNum, modNum), calidx, runLimit=10)
                    for run in calList:
                        skimTreeCal.Add("{}/latSkimDS{}_run{}_*.root".format(calDir, dsNum, run))

                    for idx, ch in enumerate(chList):
                        # Reset Tree every calidx + ch
                        skimCutCal.Reset()
                        if not os.path.exists("{}/calfs_rn/calfs_rn-DS{}-{}-ch{}.root".format(calcutDir, dsNum, calidx, ch)):
                            hCutList.append(ROOT.TH1D())
                            hFullList.append(ROOT.TH1D())
                            print ("Channel {}, calidx {} doesn't exist, skipping".format(ch, calidx))
                            continue
                        skimCutCal.Add("{}/calfs_rn/calfs_rn-DS{}-{}-ch{}.root".format(calcutDir, dsNum, calidx, ch))
                        hCutList.append(ROOT.TH1D())
                        hFullList.append(ROOT.TH1D())

                        # Add additional cut here for channel
                        hCutList[idx] = wl.H1D(skimCutCal,bins,lower,upper, "trapENFCal", cuts+"&& channel=={}".format(ch), Title="hCalDS{}_Ch{}_CalIdx{}".format(dsNum,ch,calidx), Name="hDS{}_Ch{}_CalIdx{}".format(dsNum, ch,calidx))
                        hFullList[idx] = wl.H1D(skimTreeCal,bins,lower,upper, "trapENFCal", cuts+"&&channel=={}".format(ch),Title="hFullDS_{}_Ch{}_CalIdx{}".format(dsNum,ch,calidx), Name="hCalFullDS{}_Ch{}_CalIdx{}".format(dsNum,ch,calidx))

                        # Write all histograms (even if they're empty -- for debugging purposes)
                        hCutList[idx].Write()
                        hFullList[idx].Write()
            if dType == 'Bkg':
                for bkgidx in range(nRangesBkg[0], nRangesBkg[1]+1):
                    print ("Drawing DS{} bkgidx{} mod{}".format(dsNum, bkgidx, modNum))
                    hCutList, hFullList = [], []
                    skimTreeBkg.Reset()
                    skimTreeBkg.Add("{}/latSkimDS{}_{}_*.root".format(bkgDir, dsNum, bkgidx))

                    for idx, ch in enumerate(chList):
                        # Reset Tree every bkgidx + ch
                        skimCutBkg.Reset()
                        if not os.path.exists("{}/fs_rn/fs_rn-DS{}-{}-ch{}.root".format(bkgcutDir, dsNum, bkgidx, ch)):
                            hCutList.append(ROOT.TH1D())
                            hFullList.append(ROOT.TH1D())
                            print ("Channel {}, bkgidx {} has no entries, skipping".format(ch, bkgidx))
                            continue
                        skimCutBkg.Add("{}/fs_rn/fs_rn-DS{}-{}-ch{}.root".format(bkgcutDir, dsNum, bkgidx, ch))

                        hCutList.append(ROOT.TH1D())
                        hFullList.append(ROOT.TH1D())

                        # Add additional cut here for channel
                        hCutList[idx] = wl.H1D(skimCutBkg,bins,lower,upper, "trapENFCal", cuts+"&& channel=={}".format(ch), Title="hBkgDS{}_Ch{}_BkgIdx{}".format(dsNum,ch,bkgidx), Name="hDS{}_Ch{}_BkgIdx{}".format(dsNum, ch,bkgidx))
                        hFullList[idx] = wl.H1D(skimTreeBkg,bins,lower,upper, "trapENFCal", cuts+"&&channel=={}".format(ch),Title="hFullDS_{}_Ch{}_BkgIdx{}".format(dsNum,ch,bkgidx), Name="hBkgFullDS{}_Ch{}_BkgIdx{}".format(dsNum,ch,bkgidx))

                        # Write all histograms -- for debugging
                        hCutList[idx].Write()
                        hFullList[idx].Write()

    # Write total histogram and close
    outFile.Close()
    return 0

def GetAnaylsisThreshold(dsNum=2, savePlot=False):
    """
        Returns a dictionary for analysis threshold, organized by: athreshDict[bkgidx][ch]
    """
    inDir = '/Users/brianzhu/macros/code/LAT/plots/spectra'
    f1 = ROOT.TFile("{}/AThresh_mHL1.root".format(inDir))

    if savePlot: fig, ax = plt.subplots(figsize=(10,6))

    bkgDict, bkgFullDict, athreshDict = {}, {}, {}
    for key in f1.GetListOfKeys():
        histName = key.GetName()
        if 'DS{}'.format(dsNum) not in histName: continue
        name = ''.join(c for c in histName.split('_')[1] if c.isdigit())
        idx = ''.join(c for c in histName.split('_')[2] if c.isdigit())
        ch = int(name)
        bkgidx = int(idx)
        if 'BkgIdx' in histName and 'BkgFull' not in histName:
            for xbin in range(f1.Get(histName).GetNbinsX()):
                # If greater than 50 keV -- skip
                if f1.Get(histName).GetBinCenter(xbin) > 50.: continue
                bkgFullDict.setdefault(bkgidx, {}).setdefault(ch, []).append(f1.Get(histName).GetBinContent(xbin))
                # Skip first 7 bins, -0.05 to 0.65 keV
                if f1.Get(histName).GetBinCenter(xbin) < 0.70: continue
                bkgDict.setdefault(bkgidx, {}).setdefault(ch, []).append(f1.Get(histName).GetBinContent(xbin))

    for bkgidx in bkgDict:
        print ("Accessing thresh_ds{}_bkgidx{}".format(dsNum, bkgidx) )
        tD = ds.getDBRecord("thresh_ds{}_bkgidx{}".format(dsNum, bkgidx))
        print tD
        for ch in bkgDict[bkgidx]:
            if ch not in tD:
                print ("Ch{} bkgidx{} not in thresholds database".format(ch, bkgidx))
                continue
            x = np.array(bkgDict[bkgidx][ch])
            xFull = np.array(bkgFullDict[bkgidx][ch])
            noWall = True
            if np.cumsum(x[:50])[-1] >= 4: noWall = False
            # Method 1: Find first 0 walking up in energy from threshold
            # Works if no noise wall
            thresh1 = math.ceil(tD[ch][0]*10)/10 # Round up to the 1st decimal
            sigma = tD[ch][1]
            aThresh = 0.
            if noWall:
                thresh2 = 0.1*np.where(x==0)[0][0]+0.7 # Add 0.65 or 0.7 here?
                aThresh = max(thresh1+3*sigma, thresh2)
                athreshDict.setdefault(bkgidx, {}).setdefault(ch, aThresh)
                # print "Ch%d, idx%d -- Thresh: %.1f -- AThresh: %.1f"%(ch, bkgidx, thresh1+3*sigma, thresh2)
                print ("Ch{}, idx{} -- Thresh: {:.1f} -- AThresh: {:.1f}".format(ch, bkgidx, thresh1+3*sigma, thresh2))

                if savePlot:
                    ax.cla()
                    xbins = np.linspace(0, len(xFull)/10., len(xFull))
                    ax.plot(xbins[:100], xFull[:100], ls="steps") # vals, bins instead of x-value, y-value
                    ax.plot([aThresh, aThresh],[0, np.amax(x[:100])],'k-', color='r') # Draw Line at analysis thresh
                    fig.savefig('/Users/brianzhu/macros/code/LAT/plots/AThresh/DS{}/NoWall_DS{}_ch{}_idx{}.png'.format(dsNum, dsNum, ch, bkgidx))

            # Method 2: If there's a noise wall, find the maximum index and then start the walk from there
            elif not noWall:
                amax = np.argmax(x[:50])
                thresh3 = 0.1*(np.where(x[amax:]==0)[0][0]+amax) + 0.7
                aThresh = max(thresh1, thresh3)
                athreshDict.setdefault(bkgidx, {}).setdefault(ch, aThresh)
                print ("Noise Wall Exists -- Ch{}, idx{} -- AThresh: {:.1f}".format(ch, bkgidx, thresh3))

                if savePlot:
                    ax.cla()
                    xbins = np.linspace(0, len(xFull)/10., len(xFull))
                    ax.plot(xbins[:100], xFull[:100], ls="steps") # vals, bins instead of x-value, y-value
                    ax.plot([aThresh, aThresh],[0, np.amax(x[:100])],'k-', color='r') # Draw Line at analysis thresh
                    fig.savefig('/Users/brianzhu/macros/code/LAT/plots/AThresh/DS{}/Wall_DS{}_ch{}_idx{}.png'.format(dsNum, dsNum, ch, bkgidx))


    return athreshDict

def GenerateCorrectedSpectra(dsNum = 1, dType = 'isNat', binsize = 0.1, binsize2 = 0.001, lower = 0, upper = 250):
    """
        Calculates analysis threshold and applies analysis threshold to exposure calculation
        Saves histograms into ROOT file
    """

    ROOT.gStyle.SetOptStat(0)
    bgDir, outDir = '/Users/brianzhu/project/cuts/fs_rn','/Users/brianzhu/macros/code/LAT/plots/AThresh'
    bins = int((upper-lower)/binsize)
    bins2 = int((upper-lower)/binsize2)
    cuts = "{} && gain==0 && mHL==1 && isGood && !muVeto && !(C==1&&isLNFill1) && !(C==2&&isLNFill2) && C!=0&&P!=0&&D!=0".format(dType)

    chList = ds.GetGoodChanList(dsNum, dType[2:])
    if dsNum == 5: # remove 692 and 1232 (both beges, so who cares)
        if 692 in chList: chList.remove(692)
        if 1232 in chList: chList.remove(1232)

    nRanges = [0, ds.dsMap[dsNum]]
    if dsNum == 5: nRanges[0] = 80 # exclude DS-5A

    athresh = GetAnaylsisThreshold(dsNum, True)
    threshDict = {}
    threshDictSave = {}
    specIDXDict = {}
    specDict = {}
    outFile = ROOT.TFile(outDir + '/Bkg_{}_DS{}_Test.root'.format(dType,dsNum), "RECREATE")
    cutTree = ROOT.TChain("skimTree")
    TotalSpec = ROOT.TH1D('DS{}_{}_Corr'.format(dsNum, dType),'',bins,lower,upper)
    UncorrTotalSpec = ROOT.TH1D('DS{}_{}_UnCorr'.format(dsNum, dType),'',bins,lower,upper)
    for bkgidx in range(nRanges[0], nRanges[1]+1):
        # Get Threshold dictionary
        tD = ds.getDBRecord("thresh_ds{}_bkgidx{}".format(dsNum, bkgidx))

        # Get Analysis Threshold Dictionary and Exposure Dictionary here
        for idx, ch in enumerate(chList):
            if ch in tD:
                # Reset Tree
                cutTree.Reset()

                # Check dictionaries to see if variables exist
                if ch not in athresh[bkgidx]:
                    print ("Warning: Analysis threshold doesn't exist for ch{} bkgidx{}".format(ch,bkgidx))
                    continue

                if ch not in ex.Exposure[dsNum][bkgidx]:
                    print ("Warning: Exposure doesn't exist for ch{} bkgidx{}".format(ch,bkgidx))
                    continue

                if not os.path.exists(bgDir + "/fs_rn-DS{}-{}-ch{}.root".format(dsNum,bkgidx,ch)):
                    print ("Warning: Background data doesn't exist for ch{} bkgidx{}".format(ch,bkgidx))
                    continue

                # Load Background Data with cuts applied
                cutTree.Add(bgDir + "/fs_rn-DS{}-{}-ch{}.root".format(dsNum,bkgidx,ch))

                # Create Trigger Efficiency function
                mu, sigma = tD[ch][0], tD[ch][1]
                threshFnc = ROOT.TF1("fEff_{}_{}_{}".format(dsNum, ch,bkgidx),"0.5*(1+TMath::Erf((x-[0])/(TMath::Sqrt(2)*[1])))", 0, 250)
                threshFnc.SetParameters(mu, abs(sigma))

                # Create Analysis Threshold + Exposure function
                dExp, dAThresh = ex.Exposure[dsNum][bkgidx][ch], athresh[bkgidx][ch]
                expFnc = ROOT.TF1("fExp_{}_{}_{}".format(dsNum,ch,bkgidx), "[0]*(x>[1])",0,250)
                expFnc.SetParameters(dExp, dAThresh)

                # Print out info
                print ("DS{} Ch{} BkgIdx{}  Exposure {}  Thresh {}  Analysis Threshold {} ".format(dsNum, ch, bkgidx, dExp, mu, dAThresh))

                specIDXDict.setdefault(bkgidx, {}).setdefault(ch, ROOT.TH1D())
                specIDXDict[bkgidx][ch] = wl.H1D(cutTree,bins,lower,upper, "trapENFCal", cuts+"&& trapENFCal>{:.2f} && channel=={}".format(dAThresh, ch), Title="hBkg_Ch{}_Bkgidx{}".format(ch,bkgidx), Name="hBkg_Ch{}_Bkgidx{}".format(ch,bkgidx))

                # Exposure function for scaling
                h3 = ROOT.TH1D("hCh{}_Bkgidx{}_Scale".format(ch, bkgidx), "", bins,lower,upper)
                for i in range(h3.GetNbinsX()+1):
                    if i < dAThresh*10: continue # Round up to set analysis threshold
                    h3.SetBinContent(i, dExp)
                # Scale here -- trying to make it legible instead of messy
                h3.Multiply(threshFnc)
                threshDict.setdefault(ch, h3.Clone("Efficiency_ch{}".format(ch)))
                threshDict[ch].Add(h3)

                # Save exposure function to histogram
                h4 = ROOT.TH1D("hCh{}_Bkgidx{}".format(ch, bkgidx), "", bins2,lower,upper)
                for i in range(h4.GetNbinsX()+1):
                    if i < dAThresh*1000: continue # Round up to set analysis threshold
                    h4.SetBinContent(i, dExp)
                # Scale here
                h4.Multiply(threshFnc)
                threshDictSave.setdefault(ch, h4.Clone("hEff_DS{}_ch{}".format(dsNum, ch)))
                threshDictSave[ch].Add(h4)

                if specIDXDict[bkgidx][ch].Integral() == 0:
                    print ("Ch {} has 0 counts, not saving".format(ch))
                    continue

                specDict.setdefault(ch, specIDXDict[bkgidx][ch].Clone('hBkg_Ch{}'.format(ch)))
                specDict[ch].Add(specIDXDict[bkgidx][ch])
                specIDXDict[bkgidx][ch].Write()

    EffTot = ROOT.TH1D("DS{}_EffTot_Divide".format(dsNum), "DS{} {}".format(dsNum, dType), bins,lower,upper)
    EffTotSave = ROOT.TH1D("DS{}_{}_EffTot".format(dsNum, dType), "DS{} {}".format(dsNum, dType), bins2,lower,upper)

    for ch in specDict:
        if ch not in threshDict:
            print ("Ch {} not in threshDict... you should've fixed this!".format(ch))
            continue
        TotalSpec.Add(specDict[ch])
        UncorrTotalSpec.Add(specDict[ch])
        EffTot.Add(threshDict[ch])
        EffTotSave.Add(threshDictSave[ch])

        # Save Channel specific
        h1 = specDict[ch].Clone('hDS{}_ChTot_Ch{}'.format(dsNum,ch))
        h1.SetTitle('hDS{}_ChTot_Ch{}'.format(dsNum,ch))
        threshDictSave[ch].Write()
        h1.Write()

    print ('Channels in DS{} -- {}'.format(dsNum, dType), specDict.keys())

    # Save all histograms with fancy names and axes
    EffTotSave.GetYaxis().SetTitle('Exposure (kg-day)')
    EffTotSave.GetXaxis().SetTitle('Energy (keV)')
    EffTotSave.Write()
    UncorrTotalSpec.GetXaxis().SetTitle('Energy (keV)')
    UncorrTotalSpec.GetYaxis().SetTitle('Counts')
    UncorrTotalSpec.Write()
    TotalSpec.Divide(EffTot)
    TotalSpec.GetXaxis().SetTitle('Energy (keV)')
    TotalSpec.GetYaxis().SetTitle('Counts/(0.1 keV)/kg/day')
    TotalSpec.Write()
    outFile.Close()


def MovingAve(dsNum = 1, dType = 'isNat'):
    inDir = '/Users/brianzhu/macros/code/LAT/plots/AThresh'

    chList = ds.GetGoodChanList(dsNum, dType[2:])
    if dsNum==5: # remove 692 and 1232 (both beges, so who cares)
        if 692 in chList: chList.remove(692)
        if 1232 in chList: chList.remove(1232)

    f1 = ROOT.TFile("{}/Bkg_{}_DS{}.root".format(inDir,dType,dsNum))
    bkgDict = {}
    bkgTot = []
    for key in f1.GetListOfKeys():
        histName = key.GetName()
        if "UnCorr" in histName:
            for xbin in range(f1.Get(histName).GetNbinsX()): bkgTot.append(f1.Get(histName).GetBinContent(xbin))
        if 'ChTot' not in histName: continue
        name = ''.join(c for c in histName.split('_')[2] if c.isdigit())
        ch = int(name)
        if ch not in bkgDict.keys(): bkgDict[ch] = []
        for xbin in range(f1.Get(histName).GetNbinsX()):
            bkgDict[ch].append(f1.Get(histName).GetBinContent(xbin))

    # Moving average over 1 keV
    fig, ax = plt.subplots(figsize=(10,7))
    for ch in bkgDict.keys():
        ax.cla()
        test = np.cumsum(bkgDict[ch][:200])
        moveAve = (test[5:] - test[:-5])
        moveAve2 = (test[3:] - test[:-3])
        moveAve3 = (test[7:] - test[:-7])
        moveAve4 = (test[10:] - test[:-10])
        ax.plot(np.linspace(0,len(bkgDict[ch][:200])/10.,len(bkgDict[ch][:200])), np.array(bkgDict[ch][:200]), ls='steps', label='Spectrum')
        ax.plot(np.linspace(0.3, len(moveAve2)/10., len(moveAve2)), np.array(moveAve2), label='Moving Average (0.3 keV ave)', linestyle="-")
        ax.plot(np.linspace(0.5, len(moveAve)/10., len(moveAve)), np.array(moveAve), label='Moving Average (0.5 keV ave)', linestyle="-.")
        ax.plot(np.linspace(0.7, len(moveAve3)/10., len(moveAve3)), np.array(moveAve3), label='Moving Average (0.7 keV ave)', linestyle="--")
        ax.legend()
        fig.savefig("{}/MovingAve_DS{}_{}_Ch{}.png".format(inDir,dsNum,dType[2:],ch))

    ax.cla()
    test = np.cumsum(bkgTot[:200])
    moveAve = (test[5:] - test[:-5])
    moveAve2 = (test[3:] - test[:-3])
    moveAve3 = (test[7:] - test[:-7])
    ax.plot(np.linspace(0,len(bkgTot[:200])/10.,len(bkgTot[:200])), np.array(bkgTot[:200]), ls='steps', label='Spectrum')
    ax.plot(np.linspace(0.3, len(moveAve2)/10., len(moveAve2)), np.array(moveAve2), label='Moving Average (0.3 keV ave)', linestyle="-")
    ax.plot(np.linspace(0.5, len(moveAve)/10., len(moveAve)), np.array(moveAve), label='Moving Average (0.5 keV ave)', linestyle="-.")
    ax.plot(np.linspace(0.7, len(moveAve3)/10., len(moveAve3)), np.array(moveAve3), label='Moving Average (0.7 keV ave)', linestyle="--")
    ax.legend()
    fig.savefig("{}/MovingAve_DS{}_{}_Tot.png".format(inDir,dsNum,dType[2:]))


if __name__ == "__main__":
    main()
