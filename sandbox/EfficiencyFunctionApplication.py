#!/usr/bin/env python
from __future__ import print_function
import ROOT, os
import numpy as np
from matplotlib import pyplot as plt
import waveLibs as wl
import DataSetInfo as ds
import Exposure as ex
import math

"""
    Example code for applying trigger efficiency to threshold
    Also saves exposure efficiency function to ROOT file
    Also makes spectra with analysis thresholds and saves to ROOT file
    Also makes ROOT files with analysis threshold applied

"""

ROOT.gROOT.ProcessLine("gErrorIgnoreLevel = 3001;") # suppress ROOT messages

def GetAnaylsisThreshold(dsNum=3):
    inDir = '/Users/brianzhu/macros/code/LAT/plots/spectra'
    f1 = ROOT.TFile("{}/AThresh_mHL1.root".format(inDir))

    bkgDict, athreshDict = {}, {}
    for key in f1.GetListOfKeys():
        histName = key.GetName()
        if 'DS{}'.format(dsNum) not in histName: continue
        name = ''.join(c for c in histName.split('_')[1] if c.isdigit())
        idx = ''.join(c for c in histName.split('_')[2] if c.isdigit())
        ch = int(name)
        bkgidx = int(idx)
        if 'BkgIdx' in histName and 'BkgFull' not in histName:
            if bkgidx not in bkgDict.keys(): bkgDict[bkgidx] = {}
            if ch not in bkgDict.keys(): bkgDict[bkgidx][ch] = []
            for xbin in range(f1.Get(histName).GetNbinsX()):
                # Skip first 7 bins, -0.05 to 0.65 keV
                if f1.Get(histName).GetBinCenter(xbin) < 0.7: continue
                if f1.Get(histName).GetBinCenter(xbin) > 50.: continue
                bkgDict[bkgidx][ch].append(f1.Get(histName).GetBinContent(xbin))

    for bkgidx in bkgDict.keys():
        tD = wl.getDBCalRecord("thresh_ds{}_bkgidx{}".format(dsNum, bkgidx))
        athreshDict[bkgidx] = {}
        for ch in bkgDict[bkgidx].keys():
            if ch not in tD.keys():
                print ("Ch{} bkgidx{} not in thresholds database".format(ch, bkgidx))
                continue
            x = np.array(bkgDict[bkgidx][ch])
            noWall = True
            if np.cumsum(x[:50])[-1] > 5: noWall = False
            thresh1 = tD[ch][0]

            if noWall:
                thresh2 = 0.1*np.where(x==0)[0][0]+0.7 # Add 0.65 or 0.7 here?
                athreshDict[bkgidx][ch] = max(thresh1, thresh2)

            elif not noWall:
                # print ("Noise Wall Exists, for DS{} ch{} bkgidx{}".format(dsNum, ch, bkgidx))
                amax = np.argmax(x[:50])
                thresh3 = 0.1*(np.where(x[amax:]==0)[0][0]+amax) + 0.7
                athreshDict[bkgidx][ch] = max(thresh1, thresh3)

    return athreshDict

def GenerateCorrectedSpectra(dsNum = 1, dType = 'isNat'):
        ROOT.gStyle.SetOptStat(0)
        bgDir, outDir = '/Users/brianzhu/project/cuts/fs_rn','/Users/brianzhu/macros/code/LAT/plots/spectra/PrelimSpectra'
        bins, bins2, lower, upper = 2500, 250000, 0, 250
        cuts = "{} && gain==0 && mHL==1 && isGood && !muVeto && !(C==1&&isLNFill1) && !(C==2&&isLNFill2) && C!=0&&P!=0&&D!=0".format(dType)

        chList = ds.GetGoodChanList(dsNum, dType[2:])
        if dsNum==5: # remove 692 and 1232 (both beges, so who cares)
            if 692 in chList: chList.remove(692)
            if 1232 in chList: chList.remove(1232)

        nRanges = [0, ds.dsMap[dsNum]]
        if dsNum==5: nRanges[0] = 80 # exclude DS-5A

        athresh = GetAnaylsisThreshold(dsNum)
        threshDict = {}
        threshDictSave = {}
        specIDXDict = {}
        specDict = {}
        outFile = ROOT.TFile(outDir + '/Bkg_{}_DS{}.root'.format(dType,dsNum), "RECREATE")
        cutTree = ROOT.TChain("skimTree")
        TotalSpec = ROOT.TH1D('DS{}_{}_Corr'.format(dsNum, dType),'',bins,lower,upper)
        UncorrTotalSpec = ROOT.TH1D('DS{}_{}_UnCorr'.format(dsNum, dType),'',bins,lower,upper)
        for bkgidx in range(nRanges[0], nRanges[1]+1):
            # Get Threshold dictionary
            tD = wl.getDBCalRecord("thresh_ds{}_bkgidx{}".format(dsNum, bkgidx))
            specIDXDict[bkgidx] = {}

            # Get Analysis Threshold Dictionary and Exposure Dictionary here
            for idx, ch in enumerate(chList):
                if ch in tD.keys():
                    # Reset Tree
                    cutTree.Reset()

                    if ch not in athresh[bkgidx].keys():
                        print ("Warning: Analysis threshold doesn't exist for ch{} bkgidx{}".format(ch,bkgidx))
                        continue

                    if ch not in ex.Exposure[dsNum][bkgidx].keys():
                        print ("Warning: Exposure doesn't exist for ch{} bkgidx{}".format(ch,bkgidx))
                        continue

                    if not os.path.exists(bgDir + "/fs_rn-DS%d-%d-ch%d.root"%(dsNum,bkgidx,ch)):
                        print ("Warning: Background data doesn't exist for ch%d bkgidx%d"%(ch,bkgidx))
                        continue

                    # Load Background Data
                    cutTree.Add(bgDir + "/fs_rn-DS%d-%d-ch%d.root"%(dsNum,bkgidx,ch))

                    # Create Trigger Efficiency function
                    mu, sigma = tD[ch][0], tD[ch][1]
                    threshFnc = ROOT.TF1("fEff_%d_%d_%d"%(dsNum, ch,bkgidx),"0.5*(1+TMath::Erf((x-[0])/(TMath::Sqrt(2)*[1])))", 0, 250)
                    threshFnc.SetParameters(mu, abs(sigma))

                    # Create Analysis Threshold + Exposure function
                    dExp, dAThresh = ex.Exposure[dsNum][bkgidx][ch], athresh[bkgidx][ch]
                    expFnc = ROOT.TF1("fExp_%d_%d_%d"%(dsNum,ch,bkgidx), "[0]*(x>[1])",0,250)
                    expFnc.SetParameters(dExp, dAThresh)

                    # Print out info
                    print ("DS{} Ch{} BkgIdx{}  Exposure {}  Thresh {}  Analysis Threshold {} ".format(dsNum, ch, bkgidx, dExp, mu, dAThresh))

                    if ch in specIDXDict[bkgidx].keys():
                        print ("Channel {} already exists in bkgidx{}".format(ch,bkgidx))
                        continue
                    else:
                        specIDXDict[bkgidx][ch] = ROOT.TH1D()

                    specIDXDict[bkgidx][ch] = wl.H1D(cutTree,bins,lower,upper, "trapENFCal", cuts+"&& trapENFCal>%.2f && channel==%d"%(dAThresh, ch), Title="hBkg_Ch%d_Bkgidx%d"%(ch,bkgidx), Name="hBkg_Ch%d_Bkgidx%d"%(ch,bkgidx))

                    # Exposure function for scaling
                    h3 = ROOT.TH1D("hCh%d_Bkgidx%d"%(ch, bkgidx), "", bins,lower,upper)
                    for i in range(h3.GetNbinsX()+1):
                        if i < dAThresh*10: continue # Round up to set analysis threshold
                        h3.SetBinContent(i, dExp)
                    # Scale here
                    h3.Multiply(threshFnc)
                    if ch in threshDict.keys():
                        threshDict[ch].Add(h3)
                    else:
                        threshDict[ch] = h3.Clone("Efficiency_ch{}".format(ch))

                    # Save exposure function to histogram
                    h4 = ROOT.TH1D("hCh%d_Bkgidx%d"%(ch, bkgidx), "", bins2,lower,upper)
                    for i in range(h4.GetNbinsX()+1):
                        if i < dAThresh*1000: continue # Round up to set analysis threshold
                        h4.SetBinContent(i, dExp)
                    # Scale here
                    h4.Multiply(threshFnc)
                    if ch in threshDictSave.keys():
                        threshDictSave[ch].Add(h4)
                    else:
                        threshDictSave[ch] = h4.Clone("Efficiency_ch{}".format(ch))

                    if specIDXDict[bkgidx][ch].Integral() == 0:
                        print ("Ch %d has 0 counts, not saving"%(ch))
                        continue
                    if ch in specDict.keys():
                        specDict[ch].Add(specIDXDict[bkgidx][ch])
                    else:
                        specDict[ch] = specIDXDict[bkgidx][ch].Clone('hBkg_Ch%d'%(ch))
                    specIDXDict[bkgidx][ch].Write()

        EffTot = ROOT.TH1D("DS{}_EffTot_Divide".format(dsNum), "DS{} {}".format(dsNum, dType), bins,lower,upper)
        EffTotSave = ROOT.TH1D("DS{}_{}_EffTot".format(dsNum, dType), "DS{} {}".format(dsNum, dType), bins2,lower,upper)

        for ch in specDict.keys():
            if ch not in threshDict.keys():
                print ("Ch {} not in threshDict... ".format(ch))
            TotalSpec.Add(specDict[ch])
            UncorrTotalSpec.Add(specDict[ch])
            EffTot.Add(threshDict[ch])
            EffTotSave.Add(threshDictSave[ch])
            h1 = specDict[ch].Clone('hDS{}_ch{}'.format(dsNum,ch))
            h1.Divide(threshDict[ch])
            h1.Write()

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


def GenerateCorrectedFiles(dsNum = 1):
    print ('Generating Corrected Files for DS{}'.format(dsNum))
    bgDir, outDir = '/Users/brianzhu/project/cuts/fs_rn','/Users/brianzhu/project/cuts/corrfs_rn'
    chList = ds.GetGoodChanList(dsNum)
    if dsNum==5: # remove 692 and 1232 (both beges, so who cares)
        chList.remove(692)
        chList.remove(1232)

    nRanges = [0, ds.dsMap[dsNum]]
    if dsNum==5: nRanges[0] = 80 # exclude DS-5A
    athresh = GetAnaylsisThreshold(dsNum)

    cutTree = ROOT.TChain("skimTree")
    for bkgidx in range(nRanges[0], nRanges[1]+1):
        for ch in chList:
            cutTree.Reset()
            if ch not in athresh[bkgidx].keys():
                print ("Warning: Analysis threshold doesn't exist for ch{} bkgidx{}".format(ch,bkgidx))
                continue
            dAThresh = athresh[bkgidx][ch]
            if not os.path.exists(bgDir+'/fs_rn-DS{}-{}-ch{}.root'.format(dsNum, bkgidx, ch)): continue
            cutTree.Add(bgDir+'/fs_rn-DS{}-{}-ch{}.root'.format(dsNum, bkgidx, ch))

            outFile = ROOT.TFile(outDir+'/corrfs_rn-DS{}-{}-ch{}.root'.format(dsNum, bkgidx, ch),'RECREATE')
            outTree = ROOT.TTree()
            if cutTree.GetEntries() == 0:
                print ("fs_rn-DS{}-{}-ch{}.root has no entries, skipping".format(dsNum, bkgidx, ch))
            outTree = cutTree.CopyTree("trapENFCal>{0:.3f} && channel=={chan}".format(dAThresh, chan=ch))
            print ("Copying corrfs_rn-DS{}-{}-ch{}.root: Analysis Thresh: {}".format(dsNum, bkgidx, ch, dAThresh))

            outTree.Write()
            cutUsed = ROOT.TNamed('AThreshCut',"trapENFCal>{0:.3f} && channel=={chan}".format(dAThresh, chan=ch))
            cutUsed.Write()
            outFile.Close()

if __name__ == "__main__":
    # athresh = GetAnaylsisThreshold(1)
    typeList = ['isEnr', 'isNat']
    for dType in typeList:
        GenerateCorrectedSpectra(5, dType)
    # GenerateCorrectedFiles(5)
