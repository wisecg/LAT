#!/usr/bin/env python
import sys, imp, time, glob
sys.argv.append("-b") # kill all interactive crap
ds = imp.load_source('DataSetInfo','../DataSetInfo.py')
wl = imp.load_source('waveLibs','../waveLibs.py')
import ROOT
from ROOT import TFile, TChain, TTree, TCanvas, TH1D, MGTWaveform
from ROOT import gDirectory, gStyle, std
import numpy as np
import matplotlib.pyplot as plt

def main():

    # updateFile()
    makePlots()


def updateFile():

    intMode=False

    files = glob.glob("./data/*.root")
    for fileName in files:

        start = time.clock()
        print "Now scanning",fileName

        f = TFile(fileName,"UPDATE")
        tree = f.Get("skimTree")
        nEnt = tree.GetEntries()

        # check if the wfstd branch already exists, or create a new one
        b0 = tree.GetListOfBranches().FindObject("wfstd")
        if isinstance(b0, ROOT.TBranchElement): tree.GetListOfBranches().Remove(b0)

        wfstd = std.vector("double")()
        b1 = tree.Branch("wfstd",wfstd)

        theCut = f.Get("theCut").GetTitle()
        theCut += " && Entry$ < 10"

        # fig = plt.figure(figsize=(10,7), facecolor='w')
        # p1 = plt.subplot(211)
        # p2 = plt.subplot(212)
        # plt.show(block=False)

        iList = -1
        while(True):
            iList += 1
            if intMode==True and iList !=0:
                value = raw_input()
                if value=='q': break
                if value=='p': iList -= 2  # previous
                if (value.isdigit()): iList = int(value) # go to entry
            if iList >= nEnt: break

            tree.GetEntry(iList)
            nChans = tree.channel.size()
            nWFs = tree.MGTWaveforms.size()
            if (nChans != nWFs):
                print "Wrong num entries.  Bailing!"
                exit(1)

            wfstd.assign(nChans,-88888)

            # loop over hits
            for iH in range(nWFs):
                run = tree.run
                chan = tree.channel.at(iH)
                energy = tree.trapENFCal.at(iH)
                wf = tree.MGTWaveforms.at(iH)
                signal = wl.processWaveform(wf,0,0)
                waveRaw = signal.GetWaveRaw()
                waveTS = signal.GetTS()
                # print "%d / %d  Run %d  nCh %d  chan %d  trapENF %.1f" % (iList,nList,run,nChans,chan,energy)

                wfstd[iH] = np.std(waveRaw[5:-5])
                maxAdc = max(waveRaw[5:-5])
                minAdc = min(waveRaw[5:-5])
                nBins = int(maxAdc-minAdc)

                # if intMode:
                #     p1.cla()
                #     p1.plot(waveTS[5:-5],waveRaw[5:-5])
                #     p2.cla()
                #     h2, edges = np.histogram(waveRaw[5:-5],bins=nBins)
                #     p2.bar(edges[:-1], h2, width=np.diff(edges), ec="k", align="edge")
                #     p1.set_title("wfstd: %.3f" % wfstd[iH])
                #     plt.tight_layout()
                #     plt.pause(0.00000001)
                    # plt.savefig("../plots/ralphWidth.png")

            # End loop over hits, fill branches
            if not intMode:
                b1.Fill()

        # End loop over events
        if not intMode:
            tree.Write("",ROOT.TObject.kOverwrite)
            print "  Tree entries: %d  Branch entries: %d  Time (min): %.2f" % (tree.GetEntries(),b1.GetEntries(),(time.clock()-start)/60.)

        f.Close()


def makePlots():

    gStyle.SetOptStat(0)

    tree = TChain("skimTree")
    tree.Add("./data/latSkimDS0_run4*.root")
    nEnt = tree.GetEntries()

    chan = 692

    # check if the wfstd branch already exists, or create a new one
    b0 = tree.GetListOfBranches().FindObject("wfstd")
    if isinstance(b0, ROOT.TBranchElement):
        print "Found it!"

    c = TCanvas("c","c",800,600)
    c.SetLogz(1)
    c.SetLogy(0)

    eLo, eHi = 0., 14
    kevPerBin = 0.1
    nBins = int((eHi-eLo)/kevPerBin)

    h1 = wl.H2D(tree,nBins,eLo,eHi,nBins,eLo,eHi,"wfstd:trapENFCalC","gain==0 && channel==692","Energy (keV)","wfstd"," ")
    h1.Draw("COLZ")
    c.Print("../plots/wfstd.pdf")

    # c.SetLogy(1)
    # h2 = wl.H1D(tree,nBins,eLo,eHi,"trapENFCalC","","Energy (keV)","cts"," ")
    # h2.Draw("hist")
    # c.Print("../plots/wfstdEnergy.pdf")

    # c.SetLogy(0)
    # h4 = wl.H2D(tree,200,0,1000,200,0.,250.,"wfstd:fitSlo","gain==0 && trapENFCalC < 250","fitSlo","wfstd"," ")
    # h4.SetMinimum(1)
    # h4.Draw("COLZ")
    # c.Print("../plots/wfstdVSfitslo.pdf")


if __name__=="__main__":
    main()