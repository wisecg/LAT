#!/usr/bin/env python
from ROOT import TFile, TTree, TCanvas, TH1D, MGTWaveform
from ROOT import gDirectory
import numpy as np
import matplotlib.pyplot as plt

def main():

    intMode=True

    f = TFile("./data/latSkimDS1_run10247_0.root","READ")
    tree = f.Get("skimTree")
    nEnt = tree.GetEntries()

    theCut = f.Get("theCut").GetTitle()
    # theCut += " && trapENFCal  1.1"

    print "Using cut:\n",theCut,"\n"
    tree.Draw(">>elist", theCut, "entrylist")
    elist = gDirectory.Get("elist")
    tree.SetEntryList(elist)
    nList = elist.GetN()
    print "Found",nList,"entries passing cuts."

    fig = plt.figure(figsize=(10,7), facecolor='w')
    p1 = plt.subplot(211)
    p2 = plt.subplot(212)
    plt.show(block=False)

    iList = -1
    while(True):
        iList += 1
        if intMode==True and iList !=0:
            value = raw_input()
            if value=='q': break
            if value=='p': iList -= 2  # previous
            if (value.isdigit()): iList = int(value) # go to entry
        if iList >= elist.GetN(): break

        entry = tree.GetEntryNumber(iList);
        tree.LoadTree(entry)
        tree.GetEntry(entry)
        nChans = tree.channel.size()
        nWFs = tree.MGTWaveforms.size()

        if (nWFs==0):
            print "Error - nWFs:",nWFs,"nChans",nChans
            continue

        numPass = tree.Draw("channel",theCut,"GOFF",1,iList)
        chans = tree.GetV1()
        chanList = list(set(int(chans[n]) for n in xrange(numPass)))

        # Loop over hits passing cuts
        hitList = (iH for iH in xrange(nChans) if tree.channel.at(iH) in chanList)  # a 'generator expression'
        for iH in hitList:
            run = tree.run
            chan = tree.channel.at(iH)
            energy = tree.trapENFCal.at(iH)
            wf = tree.MGTWaveforms.at(iH)
            signal = wl.processWaveform(wf,0,0)
            waveRaw = signal.GetWaveRaw()
            waveTS = signal.GetTS()
            print "%d / %d  Run %d  nCh %d  chan %d  trapENF %.1f" % (iList,nList,run,nChans,chan,energy)

            p1.plot(waveTS,waveRaw)
            plt.pause(scanSpeed)






if __name__=="__main__":
    main()