#!/usr/local/bin/python
#!/usr/common/usg/software/python/2.7.9/bin/python
import sys, imp
import numpy as np
import matplotlib.pyplot as plt
from ROOT import TFile,TTree,TEntryList,MGTWaveform,gDirectory
from ROOT import MGWFNonLinearityCorrector,
# from ROOT import GATNonLinearityCorrector
wl = imp.load_source('waveLibs', '../waveLibs.py')

def main(argv):
    """ GATNonLinearityCorrector is not accessible from pyroot.
    It makes more sense to correct for NL in wave-skim anyway.
    """

    scanSpeed = 0.2
    opt1, opt2 = "", ""
    intMode, batMode = False, False
    if (len(argv) >= 1): opt1 = argv[0]
    if (len(argv) >= 2): opt2 = argv[1]
    if "-i" in (opt1, opt2):
        intMode = True
        print "Interactive mode selected."
    if "-b" in (opt1, opt2):
        batMode = True
        print "Batch mode selected."

    inputFile = TFile("../waveSkimDS5_run21975.root")
    waveTree = inputFile.Get("skimTree")
    print "Found",waveTree.GetEntries(),"input entries."

    theCut = inputFile.Get("theCut").GetTitle()
    # theCut += " && Entry$ < 100"
    theCut += " && trapENFCal < 6"

    waveTree.Draw(">>elist", theCut, "entrylist")
    elist = gDirectory.Get("elist")
    waveTree.SetEntryList(elist)
    nList = elist.GetN()

    print "Using cut:\n",theCut,"\nFound",nList,"entries passing cuts."

    fig = plt.figure(figsize=(13,7), facecolor='w')
    p0 = plt.subplot(111)

    if not batMode: plt.show(block=False)

    # Loop over events
    iList = -1
    while(True):
        iList += 1
        if intMode==True and iList !=0:
            value = raw_input()
            if value=='q': break
            if value=='p': iList -= 2  # previous
            if (value.isdigit()): iList = int(value) # go to entry
        if iList >= elist.GetN(): break

        entry = waveTree.GetEntryNumber(iList);
        waveTree.LoadTree(entry)
        waveTree.GetEntry(entry)
        nChans = waveTree.channel.size()
        nWFs = waveTree.MGTWaveforms.size()
        if (nWFs==0):
            print "Error - nWFs:",nWFs,"nChans",nChans
            continue
        numPass = waveTree.Draw("channel",theCut,"GOFF",1,iList)
        chans = waveTree.GetV1()
        chanList = list(set(int(chans[n]) for n in xrange(numPass)))

        # Loop over hits passing cuts
        hitList = (iH for iH in xrange(nChans) if waveTree.channel.at(iH) in chanList)  # a 'generator expression'
        for iH in hitList:
            run = waveTree.run
            chan = waveTree.channel.at(iH)
            trapENFCal = waveTree.trapENFCal.at(iH)
            wf = waveTree.MGTWaveforms.at(iH)
            signal = wl.processWaveform(wf)
            data = signal.GetWaveBLSub()
            dataTS = signal.GetTS()
            _,dataNoise = signal.GetBaseNoise()
            dataTSMax = waveTree.trapENMSample.at(iH)*10. - 4000
            dataENM = waveTree.trapENM.at(iH)




            # plots
            p0.cla()
            p0.plot(dataTS,data,color='blue')

            plt.tight_layout()
            plt.pause(scanSpeed)


if __name__ == "__main__":
    main(sys.argv[1:])