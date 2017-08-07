#!/usr/local/bin/python
#!/usr/common/usg/software/python/2.7.9/bin/python
import sys, imp
import numpy as np
import matplotlib.pyplot as plt
from ROOT import TFile,TTree,MGTWaveform
wl = imp.load_source('waveLibs', '../waveLibs.py')

def main(argv):
    """ Just compare two files, one before and one after the NL correction,
    to verify that it was actually applied.
    """
    path1 = "waveSkimDS5_run21975_noNLC.root"
    # path2 = "waveSkimDS5_run21975" # pass 1 NLC
    path2 = "waveSkimDS5_run21975_NLC2.root" # pass 2 NLC


    scanSpeed = 0.2
    opt1, opt2 = "", ""
    intMode, printWF, warpMode = False, False, False
    if (len(argv) >= 1): opt1 = argv[0]
    if (len(argv) >= 2): opt2 = argv[1]
    if "-i" in (opt1, opt2):
        intMode = True
        print "Interactive mode selected."

    f1 = TFile("../waveSkimDS5_run21975_noNLC.root")
    t1 = f1.Get("skimTree")

    f2 = TFile("../waveSkimDS5_run21975.root")
    t2 = f2.Get("skimTree")

    fig = plt.figure(figsize=(9,6),facecolor='w')
    plt.show(block=False)

    # Loop over events
    iList = -1
    while(True):
        iList += 1
        if intMode==True and iList !=0:
            value = raw_input()
            if value=='q': break
            if value=='p': iList -= 2  # previous
            if (value.isdigit()): iList = int(value) # go to entry
        if iList >= t1.GetEntries(): break

        t1.GetEntry(iList)
        n1 = t1.channel.size()

        t2.GetEntry(iList)
        n2 = t2.channel.size()

        if n1!=n2:
            print "iList",iList,"entries don't match"
            return

        for iH in range(n1):

            chan1 = t1.channel.at(iH)
            wf1 = t1.MGTWaveforms.at(iH)
            s1 = wl.processWaveform(wf1)
            data1 = s1.GetWaveRaw()
            dataTS1 = s1.GetTS()

            chan2 = t2.channel.at(iH)
            wf2 = t2.MGTWaveforms.at(iH)
            s2 = wl.processWaveform(wf2)
            data2 = s2.GetWaveRaw()
            dataTS2 = s2.GetTS()

            print "chan1 %d  chan2 %d" % (chan1, chan2)

            print "arrays equal? ",np.array_equal(data1[3:],data2[3:])

            # diffs = data1 - data2
            # idx = np.where(abs(diffs)>0)
            # np.set_printoptions(threshold='nan')
            # print diffs

            plt.cla()
            plt.plot(dataTS1,data1,color='blue',label='no NLC')
            plt.plot(dataTS2,data2,color='red',label='with NLC')
            plt.legend(loc='best')
            plt.tight_layout()

            plt.pause(scanSpeed)


if __name__ == "__main__":
    main(sys.argv[1:])