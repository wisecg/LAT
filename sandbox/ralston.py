#!/usr/local/bin/python

import sys
import numpy as np
import scipy.linalg as lin
import matplotlib.pyplot as plt
from ROOT import TFile, TTree, TEntryList, gDirectory
from ROOT import GATDataSet, MGTEvent, MGTWaveform
import waveLibs as wl

def main(argv):

    intMode, batMode = True, False
    for i,opt in enumerate(argv):
        if opt == "-b":
            intMode, batMode = False, True

    # set D and R
    nBL = 512
    D = 100
    R = 80
    # lo = 80

    gFile = TFile("~/project/v2-waveskim/waveSkimDS5_90.root")
    gatTree = gFile.Get("skimTree")

    # Select a NOISE POPULATION
    theCut = gFile.Get("theCut").GetTitle()
    theCut += " && trapENFCal < 2 && channel!=692 && channel!=1232"
    theCut += " && Entry$ < 10000"

    gatTree.Draw(">>elist", theCut, "entrylist")
    elist = gDirectory.Get("elist")
    gatTree.SetEntryList(elist)
    nList = elist.GetN()
    print "Using cut:\n",theCut
    print "Found",gatTree.GetEntries(),"input entries."
    print "Found",nList,"entries passing cuts."

    nWF = 0
    rhoNoise = np.zeros([D,D])
    for iList in range(nList):
        entry = gatTree.GetEntryNumber(iList);
        gatTree.LoadTree(entry)
        gatTree.GetEntry(entry)
        nChans = gatTree.channel.size()
        numPass = gatTree.Draw("channel",theCut,"GOFF",1,iList)
        chans = gatTree.GetV1()
        chanList = list(set(int(chans[n]) for n in xrange(numPass)))
        hitList = (iH for iH in xrange(nChans) if gatTree.channel.at(iH) in chanList)
        for iH in hitList:
            run = gatTree.run
            iEvent = gatTree.iEvent
            chan = gatTree.channel.at(iH)
            dataENF = gatTree.trapENFCal.at(iH)
            wf = gatTree.MGTWaveforms.at(iH)
            # print "%d:  run %d  chan %d  trapENFCal %.2f" % (iList, run, chan, dataENF)

            signal = wl.processWaveform(wf)
            data = signal.GetWaveBLSub()
            data = data[:nBL]
            binned, nBin = BinData(data,D)
            for vec in binned:
                out = np.outer(vec,vec)
                rhoNoise = np.add(rhoNoise, out)
            nWF += nBin

    # fig = plt.figure(figsize=(9,7), facecolor='w')
    # plt.imshow(rhoNoise,cmap='jet')
    # plt.colorbar()
    # plt.tight_layout()
    # plt.show()


    # ==========================================================
    # ============ Construct the KLJ/Ralston filter ============

    rhoNoise /= nWF
    eigenvals, eigenvectors = lin.eig(rhoNoise)  # eigenvectors are returned normalized.

    eigenDict = {}
    for i in range(len(eigenvals)):
        eigenDict[eigenvals[i]] = eigenvectors[i]

    # Sort the states by the np.absolute value (hopefully that's same as sorting by noise power)
    eigenSorted = sorted(eigenvals, key=lambda eigenval: np.absolute(eigenval), reverse=True)
    # for eig in eigenSorted: print np.absolute(eig), eig
    # print len(eigenSorted)

    # Now construct filters for different ranks R:
    filterDict = {}

    for i in range(len(eigenSorted)):
    # for i in range(R+1): # don't calculate too many levels
        evec = eigenDict[eigenSorted[i]]
        pia = np.outer(evec,evec)
        filterDict[i] = pia

    piNoise = filterDict[R]
    for i in range(R+1,D): piNoise += filterDict[i]

    piSignal = np.identity(D) - piNoise
    # piSignal = piNoise # crazy test

    # fig = plt.figure(figsize=(9,7), facecolor='w')
    # plt.imshow(piSignal,cmap='jet')
    # plt.colorbar()
    # plt.tight_layout()
    # plt.show()


    # ==========================================================


    # Now that we've got the filter matrix, loop back over the data.
    gFile.Close() # it's fucking stupid I have to close and reopen.

    gFile = TFile("~/project/v2-waveskim/waveSkimDS5_90.root")
    gatTree = gFile.Get("skimTree")

    # Select a SIGNAL POPULATION
    theCut = gFile.Get("theCut").GetTitle()
    theCut += " && trapENFCal > 2 && trapENFCal < 10"

    gatTree.Draw(">>elist", theCut, "entrylist")
    elist = gDirectory.Get("elist")
    gatTree.SetEntryList(elist)
    nList = elist.GetN()
    print "Second pass using cut:\n",theCut
    print "Found",gatTree.GetEntries(),"input entries."
    print "Found",nList,"entries passing cuts."

    fig = plt.figure(figsize=(8,6),facecolor='w')
    p1 = plt.subplot(111)
    # p2 = plt.subplot(212)
    plt.show(block=False)

    # Begin loop over events
    iList = -1
    while True:
        iList += 1
        if intMode and iList != 0:
            value = raw_input()
            if value=='q': break        # quit
            if value=='p': iList -= 2   # go to previous
            if (value.isdigit()):
                iList = int(value)      # go to entry number
        elif not intMode and batMode:
            plt.pause(0.001)            # rapid-draw mode
        if iList >= nList: break        # bail out, goose!

        entry = gatTree.GetEntryNumber(iList);
        gatTree.LoadTree(entry)
        gatTree.GetEntry(entry)
        nChans = gatTree.channel.size()

        numPass = gatTree.Draw("channel",theCut,"GOFF",1,iList)
        chans = gatTree.GetV1()
        chanList = list(set(int(chans[n]) for n in xrange(numPass)))
        hitList = (iH for iH in xrange(nChans) if gatTree.channel.at(iH) in chanList)  # a 'generator expression'
        for iH in hitList:

            run = gatTree.run
            iEvent = gatTree.iEvent
            chan = gatTree.channel.at(iH)
            dataENF = gatTree.trapENFCal.at(iH)
            wf = gatTree.MGTWaveforms.at(iH)
            print "%d:  run %d  chan %d  trapENFCal %.2f" % (iList, run, chan, dataENF)

            signal = wl.processWaveform(wf)
            data = signal.GetWaveBLSub()
            dataTS = signal.GetTS()

            # partition the data into 'bins' - filter is applied to each, individually
            # filter matrix is D x D.
            # inside bins, timing resolution of order D * (sampling frequency) is expected.
            # apply the filter to the subspaces
            # transform back

            binned, nBin = BinData(data,D)
            filtBins = []
            for vec in binned:
                filtBins.append( piSignal.dot(vec) ) # dot product

            filtData = np.zeros(len(data))
            for i in range(len(filtData)):
                j, k = int(i/D), i%D
                filtData[i] = filtBins[j][k]


            # show the data
            p1.cla()
            p1.plot(dataTS,filtData,color='red',label='filt')
            p1.plot(dataTS,data,color='blue',alpha=0.7,label='data')
            p1.legend(loc=4)


            plt.tight_layout()
            plt.pause(0.0000001)


def BinData(data,D):
    # 'Bin' (divide + relabel) the first 512 samples into 512/D segments.
    # If 512/D is not an integer, that's OK.  The last vector will just be padded with some zeros.
    binned = []
    nWF, pj = 0, 0
    noiseVec = np.zeros(D)
    for i in range(len(data)):
        j, k = int(i/D), i%D
        if j!=pj:
            # print "nWF",nWF,"noiseVec",noiseVec
            binned.append(noiseVec)
            noiseVec = np.zeros(D)
            nWF += 1
        noiseVec[k] = data[i]
        # print "j %d k %d noiseVec %d" % (j, k, noiseVec[k])
        pj = j # save for next iteration
    # print "nWF",nWF,"noiseVec",noiseVec
    binned.append(noiseVec)
    nWF += 1
    return binned, nWF


if __name__ == "__main__":
    main(sys.argv[1:])