#!/usr/common/usg/software/python/2.7.9/bin/python
import sys, time, os

from ROOT import TFile, TTree, TEntryList, gDirectory, TNamed, std, TObject
import numpy as np
import pandas as pd

"""
	This script converts a ROOT TTree with a skimTree into a Pandas dataframe (HDF5 format)
	Load the dataframe with: df = pd.read_hdf(fileName, 'skimTree')
	Events become unpacked, as in multi-hit events are just separated as single events
	Using the Draw command for building event lists makes the processing significantly slower

	Usage: ./PandifySkim.py -s [DS #] [SubDS #]
	Usage: ./PandifySkim.py -f [DS #] [Run #]

"""
def main(argv):

	inDir, outDir = ".", "."
	inFile, outFile = TFile(), TFile()
	gatTree, outTree = TTree(), TTree()
	dsNum, runNum, theCut, inFileName = -1, -1, "", ""
	if len(argv) == 0: return

	for i,opt in enumerate(argv):
		if opt == "-f":
			dsNum, runNum = int(argv[i+1]), int(argv[i+2])
			inFileName = 'latSkimDS%d_run%d.root' % (dsNum, runNum)
			print "Scanning DS-%d, run %d" % (dsNum, runNum)
			print "Scanning file %s" % inFileName
			inFile = TFile("%s/%s" % (inDir, inFileName))
			gatTree = inFile.Get("skimTree")
			theCut = inFile.Get("theCut").GetTitle()
		if opt == "-s":
			dsNum, subNum, subsubNum = int(argv[i+1]), int(argv[i+2]), int(argv[i+3])
			print "Scanning DS-%d sub-range %d, subsub-range %d" % (dsNum, subNum, subsubNum)
			inFileName = 'latSkimDS%d_%d_%d.root' % (dsNum, subNum, subsubNum)
			inFile = TFile("%s/%s" % (inDir, inFileName))
			gatTree = inFile.Get("skimTree")
			theCut = inFile.Get("theCut").GetTitle()
		if opt == "-d":
			inDir, outDir = argv[i+1], argv[i+2]
			print "Custom paths: Input %s,  Output %s" % (inDir,outDir)

	print "Using cut:\n",theCut
	gatTree.Draw(">>elist", theCut, "entrylist")
	elist = gDirectory.Get("elist")
	gatTree.SetEntryList(elist)
	nList = elist.GetN()
	print "Found",gatTree.GetEntries(),"input entries."
	print "Found",nList,"entries passing cuts."

	# Gimmicky but works... this bypasses creating the branches...
	gatTree.GetEntry(0)

	channel, fails, C = std.vector("int")(), std.vector("int")(), std.vector("int")()
	threshKeV, threshSigma = std.vector("double")(), std.vector("double")()
	kvorrT = std.vector("double")()
	trapENFCal, trapENMCal, trapENMSample = std.vector("double")(), std.vector("double")(), std.vector("int")()
	trapENF, trapENM, triggerTrapt0 = std.vector("double")(), std.vector("double")(), std.vector("double")()
	waveS1, waveS2, waveS3, waveS4, waveS5 = std.vector("double")(), std.vector("double")(), std.vector("double")(), std.vector("double")(), std.vector("double")()
	bcMax, bcMin, wpar4 = std.vector("double")(), std.vector("double")(), std.vector("double")()
	tOffset, derivMax, derivTime = std.vector("double")(), std.vector("double")(), std.vector("double")()
	waveS0, bandMax, bandTime = std.vector("double")(), std.vector("double")(), std.vector("double")()
	den10, den50, den90 = std.vector("double")(), std.vector("double")(), std.vector("double")()
	raw10, raw50, raw90 = std.vector("double")(), std.vector("double")(), std.vector("double")()
	oppie, fitMatch, fitE = std.vector("double")(), std.vector("double")(), std.vector("double")()
	fitSlo, matchMax, matchWidth = std.vector("double")(), std.vector("double")(), std.vector("double")()
	matchTime, fitMax, pol0 = std.vector("double")(), std.vector("double")(), std.vector("double")()
	pol1, pol2, exp0, exp1 = std.vector("double")(), std.vector("double")(), std.vector("double")(), std.vector("double")()

	# Create map of branches to put into dataframe
	# This map is only for branches that we want to keep!
	keepMap = {
		'trapENFCal':gatTree.trapENFCal, 'trapENMCal':gatTree.trapENMCal, 'trapENMSample':gatTree.trapENMSample,
		'trapENF':gatTree.trapENF, 'trapENM':gatTree.trapENM, 'triggerTrapt0':gatTree.triggerTrapt0, "waveS0":gatTree.waveS0,
		"waveS1":gatTree.waveS1, "waveS2":gatTree.waveS2, "waveS3":gatTree.waveS3, "waveS4":gatTree.waveS4,
		"waveS5":gatTree.waveS5, "bcMax":gatTree.bcMax, "bcMin":gatTree.bcMin, "wpar4":gatTree.wpar4,
      	"tOffset":gatTree.tOffset, "derivMax":gatTree.derivMax, "derivTime":gatTree.derivTime,
      	"bandMax":gatTree.bandMax, "bandTime":gatTree.bandTime, "den10":gatTree.den10, "den50":gatTree.den50,
      	"den90":gatTree.den90, "oppie":gatTree.oppie, "fitMatch":gatTree.fitMatch, "fitE":gatTree.fitE,
      	"fitSlo":gatTree.fitSlo, "matchMax":gatTree.matchMax, "matchWidth":gatTree.matchWidth,
      	"matchTime":gatTree.matchTime, "fitMax":gatTree.fitMax, "pol0":gatTree.pol0, "pol1":gatTree.pol1,
      	"pol2":gatTree.pol2, "exp0":gatTree.exp0, "exp1":gatTree.exp1, "channel":gatTree.channel,
		"fails":gatTree.fails, "threshKeV":gatTree.threshKeV, "threshSigma":gatTree.threshSigma, "run":gatTree.run,
		"raw10":gatTree.raw10, "raw50":gatTree.raw50, "raw90":gatTree.raw90,"C":gatTree.C,"kvorrT":gatTree.kvorrT
		}

	# Create empty dictionary of branch names and arrays of values
	branchMap = {}
	for branch in gatTree.GetListOfBranches():
		if branch.GetName() in keepMap: branchMap[branch.GetName()] = np.zeros(nList)

    # Loop over events
	print "Starting event loop ..."
	iList = -1
	while True:
		iList += 1
		if iList >= nList: break
		entry = gatTree.GetEntryNumber(iList)
		gatTree.LoadTree(entry)
		gatTree.GetEntry(entry)
		nChans = gatTree.channel.size()

		numPass = gatTree.Draw("channel",theCut,"GOFF",1,iList)
		chans = gatTree.GetV1()
		chanList = list(set(int(chans[n]) for n in xrange(numPass)))
		hitList = (iH for iH in xrange(nChans) if gatTree.channel.at(iH) in chanList)  # a 'generator expression'
		for iH in hitList:
			for key, branch in keepMap.items():
			# Branches that aren't vector<T>
				if key == 'run' or key == 'mHL': branchMap[key][iList] = branch
				else: branchMap[key][iList] = branch.at(iH)

	# Convert dictionary to Pandas DataFrame
	df = pd.DataFrame(branchMap)
	print df.head()
	print df.shape
	print 'Writing to: ', '%s/%s.h5' % (outDir,inFileName.split('.')[0])
	store = pd.HDFStore('%s/%s.h5' % (outDir,inFileName.split('.')[0]))
	store.put('skimTree', df)
	store.close()

if __name__ == "__main__":
	main(sys.argv[1:])
