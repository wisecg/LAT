#!/usr/common/usg/software/python/2.7.9/bin/python
import sys, time, os, pywt
from ROOT import TFile, TTree, TEntryList, gDirectory, TNamed, std, TObject, MGTWaveform
import numpy as np
import pandas as pd
import waveLibs as wl

"""
	This script converts a ROOT TTree with a skimTree into a Pandas dataframe (HDF5 format)
	Load the dataframe with: df = pd.read_hdf(fileName, 'skimTree')
	Events become unpacked, as in multi-hit events are just separated as single events
	Using the Draw command for building event lists makes the processing significantly slower

	Usage: ./ROOTtoPandas.py -ws [DS #] [SubDS #]
	Usage: ./ROOTtoPandas.py -f [DS #] [Run #]
	Option "-l": Saves LAT parameters
	Option "-p": Saves wavelet packet parameters
	Option "-w": Saves waveforms

"""
def main(argv):

	startT = time.clock()
	inDir, outDir = ".", "."
	inFile, outFile = TFile(), TFile()
	gatTree, outTree = TTree(), TTree()
	dsNum, runNum, theCut, inFileName = -1, -1, "", ""
	saveLAT, savePacket, saveWave = False, False, False

	if len(argv) == 0: return

	for i,opt in enumerate(argv):
		if opt == "-f":
			dsNum, runNum = int(argv[i+1]), int(argv[i+2])
			inFileName = 'waveSkimDS%d_run%d.root' % (dsNum, runNum)
			print "Scanning DS-%d, run %d" % (dsNum, runNum)
		if opt == "-ws":
			dsNum, subNum = int(argv[i+1]), int(argv[i+2])
			print "Scanning DS-%d sub-range %d" % (dsNum, subNum)
			inFileName = 'waveSkimDS%d_%d.root' % (dsNum, subNum)
		if opt == "-ls":
			dsNum, subNum, subsubNum = int(argv[i+1]), int(argv[i+2], int(argv[i+3]))
			print "Scanning DS-%d sub-range %d (sub %d)" % (dsNum, subNum, subsubNum)
			inFileName = 'latSkimDS%d_%d_%d.root' % (dsNum, subNum, subsubNum)
		if opt == "-p":
			savePacket = True
			print "Saving wavelet packet"
		if opt == "-w":
			saveWave = True
			print "Saving waveform"
		if opt == "-l":
			saveLAT = True
			print "Saving LAT parameters"
		if opt == "-d":
			inDir, outDir = argv[i+1], argv[i+2]
			print "Custom paths: Input %s,  Output %s" % (inDir,outDir)

	inFile = TFile("%s/%s" % (inDir, inFileName))
	gatTree = inFile.Get("skimTree")
	theCut = inFile.Get("theCut").GetTitle()

	# Make files smaller for tests
	theCut += "&& trapETailMin < 0.5 && trapENFCal > 100.0"

	print "Using cut:\n",theCut
	gatTree.Draw(">>elist", theCut, "entrylist")
	elist = gDirectory.Get("elist")
	gatTree.SetEntryList(elist)
	nList = elist.GetN()
	print "Found",gatTree.GetEntries(),"input entries."
	print "Found",nList,"entries passing cuts."

	# Test purposes
	# nList = 50000
	# Divide up run for chunk-writing
	nDivis, nChunk = nList//2000, 2000

	# Gimmicky but works... this bypasses creating the branches...
	gatTree.GetEntry(0)

	# Mess of various branches
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
	avse, dcr99 = std.vector("double")(), std.vector("double")()

	# Create map of branches to put into dataframe
	# This map is only for branches that we want to keep!
	keepMapBase = {
		'trapENFCal':gatTree.trapENFCal, 'triggerTrapt0':gatTree.triggerTrapt0, "channel":gatTree.channel, "run":gatTree.run, "mHL":gatTree.mHL, "C":gatTree.C, "avse":gatTree.avse, "dcr99":gatTree.dcr99
		}

	keepMapLAT = {}
	if saveLAT:
		keepMapLAT = { "waveS1":gatTree.waveS1, "waveS2":gatTree.waveS2, "waveS3":gatTree.waveS3, "waveS4":gatTree.waveS4, "waveS5":gatTree.waveS5, "bcMax":gatTree.bcMax, "bcMin":gatTree.bcMin, "tOffset":gatTree.tOffset, "fitMatch":gatTree.fitMatch, "fitE":gatTree.fitE, "fitSlo":gatTree.fitSlo, "matchMax":gatTree.matchMax, "matchWidth":gatTree.matchWidth, "matchTime":gatTree.matchTime, "fitMax":gatTree.fitMax, "pol0":gatTree.pol0, "pol1":gatTree.pol1, "pol2":gatTree.pol2, "kvorrT":gatTree.kvorrT}

	# Combine dictionaries, if keepMapLAT is empty it won't add any branches
	keepMap = dict(keepMapBase)
	keepMap.update(keepMapLAT)

	dfList = []
	print 'Writing to: ', '%s/%s.h5' % (outDir,inFileName.split('.')[0])
	store = pd.HDFStore('%s/%s.h5' % (outDir,inFileName.split('.')[0]), 'w')
	iList, iChunk = -1, -1
	# Loop through number of chunks
	for chunk in xrange(0,nDivis):
		# Select size to save, depending on remaining events
		chunkSize = np.amin([nList-chunk*nChunk, nChunk])
		# Create empty dictionary of branch names and arrays of values
		branchMap = {}
		removeNBeg, removeNEnd = 0, 2
		for branch in gatTree.GetListOfBranches():
			if branch.GetName() in keepMap: branchMap[branch.GetName()] = np.zeros(chunkSize)
		if saveWave:
			wf = gatTree.MGTWaveforms.at(0)
			signal = wl.processWaveform(wf,removeNBeg,removeNEnd)
			data = signal.GetWaveBLSub()
			# Save individual samples as column
			for idx,sample in enumerate(data):
				if idx < 400 or idx > 1600: continue
				branchMap["wave%d"%(idx)] = np.zeros(chunkSize)
		# Save wavelet packet values as columns
		if savePacket:
			for xsample in range(0, 16):
				for ysample in range(0, 128):
					branchMap["wp%d_%d"%(xsample, ysample)] = np.zeros(chunkSize)

    	# Loop over events
		print "Looping chunk ", chunk
		while True:
			iList += 1
			iChunk += 1
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
				wp = []
				wave = []
				if savePacket or saveWave:
					wf = gatTree.MGTWaveforms.at(iH)
					signal = wl.processWaveform(wf,removeNBeg,removeNEnd)
					wave = signal.GetWaveBLSub()
				for key, branch in keepMap.items():
					# Save branches that aren't vector<Template> (so far only run and mHL)
					if key == 'run' or key == 'mHL': branchMap[key][iChunk] = branch
					else: branchMap[key][iChunk] = branch.at(iH)
				if saveWave:
					# branchMap['wave'][iChunk] = wave
					for idx, val in enumerate(wave):
						# Skip some samples because they're useless
						if idx < 400 or idx > 1600: continue
						branchMap['wave%d'%(idx)][iChunk] = val
				if savePacket:
					packet = pywt.WaveletPacket(wave, 'db2', 'symmetric', maxlevel=4)
					nodes = packet.get_level(4, order='freq')
					wp = np.array([n.data for n in nodes],'d')
					wp = abs(wp)
					for (x, y), val in np.ndenumerate(wp):
						branchMap['wp%d_%d'%(x, y)][iChunk] = val
			if iList == (chunk+1)*nChunk-1:
				iChunk = -1
				break
		if iList%5000 == 0 and iList!=0:
			print "%d / %d entries saved (%.2f %% done), time: %s" % (iList,nList,100*(float(iList)/nList),time.strftime('%X %x %Z'))
		# Convert dictionary to Pandas DataFrame -- make sure to set the index correctly
		df = pd.DataFrame(branchMap)
		df.index = pd.Series(df.index) + chunk*chunkSize
		store.append('skimTree', df)

	store.close()
	stopT = time.clock()
	print "Stopped:",time.strftime('%X %x %Z'),"\nProcess time (min):",(stopT - startT)/60
	print float(nList)/((stopT-startT)/60.),"entries per minute."

if __name__ == "__main__":
	main(sys.argv[1:])
