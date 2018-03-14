#!/usr/bin/env python
import sys, time, os, pywt
from ROOT import TFile, TTree, TEntryList, gDirectory, TNamed, std, TObject, MGTWaveform
import numpy as np
import pandas as pd
import waveLibs as wl
import warnings

"""
	This script converts a ROOT TTree with a skimTree into a Pandas dataframe (HDF5 format)
	Load the dataframe with: df = pd.read_hdf(fileName, 'skimTree')
	Events become unpacked, as in multi-hit events are just separated as single events
	Using the Draw command for building event lists makes the processing significantly slower

	Usage: ./ROOTtoPandas.py -f [inDir] [inFile] [outDir]

"""
def main(argv):

	startT = time.clock()
	inDir, outDir = ".", "."
	inFile, inFileName, outFileName = TFile(), "", ""
	gatTree = TTree()
	dsNum, subNum, runNum, theCut, inFileName = -1, -1, -1, "", ""
	saveLAT, savePacket, saveWave = False, False, False

	if len(argv) == 0: return
	for i,opt in enumerate(argv):
		if opt == "-f":
			# dsNum, subNum, inDir, outDir = int(argv[i+1]), int(argv[i+2]), str(argv[i+3]), str(argv[i+4])
			inDir, inFileName, outDir = str(argv[i+1]), str(argv[i+2]), str(argv[i+3])
			# inFileName = "waveSkimDS{}_{}.root".format(dsNum, subNum)
		if opt == "-r":
			inDir, inFileName, outDir = str(argv[i+1]), str(argv[i+2]), str(argv[i+3])
			# dsNum, runNum, inDir, outDir = int(argv[i+1]), int(argv[i+2]), str(argv[i+3]), str(argv[i+4])
			# inFileName = "waveSkimDS{}_run{}.root".format(dsNum, runNum)

	# inFileName = inDir + inFile
	print ("Scanning File: {}".format(inFileName))
	inFile = TFile("%s/%s" % (inDir, inFileName))
	gatTree = inFile.Get("skimTree")
	theCut = inFile.Get("theCut").GetTitle()

	# Make files smaller for tests
	# theCut += " && sumEHL > 236 && sumEHL < 240 && mHL==2 && trapENFCal < 5"
	# Select only pulsers
	# theCut += " && EventDC1Bits > 0"

	print "Using cut:\n",theCut
	gatTree.Draw(">>elist", theCut, "entrylist")
	elist = gDirectory.Get("elist")
	gatTree.SetEntryList(elist)
	nList = elist.GetN()
	print "Found",gatTree.GetEntries(),"input entries."
	print "Found",nList,"entries passing cuts."
	# Gimmicky but works... this bypasses creating the branches...
	gatTree.GetEntry(0)

	# Mess of various branches
	channel = std.vector("int")()
	trapENFCal = std.vector("double")()
	trapENM = std.vector("double")()

	# Create map of branches to put into dataframe
	# This map is only for branches that we want to keep!
	keepMapBase = {
		'trapENFCal':gatTree.trapENFCal, 'trapENM':gatTree.trapENM, "channel":gatTree.channel, "run":gatTree.run, "mHL":gatTree.mHL
		}

	# Combine dictionaries, if keepMapLAT is empty it won't add any branches
	keepMap = dict(keepMapBase)
	# keepMap.update(keepMapLAT)

	dataList = []
	print 'Writing to: ', '%s/proc%s.h5' % (outDir,inFileName.split('.')[0])

	iList, removeNBeg, removeNEnd = -1, 500, 500
	# Loop over events
	while True:
		iList += 1
		if iList >= nList: break
		# if iList >= 5000: break
		entry = gatTree.GetEntryNumber(iList)
		gatTree.LoadTree(entry)
		gatTree.GetEntry(entry)
		nChans = gatTree.channel.size()
		numPass = gatTree.Draw("channel",theCut,"GOFF",1,iList)
		chans = gatTree.GetV1()
		chanList = list(set(int(chans[n]) for n in xrange(numPass)))
		hitList = (iH for iH in xrange(nChans) if gatTree.channel.at(iH) in chanList)  # a 'generator expression'
		for iH in hitList:
			dataMap = {}
			wf = gatTree.MGTWaveforms.at(iH)
			signal = wl.processWaveform(wf,removeNBeg,removeNEnd)
			wave = np.array(signal.GetWaveRaw(), dtype=np.int16)
			for key, branch in keepMap.items():
				# Save branches that aren't vector<Template> (so far only run and mHL)
				if key == 'run' or key == 'mHL': dataMap[key] = int(branch)
				elif key == 'channel': dataMap[key] = int(branch.at(iH))
				else: dataMap[key] = float(branch.at(iH))
			dataMap['waveform'] = wave.tolist()
			dataList.append(dataMap)

		if iList%5000 == 0 and iList!=0:
			print "%d / %d entries saved (%.2f %% done), time: %s" % (iList,nList,100*(float(iList)/nList),time.strftime('%X %x %Z'))


	df = pd.DataFrame.from_dict(dataList)
	print(df.head())
	print(df.info())
	print(np.unique(df.channel))
	print(np.unique(df.mHL))
	print(np.unique(df.run))
	# Suppress stupid warning
	warnings.filterwarnings(action="ignore", module="pandas", message="^\nyour performance")

	# Chunk write like a sucker
	chunksize = 50000
	start = 0
	end = chunksize-1
	i = 0
	# for i in len(df):
	while end < df.shape[0]:
		# chunk = df.iloc[(i*chunksize):min((i+1)*chunksize,len(df))]
		chunk = df.iloc[start:end]
		try:
			chunk.to_hdf('{}/proc{}_{}.h5'.format(outDir,inFileName.split('.')[0],i), key='skimTree', data_columns=['trapENFCal', 'trapENM','channel','mHL','waveform'], format='fixed', mode='w', complevel=9)
		except (Exception) as e:
			print e
			print chunk
			print chunk.info()
		start += chunksize
		end += chunksize
		i += 1

	# df.to_hdf('%s/proc%s.h5' % (outDir,inFileName.split('.')[0]), key="skimTree", data_columns=['trapENFCal', 'trapENM','channel','mHL','waveform'], format='fixed', mode='w', complevel=9)

	stopT = time.clock()
	print("Stopped:",time.strftime('%X %x %Z'),"\nProcess time (min):",(stopT - startT)/60)
	print(float(nList)/((stopT-startT)/60.),"entries per minute.")

if __name__ == "__main__":
	main(sys.argv[1:])
