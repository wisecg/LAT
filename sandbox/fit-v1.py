#!/usr/local/bin/python
import numpy as np
import matplotlib.text as text
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
from ROOT import *
from scipy.optimize import curve_fit
from scipy import signal as sg
from scipy import integrate
import scipy.fftpack
import waveModel as wm
import sys
import pymc
from scipy.ndimage.filters import gaussian_filter
from subprocess import call

# 16 High-gain Module 1 detectors active in DS1
# BEGEs: 600 p7d1, 692 p2d1.  All others are enriched.
chanDict = {610:"P3D2", 600:"P7D1", 692:"P2D1", 598:"P7D2", 626:"P6D3", 648:"P2D2", 582:"P1D2", 640:"P2D3", 578:"P1D4", 580:"P1D3", 690:"P6D4", 592:"P7D4", 672:"P5D3", 608:"P3D3", 664:"P3D4", 632:"P6D1"}

def main():
	# take care of all the goddamn root bullshit file i/o in the main function

	waveTemplates = TChain("waveTree")
	waveTemplates.Add("./data/wave-m1cal.root")
	print "Found",waveTemplates.GetEntries(),"template waveforms."
	# PlotTemplateWaveforms(waveTemplates,calDict)

	f1 = TFile("./data/wave-standardCutDS1-newBranch.root")
	waveTree = f1.Get("waveTree")
	eventTree = f1.Get("eventTree")
	print "Found",waveTree.GetEntries(),"input waveforms."
	print "Found",eventTree.GetEntries(),"input entries."

	# control the loop with an entry list
	# theCut = "trapENFCal > 2 && trapENFCal < 5 && channel!=600 && channel!=692"
	# fft + extended run cut
	theCut = "trapENFCal > 5 && trapENFCal < 5.02 && channel!=600 && channel!=692 && d2wf0MHzTo50MHzPower < 4000"
	theCut += " && run != 9648 && run != 9663 && run != 10185 && run != 10663 && run != 10745 && run != 11175 && run != 12445 && run != 12723 && run != 12735 && run != 13004 && channel != 580"

	# get that money and take that shit to the bank
	waveTreeNew, mcmc_params = mcmc_wfs(waveTemplates,waveTree,eventTree,theCut)

	# get the results back into a ROOT file.
	print "Cloning output tree ..."

	call(['rm',"./data/wave-standardCutDS1-mcmc.root"])
	f2 = TFile("./data/wave-standardCutDS1-mcmc.root","new")

	waveTreeNew.Draw(">>elist2", theCut , "entrylist")
	elist2 = gDirectory.Get("elist2")
	waveTreeNew.SetEntryList(elist2)
	waveTreeNewSmall = waveTreeNew.CopyTree("")
	waveTreeNewSmall.SetName("waveTreeMCMC")
	f1.Close()
	f2.WriteTObject(waveTreeNewSmall)
	waveTreeMCMC = f2.Get("waveTreeMCMC")
	print waveTreeMCMC.GetEntries()

	slowness = np.zeros(1,dtype=float)
	t0 = np.zeros(1,dtype=float)
	scale = np.zeros(1,dtype=float)
	slBranch = waveTreeMCMC.Branch("slowness",slowness,"slowness/D")
	t0Branch = waveTreeMCMC.Branch("t0",t0,"t0/D")
	scBranch = waveTreeMCMC.Branch("scale",scale,"scale/D")
	for i in xrange(waveTreeMCMC.GetEntries()):
		slowness[0] = mcmc_params[i][1]
		t0[0] = mcmc_params[i][2]
		scale[0] = mcmc_params[i][3]
		slBranch.Fill()
		t0Branch.Fill()
		scBranch.Fill()
	waveTreeMCMC.Write("",TObject.kOverwrite)


def mcmc_wfs(waveTemplates, waveTree, eventTree, theCut, run_mcmc=True):

	calDict = {640:28, 672:18, 610:16, 580:19, 582:34, 648:38, 600:21, 578:39, 592:27, 664:55, 626:62, 692:8, 598:22, 690:52, 632:9, 608:7}
	calList = [[key,calDict[key]] for key in calDict]

	waveTree.Draw(">>elist", theCut , "entrylist")
	elist = gDirectory.Get("elist")
	waveTree.SetEntryList(elist)
	print "Found",elist.GetN(),"entries passing cut."

	mcmc_params = []

	# interactive mode
	# print "Press enter to start ..."
	# iList = -1
	# plt.ion()
	# fig = plt.figure(figsize=(13, 8), facecolor='w')
	# while(True):
	# 	iList += 1
	# 	if iList > elist.GetN():
	# 		exit(1)
	# 	value = raw_input()
	# 	if value=='q':		# quit
	# 		exit(1)
	# 	if value=='p': 		# previous entry
	# 		iList -= 2
	# 	if value.isdigit():	# go to entry number
	# 		iList = int(value)

	# batch mode
	print "Running in batch mode ..."
	plt.ioff()
	fig = plt.figure(figsize=(13, 8), facecolor='w')
	for iList in xrange(elist.GetN()):

		# if you're not going to space today, pull the abort switch
		if (run_mcmc == False):
			slowness = np.random.uniform()
			t0 = np.random.uniform()
			scale = np.random.uniform()
			mcmc_params.append([iList,slowness,t0,scale])
			continue

		# ok, now you're going to space today

		entry = waveTree.GetEntryNumber(iList);
		waveTree.LoadTree(entry)
		waveTree.GetEntry(entry)

		trapENFCal = waveTree.trapENFCal	# waveTree parameters
		riseTime = waveTree.blrwfFMR50
		run = waveTree.run
		chan = waveTree.channel
		runTime = waveTree.runTime

		iEvent = waveTree.iEventTree		# eventTree parameters
		eventTree.GetEntry(iEvent)
		power = eventTree.d2wf0MHzTo50MHzPower.at(waveTree.itr)

		signal = wm.CGWave(waveTree.event.GetWaveform(waveTree.itr))
		waveBLSub = signal.GetBLSub()
		waveFilt = signal.GetFilt()
		waveTS = signal.GetTS()
		waveEdge = signal.GetEdge()
		waveEdgeTS = signal.GetEdgeTS()
		loWin, hiWin, lastZero, baseAvg, noiseAvg, timeOff = signal.GetParams()

		# get template waveform

		waveTemplates.GetEntry(calDict[chan])
		tempENF = waveTemplates.trapENFCal
		tempT50 = waveTemplates.blrwfFMR50
		template = wm.CGWave(waveTemplates.event.GetWaveform(waveTemplates.itr))
		tempBLSub = template.GetBLSub()
		tempTS = template.GetTS()
		t_loWin, t_hiWin, t_lastZero, t_baseAvg, t_noiseAvg, t_timeOff = template.GetParams()

		# find an approx 50% rise time of the smoothed edge

		b, a = sg.butter(1, 0.08)
		waveEdgeFilt = sg.filtfilt(b, a, waveEdge)
		edgeMax = np.amax(waveEdgeFilt)
		# edge50idx = np.argmax(edgeMax * 0.5)
		occurences = np.where(waveEdgeFilt > edgeMax * 0.5)
		edge50idx = occurences[0][0]
		t50Time = waveEdgeTS[edge50idx]

		# make the initial guess

		amp = trapENFCal/tempENF
		diff = t50Time - tempT50
		idx2 = np.where((tempTS >= loWin*10) & (tempTS <= hiWin*10))
		tempGuess = tempBLSub[idx2] * amp
		tempTSGuess = tempTS[idx2] + diff

		# run mcmc on rising edge and save the results

		waveModel = pymc.Model( wm.WaveformModel(waveEdge, waveEdgeTS, tempBLSub, tempTS, diff, amp, noiseAvg) )
		M = pymc.MCMC(waveModel)
		M.use_step_method(pymc.Metropolis, M.switchpoint, proposal_sd=4., proposal_distribution='Normal')
		M.use_step_method(pymc.Metropolis, M.slowness, proposal_sd=0.1, proposal_distribution='Normal')
		M.use_step_method(pymc.Metropolis, M.scale, proposal_sd=0.1, proposal_distribution='Normal')
		M.sample(iter=10000, verbose=0)
		burnin = 4000
		t0 = np.median(M.trace('switchpoint')[burnin:])
		slowness = np.median(M.trace('slowness')[burnin:])
		scale = np.median(M.trace('scale')[burnin:])
		mcmc_params.append([iList,slowness,t0,scale])

		# recreate the best-fit signal

		lo = waveEdgeTS[0]
		hi = waveEdgeTS[-1]+10
		idx = np.where((tempTS + t0 >= lo) & (tempTS + t0 <= hi))
		tempStretched = gaussian_filter(tempBLSub[idx], sigma=float(slowness))

		# set up plots

		plt.clf()
		ax1 = plt.subplot2grid((4,5), (0,0), colspan=3,rowspan=2)	# original wf
		ax2 = plt.subplot2grid((4,5), (2,0), colspan=3,rowspan=2) 	# edge + mcmc
		ax3 = plt.subplot2grid((4,5), (0,3), colspan=2) 			# switchpoint trace
		ax4 = plt.subplot2grid((4,5), (1,3), colspan=2, sharex=ax3)	# scale trace
		ax5 = plt.subplot2grid((4,5), (2,3), colspan=2, sharex=ax3)	# slowness trace
		ax6 = plt.subplot2grid((4,5), (3,3), colspan=2)
		ax6.axis('off')

		# plot the original wf, smoothed wf, and time windows
		ax1.plot(waveTS, waveBLSub, color="blue", alpha=0.5, label="full wf")
		ax1.plot(waveTS, waveFilt, color="green", linewidth=3, alpha=0.9, label="filt wf")
		ax1.axvline(x=lastZero*10,color="red",linewidth=2,alpha=0.7, label="last zero")
		ax1.axvline(x=riseTime,color="orange",linewidth=2,alpha=0.7, label="50% rise")
		ax1.axvline(x=loWin*10,color="purple",linewidth=2,alpha=0.7)
		ax1.axvline(x=hiWin*10,color="purple",linewidth=2,alpha=0.7)
		ax1.axis('tight')
		handles,labels = ax1.get_legend_handles_labels()
		ax1.legend(handles, labels, loc='best')
		ax1.set_xlabel('time (ns)')

		# plot the edge and the mcmc fit

		ax2.plot(waveEdgeTS,waveEdge, label='edge wf')
		ax2.plot(tempTSGuess,tempGuess,color='orange', linewidth=2, label="guess")
		# ax2.plot(waveEdgeTS,waveEdgeFilt,color='green', linewidth=2, label="filt")
		# ax2.axvline(x=t50Time,color="red",linewidth=2,alpha=0.7, label="50%")
		ax2.plot(tempTS[idx] + t0,tempBLSub[idx] * scale,color='green',linewidth=2,alpha=0.5, label="scaled+shifted")
		ax2.plot(tempTS[idx] + t0,tempStretched * scale,color='red',linewidth=3, alpha=0.8, label="smoothed")
		ax2.axis('tight')
		handles,labels = ax2.get_legend_handles_labels()
		ax2.legend(handles, labels, loc='best')
		ax2.set_xlabel('time (ns)')

		# plot the full mcmc trace

		ax3.plot(M.trace('switchpoint')[:])
		ax3.set_ylabel('time-diff')
		ax3.set_xlabel('MCMC samples')
		ax3.axis('tight')

		ax4.plot(M.trace('slowness')[:])
		ax4.set_ylabel('slowness')
		ax4.axis('tight')

		ax5.plot(M.trace('scale')[:])
		ax5.set_ylabel('scale')
		ax5.axis('tight')

		results = "%d  %d  %d  %d  %.2fkeV  %.2fns  %dsec  %.1f \n" % (iList,entry,run,chan,trapENFCal,riseTime,runTime,power)

		results += "scale %.3f (%.3f)  t0 %.2f (%.2f)  slow %.3f  (%.2f)\n" % (scale, np.std(M.trace('scale')[burnin:]), t0, np.std(M.trace('switchpoint')[burnin:]), slowness, np.std(M.trace('slowness')[burnin:]))

		print results
		t = ax6.text(-0.1, 0.5, results,size=12)

		plt.tight_layout()

	return waveTree,mcmc_params


def PlotTemplateWaveforms(waveTree,calDict,arg="edge"):
	"""
	Plot the waveform, and then the rising edge of the waveform.
	Changing calDict lets you plot one or many at a time (in one plot).
	'arg' can be "full", "edge", or "both".
	"""

	loWin, hiWin = 0, 2016
	ts = np.arange(loWin*10, hiWin*10, 10)
	plt.style.use('mjWaveforms')  # see ~/.matplotlib/stylelib
	fig = plt.figure(figsize=(10, 7), facecolor='w')
	ax = fig.add_subplot(111)
	plt.xlim(loWin*10, hiWin*10)
	plt.ylabel('ADC')
	plt.xlabel('time [ns]')

	jet = plt.get_cmap('jet')
	cNorm  = colors.Normalize(vmin=578, vmax=672)
	scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)

	# optional: redefine calDict with only one or two channels for debugging
	# calDict = {672:18}

	# plot full waveform
	if (arg == "full" or arg == "both"):
		for chan in calDict:
			waveTree.GetEntry(calDict[chan])
			waveMGT = waveTree.event.GetWaveform(waveTree.itr)
			waveRaw = getWFArray(waveMGT)
			waveOff = waveMGT.GetTOffset()/1e9 # [seconds]

			colorVal = scalarMap.to_rgba(chan)
			label = ('%d' % chan)
			ax.plot(ts, waveRaw, color=colorVal, label=label)

		handles,labels = ax.get_legend_handles_labels()
		ax.legend(handles, labels, loc='upper right')
		plt.show()

	# plot rising edge of waveforms
	lowestWin, highestWin = 5000,0
	if (arg == "edge" or arg == "both"):
		for chan in calDict:
			waveTree.GetEntry(calDict[chan])
			waveMGT = waveTree.event.GetWaveform(waveTree.itr)
			waveRaw = getWFArray(waveMGT)
			baseAvg, noiseAvg = wm.findBaseline(waveRaw)
			waveBLSub = waveRaw
			waveBLSub[:] = [x - baseAvg for x in waveRaw]
			zeros = np.asarray(np.where(waveBLSub < 0.1))
			lastZero = zeros[0,-1]
			loWin, hiWin = lastZero-50, lastZero+200 # 50 and 200
			waveEnds = np.concatenate((np.arange(0, loWin), np.arange(hiWin, 2017)), axis=0)
			tsEdge = np.arange(loWin*10, hiWin*10, 10)
			waveEdge = np.delete(waveBLSub, waveEnds)

			ax.plot(tsEdge,waveEdge)

			lowestWin = loWin if loWin < lowestWin else lowestWin
			highestWin = hiWin if hiWin > highestWin else highestWin

		plt.xlim(lowestWin*10, highestWin*10)
		plt.show()


def fft_waveform(waveBLSub, waveEdge, tsFFTEdge):
	"""look at the FFT's of the WF and the edge
	NOTE: i'm not sure i've gotten the x-scale right.
	NOTE: also, the integral seems wrong.
	"""

	adcFFT = np.fft.rfft(waveBLSub)
	adcPwr = np.abs(adcFFT)**2.
	tsFFT = np.linspace(0, 1./(len(adcFFT)*10*1e-9), len(adcFFT))

	edgeFFT = np.fft.rfft(waveEdge)
	edgePwr = np.abs(edgeFFT)**2.
	tsFFTEdge = np.linspace(0, 1./(len(adcFFT)*10*1e-9), len(edgeFFT))

	idx = np.where(tsFFTEdge > 10000)
	idx2 = np.where(tsFFTEdge[idx] < 32000)
	integral = integrate.simps(edgePwr[idx2],tsFFTEdge[idx2])


if __name__ == "__main__":
	main()
