#!/usr/local/bin/python
from ROOT import *
import numpy as np
import matplotlib.pyplot as plt
import signal_model as sm
import pymc

def main():

	bpk = 5 # bins per keV
	plt.style.use('mjWaveforms')  # see ~/.matplotlib/stylelib

	# load un-scaled histograms into np arrays

	f1 = TFile("m1DiagnosticTree4Aug2016.root")
	t1 = f1.Get("diagTree")

	hds0enr = TH1D("ds0enr","ds0enr",bpk*25,5,30)
	t1.Project("ds0enr","calENF","enr")

	hds0nat = TH1D("ds0nat","ds0nat",bpk*25,5,30)
	t1.Project("ds0nat","calENF","!enr")

	f2 = TFile("2nuDetec.root")
	hTrit = f2.Get("tritHist")
	hFlat = f2.Get("NatflatHist")

	f3 = TFile("enrVsNat_DS1.root")
	hds1nat = f3.Get("nat")
	hds1enr = f3.Get("enr")

	ds0enr = rootToArray(hds0enr,5,30)
	ds0nat = rootToArray(hds0nat,5,30)
	ds1enr = rootToArray(hds1enr,5,30)
	ds1nat = rootToArray(hds1nat,5,30)
	ene = np.arange(5, 30, 1/float(bpk))

	tritSpec = rootToArray(hTrit,5,30)
	flatSpec = rootToArray(hFlat,5,30)
	tritSpec = tritSpec*10000
	for i in range(len(tritSpec)):
		if tritSpec[i] < 0:
			tritSpec[i]=0

	""" Graham's MALBEK thesis, pg. 78
	peak	energy (keV)	half life
	68Ge	1.3, 10.37		271.0 d
	71Ge	1.3, 10.37		11.4 d
	68Ga	9.66			67.7 m
	65Zn	1.10, 8.98		243.9 d
	55Fe	6.54			2.74 y
	49V		4.97			330 d

	MJ detector resolution: ~0.25 keV at 10.3
	Float this but only a small amount.
	"""

	ge68, fe55, zn65, flat = [], [], [], []
	for i in xrange(len(ene)):
		ge68.append(gauss_function(ene[i]+0.31,200.,10.37,0.25))
		fe55.append(gauss_function(ene[i]+0.31,200.,6.54,0.25))
		zn65.append(gauss_function(ene[i]+0.31,200.,8.98,0.25))
		flat.append(1.0)
	ge68Peak = np.asarray(ge68)
	fe55Peak = np.asarray(fe55)
	zn65Peak = np.asarray(zn65)
	flatSpec = np.asarray(flat)

	mcmc_spec(ds0nat, ene, flat, tritSpec, ge68Peak, fe55Peak, zn65Peak, flatSpec)

	# check histograms

	# fig = plt.figure(figsize=(10, 7), facecolor='w')
	# ax = fig.add_subplot(111)
	#
	# plt.plot(ene,ds0nat,label='ds0nat',ls='steps-post')
	#
	# plt.plot(ene,gausArr,label='gaus',ls='steps-post')
	#
	# handles,labels = ax.get_legend_handles_labels()
	# ax.legend(handles, labels, loc='upper right')
	# plt.show()


def mcmc_spec(spec, ene, flat, tritSpec, ge68Peak, fe55Peak, zn65Peak, flatSpec):


	spec_model = pymc.Model( sm.CreateBGModel(spec, ene, flat, tritSpec, ge68Peak, fe55Peak, zn65Peak, flatSpec) )
	M = pymc.MCMC(spec_model)

	M.use_step_method(pymc.Metropolis, M.tritScale, proposal_sd=0.3, proposal_distribution='Normal')
	M.use_step_method(pymc.Metropolis, M.flatScale, proposal_sd=1., proposal_distribution='Normal')
	M.use_step_method(pymc.Metropolis, M.fe55Scale, proposal_sd=1., proposal_distribution='Normal')
	M.use_step_method(pymc.Metropolis, M.zn65Scale, proposal_sd=1., proposal_distribution='Normal')
	M.use_step_method(pymc.Metropolis, M.ge68Scale, proposal_sd=1., proposal_distribution='Normal')

	# Run the MCMC
	M.sample(iter=10000)

	# pull out the fit parameters after burnin'
	burnin = 4000
	tritScale = np.median(M.trace('tritScale')[burnin:])
	flatScale = np.median(M.trace('flatScale')[burnin:])
	fe55Scale = np.median(M.trace('fe55Scale')[burnin:])
	zn65Scale = np.median(M.trace('zn65Scale')[burnin:])
	ge68Scale = np.median(M.trace('ge68Scale')[burnin:])


	# plot the fit result on top of the original spectrum

	fitResult = tritScale * tritSpec + ge68Scale * ge68Peak

	print "trit %.3f  flat %.3f  fe55 %.3f  zn65 %.3f  ge68 %.3f" % (tritScale, flatScale, fe55Scale, zn65Scale, ge68Scale)

	fig = plt.figure(figsize=(10, 7), facecolor='w')
	ax = fig.add_subplot(111)

	plt.plot(ene,spec,label='spec',ls='steps-post')
	plt.plot(ene,fitResult,label='bg model',ls='steps-post',color='red')

	handles,labels = ax.get_legend_handles_labels()
	ax.legend(handles, labels, loc='best')
	plt.show()


def rootToArray(hist, xLow, xHigh):
	"""Take a ROOT TH1D histogram and a range, get a numpy array.

	Note on plotting:
	Can't use np.histogram, the hist has already been done.
	So make a plot that fakes it with :
		xaxis = np.arange(xlow, xhi, step)
		plt.bar(xaxis,hist,width=step)  # 1: fake it with bar
		plt.plot(xaxis,hist,ls='steps-post')  # 2: fake it with plot
	"""
	binLow = hist.FindBin(xLow)
	binHigh = hist.FindBin(xHigh)
	loopArray = range(binLow, binHigh )
	npArray = np.empty_like(loopArray, dtype=np.float)
	for (iArray, iBin) in enumerate(loopArray):
		npArray[iArray] = np.float(hist.GetBinContent(iBin))
	return npArray

def gauss_function(x, a, x0, sigma):
	return a * np.exp(-(x - x0)**2 / (2 * sigma**2))


if __name__ == "__main__":
	main()
