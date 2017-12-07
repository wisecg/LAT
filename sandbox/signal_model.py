#!/usr/local/bin/python
import numpy as np
from pymc import Uniform, Normal, deterministic, Container

def CreateBGModel(data, ene, flat, tritSpec, ge68Peak, fe55Peak, zn65Peak, flatSpec):

	# flat = Container([])
	# flatSpec = Container([flat.append(1.0) for i in xrange(len(ene))])

	# print type(flatSpec)

	# flat = []
	# for i in xrange(len(ene)):
		# flat.append(1.0)
	# flatSpec = np.asarray(flat)


	# stochastic variables (not compleletely detmined by parent values; they have a prob. distribution.)
	# need to float resolution, energy, and offset for each peak ...
	tritScale = Uniform('tritScale', lower=0, upper=10)
	flatScale = Uniform('flatScale', lower=0, upper=10)
	fe55Scale = Uniform('fe55Scale', lower=0, upper=10)
	zn65Scale = Uniform('zn65Scale', lower=0, upper=10)
	ge68Scale = Uniform('ge68Scale', lower=0, upper=10)

	# deterministic variables : given by values of parents
	# set up the model for uncertainty (ie, the noise) and the signal (ie, the spectrum)

	@deterministic(plot=False, name="uncertainty")
	def uncertainty_model(s=tritScale):
		out = s * np.ones(len(data))
		return out

	@deterministic
	def tau(eps=uncertainty_model):
		# pymc uses this tau parameter instead of sigma to model a gaussian.  its annoying.
		return np.power(eps, -2)

	@deterministic(plot=False, name="BGModel")
	def signal_model(t=tritScale, f=flatScale, fe=fe55Scale, zn=zn65Scale, ge=ge68Scale):
		theModel = t * tritSpec + fe * fe55Peak + zn * zn65Peak + ge * ge68Peak
		return theModel

	# full model
	baseline_observed = Normal("baseline_observed", mu=signal_model, tau=tau, value=data, observed=True)

	return locals()

def sigToTau(sig):
	tau = np.power(sig, -2)
	return tau

