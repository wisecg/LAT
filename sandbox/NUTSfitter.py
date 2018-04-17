#!/usr/bin/env python
"""
    Mixture model in pymc3 for unbinned fitting of low energy data for MJD
    Note: Only works in python3.X for now, pm.Dirichlet doesn't seem to work in python2.7
"""

import os
import numpy as np
import pymc3 as pm
import matplotlib.pyplot as plt
import seaborn as sns
import theano.tensor as tt
import pandas as pd
sns.set(style='whitegrid', context='talk', palette='deep', color_codes=True)

def main():
    dsNum, dType = 1, "isNat"
    inDir = os.environ['LATDIR']+'/data/MCMC'
    df = pd.read_hdf('{}/DS{}_Spectrum_{}.h5'.format(inDir, dsNum, dType))
    dfTrit = pd.read_hdf('{}/TritSpec.h5'.format(inDir))

    Y = df['Energy'].values
    # Some values are negative (because Steve) so we need to fix it
    # Also skip some values because who needs them?
    Tritium = dfTrit['Tritium'].values[1::2]
    Tritium = np.array([x if x >= 0. else 0. for x in Tritium])
    Energy = dfTrit['Energy'].values[1::2]

    # Efficiencies -- not implemented yet
    # dfEff = pd.read_hdf('{}/DS{}_{}_Efficiency.h5'.format(inDir, dsNum, dType))
    # XEff = dfEff['Energy'].values
    # YEff = dfEff['Efficiency'].values

    # Dummy tritium spectrum
    # Energy = np.array([0.5, 1.5, 2.5,  3.5, 4.5,  5.5,  6.5,  7.5,  8.5, 9.5, 10.5, 11.5, 12.5, 13.5, 14.5, 15.5, 16.5, 17.5, 18.5])
    # Tritium = np.array([1.72e-02, 1.99e-02, 2.06e-02, 2.04e-02, 1.95e-02, 1.81e-02, 1.66e-02, 1.48e-02, 1.30e-02, 1.11e-02, 9.22e-03, 7.41e-03, 5.72e-03, 4.18e-03, 2.84e-03, 1.72e-03, 8.58e-04, 2.78e-04, 1.54e-05])

    # Peak means
    meansList = [6.54, 8.98, 10.37, 46.54]
    sdList = [GetSigma(mean, dsNum) for mean in meansList]

    with pm.Model() as model:

        Fe55 = pm.Normal("Fe55", mu=meansList[0], sd=sdList[0])
        Zn65 = pm.Normal("Zn65", mu=meansList[1], sd=sdList[1])
        Ge68 = pm.Normal("Ge68", mu=meansList[2], sd=sdList[2])
        Pb210 = pm.Normal("Pb210", mu=meansList[3], sd=sdList[3])
        Bkg = pm.Uniform('Bkg',lower=2., upper=50.)

        Tritium = pm.Interpolated("Tritium", x_points=Energy, pdf_points=Tritium, testval=5.)
        weights = pm.Dirichlet('weights', a=np.ones(len(meansList)+2))

        mix = pm.Mixture('mix', w=weights, comp_dists=[Fe55.distribution, Zn65.distribution, Ge68.distribution, Pb210.distribution, Tritium.distribution, Bkg.distribution], testval=5., observed=Y)

    with model:
        trace = pm.sample(draws=10000, n_init=5000)

    pm.traceplot(trace)
    print(pm.summary(trace))
    # ppc = pm.sample_ppc(trace, samples=500, model=model, size=100)
    # print(np.asarray(ppc['weights']).shape)

    # _, ax = plt.subplots(figsize=(12, 6))
    # ax.hist([n.mean() for n in ppc['n']], bins=19, alpha=0.5)
    # ax.axvline(df['Energy'].mean())
    # ax.set(title='Posterior predictive of the mean', xlabel='mean(x)', ylabel='Frequency');
    plt.show()


# Helper functions to save data
def ReduceData(dsNum, dType):
    """ Saves unbinned energy spectrum into dataframe """
    import ROOT
    inDir, outDir = '/projecta/projectdirs/majorana/users/bxyzhu/cuts/corrfs_rn', os.environ['LATDIR']+'/data/MCMC'
    skimTree = ROOT.TChain("skimTree")
    skimTree.Add("{}/corrfs_rn-DS{}-*.root".format(inDir, dsNum))
    theCut = "trapENFCal > 2 && trapENFCal < 50 && {}".format(dType)
    nPass = skimTree.Draw('trapENFCal:channel', theCut, 'goff')
    print ("{} events passed all cuts".format(nPass))
    nEnergy = skimTree.GetV1()
    nCh = skimTree.GetV2()
    nChList = list(int(nCh[n]) for n in xrange(nPass))
    nEnergyList = list(float(nEnergy[n]) for n in xrange(nPass))

    df = pd.DataFrame({"Energy":nEnergyList, "Channel":nChList})
    print(df.head(10))
    df.to_hdf('{}/DS{}_Spectrum_{}.h5'.format(outDir, dsNum, dType), 'skimTree', mode='w', format='table')


def ExposureEfficiency(dsNum, dType):
    """ Saves efficiency into dataframe -- not done... I forget where I was going with this """
    import ROOT
    inDir, outDir = '/projecta/projectdirs/majorana/users/bxyzhu/LATv2/plots/AThresh', os.environ['LATDIR']+'/data/MCMC'
    f1 = ROOT.TFile("{}/Bkg_{}_DS{}.root".format(inDir, dType, dsNum))
    DS1_isNat_EffTot = f1.Get('DS1_isNat_EffTot')
    DS1_isNat_EffTot.Scale(1./DS1_isNat_EffTot.GetMaximum())

    # Dummy cut efficiency
    ExtPulser = ROOT.TF1('ExtPulser', '0.5*(1+TMath::Erf((x-[0])/(TMath::Sqrt(2)*[1]) ))', 1, 250)
    ExtPulser.SetParameters(0.36,1.26)

    hExtPulser = ROOT.TH1D('hExtPulser', 'hExtPulser', 250000, 0, 250)
    for i in range(DS1_isNat_EffTot.GetNbinsX()):
        if DS1_isNat_EffTot.GetBinCenter(i) < 0.99: continue
        DS1_isNat_EffTot.SetBinContent(i, TotEffEmpirical.Eval(DS1_isNat_EffTot.GetBinCenter(i)))
        hExtPulser.SetBinContent(i, ExtPulser.Eval(hExtPulser.GetBinCenter(i)))

    df = pd.DataFrame({"Energy":nEnergyList, "Channel":nChList})
    print(df.head(10))
    # df.to_hdf('{}/DS{}_Spectrum_{}.h5'.format(outDir, dsNum, dType), 'skimTree', mode='w', format='table')


class ErrorFnc(pm.Continuous):
    """
        Custom Error Function distribution for pymc3
    """
    def __init__(self, mu, sd, *args, **kwargs):
        self._mu = tt.as_tensor_variable(mu)
        self._sd = tt.as_tensor_variable(sd)
        super(ErrorFnc, self).__init__(*args, **kwargs)

    def logp(self, value):
        """ Log-likelihood of Erf efficiency """
        mu = self._mu
        sigma = self._sd
        return tt.log( 0.5 + (1.+ tt.erf(value-mu)/(sigma*tt.sqrt(2))))


def SaveTrit():
    import ROOT
    """ Saves tritium spectrum into a pandas dataframe """
    inDir, outDir = '/projecta/projectdirs/majorana/users/bxyzhu/Axion', os.environ['LATDIR']+'/data/MCMC'
    f1 = ROOT.TFile('{}/TritSpec.root'.format(inDir))
    tritHist = f1.Get('tritHist')
    tritHist.Scale(1./tritHist.Integral("w"))
    energy = [tritHist.GetBinCenter(xbin) for xbin in range(1, tritHist.GetNbinsX())]
    val = [tritHist.GetBinContent(xbin) for xbin in range(1, tritHist.GetNbinsX())]

    df = pd.DataFrame({"Energy":energy, "Tritium":val})
    print(df.head(10))
    df.to_hdf('{}/TritSpec.h5'.format(outDir), 'TritSpec', mode='w', format='table')


def GetSigma(energy, dsNum=1):
    p0, p1, p2 = 0.,0.,0.
    if dsNum==0:
        p0 = 0.147; p1=0.0173; p2=0.0003
    elif dsNum==1:
        p0 = 0.136; p1=0.0174; p2=0.00028
    elif dsNum==3:
        p0 = 0.162; p1=0.0172; p2=0.000297
    elif dsNum==4:
        p0 = 0.218; p1=0.015; p2=0.00035
    elif dsNum==5:
        p0 = 0.2121; p1=0.01838; p2=0.00031137
    else:
        p0 = 0.2121; p1=0.01838; p2=0.00031137
    return np.sqrt(p0*p0 + p1*p1*energy + p2*p2*energy*energy)


if __name__ == '__main__':
    main()
