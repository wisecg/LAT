#!/usr/bin/env python
import random
import numpy as np
import matplotlib.pyplot as plt
from ROOT import TFile, TTree, TH1D, TCanvas, gStyle, gROOT

import matplotlib.pylab as pylab
params = {'legend.fontsize': 'large',
         'axes.labelsize': 'large',
         'axes.titlesize':'large',
         'xtick.labelsize':'large',
         'ytick.labelsize':'large'}
pylab.rcParams.update(params)

gae_mb = 2.6e-11
ex_mb = 89.5 / 365.
ex_ds0 = 478. / 365.

""" MJD g_ae Sensitivity vs. Exposure estimator.

    g_ae = (N_observed / N_expected)**(1./4.)

    This compares a "ratio method" using Graham's MALBEK results,
    and a "direct method" which assumes knowledge of the convolved axion flux.

    Ratio method:
    g_ae' = g_ae ((B'/Mt')(Mt/B))**(1./8.)

    Direct method:
    g_ae = (B dE / Mt)**(1./8.) * C**(-1./4.)

    The curves are mainly sensitive to assumptions on B.
    I can't exactly reproduce Graham's 6.5e-12 number, but he's not too specific
    about what his assumptions were.  (pg. 109 just refers to Marino's thesis, and
    Marino's thesis just gives the 0.03 number for average tritium rate.

    "Direct method" seems to be the best choice (I trust my axion flux spectrum),
    so then I plot the current MJD exposures on top of the curve.

    TODO: Calculate MJD background index from data.
    Eh, maybe just wait until we have the data cleaned spectra.
"""

def main():

    eLo, eHi = 1.5, 8.
    expos = np.arange(0.1, 100, 0.1)

    B = tritBkgIndex(eLo,eHi,0.03)
    B_mb = malbekBkgIndex(eLo,eHi)
    A = axFluxConst(eLo,eHi)
    print "B_mjd: %.2f  B_mb: %.2f  [cts/(kev-kg-d)]  A: %.2e [cts/kg-d]" % (B*365.25*(eHi-eLo), B_mb, A)
    # B = 0.03

    fig = plt.figure(figsize=(8,7),facecolor='w')

    # plot ratio method and direct method
    # plt.plot(expos, ratioMethod(expos,0.02,eLo,eHi), c='r', label='ratio, B=0.02 cts/kev/kg/d')
    # plt.plot(expos, ratioMethod(expos,0.08,eLo,eHi), c='r', label='ratio, B=0.08 cts/kev/kg/d')
    plt.plot(expos, directMethod(expos,B,eLo,eHi), c='b', label='$g_{ae}$ (proj), B=%.2f cts/kev/kg/d' % B)

    # brian's ratio method, for comparison.
    # B_bri = 1534 # cts/kg-y of malbek, depends on energy region.  Should be 1360 for 1.5--8 kev.
    # brian_gae = [np.power(np.sqrt(0.03 * 6.5 * 365.25)/(np.sqrt(x)) * ex_mb/np.sqrt(B_bri), 0.25) * 2.6e-11 for x in expos]
    # plt.plot(expos, brian_gae, c='g',label = 'brian, B=0.03')

    # plt.axhline(2.59e-11,c='g',label="edelweiss 2013") # we're way better than this
    plt.axhline(directMethod(5.24,B,eLo,eHi),c='g',label='Current MJD analysis')
    plt.axhline(6.5e-12,c='orange',label="Graham's 100kg-y MJD proj")
    plt.axhline(4.35e-12,c='r',label="PandaX PRL, Nov. 2017")

    dsExpos = {"DS 1":1.81, "DS 1+5b":3.66, "DS 1-4+5b":5.24, "DS 1-4+5b+6o":7.42, "DS 1-4+5b+6":11.82}
    ctr = 0
    cmap = plt.cm.get_cmap('hsv',len(dsExpos)+1)
    for key, val in sorted(dsExpos.iteritems(), key=lambda (k,v): (v,k)): # trick to sort by value
        val = dsExpos[key]
        plt.plot(val,directMethod(val,B,eLo,eHi),'o',c=cmap(ctr),label="%s : %.2f kg-y" % (key,val))
        ctr += 1

    plt.xlabel("Exposure (kg-y)")
    plt.ylabel("$g_{ae}$ Projected U.L.",rotation=0)
    ax = fig.gca()
    ax.yaxis.set_label_coords(0.,1.04)

    plt.legend(loc='best')
    plt.savefig("../plots/gae-proj.pdf")

    print "100 kg-y g_ae,  ratio method:", ratioMethod(expos[-1],B,eLo,eHi)
    print "               direct method:", directMethod(expos[-1],B,eLo,eHi)
    # print "              brian's method:", brian_gae[-1]


def directMethod(ex,B,eLo,eHi):
    C = axFluxConst(eLo, eHi)
    return ((B * (eHi-eLo) * 365.25)/(ex))**0.125 * C**-0.25


def ratioMethod(ex,B,eLo,eHi):
    B_mb = malbekBkgIndex(eLo, eHi)
    B_mjd = B * 365.25 * (eHi-eLo)
    return gae_mb * ((B_mjd / ex) * (ex_mb / B_mb))**0.125


def axFluxConst(eLo, eHi):
    """ A : [cts / kg-y] from the convolved axion flux spectrum. """
    f = TFile("../data/inputHists.root")
    hConv = f.Get("h4")
    ax = hConv.GetXaxis()
    cts = hConv.Integral(ax.FindBin(eLo), ax.FindBin(eHi), "width") * 365.25 # [cts / kg-y]
    return cts


def malbekBkgIndex(eLo, eHi):
    """ B : average cts/kg-y in an energy region. """
    f = TFile("../data/malbek_data.root")
    t = f.Get("malbek_wrt")
    kpb, eLo, eHi = 0.1, eLo, eHi
    bins = int((eHi-eLo) / kpb + 0.5)
    h = H1D(t,bins,eLo,eHi,"energy_keV","weight","Energy (keV)","Counts"," ")
    ax = h.GetXaxis()
    cts = h.Integral(ax.FindBin(eLo), ax.FindBin(eHi),"width") # "width" is same as multiplying by 'kpb'
    return cts / ex_mb # [cts / kg-y]


def mjdBkgIndex():
    """ B : average cts/kg-y in the DS0 5-8 kev energy region, divided by 4. """
    f = TFile("../data/m1DiagnosticTree4Aug2016.root")
    t = f.Get("diagTree")
    h = H1D(t,bins,eLo,eHi,"calENF","enr","Energy (keV)","Counts/kg-d"," ")
    ax = h.GetXaxis()
    ds0int = h.Integral(ax.FindBin(5.), ax.FindBin(8.))


def tritBkgIndex(eLo, eHi, B):
    """ B [cts / (keV-kg-d)].  Proportional to the 3H activation rate.
    B = 0.03 from Marino's thesis is an average from 0-18 keV for 15 d of surface exposure.
    """
    f = TFile("../data/inputHists.root")
    hTrit = f.Get("h5")
    ax = hTrit.GetXaxis()
    # normCts = hTrit.Integral("width")
    # print normCts # check that it's normalized to 1.

    # Scale the full histo s/t the average counts match 'B'.
    avg = 0.
    for idx in range(ax.GetNbins()):
        avg += hTrit.GetBinContent(idx)
    avg /= hTrit.GetXaxis().GetNbins()
    hTrit.Scale(B/avg)

    # compute the average between eLo and eHi
    avg2 = 0.
    for idx in range(ax.FindBin(eLo), ax.FindBin(eHi)):
        avg2 += hTrit.GetBinContent(idx)
    avg2 /= (ax.FindBin(eHi)-ax.FindBin(eLo))

    # gStyle.SetOptStat(0)
    # gROOT.ProcessLine("gErrorIgnoreLevel = 3001;") # suppress ROOT messages
    # c = TCanvas("c","c",800,600)
    # hTrit.GetXaxis().SetRangeUser(0,20)
    # hTrit.Draw('hist')
    # hTrit.SetTitle(" ")
    # hTrit.GetXaxis().SetTitle("Energy (keV)")
    # hTrit.GetYaxis().SetTitle("cts / keV / kg / d")
    # c.Print("../plots/tritPlot.pdf")

    return avg2


def H1D(tree,bins,xlo,xhi,drawStr,cutStr,xTitle="",yTitle="",Title=None, Name=None):
    nameStr, titleStr = "", ""
    if Name == None: nameStr = str(random.uniform(1.,2.))
    else: nameStr = str(Name)
    if Title == None: titleStr = str(random.uniform(1.,2.))
    else: titleStr = str(Title)
    h1 = TH1D(nameStr,titleStr,bins,xlo,xhi)
    tree.Project(nameStr,drawStr,cutStr)
    if xTitle!="": h1.GetXaxis().SetTitle(xTitle)
    if yTitle!="": h1.GetYaxis().SetTitle(yTitle)
    return h1



if __name__=="__main__":
    main()