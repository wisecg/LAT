#!/usr/bin/env python

import random
import matplotlib.pyplot as plt
import numpy as np
from ROOT import TFile, TTree, TCanvas, TH1D
from ROOT import gStyle

gStyle.SetOptStat(0)

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


def getBkgIndexes():

    mbexpo = 89.5

    f1 = TFile("../data/malbek_data.root")
    t1 = f1.Get("malbek_wrt")

    c = TCanvas("c","Bob Ross's Canvas",800,600)

    kpb, eLo, eHi = 0.1, 1., 30.
    bins = int((eHi-eLo) / kpb + 0.5)

    h1 = H1D(t1,bins,eLo,eHi,"energy_keV","weight","Energy (keV)","Counts"," ")

    ax1 = h1.GetXaxis()
    print "malbek:",h1.Integral(ax1.FindBin(1.5), ax1.FindBin(8.))

    # h1.Scale(1./mbexpo)
    h1.Draw("hist")

    c.Print("../plots/malbek-raw.pdf")


    ds0exp = 478.
    f2 = TFile("../data/m1DiagnosticTree4Aug2016.root")
    t2 = f2.Get("diagTree")

    h2 = H1D(t2,bins,eLo,eHi,"calENF","enr","Energy (keV)","Counts/kg-d"," ")

    ax2 = h2.GetXaxis()
    ds0int = h2.Integral(ax2.FindBin(5.), ax2.FindBin(8.),"width")
    print "vorren:",ds0int," B:",6.5*ds0int/4.

    h2.Scale(1./ds0exp)
    h2.Draw("hist")
    c.Print("../plots/vorren-raw.pdf")


def getExpectedCts():

        f = TFile("../data/inputHists.root")
        hConv = f.Get("h4")
        ax = hConv.GetXaxis()

        malbekExpo = 89.5 # kg-d
        nExp = hConv.Integral(ax.FindBin(1.5), ax.FindBin(8.), "width") * malbekExpo # [cts / (keV d kg)] * [kg d]
        print "nExp",nExp

        nObs = 40.8
        print "g_ae U.L.:", np.power(nObs/nExp, 1./4.)


def ratioMethod():

    dsExpos = {"malbek w/ds1bg":0.245, "1":1.81, "1+5b":3.66, "1-4,5b":5.24, "1-5":8.70, "0nbb 0-5":9.95, "1-4,5b,6o":7.42, "1-4,5b,6o,6c":11.82, "thru oct 2017":20.58}

    # bgIdx = {"malbek":0.834, "ds0":0.025, "ds1+":0.00625} # integral from 5-8 kev, div. by exposure. ds1+ is just x4 lower than ds0
    bgIdx = {"malbek":50., "ds1+":1.99}

    exps = np.arange(0.2,30,0.1)

    def sensFunc(expo):
        return ((bgIdx["malbek"] * expo) / (bgIdx["ds1+"] * dsExpos["malbek w/ds1bg"]))**0.25

    fig = plt.figure(figsize=(8,6),facecolor='w')

    plt.plot(exps, 2.6e-11/sensFunc(exps), "b", label="proj assuming ds1 bgIdx")

    ctr = 0
    cmap = plt.cm.get_cmap('hsv',len(dsExpos))
    for key in dsExpos:
        val = dsExpos[key]
        plt.plot(val,2.6e-11/sensFunc(val),'o',c=cmap(ctr),label=key)
        ctr += 1

    plt.axvline(dsExpos["1-4,5b"],color="green",label="recommended cutoff")

    plt.axhline(4.35e-12,color="red",label="most recent pandaX limit")

    plt.legend(loc="best")
    plt.xlabel("Exposure (kg-y)")
    plt.ylabel("g_ae Projected U.L.")

    plt.show(block=False)
    plt.savefig("../plots/gae-sensitivity.png")

    print "g_ae at cutoff",2.6e-11/sensFunc(dsExpos["1-4,5b"])


def directMethod():

    def gaeFunc(expo, f):
        # f is the fraction that axions contribute to the total
        # background, using a background index proportional to the # counts in DS1.
        # b = (E2-E1) * (\int_5^8 (DS-0 bkg) dE)/4
        b = 19.98
        axCts = 2.3119e44
        nObs = f * b
        return (nObs / (expo * axCts))**0.25

    dsExpos = {"malbek w/ds1bg":0.245, "1":1.81, "1+5b":3.66, "1-4,5b":5.24, "1-5":8.70, "0nbb 0-5":9.95, "1-4,5b,6o":7.42, "1-4,5b,6o,6c":11.82, "thru oct 2017":20.58}

    exps = np.arange(0.2,30,0.1)

    fig = plt.figure(figsize=(10,8),facecolor='w')

    cmap = plt.cm.get_cmap('viridis',4)
    for idx, f in enumerate([0.1, 0.03, 0.01]):
        plt.plot(exps, gaeFunc(exps,f), c=cmap(idx+1), label="proj, ax/bkg=%.2f" % f)

        cmap2 = plt.cm.get_cmap('Paired',len(dsExpos)+1)
        for idx2, key in enumerate(dsExpos):
            val = dsExpos[key]
            if idx==0:
                plt.plot(val,gaeFunc(val,f),'o',c=cmap2(idx2),label=key)
            else:
                plt.plot(val,gaeFunc(val,f),'o',c=cmap2(idx2))

    plt.axvline(dsExpos["1-4,5b"],color="green",label="recommended cutoff")

    plt.axhline(4.35e-12,color="red",label="most recent pandaX limit")

    plt.legend(loc="best")
    plt.xlabel("Exposure (kg-y)")
    plt.ylabel("g_ae Projected U.L.")

    plt.show(block=False)
    plt.savefig("../plots/gae-sensitivity.pdf")

    print "g_ae at cutoff, w/ axion cts < 0.03: ",gaeFunc(dsExpos["1-4,5b"],0.03)





if __name__=="__main__":

    # getBkgIndexes()
    # getExpectedCts()
    # ratioMethod()
    directMethod()

