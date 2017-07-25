#!/usr/common/usg/software/python/2.7.9/bin/python
#!/usr/local/bin/python
import sys
sys.argv += [ '-b' ]  # force ROOT to be loaded in batch mode
import ROOT
from ROOT import gROOT, gStyle, gPad
from ROOT import TFile, TTree, TChain, TCanvas, TH1D, TH2D, TF1, TLegend, TLine, TGraph
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import waveLibs as wl

# chBkgDS0 = [640, 608, 674, 592, 644, 662, 646, 696, 626, 688, 664, 594, 692, 610, 598, 690, 600, 624, 656, 576]
# chBkgDS1 = [608, 640, 610, 580, 582, 648, 600, 578, 592, 632, 626, 692, 598, 690, 664, 672]
# # chBkgDS2 =
# chBkgDS3 = [672, 640, 610, 580, 678, 614, 648, 632, 598, 578, 592, 664, 626, 692, 694, 690, 600, 624, 582, 608]
# chBkgDS4 = [1296, 1172, 1204, 1106, 1170, 1232, 1298, 1136, 1176, 1330, 1332, 1174, 1144, 1236]
# chBkgDS5 = [640, 1302, 648, 1298, 660, 662, 1176, 1330, 1332, 658, 1208, 1204, 1174, 1106, 1110, 1120, 610, 1124, 614, 1128, 1170, 626, 628, 1172]

def main(argv):

    # -- ROOT options --
    gROOT.ProcessLine(".x ~/env/MJDClintPlotStyle.C")
    gROOT.ProcessLine("gErrorIgnoreLevel = 3001;") # suppress ROOT messages
    # gStyle.SetOptStat(0)
    # gROOT.ForceStyle()


    # -- Load data --
    # tree = TChain("skimTree")
    # tree.Add("~/project/cal-lat/latSkimDS1_run14149.root")
    # tree.Add("~/project/lat/latSkimDS1_*.root")
    # cutFile = TFile("~/project/lat/latSkimDS1_0_0.root")
    # cutFile = TFile("~/project/cal-lat/latSkimDS1_run14149.root")
    # theCut = cutFile.Get("theCut").GetTitle()


    # -- Make plots --
    # fitSlo238(tree,theCut) # cal
    # fitSlo10(tree,theCut) # bg
    # fitArb10(tree,theCut)
    # fitBCMax()
    # fitSloChannels()
    fitSloChannelsFast( int(argv[0]) )


def H1D(tree,nameStr,bins,xlo,xhi,drawStr,cutStr):
    h1 = ROOT.TH1D(nameStr,nameStr,bins,xlo,xhi)
    tree.Project(nameStr,drawStr,cutStr)
    return h1


def npTH1D(hist):
    # this is the CORRECTED function
    bins = hist.GetNbinsX()
    xArr, yArr = np.zeros(bins),np.zeros(bins)
    for i in range(bins):
        xArr[i] = hist.GetXaxis().GetBinCenter(i)
        yArr[i] = hist.GetBinContent(i)
    return xArr,yArr


def integFunc(arr):
    integ = np.zeros(len(arr))
    sum = 0
    for i in range(0,len(arr)):
        sum+=arr[i]
        integ[i] = sum
    return integ


def GetIntegralPoints(hist):
    x_h0, y_h0 = npTH1D(hist)
    int_h0 = integFunc(y_h0)

    idx99 = np.where(int_h0 > 0.99)
    idx95 = np.where(int_h0 > 0.95)
    idx01 = np.where(int_h0 > 0.01)
    idx05 = np.where(int_h0 > 0.05)

    val99 = x_h0[idx99][0]
    val95 = x_h0[idx95][0]
    val01 = x_h0[idx01][0]
    val05 = x_h0[idx05][0]

    return val99,val95,val01,val05


def fitSlo238(tree,theCut):

    # for calibration data, keep 99% of the fast events in the 238 peak.
    # this is practice for doing this with the 10 kev peak in natural bg data.

    # scatter plot
    h1 = wl.H2D(tree,"h1",500,0,250,500,0,200,"fitSlo:trapENFCal",theCut)
    cts1 = h1.Integral( *wl.Get2DBins(h1,0,250,0,200) )

    # -- 238 peak --

    h2 = wl.H1D(tree,"h2",50,237,240,"trapENFCal",theCut)
    cts2 = h2.Integral( *wl.Get1DBins(h2,237,240) )
    h2.SetLineColor(4)

    f0 = TF1("f0","gaus",237,240)
    h2.Fit("f0","Lq") # use to set initial guesses

    f1 = TF1("f1", "[0] * x  + [1] * exp(-1.0 * (TMath::Power((x-[2]),2) / (2 * TMath::Power([3],2)) ))", 237, 240)
    wl.SetPars(f1,[1,418,238,0.381])
    f1.SetLineColor(4)
    h2.Fit("f1","q")
    [flat, norm1, mu, sig] = wl.GetPars(f1)

    h3 = wl.H1D(tree,"h3",50,237,240,"trapENFCal",theCut+" && fitSlo < 20")
    cts3 = h3.Integral( *wl.Get1DBins(h3,237,240) )

    f2 = TF1("f2", "[0] * x  + [1] * exp(-1.0 * (TMath::Power((x-[2]),2) / (2 * TMath::Power([3],2)) ))", 237, 240)
    wl.SetPars(f2,[1,418,238,0.381])
    f2.SetLineColor(2)
    h3.Fit("f2","q")
    [flat, norm2, mu, sig] = wl.GetPars(f2)

    retention = 100*norm2/norm1
    print "norm1 %.2f  norm2 %.2f  retention: %.2f" % (norm1, norm2, retention)


    # -- Plots --
    c1 = TCanvas("c1","Bob Ross's Canvas",1200,600)

    c1.Divide(2,1) # x, y
    c1.cd(1)
    gPad.SetLogz()
    h1.Draw("COLZ")

    c1.cd(2)
    h2.Draw()
    h3.SetLineColor(2) # red
    h3.Draw("same")

    l1 = TLegend(0.6,0.7,0.87,0.92)
    l1.AddEntry(h2,"basic","l")
    l1.AddEntry(f1,"basic fit:","l")
    # l1.AddEntry(f1,"flat %.2f  norm %.2f" % (flat,norm1) ,"")
    # l1.AddEntry(f1,"mu %.2f  sig %.2f" % (mu,sig) ,"")
    l1.AddEntry(h3,"+fitSlo < 20","l")
    l1.AddEntry(h3,"Retention: %.3f" % retention,"")

    l1.Draw("same")
    c1.Update()

    c1.Print("./plots/cal-fitSlo.pdf")


def fitSlo10(tree,theCut):

    theCut += " && !isEnr"

    # scatter plot
    h0 = wl.H2D(tree,"h0",100,0,20,500,0,200,"fitSlo:trapENFCal",theCut)

    # -- 10 peak --
    fitModel10 = "[0] + [1]*x + [2]*x**2 + [3] * exp(-1.0 * ((x-[4])**2 / (2 * [5]**2)) )"

    h1 = wl.H1D(tree,"h1",20,9.5,11.5,"trapENFCal",theCut)
    h1.SetLineColor(4)

    # f1 = TF1("f1","gaus",10,11)
    # h1.Fit("f1","L") # use to set initial guesses

    f1 = TF1("f1", fitModel10, 9.5,11)
    wl.SetPars(f1,[-209, 41, -2, 28, 10.3, 0.135])
    f1.SetLineColor(4)
    h1.Fit("f1","Lq")
    [p0,p1,p2,norm1,mu,sig] = wl.GetPars(f1)
    # print wl.GetPars(f1)

    h2 = wl.H1D(tree,"h2",20,9.5,11.5,"trapENFCal",theCut+" && fitSlo < 18")
    f2 = TF1("f2", fitModel10, 9.5,11)
    wl.SetPars(f2,[-209, 41, -2, 28, 10.3, 0.135])
    f2.SetLineColor(2)
    h2.Fit("f2","Lq")
    [p0,p1,p2,norm2,mu,sig] = wl.GetPars(f2)

    retention1 = 100*norm2/norm1


    # -- 6.54 (Fe55) peak --
    # TODO: This doesn't work well.  wait till you're using RooFit.

    fitModel6 = "[0] * x + [1] * exp( -1.0 * ((x - [2])**2 / (2 * [3]**2)) )"
    # fitModel6 = "[0] * exp( -1.0 * ((x - [1])**2 / (2 * [2]**2)) )"

    h3 = wl.H1D(tree,"h3",20,5.5,7.5,"trapENFCal",theCut)

    # f3 = TF1("f3","gaus",6,7)
    # h3.Fit("f3","L") # use to set initial guesses

    # f3 = TF1("f3", fitModel6, 6,7)
    # wl.SetPars(f3,[-1, 19.4, 6.47, 0.631])
    # f3.SetLineColor(4)
    # h3.Fit("f3","Lq")
    # [lin,norm3,mu,sig] = wl.GetPars(f3)
    #
    # h4 = wl.H1D(tree,"h4",20,5.5,7.5,"trapENFCal",theCut+" && fitSlo < 18")
    # f4 = TF1("f4", fitModel6, 6,7)
    # wl.SetPars(f4,[-1, 19.4, 6.47, 0.631])
    # f4.SetLineColor(2)
    # h4.Fit("f4","Lq")
    # [lin,norm4,mu,sig] = wl.GetPars(f4)

    # retention2 = 100*norm4/norm3
    # print "norm3 %.2f  norm4 %.2f  retention2: %.2f" % (norm3, norm4, retention2)


    # -- Plots --
    c1 = TCanvas("c1","Bob Ross's Canvas",1200,600)

    # c1.Divide(3,1,0.00001) # TPad::Divide
    c1.Divide(2,1,0.00001)
    c1.cd(1)
    gPad.SetLogz()
    h0.SetMinimum(1)
    h0.Draw("COLZ")

    c1.cd(2)
    h1.SetLineColor(4) # blue
    h1.Draw()
    h2.SetLineColor(2) # red
    h2.Draw("same")

    l1 = TLegend(0.6,0.6,0.87,0.92)
    l1.AddEntry(h1,"basic","l")
    l1.AddEntry(h2,"+fitSlo < 20","l")
    l1.AddEntry(h2,"Retention: %.3f" % retention1,"")
    l1.Draw("same")

    # c1.cd(3)
    # h3.Draw()
    # h4.SetLineColor(2) # red
    # h4.Draw("same")
    #
    # l2 = TLegend(0.6,0.6,0.87,0.92)
    # l2.AddEntry(h3,"basic","l")
    # l2.AddEntry(h4,"+fitSlo < 20","l")
    # l2.AddEntry(h4,"Retention: %.3f" % retention2,"")
    # l2.Draw("same")

    c1.Print("./plots/ds1-fitSlo.pdf")


def fitArb10(tree,theCut):
    # fit an arbitrary low-energy cut

    print theCut

    theCut += " && !isEnr"
    checkCut = ""
    thisVar = "bandTime"
    fileName = "./plots/lat-max.pdf"
    parLim = [0,20000]

    # scatter plot
    h0 = wl.H2D(tree,"h0",100,0,20,500,parLim[0],parLim[1],thisVar+":trapENFCal",theCut)

    # -- 10 peak --
    fitModel10 = "[0] + [1]*x + [2]*x**2 + [3] * exp(-1.0 * ((x-[4])**2 / (2 * [5]**2)) )"

    h1 = wl.H1D(tree,"h1",20,9.5,11.5,"trapENFCal",theCut)
    h2 = wl.H1D(tree,"h2",20,9.5,11.5,"trapENFCal",theCut+checkCut)

    f0 = TF1("f0","gaus",10,11)
    h1.Fit("f0","Lq")
    [norm,mu,sig] = wl.GetPars(f0)

    f1 = TF1("f1", fitModel10, 9.5,11)
    wl.SetPars(f1,[-209, 41, -2, norm, mu, sig])
    f1.SetLineColor(4)
    h1.Fit("f1","Lq")
    [p0,p1,p2,norm1,mu,sig] = wl.GetPars(f1)

    f2 = TF1("f2", fitModel10, 9.5,11)
    wl.SetPars(f2,[-209, 41, -2, norm1, mu, sig])
    f2.SetLineColor(2)
    h2.Fit("f2","Lq")
    [p0,p1,p2,norm2,mu,sig] = wl.GetPars(f2)

    retention = 100*norm2/norm1
    print "Retention:",retention," Cut used:",checkCut

    # -- Energy Spectrum --
    h3 = wl.H1D(tree,"h3",250,0,50,"trapENFCal",theCut)
    h4 = wl.H1D(tree,"h4",250,0,50,"trapENFCal",theCut+checkCut)

    cts3 = h3.Integral( *wl.Get1DBins(h3,0,20))
    cts4 = h4.Integral( *wl.Get1DBins(h4,0,20))
    cts5 = h3.Integral( *wl.Get1DBins(h3,20,40))
    cts6 = h4.Integral( *wl.Get1DBins(h4,20,40))

    print "Cts. 0-20 %d before  %d after.  20-40 %d before %d after" % (cts3,cts4,cts5,cts6)

    h3.GetXaxis().SetRangeUser(0,20)
    h4.GetXaxis().SetRangeUser(0,20)


    # -- Plots --
    c1 = TCanvas("c1","Bob Ross's Canvas",1800,600)

    c1.Divide(3,1,0.00001) # TPad::Divide
    c1.cd(1)
    gPad.SetLogz()
    h0.SetMinimum(1)
    h0.Draw("COLZ")

    c1.cd(2)
    gPad.SetLogy()
    h3.SetLineColor(4)
    h3.Draw()
    h4.SetLineColor(2)
    h4.Draw("same")
    l1 = TLegend(0.4,0.7,0.87,0.92)
    l1.AddEntry(h1,"basic","l")
    l1.AddEntry(h2,thisVar,"l")
    l1.Draw("same")

    c1.cd(3)
    h1.SetLineColor(4) # blue
    h1.Draw()
    h2.SetLineColor(2) # red
    h2.Draw("same")

    l2 = TLegend(0.6,0.6,0.87,0.92)
    l2.AddEntry(h1,"basic","l")
    l2.AddEntry(h2,thisVar,"l")
    l2.AddEntry(h2,"Retention: %.3f" % retention,"")
    l2.Draw("same")

    c1.Print(fileName)


def fitBCMax():

    tree = TChain("skimTree")
    tree.Add("~/project/cal-lat/latSkimDS3*")

    cutFile = TFile("~/project/lat/latSkimDS3_0_0.root")
    theCut = cutFile.Get("theCut").GetTitle()

    # scatter plot
    h0 = wl.H2D(tree,"h0",100,0,20,500,parLim[0],parLim[1],thisVar+":trapENFCal",theCut)

    # fitModel10 = "[0] + [1]*x + [2]*x**2 + [3] * exp(-1.0 * ((x-[4])**2 / (2 * [5]**2)) )"

    h1 = wl.H1D(tree,"h1",20,9.5,11.5,"trapENFCal",theCut)
    h2 = wl.H1D(tree,"h2",20,9.5,11.5,"trapENFCal",theCut+checkCut)

    f0 = TF1("f0","gaus",10,11)
    h1.Fit("f0","Lq")
    [norm,mu,sig] = wl.GetPars(f0)

    f1 = TF1("f1", fitModel10, 9.5,11)
    wl.SetPars(f1,[-209, 41, -2, norm, mu, sig])
    f1.SetLineColor(4)
    h1.Fit("f1","Lq")
    [p0,p1,p2,norm1,mu,sig] = wl.GetPars(f1)

    f2 = TF1("f2", fitModel10, 9.5,11)
    wl.SetPars(f2,[-209, 41, -2, norm1, mu, sig])
    f2.SetLineColor(2)
    h2.Fit("f2","Lq")
    [p0,p1,p2,norm2,mu,sig] = wl.GetPars(f2)

    retention = 100*norm2/norm1
    print "Retention:",retention," Cut used:",checkCut


def fitSloChannels():
    ds = 1

    f1 = TFile("~/project/lat/latSkimDS%d_0_0.root"%ds)
    theCut = f1.Get("theCut").GetTitle()
    calib = TChain("skimTree")
    calib.Add("~/project/cal-lat/latSkimDS%d*.root"%ds)

    c = TCanvas("c","c",1600,1200)
    c.Divide(3,2,0.00001,0.00001)

    # ch = 626 # the golden detector

    for ch in chCalDS1:

        # 10-200 spectrum
        var = "fitSlo"
        eb, elo, ehi = 10,10,200
        vb, vlo, vhi = 200,0,30
        d2Draw = "fitSlo:trapENFCal>>a(%d,%d,%d,%d,%d,%d)"%(eb,elo,ehi,vb,vlo,vhi)
        d2Cut = theCut + " && channel==%d" % ch
        d1Cut = theCut + " && trapENFCal > 5 && channel==%d" % ch

        c.cd(2)
        gPad.SetLogy(0)
        h0 = TH1D("h0","h0",vb,vlo,vhi)
        calib.Project("h0",var,d1Cut)
        h0.Scale(1/h0.Integral())
        h0.Draw("hist")

        c.cd(1)
        gPad.SetLogy(1)
        calib.Draw(d2Draw,d2Cut)
        sp99,sp95,sp01,sp05 = GetIntegralPoints(h0)
        line99 = TLine()
        line99.SetLineColor(ROOT.kGreen)
        line95 = TLine()
        line95.SetLineColor(ROOT.kRed)
        line99.DrawLine(elo, sp99, ehi, sp99)
        line95.DrawLine(elo, sp95, ehi, sp95)
        line95.DrawLine(elo, sp05, ehi, sp05)
        line99.DrawLine(elo, sp01, ehi, sp01)

        c.cd(3)
        x_h0, y_h0 = npTH1D(h0)
        int_h0 = integFunc(y_h0)
        g1 = TGraph(len(x_h0), x_h0, int_h0)
        g1.Draw("ACP")
        line99.DrawLine(sp99, 0, sp99, 1)
        line95.DrawLine(sp95, 0, sp95, 1)
        line99.DrawLine(sp01, 0, sp01, 1)
        line95.DrawLine(sp05, 0, sp05, 1)


        # 238 peak spectrum
        eb, elo, ehi = 10,237,240
        d2Draw = "fitSlo:trapENFCal>>b(%d,%d,%d,%d,%d,%d)"%(eb,elo,ehi,vb,vlo,vhi)
        d2Cut = theCut + " && channel==%d" % ch
        d1Cut = theCut + " && trapENFCal > %d && trapENFCal < %d && channel==%d" % (elo,ehi,ch)

        c.cd(5)
        gPad.SetLogy(0)
        h1 = TH1D("h1","h1",vb,vlo,vhi)
        calib.Project("h1",var,d1Cut)
        h1.Scale(1/h1.Integral())
        h1.Draw("hist")

        c.cd(4)
        gPad.SetLogy(0)
        calib.Draw(d2Draw,d2Cut)
        pk99,pk95,pk01,pk05 = GetIntegralPoints(h1)
        line99 = TLine()
        line99.SetLineColor(ROOT.kGreen)
        line95 = TLine()
        line95.SetLineColor(ROOT.kRed)
        line99.DrawLine(elo, pk99, ehi, pk99)
        line95.DrawLine(elo, pk95, ehi, pk95)
        line95.DrawLine(elo, pk05, ehi, pk05)
        line99.DrawLine(elo, pk01, ehi, pk01)

        c.cd(6)
        x_h1, y_h1 = npTH1D(h1)
        int_h1 = integFunc(y_h1)
        g2 = TGraph(len(x_h1), x_h1, int_h1)
        g2.Draw("ACP")
        line99.DrawLine(pk99, 0, pk99, 1)
        line95.DrawLine(pk95, 0, pk95, 1)
        line99.DrawLine(pk01, 0, pk01, 1)
        line95.DrawLine(pk05, 0, pk05, 1)

        c.Print("./plots/fitSlo/ds%d_ch%d.png" % (ds,ch))

        print "DS-%d  Channel %d:" % (ds,ch)
        print "  Peak - 99: %.2f  95: %.2f  5: %.2f  1: %.2f" % (sp99,sp95,sp05,sp01)
        print "  Spec - 99: %.2f  95: %.2f  5: %.2f  1: %.2f" % (pk99,pk95,pk05,pk01)



def fitSloChannelsFast(dsNum):

    # Just get the 95% upper value from the peak.
    # NOTE: the slowness model might not be good at the 238 peak (S/N too high ??)
    # Ben is shaking his head in disapproval.

    ds = dsNum
    chCal = {}
    chCal[0] = [640, 644, 646, 656, 662, 664, 674, 688, 690, 692, 696, 576, 592, 594, 598, 600, 608, 610, 624, 626]
    chCal[1] = [672, 640, 610, 580, 578, 582, 648, 600, 626, 592, 664, 594, 692, 598, 690, 632, 608]
    chCal[3] = [608, 624, 610, 580, 614, 582, 648, 664, 694, 578, 678, 592, 600, 626, 692, 598, 690, 632, 640, 672]
    chCal[4] = [1232, 1236, 1332, 1298, 1330, 1170, 1296, 1136, 1176, 1106, 1204, 1174, 1144, 1172]
    chCal[5] = [640, 1302, 648, 1170, 1172, 662, 1176, 672, 678, 680, 1330, 688, 690, 1332, 694, 1208, 1204, 1110, 1174, 608, 584, 1298, 592, 1106, 598, 1120, 610, 1124, 614, 1128, 658, 624, 626, 628, 632, 660]

    f1 = TFile("~/project/lat/latSkimDS%d_0_0.root"%ds)
    theCut = f1.Get("theCut").GetTitle()
    calib = TChain("skimTree")
    calib.Add("~/project/cal-lat/latSkimDS%d*.root"%ds)

    c = TCanvas("c","c",1600,600)
    c.Divide(3,1,0.00001,0.00001)

    # ch = 626 # the golden detector
    for ch in chCal[ds]:

        var = "fitSlo"
        eb, elo, ehi = 10,237,240
        vb, vlo, vhi = 200,0,30
        d2Draw = "fitSlo:trapENFCal>>b(%d,%d,%d,%d,%d,%d)"%(eb,elo,ehi,vb,vlo,vhi)
        d2Cut = theCut + " && channel==%d" % ch
        d1Cut = theCut + " && trapENFCal > %d && trapENFCal < %d && channel==%d" % (elo,ehi,ch)

        c.cd(2)
        gPad.SetLogy(0)
        h1 = TH1D("h1","h1",vb,vlo,vhi)
        calib.Project("h1",var,d1Cut)
        h1.Scale(1/h1.Integral())
        h1.Draw("hist")

        c.cd(1)
        gPad.SetLogy(0)
        calib.Draw(d2Draw,d2Cut)
        pk99,pk95,pk01,pk05 = GetIntegralPoints(h1)
        line99 = TLine()
        line99.SetLineColor(ROOT.kGreen)
        line95 = TLine()
        line95.SetLineColor(ROOT.kRed)
        line99.DrawLine(elo, pk99, ehi, pk99)
        line95.DrawLine(elo, pk95, ehi, pk95)
        line95.DrawLine(elo, pk05, ehi, pk05)
        line99.DrawLine(elo, pk01, ehi, pk01)

        c.cd(3)
        x_h1, y_h1 = npTH1D(h1)
        int_h1 = integFunc(y_h1)
        g2 = TGraph(len(x_h1), x_h1, int_h1)
        g2.Draw("ACP")
        line99.DrawLine(pk99, 0, pk99, 1)
        line95.DrawLine(pk95, 0, pk95, 1)
        line99.DrawLine(pk01, 0, pk01, 1)
        line95.DrawLine(pk05, 0, pk05, 1)

        c.Print("./plots/fitSlo/fast_ds%d_ch%d.png" % (ds,ch))
        print "DS-%d  Ch %d  fitSlo95 %.2f " % (ds,ch,pk95)





if __name__ == "__main__":
    main(sys.argv[1:])
