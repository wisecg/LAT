#!/usr/local/bin/python
#!/usr/common/usg/software/python/2.7.9/bin/python
""" feed-the-beast.py
    Low energy data cleaning plot factory.
    C. Wiseman, 3/25/17
"""
import ROOT
from ROOT import TFile,TChain,TCanvas,TH1D,TH2D,TH3D,TLegend,TEntryList,gDirectory,TNamed,TObject,gROOT,gStyle

def main():
    # gROOT.ProcessLine("gErrorIgnoreLevel = 3001;") # suppress ROOT messages
    # gStyle.SetOptStat(0)
    gROOT.ProcessLine(".x ~/env/MJDClintPlotStyle.C")
    # gROOT.ForceStyle()

    dsList = [3]

    # checkSkims()
    # validateFiles()
    # mLVsTOffset(dsList)
    # waveletPlots2D(dsList)
    # butterPlots2D(dsList)
    # cutSpectra(dsList)
    # cutVsRun(dsList)
    # noisyRuns(dsList)
    # allCuts(dsList)
    # thresholds(dsList)
    checkEnrNat([1])

# =======================================================

def checkSkims():
    """ Due diligence: make sure nothing surprising got cut. """

    theCut = "trapENFCal > 0.8 && gain==0 && mH==1 && isGood && !muVeto && !wfDCBits && !isLNFill1 && !isLNFill2 && trapETailMin < 0 && channel!=596 && channel!=676 && channel!=676 && channel!=612 && channel!=1104 && channel!=1200 && channel!=1334 && channel!=1336"

    for dsNum in [0,1,2,4]:
        print "Plotting DS-%d ..." % dsNum

        ch=ch2=ch3=ch4 = TChain("skimTree")
        ch.Add(str("~/project/skim/skimDS%d*" % dsNum))
        ch2.Add(str("~/project/skim-basic/skim-basicDS%d*" % dsNum))
        ch3.Add(str("~/project/wave-skim/waveSkimDS%d*" % dsNum))
        ch4.Add(str("~/project/wavelet-skim/waveletSkimDS%d*" % dsNum))

        h1 = TH1D("h1","h1",100,0,100)
        h2 = TH1D("h2","h2",100,0,100)
        h3 = TH1D("h3","h3",100,0,100)
        h4 = TH1D("h4","h4",100,0,100)
        ch.Project("h1","trapENFCal",theCut)
        ch2.Project("h2","trapENFCal",theCut)
        ch3.Project("h3","trapENFCal",theCut)
        ch4.Project("h4","trapENFCal",theCut)

        c = TCanvas("","",800,600)
        c.SetLogy()
        h1.SetLineColor(ROOT.kBlack)
        h2.SetLineColor(ROOT.kBlue)
        h3.SetLineColor(ROOT.kGreen)
        h4.SetLineColor(ROOT.kRed)
        h1.Draw("hist")
        h2.Draw("hist same")
        h3.Draw("hist same")
        h4.Draw("hist same")

        l = TLegend(0.7,0.65,0.87,0.92)
        l.AddEntry(h1,"skim","l")
        l.AddEntry(h2,"basic","l")
        l.AddEntry(h3,"+wave","l")
        l.AddEntry(h4,"wavelet","l")
        l.Draw("same")
        c.Update()
        c.Print("./plots/checkSkimDS%d.pdf" % dsNum)


def validateFiles():
    """ wave-view is saying it can't find an MGTWaveforms branch sometimes. """

    import glob
    files = glob.glob("/Users/wisecg/project/wavelet-skim/waveletSkimDS3*")

    for inFile in files:
        f = TFile(inFile)
        if (f.IsZombie()):
            print "file dead!  continuing ..."
            continue

        tree = f.Get("skimTree")

        found = False
        for branch in tree.GetListOfBranches():
            if( branch.GetName() == "MGTWaveforms" ):
                found = True

        if found: print inFile, " .... found!"
        else: print inFile, " .... branch found!"


def mLVsTOffset(dsList):
    """ is mL==1 the same as the tOffset cut? """

    theCut = "trapENFCal > 0.8 && gain==0 && mH==1 && isGood && !muVeto && !wfDCBits && !isLNFill1 && !isLNFill2 && trapETailMin < 0 && channel!=596 && channel!=676 && channel!=676 && channel!=612 && channel!=1104 && channel!=1200 && channel!=1334 && channel!=1336"

    bn, lo, hi = 100,0,20
    c = TCanvas("","",800,600)
    for dsNum in dsList:
        print "Plotting DS-%d ..." % dsNum
        ch = TChain("skimTree")
        ch.Add(str("~/project/wavelet-skim/waveletSkimDS%d*" % dsNum))

        h1 = TH1D("h1","h1",bn,lo,hi)
        ch.Project("h1","trapENFCal",theCut)
        h1.GetXaxis().SetTitle("Energy (keV)")
        h1.SetLineColor(ROOT.kBlack)

        h2 = TH1D("h2","h2",bn,lo,hi)
        ch.Project("h2","trapENFCal",theCut + " && tOffset > 10")
        h2.SetLineColor(ROOT.kRed)

        h3 = TH1D("h3","h3",bn,lo,hi)
        ch.Project("h3",theCut + " && mL==1")
        h3.SetLineColor(ROOT.kBlue)

        c.SetLogy()
        h1.Draw("hist")
        h2.Draw("hist same")
        h3.Draw("hist same")
        l1 = TLegend(0.7,0.60,0.87,0.87)
        l1.AddEntry(h1,"basic","l")
        l1.AddEntry(h2,"tOffset","l")
        l1.AddEntry(h3,"mL==1","l")
        l1.Draw("same")
        c.Print("./plots/mL-tOffset_DS%d.pdf" % dsNum)


def waveletPlots2D(dsList):

    c = TCanvas("","",800,600)
    for dsNum in dsList:
        print "Plotting DS-%d ..." % dsNum
        ch = TChain("skimTree")
        ch.Add(str("~/project/wavelet-skim/waveletSkimDS%d*" % dsNum))

        h1 = TH2D("h1","h1",100,0,3000,1000,0,7000)
        ch.Project("h1","waveS5/TMath::Power(trapENFCal,1/4):trapENFCal",theCut)
        h1.GetXaxis().SetTitle("Energy (keV)")
        h1.GetYaxis().SetTitle("waveS5")
        h1.Draw()
        c.Print("./plots/waveS5CutDS%d_3000.pdf" % dsNum)

        h2 = TH2D("h2","h2",100,0,20,1000,0,7000)
        ch.Project("h2","waveS5/TMath::Power(trapENFCal,1/4):trapENFCal",theCut)
        h2.GetXaxis().SetTitle("Energy (keV)")
        h2.GetYaxis().SetTitle("waveS5")
        h2.Draw()
        c.Print("./plots/waveS5CutDS%d_20.pdf" % dsNum)


def butterPlots2D(dsList):

    basicCut = "trapENFCal > 0.8 && gain==0 && mH==1 && isGood && !muVeto && !wfDCBits && !isLNFill1 && !isLNFill2 && trapETailMin < 0 && channel!=596 && channel!=676 && channel!=676 && channel!=612 && channel!=1104 && channel!=1200 && channel!=1334 && channel!=1336"

    bcrCut = " && tOffset < 10 && waveS5/TMath::Power(trapENFCal,1/4) < 1200 && (waveS3-waveS2)/trapENFCal > 100 && (waveS3-waveS2)/trapENFCal < 300"

    # bcrCut += " && (waveS3-waveS4)/trapENFCal < 100 && (waveS3-waveS4)/trapENFCal > 0"  # this doesn't do much

    ds3noisyRunCut = " && !(channel==692 && (run==16974 || run==16975 || run==16976 || run==16977 || run==16979))"

    bn, lo, hi = 100,0,20
    c = TCanvas("","",800,600)

    for dsNum in dsList:
        print "Plotting DS-%d ..." % dsNum
        ch = TChain("skimTree")
        ch.Add(str("~/project/wavelet-skim/waveletSkimDS%d*" % dsNum))

        h1 = TH2D("h1","h1",100,0,20,100,0,22000)
        ch.Project("h1","butterTime:trapENFCal",basicCut + bcrCut + ds3noisyRunCut)
        h1.GetXaxis().SetTitle("Energy (keV)")
        h1.GetYaxis().SetTitle("Deriv.MaxTime")
        h1.GetYaxis().SetTitleOffset(1.5)
        h1.Draw()
        c.SetLeftMargin(0.15)
        c.Print("./plots/butterTime_DS%d.pdf" % dsNum)


def cutSpectra(dsList):

    basicCut = "trapENFCal > 0.8 && gain==0 && mH==1 && isGood && !muVeto && !wfDCBits && !isLNFill1 && !isLNFill2 && trapETailMin < 0 && channel!=596 && channel!=676 && channel!=676 && channel!=612 && channel!=1104 && channel!=1200 && channel!=1334 && channel!=1336"

    bcrCut = " && tOffset < 10 && waveS5/TMath::Power(trapENFCal,1/4) < 1200 && (waveS3-waveS2)/trapENFCal > 100 && (waveS3-waveS2)/trapENFCal < 300"

    # bcrCut += " && (waveS3-waveS4)/trapENFCal < 100 && (waveS3-waveS4)/trapENFCal > 0"  # this doesn't do much

    ds3noisyRunCut = " && !(channel==692 && (run==16974 || run==16975 || run==16976 || run==16977 || run==16979))"

    butterCut = " && butterTime > 4000 && butterTime < 10500"

    bn, lo, hi = 100,0,20
    c = TCanvas("","",800,600)
    c.SetLogy()

    for dsNum in dsList:
        print "Plotting DS-%d ..." % dsNum
        ch = TChain("skimTree")
        ch.Add(str("~/project/wavelet-skim/waveletSkimDS%d*" % dsNum))

        h1 = TH1D("h1","h1",bn,lo,hi)
        ch.Project("h1","trapENFCal",basicCut)
        h1.GetXaxis().SetTitle("Energy (keV)")
        h1.SetLineColor(ROOT.kBlack)

        h2 = TH1D("h2","h2",bn,lo,hi)
        ch.Project("h2","trapENFCal",basicCut + " && tOffset < 10")
        h2.SetLineColor(ROOT.kBlue)

        h3 = TH1D("h3","h3",bn,lo,hi)
        ch.Project("h3","trapENFCal",basicCut + " && waveS5/TMath::Power(trapENFCal,1/4) < 1200")
        h3.SetLineColor(ROOT.kRed)

        h4 = TH1D("h4","h4",bn,lo,hi)
        ch.Project("h4","trapENFCal",basicCut + " && (waveS3-waveS2)/trapENFCal > 100 && (waveS3-waveS2)/trapENFCal < 300")
        h4.SetLineColor(ROOT.kGreen+3)

        h5 = TH1D("h5","h5",bn,lo,hi)
        ch.Project("h5","trapENFCal",basicCut + " && (waveS3-waveS4)/trapENFCal < 100 && (waveS3-waveS4)/trapENFCal > 0")
        h5.SetLineColor(ROOT.kCyan)

        h6 = TH1D("h6","h6",bn,lo,hi)
        ch.Project("h6","trapENFCal",basicCut + " && butterTime > 4000 && butterTime < 10500")
        h6.SetLineColor(ROOT.kMagenta)

        h7 = TH1D("h7","h7",bn,lo,hi)
        ch.Project("h7","trapENFCal",basicCut + bcrCut + ds3noisyRunCut + butterCut)
        h7.SetLineColor(ROOT.kOrange+7)

        h1.Draw("hist")
        h2.Draw("hist same")
        h3.Draw("hist same")
        h4.Draw("hist same")
        h5.Draw("hist same")
        h6.Draw("hist same")
        h7.Draw("hist same")
        l1 = TLegend(0.67,0.62,0.87,0.9)
        l1.AddEntry(h1,"basic","l")
        l1.AddEntry(h2,"tOffset","l")
        l1.AddEntry(h3,"S5/E^(1/4)","l")
        l1.AddEntry(h4,"(S3-S2)/E","l")
        l1.AddEntry(h5,"(S3-S4)/E","l")
        l1.AddEntry(h6,"butterTime","l")
        l1.AddEntry(h7,"+noisy runs","l")
        l1.Draw("same")
        c.Print("./plots/cutSpectra_DS%d.pdf" % dsNum)


def cutVsRun(dsList):

    theCut = "trapENFCal > 0.8 && gain==0 && mH==1 && isGood && !muVeto && !wfDCBits && !isLNFill1 && !isLNFill2 && trapETailMin < 0 && channel!=596 && channel!=676 && channel!=676 && channel!=612 && channel!=1104 && channel!=1200 && channel!=1334 && channel!=1336"
    bn, lo, hi = 1300,16700,18000
    c = TCanvas("","",800,600)

    for dsNum in dsList:
        print "Plotting DS-%d ..." % dsNum
        ch = TChain("skimTree")
        ch.Add(str("~/project/wavelet-skim/waveletSkimDS%d*" % dsNum))

        h1 = TH1D("h1","h1",bn,lo,hi)
        ch.Project("h1","run",theCut)
        h1.GetXaxis().SetTitle("Run")
        h1.SetLineColor(ROOT.kBlack)

        h2 = TH1D("h2","h2",bn,lo,hi)
        ch.Project("h2","run",theCut + " && tOffset < 10")
        h2.SetLineColor(ROOT.kBlue)

        h3 = TH1D("h3","h3",bn,lo,hi)
        ch.Project("h3","run",theCut + " && waveS5/TMath::Power(trapENFCal,1/4) < 1200")
        h3.SetLineColor(ROOT.kRed)

        h4 = TH1D("h4","h4",bn,lo,hi)
        ch.Project("h4","run",theCut + " && (waveS3-waveS2)/trapENFCal > 100 && (waveS3-waveS2)/trapENFCal < 300")
        h4.SetLineColor(ROOT.kGreen+3)

        # 1. cut comparison
        h1.Draw("hist")
        h2.Draw("hist same")
        h3.Draw("hist same")
        h4.Draw("hist same")
        l1 = TLegend(0.7,0.6,0.87,0.87)
        l1.AddEntry(h1,"basic","l")
        l1.AddEntry(h2,"tOffset","l")
        l1.AddEntry(h3,"S5/E^(1/4)","l")
        l1.AddEntry(h4,"(S3-S2)/E","l")
        l1.Draw("same")
        c.Print("./plots/cutVsRun_DS%d.pdf" % dsNum)

        # 2. total cut + noisy runs cut
        c.SetLogy()
        h5 = TH1D("h5","h5",bn,lo,hi)
        ch.Project("h5","run",theCut + " && tOffset < 10 && waveS5/TMath::Power(trapENFCal,1/4) < 1200 && (waveS3-waveS2)/trapENFCal > 100 && (waveS3-waveS2)/trapENFCal < 300")
        h5.SetLineColor(ROOT.kOrange+7)


        h1.SetMinimum(1)
        h1.Draw("hist")
        h5.Draw("hist same")
        h6.Draw("hist same")
        c.Print("./plots/cutTotalVsRun_DS%d.pdf" % dsNum)


def noisyRuns(dsList):

    theCut = "trapENFCal > 0.8 && gain==0 && mH==1 && isGood && !muVeto && !wfDCBits && !isLNFill1 && !isLNFill2 && trapETailMin < 0 && channel!=596 && channel!=676 && channel!=676 && channel!=612 && channel!=1104 && channel!=1200 && channel!=1334 && channel!=1336"

    theCut += " && tOffset < 10 && waveS5/TMath::Power(trapENFCal,1/4) < 1200 && (waveS3-waveS2)/trapENFCal > 100 && (waveS3-waveS2)/trapENFCal < 300"

    # theCut += " && run!=16974 && run!=16975 && run!=16976 && run!=16977 && run!=16979"  # DS3 noisy runs

    bn, lo, hi = 1300,16700,18000
    # bn, lo, hi = 20,16970,16990
    c = TCanvas("","",800,600)

    for dsNum in dsList:
        print "Plotting DS-%d ..." % dsNum
        ch = TChain("skimTree")
        ch.Add(str("~/project/wavelet-skim/waveletSkimDS%d*" % dsNum))

        h1 = TH1D("h1","h1",bn,lo,hi)
        ch.Project("h1","run",theCut)
        h1.SetLineColor(ROOT.kOrange+7)

        h1.SetMinimum(1)
        h1.Draw("hist")
        h1.GetXaxis().SetNdivisions(505)
        c.Print("./plots/noisyRuns_DS%d.pdf" % dsNum)

        # find and print noisy runs
        noisyRuns = []
        for i in range(bn):
            rate = h1.GetBinContent(i+1)
            if rate > 200: noisyRuns.append((i+lo,rate))
        for run,rate in noisyRuns:
            print run,rate
        # result: 16974, 16975, 16976, 16977, 16979 have rate over 200 cts/run

        # ... so which CHANNELS are responsible for this high rate?


def allCuts(dsList):

    basicCut = "trapENFCal > 0.8 && gain==0 && mH==1 && isGood && !muVeto && !wfDCBits && !isLNFill1 && !isLNFill2 && trapETailMin < 0 && channel!=596 && channel!=676 && channel!=676 && channel!=612 && channel!=1104 && channel!=1200 && channel!=1334 && channel!=1336"

    bcrCut = " && tOffset < 10 && waveS5/TMath::Power(trapENFCal,1/4) < 1200 && (waveS3-waveS2)/trapENFCal > 100 && (waveS3-waveS2)/trapENFCal < 300"

    # bcrCut += " && (waveS3-waveS4)/trapENFCal < 100 && (waveS3-waveS4)/trapENFCal > 0"  # this doesn't do much

    ds3noisyRunCut = " && !(channel==692 && (run==16974 || run==16975 || run==16976 || run==16977 || run==16979))"

    butterCut = " && butterTime > 5000 && butterTime < 10500"

    bn, lo, hi = 100,0,20
    c = TCanvas("","",800,600)

    for dsNum in dsList:
        print "Plotting DS-%d ..." % dsNum
        ch = TChain("skimTree")
        ch.Add(str("~/project/wavelet-skim/waveletSkimDS%d*" % dsNum))

        h1 = TH1D("h1","h1",bn,lo,hi)
        ch.Project("h1","trapENFCal",basicCut + bcrCut + ds3noisyRunCut + butterCut)
        h1.GetXaxis().SetTitle("Energy (keV)")
        h1.SetLineColor(ROOT.kBlue)
        h1.SetMaximum(100)

        h1.Draw("hist")
        c.Print("./plots/allCuts_DS%d.pdf" % dsNum)


def checkEnrNat(dsList):
    """."""
    basicCut = "trapENFCal > 0.8 && gain==0 && mH==1 && isGood && !muVeto && !wfDCBits && !isLNFill1 && !isLNFill2 && trapETailMin < 0 && channel!=596 && channel!=676 && channel!=676 && channel!=612 && channel!=1104 && channel!=1200 && channel!=1334 && channel!=1336"

    # bcrCut = " && tOffset < 10 && waveS5/TMath::Power(trapENFCal,1/4) < 1200 && (waveS3-waveS2)/trapENFCal > 100 && (waveS3-waveS2)/trapENFCal < 300"
    # bcrCut += " && (waveS3-waveS4)/trapENFCal < 100 && (waveS3-waveS4)/trapENFCal > 0"  # this doesn't do much
    # ds3noisyRunCut = " && !(channel==692 && (run==16974 || run==16975 || run==16976 || run==16977 || run==16979))"
    # butterCut = " && butterTime > 5000 && butterTime < 10500"

    # bn, lo, hi = 100,0,20
    c = TCanvas("","",800,600)

    for dsNum in dsList:
        print "Plotting DS-%d ..." % dsNum
        ch = TChain("skimTree")
        ch.Add(str("~/project/wavelet-skim/waveletSkimDS%d*" % dsNum))

        h1 = TH2D("h1","h1",2,-0.5,1.5,2,-0.5,1.5)
        ch.Project("h1","isEnr:isNat",basicCut)

        h1.Draw("COLZ")
        c.Print("./plots/enrVsNat%d.pdf" % dsNum)



if __name__ == "__main__":
    # main(sys.argv[1:])
    main()
