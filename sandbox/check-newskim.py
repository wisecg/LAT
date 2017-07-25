#!/usr/local/bin/python
#!/usr/common/usg/software/python/2.7.9/bin/python
import ROOT
from ROOT import TFile,TChain,TCanvas,TH1D,TH2D,TH3D,TLegend,TEntryList,gDirectory,TNamed,TObject,gROOT,gStyle

def main():
    # gROOT.ProcessLine("gErrorIgnoreLevel = 3001;") # suppress ROOT messages
    # gStyle.SetOptStat(0)
    gROOT.ProcessLine(".x ~/env/MJDClintPlotStyle.C")
    # gROOT.ForceStyle()

    dsNum = 3
    print "Plotting DS-%d ..." % dsNum

    # old = new = TChain("skimTree") # DON'T DO THIS, it links the chains together
    new, old = TChain("skimTree"), TChain("skimTree")

    mode = "bg"
    old.Add("~/project/newskim/skimDS3_0.root")
    new.Add("~/project/newskim/skimDS3_0_new.root")

    # mode = "cal"
    # old.Add("~/project/newskim/skimDS3_run17959.root")
    # new.Add("~/project/newskim/skimDS3_run17959_new.root")

    h1 = TH1D("h1","h1",100,100,6000)
    new.Project("h1","trapENFCal","gain==0")
    h1.SetLineColor(ROOT.kBlue)
    h1.GetXaxis().SetTitle("Energy (keV)")
    h1.GetYaxis().SetTitle("Counts")
    h1.SetLineWidth(3)

    h2 = TH1D("h2","h2",100,100,6000)
    new.Project("h2","trapENFCal","gain==1")
    h2.SetLineColor(ROOT.kRed)
    h2.SetLineWidth(3)

    h3 = TH1D("h3","h3",100,100,6000)
    old.Project("h3","trapENFCal","gain==0")
    h3.SetLineColor(ROOT.kCyan)

    h4 = TH1D("h4","h4",100,100,6000)
    old.Project("h4","trapENFCal","gain==1")
    h4.SetLineColor(ROOT.kGreen + 2)

    c = TCanvas("","",800,600)
    c.SetLogy()

    h1.Draw("hist")
    h2.Draw("hist same")
    h3.Draw("hist same")
    h4.Draw("hist same")

    l = TLegend(0.7,0.65,0.87,0.92)
    l.AddEntry(h1,"New HG","l")
    l.AddEntry(h2,"New LG","l")
    l.AddEntry(h3,"Old HG","l")
    l.AddEntry(h4,"Old LG","l")
    l.Draw("same")
    c.Print("./plots/checkSkimDS%d_%s.pdf" % (dsNum,mode))




if __name__ == "__main__":
    main()
