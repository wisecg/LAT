// Data cleaning.
// Creates plots to justify TCuts.
// C. Wiseman, 12/11/2016

#include <iostream>
#include <fstream>
#include "TFile.h"
#include "TROOT.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TEntryList.h"
#include "GATDataSet.hh"

using namespace std;

void LargeFileCuts(int dsNum);
void TCutSkimmer(int dsNum);
void GetExposure(int dsNum, double &enrExpo, double &natExpo);
void GenerateSpectra(int dsNum);
void CombineSpectra();
void ThresholdCut(int dsNum);
void LowESpectra(int dsNum);

// ====================================================================================

int main(int argc, char** argv)
{
  gStyle->SetOptStat(0);
  gROOT->ProcessLine(".x ~/env/MJDClintPlotStyle.C");
  gROOT->ForceStyle();
  // gROOT->ProcessLine( "gErrorIgnoreLevel = 2001;");

  int dsNum = 0;
  if (argc > 1) dsNum = stoi(argv[1]);
  // LargeFileCuts(dsNum);
  // TCutSkimmer(dsNum);
  // GenerateSpectra(dsNum);
  CombineSpectra();
  // LowESpectra(dsNum);
  // ThresholdCut(dsNum);
}

// ====================================================================================

char theCut[1000];

char basicCut[1000] = "trapENFCal > 0 && trapENFCal < 3000 && gain==0 && isGood && !muVeto && !wfDCBits && !isLNFill1 && !isLNFill2 && channel!=596 && channel!=676 && channel!=676 && channel!=612 && channel!=1104 && channel!=1200 && channel!=1334 && channel!=1336";

char ds1burstCut[1000] = "!(time_s > 2192e3 && time_s < 2195e3) && !(time_s > 7370e3 && time_s <7371e3) && !(time_s > 7840e3 && time_s < 7860e3) && !(time_s > 8384e3 && time_s < 8387e3) && !(time_s > 8984e3 && time_s < 8985e3) && !(time_s > 9002e3 && time_s < 9005e3) && run != 13075 && run != 13093 && run != 13116";

char ds1noisyRunsCut[1000] = "run!=9648 && run!=10663 && run!=10745 && run!=11175 && run!=12445 && run!=12723 && run!=12735 && run!=12745 && run!=12746 && run!=12765 && run!=12766 && run!=12767 && run!=13004";

char ds0_toeCut[1000] = "(kvorrT/trapENFCal) >= 1.2 && (kvorrT/trapENFCal) <= 2.1";

char ds1_toeCut[1000] = "(((kvorrT/trapENFCal) >= 1.2 && (kvorrT/trapENFCal) <= 2.1) || (channel==580 && ((kvorrT/trapENFCal) >= 0.9 && (kvorrT/trapENFCal) <= 2.1)) || (channel==664 && ((kvorrT/trapENFCal) >= 1.1 && (kvorrT/trapENFCal) <= 2.1)))";

char ds1standardCut[1000] = "trapENFCal > 0 && trapENFCal < 3000 && gain==0 && mH==1 && isGood && !muVeto && !wfDCBits && !isLNFill1 && !isLNFill2 && channel!=596 && channel!=676 && channel!=676 && channel!=612 && channel!=1104 && channel!=1200 && channel!=1334 && channel!=1336 && trapETailMin < 0";


void LargeFileCuts(int dsNum)
{
  cout << "DS-" << dsNum << endl;
  string inFile = Form("~/project/skim-files/skimDS%i*",dsNum);
  string outFile = Form("./data/skimDS%i_pt8_basicCutSpec.root",dsNum);
  double bins=299, xlo=100, xhi=3000;

  TFile *f = new TFile(outFile.c_str(),"RECREATE");
  TChain *skimTree = new TChain("skimTree");
  skimTree->Add(inFile.c_str());
  cout << "Found " << skimTree->GetEntries() << " entries.\n";

  // 1 basic cut - raw HG energy spectrum
  TH1D *h1 = new TH1D("basicCut","h1",bins,xlo,xhi);
  sprintf(theCut,"%s",basicCut);
  int cts = skimTree->Project("basicCut","trapENFCal",theCut);
  h1->Write();

  // 2 trapETailMin cut - significantly cut down on file size for DS-1.
  TH1D *h2 = new TH1D("basicTETMCut","h2",bins,xlo,xhi);
  sprintf(theCut,"%s && trapETailMin < 0",basicCut);
  int cts2 = skimTree->Project("basicTETMCut","trapENFCal",theCut);
  h2->Write();

  cout << "basic cut: " << cts << "  basic+TETM cut: " << cts2 << endl;
  f->Close();
}

void TCutSkimmer(int dsNum)
{
  cout << "DS-" << dsNum << endl;
  string inFile = Form("~/project/skim-files/skimDS%i*",dsNum);
  string outFile = Form("./data/skimDS%i_pt8_basicCuts.root",dsNum);
  sprintf(theCut,"%s && trapETailMin < 0",basicCut);
  cout << "Skimming file " << inFile << " using this cut: " << theCut << "\n\n";
  TChain *skim = new TChain("skimTree");
  skim->Add(inFile.c_str());
  skim->Draw(">>entList",theCut,"entrylist GOFF");
  TEntryList *entList = (TEntryList*)gDirectory->Get("entList");
  skim->SetEntryList(entList);
  TFile *f2 = new TFile(outFile.c_str(),"recreate");
  TTree *small = skim->CopyTree("");
  small->Write();
  cout << "Wrote " << small->GetEntries() << " entries.\n";
  TNamed thisCut("cutUsedHere",theCut);	// save the cut used into the file.
  thisCut.Write();
  f2->Close();
}

void GetExposure(int dsNum, double &enrExpo, double &natExpo)
{
  // NOTE: This is CALCULATED from the code in ds_livetime.cc,
  // NOT using the official numbers !
  if (dsNum == 0)      { enrExpo = 508.6206;  natExpo = 186.0343; }
  else if (dsNum == 1) { enrExpo = 679.9394;  natExpo = 67.3326;  }
  else if (dsNum == 2) { enrExpo = 114.8518;  natExpo = 10.8183;  }
  else if (dsNum == 3) { enrExpo = 377.6657;  natExpo = 83.1516;  }
  else if (dsNum == 4) { enrExpo = 128.9845;  natExpo = 93.1219;  }
  else if (dsNum == 5) { enrExpo = 1513.5387; natExpo = 603.2865; }
}

void GenerateSpectra(int dsNum)
{
  string inFile = Form("./data/skimDS%i_pt8_basicCuts.root",dsNum);
  double bins=290, xlo=100, xhi=3000;
  double enrExp=0, natExp=0;
  GetExposure(dsNum,enrExp,natExp);

  TCanvas *c1 = new TCanvas("c1","Bob Ross's Canvas",800,600);
  TFile *f1 = new TFile(inFile.c_str());
  TTree *skim = (TTree*)f1->Get("skimTree");

  // ======= "Raw" (basic cut) Spectrum  =======
  sprintf(theCut,"%s && trapETailMin < 0 && isEnr",basicCut);

  TH1D *h1 = new TH1D("h1","h1",bins,xlo,xhi);
  skim->Project("h1","trapENFCal",theCut);
  h1->GetXaxis()->SetTitle("Energy (trapENFCal)");
  h1->GetYaxis()->SetTitle(Form("Counts/%.1fkev-kg-days",(xhi-xlo)/bins));
  h1->GetYaxis()->SetTitleOffset(1.2);
  h1->SetLineColor(kBlue);
  h1->Scale(1/enrExp);
  h1->SetMaximum(0.25);  // 0.25 for 10 kev bins, 0.12 for 2 kev bins
  h1->Draw("hist");
  c1->Print(Form("./plots/ds%ienr_basicCut.pdf",dsNum));

  // ======= R + AvE =======
  sprintf(theCut,"%s && trapETailMin < 0 && avse>-1 && isEnr",basicCut);

  TH1D *h2 = new TH1D("h2","h2",bins,xlo,xhi);
  skim->Project("h2","trapENFCal",theCut);
  h2->GetXaxis()->SetTitle("Energy (trapENFCal)");
  h2->GetYaxis()->SetTitle(Form("Counts/%.1fkev-kg-days",(xhi-xlo)/bins));
  h2->GetYaxis()->SetTitleOffset(1.2);
  h2->SetLineColor(kRed);
  h2->Scale(1/enrExp);

  h1->Draw("hist");
  h2->Draw("hist same");

  TLegend* l1 = new TLegend(0.5,0.65,0.87,0.92);
  l1->AddEntry(h1,"basic cut","l");
  l1->AddEntry(h2,"basic + avse","l");
  l1->Draw("SAME");
  c1->Update();

  c1->Print(Form("./plots/ds%ienr_avse.pdf",dsNum));

  // ======= R + DCR =======
  if (dsNum !=2 && dsNum !=5)
  {
    sprintf(theCut,"%s && trapETailMin < 0 && dcrctc90 < 0 && isEnr",basicCut);

    TH1D *h3 = new TH1D("h3","h3",bins,xlo,xhi);
    skim->Project("h3","trapENFCal",theCut);
    h3->GetXaxis()->SetTitle("Energy (trapENFCal)");
    h3->GetYaxis()->SetTitle(Form("Counts/%.1fkev-kg-days",(xhi-xlo)/bins));
    h3->SetLineColor(kRed);
    h3->Scale(1/enrExp);

    h1->SetMaximum(0.25);  // 0.1 for 2 kev bins
    h1->Draw("hist");
    h3->Draw("hist same");

    TLegend* l2 = new TLegend(0.5,0.65,0.87,0.92);
    l2->AddEntry(h1,"basic cut","l");
    l2->AddEntry(h3,"basic + dcr","l");
    l2->Draw("SAME");
    c1->Update();

    c1->Print(Form("./plots/ds%ienr_dcr.pdf",dsNum));
  }

  // ======= R + AvE + DCR =======

  if (dsNum !=2 && dsNum !=5)
  {
    sprintf(theCut,"%s && trapETailMin < 0 && dcrctc90 < 0 && avse>-1 && isEnr",basicCut);

    TH1D *h4 = new TH1D("h4","h4",bins,xlo,xhi);
    skim->Project("h4","trapENFCal",theCut);
    h4->GetXaxis()->SetTitle("Energy (trapENFCal)");
    h4->GetYaxis()->SetTitle(Form("Counts/%.1fkev-kg-days",(xhi-xlo)/bins));
    h4->GetYaxis()->SetTitleOffset(1.2);
    h4->SetLineColor(kRed);
    h4->Scale(1/enrExp);

    h1->SetMaximum(0.25);
    h1->Draw("hist");
    h4->Draw("hist same");

    TLegend* l3 = new TLegend(0.35,0.65,0.87,0.92);
    l3->AddEntry(h1,"basic cut","l");
    l3->AddEntry(h4,"basic + avse + dcr","l");
    l3->Draw("SAME");
    c1->Update();

    c1->Print(Form("./plots/ds%ienr_avse_dcr.pdf",dsNum));
  }

  // ======= Enriched vs. natural, normalized =======

  TH1D *h5 = new TH1D("h5","h5",bins,xlo,xhi);
  sprintf(theCut,"%s && !isEnr && avse>-1",basicCut);
  skim->Project("h5","trapENFCal",theCut);
  h5->GetXaxis()->SetTitle("Energy (trapENFCal)");
  h5->GetYaxis()->SetTitle(Form("Counts/%.1fkev-kg-days",(xhi-xlo)/bins));
  h5->GetYaxis()->SetTitleOffset(1.2);
  h5->SetLineColor(kBlue);
  h5->Scale(1/natExp);
  h5->Draw("hist");

  TH1D *h6 = new TH1D("h6","h6",bins,xlo,xhi);
  sprintf(theCut,"%s && isEnr && avse>-1",basicCut);
  skim->Project("h6","trapENFCal",theCut);
  h6->Scale(1/enrExp);
  h6->SetLineColor(kRed);
  h6->Draw("hist same");

  TLegend* l4 = new TLegend(0.3,0.7,0.87,0.92);
  l4->AddEntry(h5,Form("DS-%i Natural+avse: %.4f kg-d",dsNum,natExp),"l");
  l4->AddEntry(h6,Form("DS-%i Enriched+avse: %.4f kg-d",dsNum,enrExp),"l");
  l4->Draw("SAME");
  c1->Update();

  c1->Print(Form("./plots/ds%i_enrVsNat_avse.pdf",dsNum));

}

void CombineSpectra()
{
  double bins=290, xlo=100, xhi=3000;
  TCanvas *c1 = new TCanvas("c1","Bob Ross's Canvas",800,600);

  double ds34enrExp=0, ds34natExp=0;
  TH1D *h1 = new TH1D("h1","h1",bins,xlo,xhi);  // === R + DCR, Enr ===
  TH1D *h2 = new TH1D("h2","h2",bins,xlo,xhi);  // === R + DCR, Nat ===
  TH1D *h3 = new TH1D("h3","h3",bins,xlo,xhi);  // === R + DCR + AvE, Enr ===
  TH1D *h4 = new TH1D("h4","h4",bins,xlo,xhi);  // === R + DCR + AvE, Nat ===
  vector<int> ds = {3,4};
  for (auto dsNum : ds)
  {
    cout << "drawing ds " << dsNum << endl;

    string inFile = Form("./data/skimDS%i_pt8_basicCuts.root",dsNum);
    TFile *f1 = new TFile(inFile.c_str());
    TTree *skim = (TTree*)f1->Get("skimTree");
    double enrExp=0, natExp=0;
    GetExposure(dsNum,enrExp,natExp);
    ds34enrExp+=enrExp;
    ds34natExp+=natExp;

    sprintf(theCut,"%s && trapETailMin < 0 && dcrctc90 < 0 && isEnr",basicCut);
    TH1D *h5 = new TH1D("h5","h5",bins,xlo,xhi);
    skim->Project("h5","trapENFCal",theCut);
    h5->Scale(1/enrExp);
    h1->Add(h5);

    sprintf(theCut,"%s && trapETailMin < 0 && dcrctc90 < 0 && !isEnr",basicCut);
    TH1D *h6 = new TH1D("h6","h6",bins,xlo,xhi);
    skim->Project("h6","trapENFCal",theCut);
    h6->Scale(1/natExp);
    h2->Add(h6);

    sprintf(theCut,"%s && trapETailMin < 0 && dcrctc90 < 0 && avse>-1 && isEnr",basicCut);
    TH1D *h7 = new TH1D("h7","h7",bins,xlo,xhi);
    skim->Project("h7","trapENFCal",theCut);
    h7->Scale(1/enrExp);
    h3->Add(h7);

    sprintf(theCut,"%s && trapETailMin < 0 && dcrctc90 < 0 && avse>-1 && !isEnr",basicCut);
    TH1D *h8 = new TH1D("h8","h8",bins,xlo,xhi);
    skim->Project("h8","trapENFCal",theCut);
    h8->Scale(1/natExp);
    h4->Add(h8);

    delete h5;
    delete h6;
    delete h7;
    delete h8;
    delete f1;
  }

  // === R + DCR, Enr ===
  h1->SetLineColor(kBlue);
  h1->GetXaxis()->SetTitle("Energy (trapENFCal)");
  h1->GetYaxis()->SetTitle(Form("Counts/%.1fkev-kg-days",(xhi-xlo)/bins));
  h1->GetYaxis()->SetTitleOffset(1.2);
  h1->Draw("hist");
  TLegend* l1 = new TLegend(0.4,0.7,0.87,0.92);
  l1->AddEntry(h1,Form("DS3,4 enr dcr: %.4f kg-d",ds34enrExp),"l");
  l1->Draw("same");
  c1->Update();
  c1->Print("./plots/sum_enr_dcr.pdf");

  // === R + DCR, Nat ===
  h2->SetLineColor(kBlue);
  h2->GetXaxis()->SetTitle("Energy (trapENFCal)");
  h2->GetYaxis()->SetTitle(Form("Counts/%.1fkev-kg-days",(xhi-xlo)/bins));
  h2->GetYaxis()->SetTitleOffset(1.2);
  h2->Draw("hist");
  TLegend* l2 = new TLegend(0.4,0.7,0.87,0.92);
  l2->AddEntry(h2,Form("DS3,4 nat dcr: %.4f kg-d",ds34natExp),"l");
  l2->Draw("same");
  c1->Update();
  c1->Print("./plots/sum_nat_dcr.pdf");

  // === R + DCR + AvE, Enr ===
  h3->SetLineColor(kBlue);
  h3->GetXaxis()->SetTitle("Energy (trapENFCal)");
  h3->GetYaxis()->SetTitle(Form("Counts/%.1fkev-kg-days",(xhi-xlo)/bins));
  h3->GetYaxis()->SetTitleOffset(1.2);
  h3->Draw("hist");
  TLegend* l3 = new TLegend(0.4,0.7,0.87,0.92);
  l3->AddEntry(h3,Form("DS3,4 enr dcr+avse: %.4f kg-d",ds34enrExp),"l");
  l3->Draw("same");
  c1->Update();
  c1->Print("./plots/sum_enr_dcr_avse.pdf");

  // === R + DCR + AvE, Nat ===
  h4->SetLineColor(kBlue);
  h4->GetXaxis()->SetTitle("Energy (trapENFCal)");
  h4->GetYaxis()->SetTitle(Form("Counts/%.1fkev-kg-days",(xhi-xlo)/bins));
  h4->GetYaxis()->SetTitleOffset(1.2);
  h4->Draw("hist");
  TLegend* l4 = new TLegend(0.4,0.7,0.87,0.92);
  l4->AddEntry(h4,Form("DS3,4 nat dcr+avse: %.4f kg-d",ds34natExp),"l");
  l4->Draw("same");
  c1->Update();
  c1->Print("./plots/sum_nat_dcr_avse.pdf");

  // ================= R + AvE (0vBB) ==================
  double ds1234enrExp=0, ds1234natExp=0;
  TH1D *h9 = new TH1D("h9","h9",bins,xlo,xhi);    // === R + AvE, Nat ===
  TH1D *h10 = new TH1D("h10","h10",bins,xlo,xhi); // === R + AvE, Enr ===
  vector<int> ds_a = {1,2,3,4};
  for (auto dsNum : ds_a)
  {
    cout << "drawing ds " << dsNum << endl;

    string inFile = Form("./data/skimDS%i_pt8_basicCuts.root",dsNum);
    TFile *f1 = new TFile(inFile.c_str());
    TTree *skim = (TTree*)f1->Get("skimTree");
    double enrExp=0, natExp=0;
    GetExposure(dsNum,enrExp,natExp);
    ds1234enrExp+=enrExp;
    ds1234natExp+=natExp;

    sprintf(theCut,"%s && trapETailMin < 0 && avse>-1 && !isEnr",basicCut);
    TH1D *h11 = new TH1D("h11","h11",bins,xlo,xhi);
    skim->Project("h11","trapENFCal",theCut);
    h11->Scale(1/natExp);
    h9->Add(h11);

    sprintf(theCut,"%s && trapETailMin < 0 && avse>-1 && isEnr",basicCut);
    TH1D *h12 = new TH1D("h12","h12",bins,xlo,xhi);
    skim->Project("h12","trapENFCal",theCut);
    h12->Scale(1/enrExp);
    h10->Add(h12);

    delete h11;
    delete h12;
    delete f1;
  }

  // === R + AvE, Nat ===
  h9->SetLineColor(kBlue);
  h9->GetXaxis()->SetTitle("Energy (trapENFCal)");
  h9->GetYaxis()->SetTitle(Form("Counts/%.1fkev-kg-days",(xhi-xlo)/bins));
  h9->GetYaxis()->SetTitleOffset(1.2);
  h9->Draw("hist");
  TLegend* l5 = new TLegend(0.4,0.7,0.87,0.92);
  l5->AddEntry(h9,Form("DS1-4 nat avse: %.4f kg-d",ds1234natExp),"l");
  l5->Draw("same");
  c1->Update();
  c1->Print("./plots/sum_nat_avse.pdf");

  // === R + AvE, Enr ===
  h10->SetLineColor(kBlue);
  h10->GetXaxis()->SetTitle("Energy (trapENFCal)");
  h10->GetYaxis()->SetTitle(Form("Counts/%.1fkev-kg-days",(xhi-xlo)/bins));
  h10->GetYaxis()->SetTitleOffset(1.2);
  h10->Draw("hist");
  TLegend* l6 = new TLegend(0.4,0.7,0.87,0.92);
  l6->AddEntry(h10,Form("DS1-4 enr avse: %.4f kg-d",ds1234enrExp),"l");
  l6->Draw("same");
  c1->Update();
  c1->Print("./plots/sum_enr_avse.pdf");

}

void ThresholdCut(int dsNum)
{
  string inFile = Form("./data/skimDS%i_pt8_basicCuts.root",dsNum);
  TFile *f1 = new TFile(inFile.c_str());
  TTree *skim = (TTree*)f1->Get("skimTree");
  double bins=100, xlo=0, xhi=20;
  double enrExp=0, natExp=0;
  GetExposure(dsNum,enrExp,natExp);

  TCanvas *c1 = new TCanvas("c1","Bob Ross's Canvas",800,600);

  // skim->Draw("threshkeV:run>>h3",ds1standardCut,"GOFF");
  // TH2D *h3 = (TH2D*)gDirectory->Get("h3");
  // h3->Draw("COLZ");

  // m1 channel boundaries: 575 - 700, m2 channel boundaries:
  // keep a list of run boundaries
  map<int,pair<int,int>> runMap = { {0,make_pair(2580,6963)}, {1,make_pair(9422,14387)}, {2,make_pair(14775,15803)}, {3,make_pair(16797,17980)}, {4,make_pair(60000802,60001888)}, {5,make_pair(18623,22142)} };

  double firstRun = runMap[dsNum].first;
  double lastRun = runMap[dsNum].second;
  int runBins = lastRun-firstRun;

  TH2D *h3 = new TH2D("h3","h3",280,0,70,runBins,firstRun,lastRun);
  skim->Project("h3","run:threshkeV",ds1standardCut);
  h3->Draw("");
  c1->Print(Form("./plots/ds%iLow_thresh.pdf",dsNum));
}

// Develop the "standard cut" for each dataset:
// Basic Cut, Burst Cut, Noisy Run Cut, Threshold Cut
// Will the wavelet denoising make the noisy run / burst cuts unnecessary?
void LowESpectra(int dsNum)
{
  string inFile = Form("./data/skimDS%i_pt8_basicCuts.root",dsNum);
  double bins=100, xlo=0, xhi=20;
  double enrExp=0, natExp=0;
  GetExposure(dsNum,enrExp,natExp);

  // ==== Standard DC Cut spectrum ====

  TCanvas *c1 = new TCanvas("c1","Bob Ross's Canvas",800,600);

  TFile *f1 = new TFile(inFile.c_str());
  TTree *skim = (TTree*)f1->Get("skimTree");

  sprintf(theCut,"%s && !isEnr",ds1standardCut);
  TH1D *h1 = new TH1D("h1","h1",bins,xlo,xhi);
  skim->Project("h1","trapENFCal",theCut);
  h1->GetXaxis()->SetTitle("Energy (trapENFCal)");
  h1->GetYaxis()->SetTitle(Form("Counts/%.1fkev-kg-days",(xhi-xlo)/bins));
  h1->GetYaxis()->SetTitleOffset(1.2);
  h1->SetLineColor(kBlue);
  h1->Scale(1/natExp);

  sprintf(theCut,"%s && isEnr",ds1standardCut);
  TH1D *h2 = new TH1D("h2","h2",bins,xlo,xhi);
  skim->Project("h2","trapENFCal",theCut);
  h2->SetLineColor(kRed);
  h2->Scale(1/enrExp);

  c1->SetLogy();
  h1->Draw("hist");
  h2->Draw("hist same");

  TLegend* l1 = new TLegend(0.4,0.7,0.87,0.92);
  l1->AddEntry(h1,Form("DS-%i Nat: %.4f kg-d",dsNum,natExp),"l");
  l1->AddEntry(h2,Form("DS-%i Enr: %.4f kg-d",dsNum,enrExp),"l");
  l1->Draw("same");
  c1->Update();

  c1->Print(Form("./plots/ds%iLow_standard.pdf",dsNum));


}

// Then finally generate your cleaned files you pull the waveforms for.
// But don't make too many low-energy cosmogenic plots here,
// because you don't have enough data cleaning power without the waveforms.