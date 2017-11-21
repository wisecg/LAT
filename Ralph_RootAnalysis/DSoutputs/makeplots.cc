#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TString.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TProof.h"
#include "TParameter.h"

#include <iostream>
#include <string>
#include <vector>
#include <cstdlib>
#include <fstream>
#include <ctime>

using namespace std;

void makeplots(){
	TSystemDirectory dir("./", "./");
  TList *files = dir.GetListOfFiles();
	files->Sort();
  TSystemFile *file;
  TString fname;
	TString cname;
	TString title;
  TIter next(files);
	TFile* f;

  TH1D* h_exposureDSNat;
  TH1D* h_energyBkgTot1Nat;
  TH1D* h_energyBkgTot2Nat;
  TH1D* h_energyBkgTot3Nat;
  TH1D* h_energyCalTot1Nat;
  TH1D* h_energyCalTot2Nat;
  TH1D* h_EffNat;
  TH1D* h_EffNat2;
  TH1D* h_exposureDSEnr;
  TH1D* h_energyBkgTot1Enr;
  TH1D* h_energyBkgTot2Enr;
  TH1D* h_energyBkgTot3Enr;
  TH1D* h_energyCalTot1Enr;
  TH1D* h_energyCalTot2Enr;
  TH1D* h_EffEnr;
  TH1D* h_EffEnr2;

	TH1D* h_totalExposureEnr = new TH1D("h_totalExposureEnr","h_totalExposureEnr",3000,0,300);
	TH1D* h_totalExposureNat = new TH1D("h_totalExposureNat","h_totalExposureNat",3000,0,300);
	TH1D* h_totalSpectrumEnr = new TH1D("h_totalSpectrumEnr","h_totalSpectrumEnr",3000,0,300);
	TH1D* h_totalSpectrumNat = new TH1D("h_totalSpectrumNat","h_totalSpectrumNat",3000,0,300);

	int dataset;
	int counter = 0;
	double value1;
	double value2;

////////////////////////////////////////////////////////////////////
	while ((file=(TSystemFile*)next())) {
   	fname = file->GetName();
    if (!file->IsDirectory() && fname.EndsWith(".root") && fname.BeginsWith("DS")) {
			cname = fname(2,1);
			dataset = cname.Atoi();
			cout << fname << " " << dataset << endl;
			f = new TFile(fname);
			//f->ls();			
			h_exposureDSNat = (TH1D*)f->Get("h_exposureDSNat");
			h_energyBkgTot1Nat = (TH1D*)f->Get("h_energyBkgTot1Nat");
			h_energyBkgTot2Nat = (TH1D*)f->Get("h_energyBkgTot2Nat");
			h_energyCalTot1Nat = (TH1D*)f->Get("h_energyCalTot1Nat");
			h_energyCalTot2Nat = (TH1D*)f->Get("h_energyCalTot2Nat");
			h_exposureDSEnr = (TH1D*)f->Get("h_exposureDSEnr");
			h_energyBkgTot1Enr = (TH1D*)f->Get("h_energyBkgTot1Enr");
			h_energyBkgTot2Enr = (TH1D*)f->Get("h_energyBkgTot2Enr");
			h_energyCalTot1Enr = (TH1D*)f->Get("h_energyCalTot1Enr");
			h_energyCalTot2Enr = (TH1D*)f->Get("h_energyCalTot2Enr");
			h_exposureDSNat->SetDirectory(0);
			h_energyBkgTot1Nat->SetDirectory(0);
			h_energyBkgTot2Nat->SetDirectory(0);
			h_energyCalTot1Nat->SetDirectory(0);
			h_energyCalTot2Nat->SetDirectory(0);
			h_exposureDSEnr->SetDirectory(0);
			h_energyBkgTot1Enr->SetDirectory(0);
			h_energyBkgTot2Enr->SetDirectory(0);
			h_energyCalTot1Enr->SetDirectory(0);
			h_energyCalTot2Enr->SetDirectory(0);

			h_EffNat = new TH1D("h_EffNat","h_EffNat",3000,0,300);
			h_EffNat2 = new TH1D("h_EffNat2","h_EffNat2",300,0,300);
			h_EffNat->SetDirectory(0);
			h_EffNat2->SetDirectory(0);
			h_EffEnr = new TH1D("h_EffEnr","h_EffEnr",3000,0,300);
			h_EffEnr2 = new TH1D("h_EffEnr2","h_EffEnr2",300,0,300);
			h_EffEnr->SetDirectory(0);
			h_EffEnr2->SetDirectory(0);
			h_energyBkgTot3Nat = new TH1D("h_energyBkgTot3Nat","h_energyBkgTot3Nat",3000,0,300);
			h_energyBkgTot3Enr = new TH1D("h_energyBkgTot3Enr","h_energyBkgTot3Enr",3000,0,300);
			h_energyBkgTot3Nat->SetDirectory(0);
			h_energyBkgTot3Enr->SetDirectory(0);


			for (int i = 0;i<3000;i++){
				//efficiency
				value1 = h_energyCalTot1Nat->GetBinContent(i);
				value2 = h_energyCalTot2Nat->GetBinContent(i);
				if (value1 > 0) {
					h_EffNat->SetBinContent(i,value2/value1);
					h_EffNat->SetBinError(i,sqrt(value2)/value1+sqrt(value1)*value2/value1/value1);
					h_EffNat2->Fill(i/10.,value2/value1*0.1);
				}
				value1 = h_energyCalTot1Enr->GetBinContent(i);
				value2 = h_energyCalTot2Enr->GetBinContent(i);
				if (value1 > 0){
					h_EffEnr->SetBinContent(i,value2/value1);
					h_EffEnr->SetBinError(i,sqrt(value2)/value1+sqrt(value1)*value2/value1/value1);	
					h_EffEnr2->Fill(i/10.,value2/value1*0.1);
				}
				//corrected spectra
				value1 = 0;
				value2 = 0;
				if ((h_EffNat->GetBinContent(i)>0.6)&&(h_exposureDSNat->GetBinContent(i)>0)){
					value1 = h_energyBkgTot2Nat->GetBinContent(i)/h_EffNat->GetBinContent(i)/h_exposureDSNat->GetBinContent(i);
					value2 = value1*( sqrt(h_energyBkgTot2Nat->GetBinContent(i)) / h_energyBkgTot2Nat->GetBinContent(i)    //error in counts
													 	+ h_EffNat->GetBinError(i) / h_EffNat->GetBinContent(i)															// error in eff
														+ 0.02);																																						// error in exposure
				}
				h_energyBkgTot3Nat->SetBinContent(i,value1);
				h_energyBkgTot3Nat->SetBinError(i,value2);
				value1 = 0;
				value2 = 0;
				if ((h_EffEnr->GetBinContent(i)>0.6)&&(h_exposureDSEnr->GetBinContent(i)>0)){
					value1 = h_energyBkgTot2Enr->GetBinContent(i)/h_EffEnr->GetBinContent(i)/h_exposureDSEnr->GetBinContent(i);
					value2 = value1*( sqrt(h_energyBkgTot2Enr->GetBinContent(i)) / h_energyBkgTot2Enr->GetBinContent(i)    //error in counts
													 	+ h_EffEnr->GetBinError(i) / h_EffEnr->GetBinContent(i)															// error in eff
														+ 0.02);																																						// error in exposure
				}
				h_energyBkgTot3Enr->SetBinContent(i,value1);
				h_energyBkgTot3Enr->SetBinError(i,value2);
			}

			if (dataset>=0){
				h_totalExposureEnr->Add(h_exposureDSEnr);
				h_totalExposureNat->Add(h_exposureDSNat);
				h_totalSpectrumNat->Add(h_energyBkgTot3Nat);
				h_totalSpectrumEnr->Add(h_energyBkgTot3Enr);
				counter++;
			}
/*
/////////////////////////////////////////////////////////
			cname = fname(0,fname.First(".")-4);
			cname += "_BeforeAfter_Nat";

			title = fname(0,fname.First("_"));
			title += " ";
			title += fname(fname.First("_")+1,2);
			title += " background in natural detectors";

			TCanvas* c1 = new TCanvas(cname,cname,800,600);
			c1->SetTopMargin(0.1);
		  c1->SetLeftMargin(0.15);
		  c1->SetRightMargin(0.15);
		  c1->SetBottomMargin(0.15);
		  c1->SetFrameLineWidth(3);

			gPad->SetLogy();
			h_energyBkgTot1Nat->SetTitleSize(0.04);
			h_energyBkgTot1Nat->SetTitle(title);

			h_energyBkgTot2Nat->GetXaxis()->SetRangeUser(0,40);
			h_energyBkgTot1Nat->GetXaxis()->SetRangeUser(0,40);
			h_energyBkgTot1Nat->GetXaxis()->SetTitle("energy (keV)");
			h_energyBkgTot1Nat->GetXaxis()->SetTitleOffset(1);
			h_energyBkgTot1Nat->GetXaxis()->SetTitleSize(0.04);
			h_energyBkgTot1Nat->GetXaxis()->SetNoExponent(kTRUE);
			h_energyBkgTot1Nat->GetXaxis()->SetLabelOffset(0.005);
			h_energyBkgTot1Nat->GetXaxis()->SetLabelSize(0.04);
			h_energyBkgTot1Nat->GetXaxis()->CenterTitle(1);

			h_energyBkgTot1Nat->GetYaxis()->SetRangeUser(0.11,9E6);
			h_energyBkgTot1Nat->GetYaxis()->SetNdivisions(404);
			h_energyBkgTot1Nat->GetYaxis()->SetTitle("counts (a.u.)");
			h_energyBkgTot1Nat->GetYaxis()->SetTitleOffset(1.4);
			h_energyBkgTot1Nat->GetYaxis()->SetTitleSize(0.04);
			h_energyBkgTot1Nat->GetYaxis()->SetLabelSize(0.04);
			h_energyBkgTot1Nat->GetYaxis()->CenterTitle(1);
			h_energyBkgTot1Nat->SetStats(0);


			h_energyBkgTot1Nat->SetLineColor(kBlack);
			h_energyBkgTot1Nat->SetLineStyle(1);						
			h_energyBkgTot1Nat->SetLineWidth(2);
			h_energyBkgTot1Nat->Draw("Hist");
			h_energyBkgTot2Nat->SetLineColor(kRed);
			h_energyBkgTot2Nat->SetLineStyle(1);						
			h_energyBkgTot2Nat->SetLineWidth(2);
			h_energyBkgTot2Nat->Draw("Hist SAME");

			TLegend* leg_c1 = new TLegend(0.40,0.55,0.70,0.70);
		  leg_c1->SetLineWidth(0);
			leg_c1->AddEntry(h_energyBkgTot1Nat,"before cuts","l");
			leg_c1->AddEntry(h_energyBkgTot2Nat,"after RMS cut","l");
			leg_c1->Draw();
			
			c1->Update();
			cname+=".png";
			c1->Print(cname);	
			
/////////////////////////////////////////////////////////
			cname = fname(0,fname.First(".")-4);
			cname += "_BeforeAfter_Enr";

			title = fname(0,fname.First("_"));
			title += " ";
			title += fname(fname.First("_")+1,2);
			title += " background in enriched detectors";

			TCanvas* c2 = new TCanvas(cname,cname,800,600);
			c2->SetTopMargin(0.1);
		  c2->SetLeftMargin(0.15);
		  c2->SetRightMargin(0.15);
		  c2->SetBottomMargin(0.15);
		  c2->SetFrameLineWidth(3);

			gPad->SetLogy();
			h_energyBkgTot1Enr->SetTitleSize(0.04);
			h_energyBkgTot1Enr->SetTitle(title);

			h_energyBkgTot2Enr->GetXaxis()->SetRangeUser(0,40);
			h_energyBkgTot1Enr->GetXaxis()->SetRangeUser(0,40);
			h_energyBkgTot1Enr->GetXaxis()->SetTitle("energy (keV)");
			h_energyBkgTot1Enr->GetXaxis()->SetTitleOffset(1);
			h_energyBkgTot1Enr->GetXaxis()->SetTitleSize(0.04);
			h_energyBkgTot1Enr->GetXaxis()->SetNoExponent(kTRUE);
			h_energyBkgTot1Enr->GetXaxis()->SetLabelOffset(0.005);
			h_energyBkgTot1Enr->GetXaxis()->SetLabelSize(0.04);
			h_energyBkgTot1Enr->GetXaxis()->CenterTitle(1);

			h_energyBkgTot1Enr->GetYaxis()->SetRangeUser(0.11,9E6);
			h_energyBkgTot1Enr->GetYaxis()->SetNdivisions(404);
			h_energyBkgTot1Enr->GetYaxis()->SetTitle("counts (a.u.)");
			h_energyBkgTot1Enr->GetYaxis()->SetTitleOffset(1.4);
			h_energyBkgTot1Enr->GetYaxis()->SetTitleSize(0.04);
			h_energyBkgTot1Enr->GetYaxis()->SetLabelSize(0.04);
			h_energyBkgTot1Enr->GetYaxis()->CenterTitle(1);
			h_energyBkgTot1Enr->SetStats(0);


			h_energyBkgTot1Enr->SetLineColor(kBlack);
			h_energyBkgTot1Enr->SetLineStyle(1);						
			h_energyBkgTot1Enr->SetLineWidth(2);
			h_energyBkgTot1Enr->Draw("Hist");
			h_energyBkgTot2Enr->SetLineColor(kRed);
			h_energyBkgTot2Enr->SetLineStyle(1);						
			h_energyBkgTot2Enr->SetLineWidth(2);
			h_energyBkgTot2Enr->Draw("Hist SAME");

			TLegend* leg_c2 = new TLegend(0.40,0.55,0.70,0.70);
		  leg_c2->SetLineWidth(0);
			leg_c2->AddEntry(h_energyBkgTot1Enr,"before cuts","l");
			leg_c2->AddEntry(h_energyBkgTot2Enr,"after RMS cut","l");
			leg_c2->Draw();
			
			c2->Update();
			cname+=".png";
			c2->Print(cname);	
			
/////////////////////////////////////////////////////////

			cname = fname(0,fname.First(".")-4);
			cname += "_Efficiency_Nat";

			title = fname(0,fname.First("_"));
			title += " ";
			title += fname(fname.First("_")+1,2);
			title += " Efficiency in natural detectors";

			TCanvas* c3 = new TCanvas(cname,cname,800,600);
			c3->SetTopMargin(0.1);
		  c3->SetLeftMargin(0.15);
		  c3->SetRightMargin(0.15);
		  c3->SetBottomMargin(0.15);
		  c3->SetFrameLineWidth(3);

			h_EffNat->SetTitleSize(0.04);
			h_EffNat->SetTitle(title);

			h_EffNat->GetXaxis()->SetRangeUser(0,40);
			h_EffNat->GetXaxis()->SetTitle("energy (keV)");
			h_EffNat->GetXaxis()->SetTitleOffset(1);
			h_EffNat->GetXaxis()->SetTitleSize(0.04);
			h_EffNat->GetXaxis()->SetNoExponent(kTRUE);
			h_EffNat->GetXaxis()->SetLabelOffset(0.005);
			h_EffNat->GetXaxis()->SetLabelSize(0.04);
			h_EffNat->GetXaxis()->CenterTitle(1);

			h_EffNat->GetYaxis()->SetRangeUser(0,2);
			h_EffNat->GetYaxis()->SetNdivisions(404);
			h_EffNat->GetYaxis()->SetTitle("efficiency");
			h_EffNat->GetYaxis()->SetTitleOffset(1.4);
			h_EffNat->GetYaxis()->SetTitleSize(0.04);
			h_EffNat->GetYaxis()->SetLabelSize(0.04);
			h_EffNat->GetYaxis()->CenterTitle(1);
			h_EffNat->SetStats(0);


			h_EffNat->SetLineColor(kGray+1);
			h_EffNat->SetLineStyle(1);						
			h_EffNat->SetLineWidth(2);
			h_EffNat->Draw("");
			h_EffNat2->SetLineColor(kBlack);
			h_EffNat2->SetLineStyle(2);						
			h_EffNat2->SetLineWidth(3);
			h_EffNat2->Draw("HIST SAME");

			TLegend* leg_c3 = new TLegend(0.55,0.70,0.83,0.83);
		  leg_c3->SetLineWidth(0);
			leg_c3->AddEntry(h_EffNat,"0.1 keV steps","l");
			leg_c3->AddEntry(h_EffNat2,"1 keV averaged","l");
			leg_c3->Draw();
			
			c3->Update();
			cname+=".png";
			c3->Print(cname);	


/////////////////////////////////////////////////////////
			cname = fname(0,fname.First(".")-4);
			cname += "_Efficiency_Enr";

			title = fname(0,fname.First("_"));
			title += " ";
			title += fname(fname.First("_")+1,2);
			title += " Efficiency in enriched detectors";

			TCanvas* c4 = new TCanvas(cname,cname,800,600);
			c4->SetTopMargin(0.1);
		  c4->SetLeftMargin(0.15);
		  c4->SetRightMargin(0.15);
		  c4->SetBottomMargin(0.15);
		  c4->SetFrameLineWidth(3);

			h_EffEnr->SetTitleSize(0.04);
			h_EffEnr->SetTitle(title);

			h_EffEnr->GetXaxis()->SetRangeUser(0,40);
			h_EffEnr->GetXaxis()->SetTitle("energy (keV)");
			h_EffEnr->GetXaxis()->SetTitleOffset(1);
			h_EffEnr->GetXaxis()->SetTitleSize(0.04);
			h_EffEnr->GetXaxis()->SetNoExponent(kTRUE);
			h_EffEnr->GetXaxis()->SetLabelOffset(0.005);
			h_EffEnr->GetXaxis()->SetLabelSize(0.04);
			h_EffEnr->GetXaxis()->CenterTitle(1);

			h_EffEnr->GetYaxis()->SetRangeUser(0,2);
			h_EffEnr->GetYaxis()->SetNdivisions(404);
			h_EffEnr->GetYaxis()->SetTitle("efficiency");
			h_EffEnr->GetYaxis()->SetTitleOffset(1.4);
			h_EffEnr->GetYaxis()->SetTitleSize(0.04);
			h_EffEnr->GetYaxis()->SetLabelSize(0.04);
			h_EffEnr->GetYaxis()->CenterTitle(1);
			h_EffEnr->SetStats(0);


			h_EffEnr->SetLineColor(kGray+1);
			h_EffEnr->SetLineStyle(1);						
			h_EffEnr->SetLineWidth(2);
			h_EffEnr->Draw("");
			h_EffEnr2->SetLineColor(kBlack);
			h_EffEnr2->SetLineStyle(2);						
			h_EffEnr2->SetLineWidth(3);
			h_EffEnr2->Draw("HIST SAME");

			TLegend* leg_c4 = new TLegend(0.55,0.70,0.83,0.83);
		  leg_c4->SetLineWidth(0);
			leg_c4->AddEntry(h_EffEnr,"0.1 keV steps","l");
			leg_c4->AddEntry(h_EffEnr2,"1 keV averaged","l");
			leg_c4->Draw();
			
			c4->Update();
			cname+=".png";
			c4->Print(cname);	
/////////////////////////////////////////////////////////

			cname = fname(0,fname.First(".")-4);
			cname += "_Exposure";

			title = fname(0,fname.First("_"));
			title += " ";
			title += fname(fname.First("_")+1,2);
			title += " Exposure";

			TCanvas* c5 = new TCanvas(cname,cname,800,600);
			c5->SetTopMargin(0.1);
		  c5->SetLeftMargin(0.15);
		  c5->SetRightMargin(0.15);
		  c5->SetBottomMargin(0.15);
		  c5->SetFrameLineWidth(3);
			//gPad->SetLogx();
			h_exposureDSEnr->SetTitleSize(0.04);
			h_exposureDSEnr->SetTitle(title);

			h_exposureDSNat->GetXaxis()->SetRangeUser(0.1,20);
			h_exposureDSEnr->GetXaxis()->SetRangeUser(0.1,20);
			h_exposureDSEnr->GetXaxis()->SetTitle("energy (keV)");
			h_exposureDSEnr->GetXaxis()->SetTitleOffset(1);
			h_exposureDSEnr->GetXaxis()->SetTitleSize(0.04);
			h_exposureDSEnr->GetXaxis()->SetNoExponent(kTRUE);
			h_exposureDSEnr->GetXaxis()->SetLabelOffset(0.005);
			h_exposureDSEnr->GetXaxis()->SetLabelSize(0.04);
			h_exposureDSEnr->GetXaxis()->CenterTitle(1);

		//	h_exposureDSEnr->GetYaxis()->SetRangeUser(0,2);
			h_exposureDSEnr->GetYaxis()->SetNdivisions(406);
			h_exposureDSEnr->GetYaxis()->SetTitle("exposure (kg yr)");
			h_exposureDSEnr->GetYaxis()->SetTitleOffset(1.4);
			h_exposureDSEnr->GetYaxis()->SetTitleSize(0.04);
			h_exposureDSEnr->GetYaxis()->SetLabelSize(0.04);
			h_exposureDSEnr->GetYaxis()->CenterTitle(1);
			h_exposureDSEnr->SetStats(0);


			h_exposureDSEnr->SetLineColor(kBlack);
			h_exposureDSEnr->SetLineStyle(1);						
			h_exposureDSEnr->SetLineWidth(2);
			h_exposureDSEnr->Draw("Hist");
			h_exposureDSNat->SetLineColor(kRed);
			h_exposureDSNat->SetLineStyle(1);						
			h_exposureDSNat->SetLineWidth(2);
			h_exposureDSNat->Draw("Hist SAME");

			TLegend* leg_c5 = new TLegend(0.40,0.65,0.70,0.80);
		  leg_c5->SetLineWidth(0);
			leg_c5->AddEntry(h_exposureDSEnr,"enriched detectors","l");
			leg_c5->AddEntry(h_exposureDSNat,"natural detectors","l");
			leg_c5->Draw();
			
			c5->Update();
			cname+=".png";
			c5->Print(cname);	

/////////////////////////////////////////////////////////

			cname = fname(0,fname.First(".")-4);
			cname += "_Spectrum";

			title = fname(0,fname.First("_"));
			title += " ";
			title += fname(fname.First("_")+1,2);
			title += " Spectrum, exposure and efficiency corrected";

			TCanvas* c6 = new TCanvas(cname,cname,800,600);
			c6->SetTopMargin(0.1);
		  c6->SetLeftMargin(0.15);
		  c6->SetRightMargin(0.15);
		  c6->SetBottomMargin(0.15);
		  c6->SetFrameLineWidth(3);
			gPad->SetLogy();
			h_energyBkgTot3Nat->SetTitleSize(0.04);
			h_energyBkgTot3Nat->SetTitle(title);


			h_energyBkgTot3Enr->Rebin(2);
			h_energyBkgTot3Nat->Rebin(2);
			h_energyBkgTot3Enr->Scale(5/365.);
			h_energyBkgTot3Nat->Scale(5/365.);
			h_energyBkgTot3Enr->GetXaxis()->SetRangeUser(0,40);
			h_energyBkgTot3Nat->GetXaxis()->SetRangeUser(0,40);
			h_energyBkgTot3Nat->GetXaxis()->SetTitle("energy (keV)");
			h_energyBkgTot3Nat->GetXaxis()->SetTitleOffset(1);
			h_energyBkgTot3Nat->GetXaxis()->SetTitleSize(0.04);
			h_energyBkgTot3Nat->GetXaxis()->SetNoExponent(kTRUE);
			h_energyBkgTot3Nat->GetXaxis()->SetLabelOffset(0.005);
			h_energyBkgTot3Nat->GetXaxis()->SetLabelSize(0.04);
			h_energyBkgTot3Nat->GetXaxis()->CenterTitle(1);

			h_energyBkgTot3Nat->GetYaxis()->SetRangeUser(0.01,10);
			h_energyBkgTot3Nat->GetYaxis()->SetNdivisions(406);
			h_energyBkgTot3Nat->GetYaxis()->SetTitle("cts (kg^{-1} day^{-1} keV^{-1})");
			h_energyBkgTot3Nat->GetYaxis()->SetTitleOffset(1.4);
			h_energyBkgTot3Nat->GetYaxis()->SetTitleSize(0.04);
			h_energyBkgTot3Nat->GetYaxis()->SetLabelSize(0.04);
			h_energyBkgTot3Nat->GetYaxis()->CenterTitle(1);
			h_energyBkgTot3Nat->SetStats(0);

			h_energyBkgTot3Nat->SetLineColor(kBlack);
			h_energyBkgTot3Nat->SetLineStyle(1);						
			h_energyBkgTot3Nat->SetLineWidth(2);
			h_energyBkgTot3Nat->Draw("Hist");
			h_energyBkgTot3Enr->SetLineColor(kRed);
			h_energyBkgTot3Enr->SetLineStyle(1);						
			h_energyBkgTot3Enr->SetLineWidth(2);
			h_energyBkgTot3Enr->Draw("Hist SAME");

			TLegend* leg_c6 = new TLegend(0.40,0.65,0.80,0.80);
		  leg_c6->SetLineWidth(0);
			leg_c6->AddEntry(h_energyBkgTot3Enr,"enriched detectors, 0.2-keV binning","l");
			leg_c6->AddEntry(h_energyBkgTot3Nat,"natural detectors, 0.2-keV binning","l");
			leg_c6->Draw();
			
			c6->Update();
			cname+=".png";
			c6->Print(cname);	
			cname = fname(0,fname.First(".")-4);
			cname += "_Spectrum";
			cname+=".C";
			c6->Print(cname);	


/////////////////////////////////////////////////////////
*/
	
	
			f->Close();
			delete f;
		}
	}
	cout << counter  << endl;
	
/////////////////////////////////////////////////////////
	cname = "DSall_Spectrum";
	title = "DSall ";
	title += "Spectrum, exposure and efficiency corrected";

	TCanvas* c7 = new TCanvas(cname,cname,800,600);
	c7->SetTopMargin(0.1);
	c7->SetLeftMargin(0.15);
	c7->SetRightMargin(0.15);
	c7->SetBottomMargin(0.15);
	c7->SetFrameLineWidth(3);
	gPad->SetLogy();
	h_totalSpectrumNat->SetTitleSize(0.04);
	h_totalSpectrumNat->SetTitle(title);

	h_totalSpectrumEnr->Scale(10./365./counter);
	h_totalSpectrumNat->Scale(10./365./counter);
	h_totalSpectrumEnr->GetXaxis()->SetRangeUser(0,40);
	h_totalSpectrumNat->GetXaxis()->SetRangeUser(0,40);
	h_totalSpectrumNat->GetXaxis()->SetTitle("energy (keV)");
	h_totalSpectrumNat->GetXaxis()->SetTitleOffset(1);
	h_totalSpectrumNat->GetXaxis()->SetTitleSize(0.04);
	h_totalSpectrumNat->GetXaxis()->SetNoExponent(kTRUE);
	h_totalSpectrumNat->GetXaxis()->SetLabelOffset(0.005);
	h_totalSpectrumNat->GetXaxis()->SetLabelSize(0.04);
	h_totalSpectrumNat->GetXaxis()->CenterTitle(1);

	h_totalSpectrumNat->GetYaxis()->SetRangeUser(0.01,10);
	h_totalSpectrumNat->GetYaxis()->SetNdivisions(406);
	h_totalSpectrumNat->GetYaxis()->SetTitle("cts (kg^{-1} day^{-1} keV^{-1})");
	h_totalSpectrumNat->GetYaxis()->SetTitleOffset(1.4);
	h_totalSpectrumNat->GetYaxis()->SetTitleSize(0.04);
	h_totalSpectrumNat->GetYaxis()->SetLabelSize(0.04);
	h_totalSpectrumNat->GetYaxis()->CenterTitle(1);
	h_totalSpectrumNat->SetStats(0);

	h_totalSpectrumNat->SetLineColor(kBlack);
	h_totalSpectrumNat->SetLineStyle(1);				
	h_totalSpectrumNat->SetLineWidth(2);
	h_totalSpectrumNat->Draw("Hist");
	h_totalSpectrumEnr->SetLineColor(kRed);
	h_totalSpectrumEnr->SetLineStyle(1);				
	h_totalSpectrumEnr->SetLineWidth(2);
	h_totalSpectrumEnr->Draw("Hist SAME");

	TLegend* leg_c7 = new TLegend(0.40,0.65,0.80,0.80);
	leg_c7->SetLineWidth(0);
	leg_c7->AddEntry(h_totalSpectrumEnr,"enriched detectors, 0.1-keV binning","l");
	leg_c7->AddEntry(h_totalSpectrumNat,"natural detectors, 0.1-keV binning","l");
	leg_c7->Draw();
	
	c7->Update();
	cname+=".png";
	c7->Print(cname);

/////////////////////////////////////////////////////////
	cname = "DSall_Exposure";
	title = "DSall ";
	title += "Exposure";

	TCanvas* c8 = new TCanvas(cname,cname,800,600);
	c8->SetTopMargin(0.1);
	c8->SetLeftMargin(0.15);
	c8->SetRightMargin(0.15);
	c8->SetBottomMargin(0.15);
	c8->SetFrameLineWidth(3);
	h_totalExposureEnr->SetTitleSize(0.04);
	h_totalExposureEnr->SetTitle(title);

	h_totalExposureNat->GetXaxis()->SetRangeUser(0,20);
	h_totalExposureEnr->GetXaxis()->SetRangeUser(0,20);
	h_totalExposureEnr->GetXaxis()->SetTitle("energy (keV)");
	h_totalExposureEnr->GetXaxis()->SetTitleOffset(1);
	h_totalExposureEnr->GetXaxis()->SetTitleSize(0.04);
	h_totalExposureEnr->GetXaxis()->SetNoExponent(kTRUE);
	h_totalExposureEnr->GetXaxis()->SetLabelOffset(0.005);
	h_totalExposureEnr->GetXaxis()->SetLabelSize(0.04);
	h_totalExposureEnr->GetXaxis()->CenterTitle(1);

	h_totalExposureEnr->GetYaxis()->SetRangeUser(0.01,10);
	h_totalExposureEnr->GetYaxis()->SetNdivisions(406);
	h_totalExposureEnr->GetYaxis()->SetTitle("exposure (kg yr)");
	h_totalExposureEnr->GetYaxis()->SetTitleOffset(1.4);
	h_totalExposureEnr->GetYaxis()->SetTitleSize(0.04);
	h_totalExposureEnr->GetYaxis()->SetLabelSize(0.04);
	h_totalExposureEnr->GetYaxis()->CenterTitle(1);
	h_totalExposureEnr->SetStats(0);

	h_totalExposureEnr->SetLineColor(kBlack);
	h_totalExposureEnr->SetLineStyle(1);				
	h_totalExposureEnr->SetLineWidth(2);
	h_totalExposureEnr->Draw("Hist");
	h_totalExposureNat->SetLineColor(kRed);
	h_totalExposureNat->SetLineStyle(1);				
	h_totalExposureNat->SetLineWidth(2);
	h_totalExposureNat->Draw("Hist SAME");

	TLegend* leg_c8 = new TLegend(0.40,0.65,0.80,0.80);
	leg_c8->SetLineWidth(0);
	leg_c8->AddEntry(h_totalExposureEnr,"enriched detectors","l");
	leg_c8->AddEntry(h_totalExposureNat,"natural detectors","l");
	leg_c8->Draw();
	
	c8->Update();
	cname+=".png";
	c8->Print(cname);		
/////////////////////////////////////////////////////////
	cname = "DSall_Spectrum2";
	title = "DSall ";
	title += "Spectrum, exposure and efficiency corrected";

	TCanvas* c9 = new TCanvas(cname,cname,800,600);
	c9->SetTopMargin(0.1);
	c9->SetLeftMargin(0.15);
	c9->SetRightMargin(0.15);
	c9->SetBottomMargin(0.15);
	c9->SetFrameLineWidth(3);
	gPad->SetLogy();
	h_totalSpectrumNat->SetTitleSize(0.04);
	h_totalSpectrumNat->SetTitle(title);

	h_totalSpectrumEnr->Rebin(10);
	h_totalSpectrumNat->Rebin(10);
	h_totalSpectrumEnr->Scale(100./365./counter);
	h_totalSpectrumNat->Scale(100./365./counter);
	h_totalSpectrumEnr->GetXaxis()->SetRangeUser(0,250);
	h_totalSpectrumNat->GetXaxis()->SetRangeUser(0,250);
	h_totalSpectrumNat->GetXaxis()->SetTitle("energy (keV)");
	h_totalSpectrumNat->GetXaxis()->SetTitleOffset(1);
	h_totalSpectrumNat->GetXaxis()->SetTitleSize(0.04);
	h_totalSpectrumNat->GetXaxis()->SetNoExponent(kTRUE);
	h_totalSpectrumNat->GetXaxis()->SetLabelOffset(0.005);
	h_totalSpectrumNat->GetXaxis()->SetLabelSize(0.04);
	h_totalSpectrumNat->GetXaxis()->CenterTitle(1);

	h_totalSpectrumNat->GetYaxis()->SetRangeUser(0.001,1);
	h_totalSpectrumNat->GetYaxis()->SetNdivisions(406);
	h_totalSpectrumNat->GetYaxis()->SetTitle("cts (kg^{-1} day^{-1} keV^{-1})");
	h_totalSpectrumNat->GetYaxis()->SetTitleOffset(1.4);
	h_totalSpectrumNat->GetYaxis()->SetTitleSize(0.04);
	h_totalSpectrumNat->GetYaxis()->SetLabelSize(0.04);
	h_totalSpectrumNat->GetYaxis()->CenterTitle(1);
	h_totalSpectrumNat->SetStats(0);

	h_totalSpectrumNat->SetLineColor(kBlack);
	h_totalSpectrumNat->SetLineStyle(1);				
	h_totalSpectrumNat->SetLineWidth(2);
	h_totalSpectrumNat->Draw("Hist");
	h_totalSpectrumEnr->SetLineColor(kRed);
	h_totalSpectrumEnr->SetLineStyle(1);				
	h_totalSpectrumEnr->SetLineWidth(2);
	h_totalSpectrumEnr->Draw("Hist SAME");

	TLegend* leg_c9 = new TLegend(0.40,0.65,0.80,0.80);
	leg_c9->SetLineWidth(0);
	leg_c9->AddEntry(h_totalSpectrumEnr,"enriched detectors, 1-keV binning","l");
	leg_c9->AddEntry(h_totalSpectrumNat,"natural detectors, 1-keV binning","l");
	leg_c9->Draw();
	
	c9->Update();
	cname+=".png";
	c9->Print(cname);



	return;



}
