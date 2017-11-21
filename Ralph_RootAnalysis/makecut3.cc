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

int makecut3(TString DS){

	TSystemDirectory dir("./", "./");
  TList *files = dir.GetListOfFiles();
	files->Sort();
  TSystemFile *file;
  TString fname;
  TString fname2;
  TIter next(files);
	TFile* f;

	int n =0;

	TString detector;
	TString run;
	int run_a, run_b;
	int detector_a, detector_b;
	double factor;
	TString hname;
	
	TH2D* h_energyRMSCal;
	TH2D* h_energyRMSBkg;
	TH1D* h_exposure;
	TH1D* h_projection;
	TGraph* g_energyFit;
	vector<double> forg_energies;
	vector<double> forg_sigmaMax;

	TH1D* h_energyCal1 = new TH1D("h_energyCal1","h_energyCal1",3000,0,300);
	TH1D* h_energyCal2 = new TH1D("h_energyCal2","h_energyCal2",3000,0,300);
	TH1D* h_energyBkg1 = new TH1D("h_energyBkg1","h_energyBkg1",3000,0,300);
	TH1D* h_energyBkg2 = new TH1D("h_energyBkg2","h_energyBkg2",3000,0,300);

  TH1D* h_exposureDSNat = new TH1D("h_exposureDSNat","h_exposureDSNat",3000,0,300);
  TH1D* h_exposureDSEnr = new TH1D("h_exposureDSEnr","h_exposureDSEnr",3000,0,300);
  TH1D* h_chi = new TH1D("h_chi","h_chi",1000,0,2);

  TH1D* h_energyBkgTot1Nat = new TH1D("h_energyBkgTot1Nat","h_energyBkgTot1Nat",3000,0,300);
  TH1D* h_energyBkgTot2Nat = new TH1D("h_energyBkgTot2Nat","h_energyBkgTot2Nat",3000,0,300);
  TH1D* h_energyBkgTot1Enr = new TH1D("h_energyBkgTot1Enr","h_energyBkgTot1Enr",3000,0,300);
  TH1D* h_energyBkgTot2Enr = new TH1D("h_energyBkgTot2Enr","h_energyBkgTot2Enr",3000,0,300);

  TH1D* h_energyCalTot1Nat = new TH1D("h_energyCalTot1Nat","h_energyCalTot1Nat",3000,0,300);
  TH1D* h_energyCalTot2Nat = new TH1D("h_energyCalTot2Nat","h_energyCalTot2Nat",3000,0,300);
  TH1D* h_energyCalTot1Enr = new TH1D("h_energyCalTot1Enr","h_energyCalTot1Enr",3000,0,300);
  TH1D* h_energyCalTot2Enr = new TH1D("h_energyCalTot2Enr","h_energyCalTot2Enr",3000,0,300);
	double exposure;
	double baseline;
	double energy;
	double sigmaMax;
	double chi2;

	double a0, a1, a2, a3, a4, a5;
	double value_up, value_low;
	double threshold_cal, threshold_bkg;

	TString answer;

	TF1* f_fit = new TF1("f_fit","pol4",0,200);
	TF1* f_draw = new TF1("f_draw","sqrt(pow([0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x,2)+[5])",0,300);


	TCanvas *c1 = new TCanvas (DS+"_c1",DS+"_c1",1000,700);
	c1->Divide(3,2);
	TCanvas *c2 = new TCanvas (DS+"_c2",DS+"_c2",1000,700);
	c2->Divide(2,2);

////////////////////////////////////////////////////////	
	while ((file=(TSystemFile*)next())) {
   	fname = file->GetName();
		//cout << fname << endl;
    if (!file->IsDirectory() && fname.EndsWith(".root") && fname.Contains("Det")) {
			run = fname(0,fname.First("_"));
			detector = fname(fname.First("_")+5,fname.Length() - fname.First(".")-2);

			if (DS.Contains("DS0_M1")){
				run_a = 1;
				run_b = 34;
				detector_a = 100;
				detector_b = 200;
				factor = 1;
			}
			else if (DS.Contains("DS1_M1")){
				run_a = 35;
				run_b = 91;
				detector_a = 100;
				detector_b = 200;
				factor = 1;
			}
			else if (DS.Contains("DS2_M1")){
				run_a = 92;
				run_b = 102;
				detector_a = 100;
				detector_b = 200;
				factor = 2.5;
			}		
			else if (DS.Contains("DS3_M1")){
				run_a = 103;
				run_b = 111;
				detector_a = 100;
				detector_b = 200;
				factor = 1;
			}		
			else if (DS.Contains("DS4_M2")){
					run_a = 112;
					run_b = 122;
					detector_a = 200;
					detector_b = 300;
					factor = 1;
			}		
			else if (DS.Contains("DS5_M1")){
				run_a = 123;
				run_b = 139;
				detector_a = 100;
				detector_b = 200;
				factor = 1;
			}	
			else if (DS.Contains("DS5_M2")){
				run_a = 123;
				run_b = 139;
				detector_a = 200;
				detector_b = 300;
				factor = 1.2;
			}	
			else{
				run_a = 1;
				run_b = 139;
				detector_a = 100;
				detector_b = 300;
				factor = 1.2;
			}	
			run_a--;
			run_b--;
			if (run.Atoi()<run_a) continue;
			if (run.Atoi()>run_b) continue;
			if (detector.Atoi()<detector_a) continue;
			if (detector.Atoi()>detector_b) continue;

			f = new TFile(fname);

			//f->ls();
			hname = "h_energyRMSCal_";
			hname+=detector;
			hname+="_";
			hname+=run.Atoi();		
			h_energyRMSCal = (TH2D*)f->Get(hname);
			h_energyRMSCal->SetDirectory(0);
			hname = "h_energyRMSBkg_";
			hname+=detector;
			hname+="_";
			hname+=run.Atoi();		
			h_energyRMSBkg = (TH2D*)f->Get(hname);
			h_energyRMSBkg->SetDirectory(0);
			hname = "h_exposure_";
			hname+=detector;
			hname+="_";
			hname+=run.Atoi();		
			h_exposure = (TH1D*)f->Get(hname);
			h_exposure->SetDirectory(0);




			forg_energies.clear();
			forg_sigmaMax.clear();
			exposure = h_exposure->GetBinContent(1) * h_exposure->GetBinContent(2);
			
			if ((exposure > 0)&&(h_energyRMSCal->GetEntries()>2000)&&(h_energyRMSBkg->GetEntries()>0)){
				cout << "-----------------" << endl;
				cout << n << " " << fname << " : "  << run << " " << detector << " " << exposure << " kg*yr" << endl;
				
				baseline = h_energyRMSCal->ProjectionY(" ",1,1)->FindFirstBinAbove(10,1)/10.;
				//baseline = h_energyRMSCal->ProjectionY(" ",1,1)->GetMaximumBin()/10.;

				for (int k=2400;k>5;k--){
					h_projection = (TH1D*)h_energyRMSCal->ProjectionY(" ",k-3,k+3);
					energy = k/10.;
					
					sigmaMax = h_projection->GetMaximumBin()/10.;
					//cout << energy << " " << sigmaMax <<  endl;
					if (sigmaMax>baseline*factor){
						if ((energy <10)||(k%10 ==0)){
							forg_energies.push_back(energy);
							forg_sigmaMax.push_back(sqrt(sigmaMax*sigmaMax-baseline*baseline*factor*factor));
						}
					}
				}
				g_energyFit = new TGraph(forg_energies.size(), &forg_energies[0], &forg_sigmaMax[0]);
				f_fit->SetParameter(0,0);
				f_fit->SetParLimits(0,-0.5*factor*factor,0.5*factor*factor);
				//f_fit->SetParLimits(2,0,10);
				f_fit->SetRange(2,300);
				g_energyFit->Fit("f_fit","WQMR");
				chi2 = f_fit->GetChisquare()/f_fit->GetNDF() ;
				h_chi->Fill(chi2);

				cout << "base : " << baseline*factor << endl;
				cout << "chi" <<  chi2 << endl;

				f_draw->SetParameter(0,f_fit->GetParameter(0));
				f_draw->SetParameter(1,f_fit->GetParameter(1));
				f_draw->SetParameter(2,f_fit->GetParameter(2));
				f_draw->SetParameter(3,f_fit->GetParameter(3));
				f_draw->SetParameter(4,f_fit->GetParameter(4));	
				//f_draw->SetParameter(3,0);
				//f_draw->SetParameter(4,0);
				f_draw->SetParameter(5,baseline*baseline*factor*factor);

				a0=f_draw->GetParameter(0);
				a1=f_draw->GetParameter(1);
				a2=f_draw->GetParameter(2);
				a3=f_draw->GetParameter(3);
				a4=f_draw->GetParameter(4);	
				a5=f_draw->GetParameter(5);	


				h_energyCal1->Reset();
				h_energyCal2->Reset();
				h_energyBkg1->Reset();
				h_energyBkg2->Reset();
				
				threshold_cal =0;
				threshold_bkg =0;
				for (int k=5;k<2500;k++){
					h_projection = (TH1D*)h_energyRMSCal->ProjectionY(" ",k,k);
					energy = (k+0.5)/10.;
					value_up = sqrt(pow(a0+a1*energy+a2*energy*energy+a3*energy*energy*energy+a4*energy*energy*energy*energy,2)+a5) + 0.25 + 0.003*energy;
					value_low = sqrt(pow(a0+a1*energy+a2*energy*energy+a3*energy*energy*energy+a4*energy*energy*energy*energy,2)+a5) - 0.25 - 0.003*energy;
					h_energyCal1->Fill(energy,h_projection->Integral());
					h_energyCal2->Fill(energy,h_projection->Integral(value_low*10,value_up*10));
					if ((threshold_cal==0)&&(h_projection->Integral(value_low*10,value_up*10)>5)) threshold_cal = energy;

					h_projection = (TH1D*)h_energyRMSBkg->ProjectionY(" ",k,k);
					h_energyBkg1->Fill(energy,h_projection->Integral());
					h_energyBkg2->Fill(energy,h_projection->Integral(value_low*10,value_up*10));
					if ((threshold_cal!=0)&&(threshold_bkg==0)&&(h_projection->Integral(value_low*10,value_up*10)<1)) threshold_bkg = energy;
					if (h_projection->Integral(value_low*10,value_up*10)>3) threshold_bkg = 0;
				
					
				}

				cout << "threshold calibration : " << threshold_cal << endl;
				cout << "threshold background : " << threshold_bkg << endl;



				n++;
	
				c1->cd(1);
				gPad->SetLogz();
				h_energyRMSCal->GetXaxis()->SetRangeUser(0,threshold_bkg+10);
				h_energyRMSCal->GetYaxis()->SetRangeUser(0,threshold_bkg+10);
				h_energyRMSCal->Draw("COLZ");
				f_draw->Draw("SAME");
				c1->cd(2);
				g_energyFit->Draw();
				f_fit->Draw("SAME");
				c1->cd(3);
				h_energyRMSBkg->GetXaxis()->SetRangeUser(0,20);
				h_energyRMSBkg->GetYaxis()->SetRangeUser(0,20);
				h_energyRMSBkg->Draw();
				f_draw->Draw("SAME");
				c1->cd(4);
				gPad->SetLogy();
				h_energyCal1->GetXaxis()->SetRangeUser(0,threshold_bkg+10);
				h_energyCal1->Draw("HIST");
				h_energyCal2->SetLineColor(kRed);
				h_energyCal2->Draw("HIST SAME");
	 			TText *t_1 = new TText(2,50,Form("%.2f keV",threshold_cal));
				t_1->Draw("SAME");
				c1->cd(5);
				gPad->SetLogy();
				h_energyBkg1->GetXaxis()->SetRangeUser(0,threshold_bkg+10);
				h_energyBkg1->Draw("HIST");
				h_energyBkg2->SetLineColor(kRed);
				h_energyBkg2->Draw("HIST SAME");
	 			TText *t_2 = new TText(2,5,Form("%.2f keV",threshold_bkg));
				t_2->Draw("SAME");
				c1->Update();
				hname = fname(0,fname.First("."));
				hname+=".png";
				c1->Print(hname);			

				if (n>0) break;

				answer = "n";
				//cin >> answer;
				if (chi2 <0.18) answer = "y";
				if ((answer.Contains("y"))&&(threshold_bkg>0)&&(threshold_bkg<20)){
					if ((detector.Contains("121"))||
						(detector.Contains("131"))||
						(detector.Contains("141"))||
						(detector.Contains("142"))||
						(detector.Contains("143"))||
						(detector.Contains("144"))||
						(detector.Contains("145"))||
						(detector.Contains("151"))||
						(detector.Contains("171"))||
						(detector.Contains("211"))||
						(detector.Contains("221"))||
						(detector.Contains("222"))||
						(detector.Contains("223"))||
						(detector.Contains("224"))||
						(detector.Contains("225"))||
						(detector.Contains("241"))||
						(detector.Contains("242"))||
						(detector.Contains("243"))||
						(detector.Contains("244"))||
						(detector.Contains("245"))||
						(detector.Contains("255"))||
						(detector.Contains("272"))||
						(detector.Contains("274"))){
						for (int k = threshold_bkg*10 +2 ; k<3000; k++){
							h_exposureDSNat->Fill((k+0.5)/10.,exposure);
							h_energyBkgTot2Nat->Fill((k+0.5)/10.-0.1,h_energyBkg2->GetBinContent(k));
							h_energyCalTot1Nat->Fill((k+0.5)/10.-0.1,h_energyCal1->GetBinContent(k));
							h_energyCalTot2Nat->Fill((k+0.5)/10.-0.1,h_energyCal2->GetBinContent(k));
						}
						h_energyBkgTot1Nat->Add(h_energyBkg1);
						//h_energyCalTot1Nat->Add(h_energyCal1);
					}	
					else{
						for (int k = threshold_bkg*10+2; k<3000; k++){
							h_exposureDSEnr->Fill((k+0.5)/10.,exposure);
							h_energyBkgTot2Enr->Fill((k+0.5)/10.-0.1,h_energyBkg2->GetBinContent(k));
							h_energyCalTot1Enr->Fill((k+0.5)/10.-0.1,h_energyCal1->GetBinContent(k));
							h_energyCalTot2Enr->Fill((k+0.5)/10.-0.1,h_energyCal2->GetBinContent(k));
						}
						h_energyBkgTot1Enr->Add(h_energyBkg1);
						//h_energyCalTot1Enr->Add(h_energyCal1);
					}

					c2->cd(1);
					gPad->SetLogy();
					h_energyBkgTot2Enr->GetXaxis()->SetRangeUser(0,20);
					h_energyBkgTot2Enr->Draw("HIST");
					h_energyBkgTot1Enr->SetLineStyle(2);
					h_energyBkgTot1Enr->Draw("HIST SAME");
					c2->cd(2);
					gPad->SetLogy();
					h_energyBkgTot2Nat->GetXaxis()->SetRangeUser(0,20);
					h_energyBkgTot2Nat->Draw("HIST");
					h_energyBkgTot1Nat->SetLineStyle(2);
					h_energyBkgTot1Nat->Draw("HIST SAME");
					c2->cd(3);
					h_exposureDSEnr->GetXaxis()->SetRangeUser(0,20);
					h_exposureDSEnr->Draw("HIST");
					h_exposureDSNat->SetLineColor(kRed);
					h_exposureDSNat->Draw("HIST SAME");
					c2->cd(4);
					gPad->SetLogy();
					h_chi->GetXaxis()->SetRangeUser(0,0.25);
					h_chi->Draw("HIST");
					c2->Update();
				}
	
				fname2 = DS;
				fname2+= ".dat";
				ofstream myfile;
			  myfile.open (fname2,ios::out | ios::app);
			  myfile << run << " " << detector 
						<< " " << answer 
						<< " exposure " << exposure 
						<< " thresh " << threshold_bkg << " : " 
						<< " fit " << a0 << " " << a1 << " " << a2 << " " << a3 << " " << a4 << " " << a5 
						<< " Band " << 0.25 << " " << 0.003 
						<< " chi2 " << chi2
						<< endl;
			  myfile.close();
			//	cin.get();

				delete g_energyFit;
			
			}
			f->Close();
			delete f;
			delete h_energyRMSCal;
			delete h_energyRMSBkg;
		}
	}

	hname = DS;
	hname+="_final.png";
	c2->Print(hname);			


	fname=DS;
	fname+="_out.root";
	TFile *MyFile = new TFile(fname,"RECREATE");
	MyFile->cd();
	h_exposureDSNat->Write();
	h_energyBkgTot2Nat->Write();
	h_energyBkgTot1Nat->Write();
	h_energyCalTot2Nat->Write();
	h_energyCalTot1Nat->Write();


	h_exposureDSEnr->Write();
	h_energyBkgTot2Enr->Write();
	h_energyBkgTot1Enr->Write();
	h_energyCalTot2Enr->Write();
	h_energyCalTot1Enr->Write();

	h_chi->Write();

	MyFile->Close();

	return 0;
}
