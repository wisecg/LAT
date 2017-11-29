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

int makecut(int subset, int det){

	vector<int> run_Module;
	vector<int> run_DS;	
	vector<int> run_cal_start;	
	vector<int> run_cal_stop;
	vector<int> run_back_start;
	vector<int> run_back_stop;
	vector<int> run_backfile_start;
	vector<int> run_backfile_stop;

	int in_data_1;
	int in_data_2;
	int in_data_3;
	int in_data_4;
	int in_data_5;
	int in_data_6;
	int in_data_7;
	int in_data_8;

	ifstream in;
	in.open("data.txt");
	while (in >> in_data_2 >> in_data_3 >> in_data_4 >> in_data_5 >> in_data_6 >> in_data_7 >> in_data_8){
		run_DS.push_back(in_data_2);
		run_cal_start.push_back(in_data_3);
		run_cal_stop.push_back(in_data_4);
		run_back_start.push_back(in_data_5);
		run_back_stop.push_back(in_data_6);
		run_backfile_start.push_back(in_data_7);
		run_backfile_stop.push_back(in_data_8);
	}
	in.close();
	//cout << "--" << run_back_stop.size() << "  range" << endl;
//////////////////////////////////////////////////////////////////////////////////////////////

	gROOT->ProcessLine(".x $MGDODIR/Majorana/LoadMGDOMJClasses.C");
	
	TString path_background = "/mnt/mjdDisk1/Majorana/data/sandbox/latv3/lat/latSkimDS";
	TString path_calibration = "/mnt/mjdDisk1/Majorana/data/sandbox/latv3/cal-lat/";
	TString path;

	TChain* chain_data;
	int nentries;
	int detector;
	double totaltime;
	double totalmass;

	int run;

	double startTime_s;
	double stopTime_s;
	vector<double> *mAct_g = 0;
	vector<int> *gain = 0;
	vector<int> *P = 0;
	vector<int> *D = 0;
	vector<int> *C = 0;
	vector<double> *energy = 0;
 	vector<MGTWaveform*> *waveform =0;

	TH1D* h_wf;
	TH1D* h_ADC;
	TH2D* h_energyRMSCal = new TH2D(Form("h_energyRMSCal_%i_%i",det,subset),Form("h_energyRMSCal_%i_%i",det,subset),3000,0,300,3000,0,300);
	TH2D* h_energyRMSBkg = new TH2D(Form("h_energyRMSBkg_%i_%i",det,subset),Form("h_energyRMSBkg_%i_%i",det,subset),3000,0,300,3000,0,300);
	double ADCRMS;
	

//////////////////////////////////////////////////////////////////////////////////////////////
	for(int range = 0; range < run_back_stop.size(); range++){
		if (range!=subset) continue;
		cout << range << " of " << run_back_stop.size()  << " : "  << det << " : " <<  run_cal_start[range] << " ... " << run_cal_stop[range] << endl;

		//////////////////////////////////////
		//get calibration histogram
		chain_data = new TChain("skimTree");
		
		for (int runnumber = run_cal_start[range]; runnumber <= run_cal_stop[range]; runnumber++){
			path = path_calibration;
			path +="*_run";
			path +=runnumber;
			path +="_*.root";			
			chain_data->Add(path);
		}
		nentries = chain_data->GetEntries();	
		
    gain = 0;
		mAct_g = 0;
	  P = 0;
	  D = 0;
	  C = 0;
	  energy = 0;
 	  waveform =0;

		chain_data->SetBranchAddress("trapENFCalC", &energy);
		chain_data->SetBranchAddress("gain", &gain);
		chain_data->SetBranchAddress("P", &P);
		chain_data->SetBranchAddress("D", &D);
		chain_data->SetBranchAddress("C", &C);
		chain_data->SetBranchAddress("MGTWaveforms", &waveform);
		chain_data->SetBranchAddress("mAct_g", &mAct_g);
		
	
		for (int event =0; event<nentries;event++){
			chain_data->GetEntry(event);
			for(int j=0;j<energy->size();j++){	
				detector = C->at(j)*100 +  P->at(j)*10 +  D->at(j)*1;
				if (detector!=det) continue;
				if (gain->at(j)!=0) continue;
				//cout << det << " " << detector << " " << C->at(j) << " " << P->at(j) << " " << D->at(j) << chain_data->GetTree()->GetCurrentFile()->GetName() << endl;
				//if (event%1000==0) cout << event << " " << nentries << endl;
				totalmass=mAct_g->at(j);
				h_wf = (TH1D*)waveform->at(j)->GimmeHist(Form("wf_%i",detector));
				h_wf->GetXaxis()->SetRangeUser(3,(h_wf->GetNbinsX()-3)*10);

				h_ADC = new TH1D("h_ADC","h_ADC",h_wf->GetMaximum()-h_wf->GetMinimum()+100,h_wf->GetMinimum()-50,h_wf->GetMaximum()+50);
				h_wf->SetDirectory(0);	
				h_ADC->SetDirectory(0);
				for (int k=1;k<h_wf->GetNbinsX()-3;k++)	h_ADC->Fill(h_wf->GetBinContent(k));
				ADCRMS = h_ADC->GetRMS();
				h_energyRMSCal->Fill(energy->at(j),ADCRMS);

				if (energy->at(j) > 100){				
					h_ADC->GetXaxis()->SetRangeUser(h_wf->GetMinimum()-50,h_wf->GetMinimum()+50);
					h_ADC->GetXaxis()->SetRangeUser(h_ADC->GetMean()-10,h_ADC->GetMean()+10);
					ADCRMS = h_ADC->GetRMS();
					if (ADCRMS>0) h_energyRMSCal->Fill(0.05,ADCRMS);
				}

				delete h_wf;
				delete h_ADC;
			}
			
		}
		delete chain_data;
		//cout << h_energyRMSCal->GetEntries() << " " <<  h_energyRMSCal->GetRMS(1) << endl;
		if (h_energyRMSCal->GetRMS(1) == 0) return 1;

		//////////////////////////////////////
		//get background
		chain_data = new TChain("skimTree");
		cout << range << " of " << run_back_stop.size()  << " : " << det << " : " <<  run_back_start[range] << " ... " << run_back_stop[range] << endl;
		for (int runnumber = run_backfile_start[range]; runnumber <= run_backfile_stop[range]; runnumber++){
			path = path_background;
			path +=run_DS[range];
			path +="_";
			path +=runnumber;
			path +="_*.root";
			chain_data->Add(path);
			cout << path << endl;
		}
		nentries = chain_data->GetEntries();	
		//cout <<"--nentries: " << nentries << endl;

		double time[run_back_stop[range]-run_back_start[range]];
		for (int i = 0; i<run_back_stop[range]-run_back_start[range];i++) time[i]=0;
		

		run = 0;

		startTime_s = 0;
		stopTime_s = 0;
    gain = 0;
	  P = 0;
	  D = 0;
	  C = 0;
	  energy = 0;
 	  waveform =0;

		chain_data->SetBranchAddress("run", &run);
		chain_data->SetBranchAddress("startTime_s", &startTime_s);
		chain_data->SetBranchAddress("stopTime_s", &stopTime_s);

		chain_data->SetBranchAddress("trapENFCalC", &energy);
		chain_data->SetBranchAddress("gain", &gain);
		chain_data->SetBranchAddress("P", &P);
		chain_data->SetBranchAddress("D", &D);
		chain_data->SetBranchAddress("C", &C);
		chain_data->SetBranchAddress("MGTWaveforms", &waveform);

		for (int event =0; event<nentries;event++){
			chain_data->GetEntry(event);
			//cout << run << " " << stopTime_s-startTime_s<< endl;
			if (run < run_back_start[range]) continue;
			if (run > run_back_stop[range]) continue;
			time[run_back_stop[range]-run]=stopTime_s-startTime_s;
			for(int j=0;j<energy->size();j++){	
				detector = C->at(j)*100 +  P->at(j)*10 +  D->at(j)*1;
				if (detector!=det) continue;
				if (gain->at(j)!=0) continue;
				//if (event%1000==0) cout << "run# "<< run << " , " << event << " " << nentries << endl;
				h_wf = (TH1D*)waveform->at(j)->GimmeHist(Form("wf_%i",detector));
				h_wf->GetXaxis()->SetRangeUser(3,(h_wf->GetNbinsX()-3)*10);

				h_ADC = new TH1D("h_ADC","h_ADC",h_wf->GetMaximum()-h_wf->GetMinimum()+100,h_wf->GetMinimum()-50,h_wf->GetMaximum()+50);
				h_wf->SetDirectory(0);	
				h_ADC->SetDirectory(0);
				for (int k=1;k<h_wf->GetNbinsX()-3;k++)	h_ADC->Fill(h_wf->GetBinContent(k));
				ADCRMS = h_ADC->GetRMS();
				h_energyRMSBkg->Fill(energy->at(j),ADCRMS);


				delete h_wf;
				delete h_ADC;
			}
		}
		delete chain_data;
		totaltime =0;
		for (int i = 0; i<run_back_stop[range]-run_back_start[range];i++){
			totaltime+=time[i];
			//cout << time[i] << endl;
		}

		totalmass=totalmass/1000.;
		totaltime=totaltime/31536000;
		//cout << totaltime << " " << totalmass << endl;
	}
	
	
/*
	TCanvas *c = new TCanvas("c","c",1200,800);
	c->Divide(2);
	c->cd(1);
	h_energyRMSCal->Draw("COLZ");	
	c->cd(2);
	h_energyRMSBkg->Draw("COLZ");
*/
	TH1D *h_exposure = new TH1D(Form("h_exposure_%i_%i",det,subset),Form("h_exposure_%i_%i",det,subset),10,0,10);
	h_exposure->SetBinContent(1,totaltime);
	h_exposure->SetBinContent(2,totalmass);
	
	
	TFile *MyFile = new TFile(Form("%i_Det_%i.root",subset,det),"NEW");
	MyFile->cd();
	h_energyRMSCal->Write();
	h_energyRMSBkg->Write();
	h_exposure->Write();
	MyFile->Close();


	delete h_energyRMSCal;
	delete h_energyRMSBkg;
	delete h_exposure;
	return 0;
}




///////////////////////
void makecut(int det){
	for (int i=0; i<139;i++){
		makecut(i,det);
	}
}





///////////////////////
void makecut(){
	makecut(0,111);
}


