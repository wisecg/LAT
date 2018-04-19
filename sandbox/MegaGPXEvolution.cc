// Run the evolved version of Wenqin

#include "GPXFitter.hh"
#include "TStyle.h"
#include "TChain.h"
#include "TFile.h"
#include "TH1D.h"
#include "RooFitResult.h"
#include "RooWorkspace.h"
#include "RooAbsReal.h"
#include "RooArgList.h"
#include "RooAbsArg.h"
#include "RooMsgService.h"
#include "RooRealVar.h"
#include "TROOT.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TMath.h"
#include "TLegend.h"
// #include <Python/Python.h>

using namespace RooFit;
using namespace std;

map<string, vector<double>> RunWenqin(vector<string> argS, int fDS, float fitMin, float fitMax, string theCut);
map<string, vector<double>> RunWenqin(vector<string> argS, int fDS, float fitMin, float fitMax, string ftype, string theCut);
void RunBasicFit(int fDS, double fitMin, double fitMax, string ftype);
void RunCutComparison(int fDS, double fitMin, double fitMax, string ftype);
void CutEfficiencyStudy(int fDS, double fitMin, double fitMax, string ftype);

int main(int argc, char** argv)
{
  	// gROOT->ProcessLine("gErrorIgnoreLevel = 3001;");
  	// gROOT->ProcessLine("RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);");
	if(argc <= 4)
	{
		cout << "Usage: " << argv[0] << " [DS] [Fit Min] [Fit Max] [Nat/Enr/Cut]" << endl;
		return 0;
	}

	int fDS = atoi(argv[1]);
  float fitMin = atof(argv[2]);
	float fitMax = atof(argv[3]);
	string ftype = argv[4];
	gStyle->SetOptStat(0);
	// gStyle->SetPalette(kRainBow);

  // Prelim exposure
  std::map<int, std::vector<double>> liveTime;
	liveTime[0] = {460.054, 171.021};
	liveTime[1] = {661.811, 63.2937};
  liveTime[2] = {106.286, 10.6791};
	liveTime[3] = {368.52, 81.7408};
	liveTime[4] = {102.858, 73.8446};
	liveTime[5] = {492.158+182.193, 138.461+197.769};
  liveTime[6] = {460.054+661.811+106.286+368.52+102.858+492.158+182.193, 171.021+63.2937+10.6791+81.7408+73.8446+138.461+197.769};

  // RunCutComparison(fDS, fitMin, fitMax, ftype);
	RunBasicFit(fDS, fitMin, fitMax, ftype);

	return 0;
}

map<string, vector<double>> RunWenqin(vector<string> argS, int fDS, float fitMin, float fitMax, string theCut, int idx)
{
	string ftype = Form("Ch%d",idx);
	map<string, vector<double>> WenqinMap = RunWenqin(argS, fDS, fitMin, fitMax, ftype, theCut);

	return WenqinMap;
}

map<string, vector<double>> RunWenqin(vector<string> argS, int fDS, float fitMin, float fitMax, string ftype, string theCut)
{
	GPXFitter *fitter = new GPXFitter(fDS, fitMin, fitMax);
	// This is just a string for output files
	fitter->SetSavePrefix(Form("DS%d_%s_%.1f_%.1f", fDS, ftype.c_str(), fitMin, fitMax));

	// Load data from TChain with a cut string
	TChain *skimTree = new TChain("skimTree");
    if(fDS == 6) {
      skimTree->Add("/mnt/mjdDisk1/Majorana/data/sandbox/latv4/cuts/corrfs_rn/corrfs_rn-DS1-*.root" );
      skimTree->Add("/mnt/mjdDisk1/Majorana/data/sandbox/latv4/cuts/corrfs_rn/corrfs_rn-DS2-*.root" );
      skimTree->Add("/mnt/mjdDisk1/Majorana/data/sandbox/latv4/cuts/corrfs_rn/corrfs_rn-DS3-*.root" );
      skimTree->Add("/mnt/mjdDisk1/Majorana/data/sandbox/latv4/cuts/corrfs_rn/corrfs_rn-DS4-*.root" );
      skimTree->Add("/mnt/mjdDisk1/Majorana/data/sandbox/latv4/cuts/corrfs_rn/corrfs_rn-DS5-*.root" );
    }
  	else skimTree->Add(Form("/mnt/mjdDisk1/Majorana/data/sandbox/latv4/cuts/corrfs_rn/corrfs_rn-DS%d-*.root", fDS) );
  	fitter->LoadChainData(skimTree, theCut);

	// Construct PDF and do fit
	fitter->ConstructPDF();
  fitter->DoFit();
	fitter->DrawBasicShit(0.2, false, false);
	fitter->GetFitResult()->Print("v");
	map<string, vector<double>> WenqinMap;
	for(auto &argN: argS){
		vector<double> WenqinParameter = fitter->GetVar( Form("%s", argN.c_str()) );
		WenqinMap[argN.c_str()] = WenqinParameter;
	}
	return WenqinMap;
}


void CutEfficiencyStudy(int fDS, double fitMin, double fitMax, string ftype)
{
  ////// First round of fits
    int bNat = 0;
    string theCut0 = "";
    theCut0 += Form("trapENFCal>=%.2f&&trapENFCal<=%.2f", fitMin, fitMax); // Energy cut for fit range

    if(ftype == "Nat" || ftype == "nat"){
      theCut0 += "&&isNat"; // Set Enriched or Natural
      bNat = 1;
    }
    else if(ftype == "Enr" || ftype == "enr") {
      theCut0 += "&&isEnr";
      bNat = 0;
    }
    else theCut0 = ftype;
    GPXFitter *fitter0 = new GPXFitter(fDS, fitMin, fitMax);
    fitter0->SetSavePrefix(Form("EffTest0_DS%d_%s_%.1f_%.1f", fDS, ftype.c_str(), fitMin, fitMax));

    // Load data from TChain with a cut string
    TChain *skimTree0 = new TChain("skimTree");
    if(fDS == 6) {
    skimTree0->Add("/mnt/mjdDisk1/Majorana/data/sandbox/latv4/cuts/corrfs_rn/corrfs_rn-DS1-*.root" );
    skimTree0->Add("/mnt/mjdDisk1/Majorana/data/sandbox/latv4/cuts/corrfs_rn/corrfs_rn-DS2-*.root" );
    skimTree0->Add("/mnt/mjdDisk1/Majorana/data/sandbox/latv4/cuts/corrfs_rn/corrfs_rn-DS3-*.root" );
    skimTree0->Add("/mnt/mjdDisk1/Majorana/data/sandbox/latv4/cuts/corrfs_rn/corrfs_rn-DS4-*.root" );
    skimTree0->Add("/mnt/mjdDisk1/Majorana/data/sandbox/latv4/cuts/corrfs_rn/corrfs_rn-DS5-*.root" );
    }
    else skimTree0->Add(Form("/mnt/mjdDisk1/Majorana/data/sandbox/latv4/cuts/corrfs_rn/corrfs_rn-DS%d-*.root", fDS) );
    fitter0->LoadChainData(skimTree0, theCut0);

    // Construct PDF and do fit
    fitter0->ConstructPDF();
    // fitter0->DrawModels(0.2);
    fitter0->DoFit();
    fitter0->DrawBasicShit(0.2, false, false);
    vector<double> valTrit0 = fitter0->GetVar("Tritium");
    double tritVal0 = valTrit0[0];
    double tritErr0 = valTrit0[1];

    vector<double> valGe0 = fitter0->GetVar("Ge68");
    double Ge68Val0 = valGe0[0];
    double Ge68Err0 = valGe0[1];

    vector<double> valFe0 = fitter0->GetVar("Fe55");
    double Fe55Val0 = valFe0[0];
    double Fe55Err0 = valFe0[1];


  	//// Set cuts here
    vector<double> minList = {1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 4.25, 4.5, 4.75, 5.0};
    vector<double> minErr = {0., 0., 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    std::vector<double> tritVal;
    std::vector<double> tritErr;
    std::vector<double> tritValEff;
    std::vector<double> tritErrEff;
    std::vector<double> Ge68Val;
    std::vector<double> Ge68Err;

    for(auto minVal : minList)
    {
          string theCut = "";
          theCut += Form("trapENFCal>=%.2f&&trapENFCal<=%.2f", minVal, fitMax); // Energy cut for fit range

        	if(ftype == "Nat" || ftype == "nat"){
        		theCut += "&&isNat"; // Set Enriched or Natural
        		bNat = 1;
        	}
        	else if(ftype == "Enr" || ftype == "enr") {
        		theCut += "&&isEnr";
        		bNat = 0;
        	}
        	else theCut = ftype;

          GPXFitter *fitter = new GPXFitter(fDS, minVal, fitMax);
      	  // This is just a string for output files
          fitter->SetSavePrefix(Form("EffTest_DS%d_%s_%.1f_%.1f", fDS, ftype.c_str(), minVal, fitMax));

      	  // Load data from TChain with a cut string
      	  TChain *skimTree = new TChain("skimTree");
          if(fDS == 6) {
          skimTree->Add("/mnt/mjdDisk1/Majorana/data/sandbox/latv4/cuts/corrfs_rn/corrfs_rn-DS1-*.root" );
          skimTree->Add("/mnt/mjdDisk1/Majorana/data/sandbox/latv4/cuts/corrfs_rn/corrfs_rn-DS2-*.root" );
          skimTree->Add("/mnt/mjdDisk1/Majorana/data/sandbox/latv4/cuts/corrfs_rn/corrfs_rn-DS3-*.root" );
          skimTree->Add("/mnt/mjdDisk1/Majorana/data/sandbox/latv4/cuts/corrfs_rn/corrfs_rn-DS4-*.root" );
          skimTree->Add("/mnt/mjdDisk1/Majorana/data/sandbox/latv4/cuts/corrfs_rn/corrfs_rn-DS5-*.root" );
          }
          else skimTree->Add(Form("/mnt/mjdDisk1/Majorana/data/sandbox/latv4/cuts/corrfs_rn/corrfs_rn-DS%d-*.root", fDS) );
          fitter->LoadChainData(skimTree, theCut);

      	  // Construct PDF and do fit
      	  fitter->ConstructPDF();
          fitter->SetParameter("Ge68", Ge68Val0, true);
          fitter->SetParameter("Fe55", Fe55Val0, true);
          fitter->DoFit();
          fitter->GetFitResult()->Print("v");
          // cout << "Drawing" << endl;
          fitter->DrawBasicShit(0.3, false, false);
          vector<double> vals = fitter->GetVar("Tritium");
          tritVal.push_back(vals[0]);
          tritErr.push_back(vals[1]);
          cout << "Fit Values: " << vals[0] << " +/- " << vals[1] << endl;
        // vector<double> vals2 = fitter->GetVar("Ge68");
        // Ge68Val.push_back(10*vals2[0]);
        // Ge68Err.push_back(10*vals2[1]);
      }

        // vector<string> argList = {"Tritium", "Ge68", "Fe55", "Bkg"};
        // vector<string> argList = {"Tritium"};
        // map<string, vector<double>> ValMap;
      	// if(ftype == "Nat" || ftype == "nat" || ftype == "Enr" || ftype == "enr") ValMap = RunWenqin(argList, fDS, fitMin, fitMax, ftype, theCut);
      	// else  ValMap = RunWenqin(argList, fDS, fitMin, fitMax, theCut);
      	// cout << "LiveTime: " << liveTime[fDS][bNat] << endl;
  /*


        // Calculate Rates
        for(auto &kv : ValMap)
      	{
      		if(kv.first == "Tritium") {
      			cout << kv.first << ":[" << kv.second[0]*tritScale << "," << kv.second[1]*tritScale << "," << kv.second[2]*tritScale << "]"<< endl;
      			cout << "Rates (c/keV/kg/day):" << kv.second[0]*tritScale/liveTime[fDS][bNat]/(fitMax-fitMin) << " +" << kv.second[1]*tritScale/liveTime[fDS][bNat]/(fitMax-fitMin) << " " << kv.second[2]*tritScale/liveTime[fDS][bNat]/(fitMax-fitMin) << endl;
      			cout << "Rates (2-4 keV) (c/keV/kg/day):" << kv.second[0]*tritScale*0.224487/liveTime[fDS][bNat]/(fitMax-fitMin) << " +" << kv.second[1]*tritScale*0.224487/liveTime[fDS][bNat]/(fitMax-fitMin) << " " << kv.second[2]*tritScale*0.224487/liveTime[fDS][bNat]/(fitMax-fitMin) << endl;
      		}
      		else {
      			cout << kv.first << ":[" << kv.second[0] << "," << kv.second[1] << "," << kv.second[2] << "]"<< endl;
      			cout << "Rates (c/keV/kg/day):" << kv.second[0]/liveTime[fDS][bNat]/(fitMax-fitMin) << " +" << kv.second[1]/liveTime[fDS][bNat]/(fitMax-fitMin) << " " << kv.second[2]/liveTime[fDS][bNat]/(fitMax-fitMin) << endl;
      		}
      	}
  */

    std::string tritDir = "/mnt/mjdDisk1/Majorana/users/psz/CUORE/MJDAnalysis/Wenqin/Data";
    TFile *tritFile = new TFile(Form("%s/TritSpec.root", tritDir.c_str()));
    TH1D *tritSpec = dynamic_cast<TH1D*>(tritFile->Get("tritHist"));
    double tritScale0 = 1./tritSpec->Integral(tritSpec->FindBin(fitMin), tritSpec->FindBin(fitMax));
    double tritValCorr0 = tritVal0*tritScale0;
    double tritErrCorr0 = tritErr0*tritScale0;
    std::vector<double> tritValCorr;
    std::vector<double> tritErrCorr;

    cout << "Initial Corrected Tritium: " << tritValCorr0 << " +/- " << tritErrCorr0 << endl;
    for(int i = 0; i < minList.size(); i++)
    {
      double tritScale = 1./tritSpec->Integral(tritSpec->FindBin(minList[i]), tritSpec->FindBin(fitMax));
      cout <<"Minval: " << minList[i] <<"  Tritium (uncorrected): " << tritVal[i]  << "  Tritium (counts): " << tritVal[i]*tritScale << " +/- " << tritErr[i]*tritScale << endl;
      tritValCorr.push_back(tritVal[i]*tritScale/tritValCorr0);
      tritErrCorr.push_back(tritErr[i]*tritScale/tritValCorr0);
    }

    double *mList = &minList[0];
    double *mErr = &minErr[0];
    double *tList = &tritValCorr[0];
    double *tErr = &tritErrCorr[0];

    TCanvas *c1 = new TCanvas("c1", "c1", 1200, 800);
  	TGraphErrors *g1 = new TGraphErrors((int)minList.size(), mList, tList, mErr, tErr);
    g1->SetTitle("Fit Range vs Tritium Efficiency");
  	g1->GetYaxis()->SetTitle("Efficiency");
  	g1->GetXaxis()->SetTitle("Minimum fit range (keV)");
  	g1->SetMarkerStyle(21);
  	g1->SetMarkerColor(kBlue);
  	g1->SetFillColor(kBlue);

    auto tmin = min_element(tritValCorr.begin(), tritValCorr.end());
    auto tmax = max_element(tritValCorr.begin(), tritValCorr.end());
    g1->SetMaximum(*tmax + *tmax*0.1);
    g1->SetMinimum(*tmin - *tmin*0.1);
  	g1->Draw("APL");
    // TLegend *leg1 = new TLegend(0.38, 0.70, 0.65, 0.88);
    // leg1->AddEntry(g1,"Tritium Corrected", "p");
    // leg1->AddEntry(g2,"Ge68 (x10)", "p");
    // leg1->Draw();
    c1->SaveAs(Form("DS%d_TritEff_%s.pdf", fDS, ftype.c_str()));
}

void RunCutComparison(int fDS, double fitMin, double fitMax, string ftype)
{
  ////// First round of fits
		string inDir = "/mnt/mjdDisk1/Majorana/data/sandbox/latv4/";
		string theCut0 = "";
    int bNat = 0;
    theCut0 += Form("trapENFCal>=%.2f&&trapENFCal<=%.2f", fitMin, fitMax); // Energy cut for fit range

    if(ftype == "Nat" || ftype == "nat"){
      theCut0 += "&&isNat"; // Set Enriched or Natural
      bNat = 1;
    }
    else if(ftype == "Enr" || ftype == "enr") {
      theCut0 += "&&isEnr";
      bNat = 0;
    }
    else theCut0 = ftype;
    GPXFitter *fitter0 = new GPXFitter(fDS, fitMin, fitMax);
    fitter0->SetSavePrefix(Form("CutData_DS%d_%s_%.1f_%.1f", fDS, ftype.c_str(), fitMin, fitMax));

    // Load data from TChain with a cut string
    TChain *skimTree0 = new TChain("skimTree");
    if(fDS == 6) {
    skimTree0->Add(Form("%scuts/corrfs_rn/corrfs_rn-DS1-*.root", inDir.c_str()) );
    skimTree0->Add(Form("%scuts/corrfs_rn/corrfs_rn-DS2-*.root", inDir.c_str()) );
    skimTree0->Add(Form("%scuts/corrfs_rn/corrfs_rn-DS3-*.root", inDir.c_str()) );
    skimTree0->Add(Form("%scuts/corrfs_rn/corrfs_rn-DS4-*.root", inDir.c_str()) );
    skimTree0->Add(Form("%scuts/corrfs_rn/corrfs_rn-DS5-*.root", inDir.c_str()) );
    }
    else skimTree0->Add(Form("%scuts/corrfs_rn/corrfs_rn-DS%d-*.root", inDir.c_str(), fDS) );
    fitter0->LoadChainData(skimTree0, theCut0);

    // Construct PDF and do fit
    fitter0->ConstructPDF();
    // fitter0->DrawModels(0.2);
    fitter0->DoFit("Minuit");
    fitter0->DrawBasicShit(0.2, false, false);
    fitter0->GetFitResult()->Print("v");

    vector<double> valTrit0 = fitter0->GetVar("Tritium");
    double tritVal0 = valTrit0[0];
    double tritErr0 = valTrit0[1];

    vector<double> valGe0 = fitter0->GetVar("Ge68");
    double Ge68Val0 = valGe0[0];
    double Ge68Err0 = valGe0[1];

    vector<double> valFe0 = fitter0->GetVar("Fe55");
    double Fe55Val0 = valFe0[0];
    double Fe55Err0 = valFe0[1];


    GPXFitter *fitter1 = new GPXFitter(fDS, fitMin, fitMax);
    fitter1->SetSavePrefix(Form("FullData_DS%d_%s_%.1f_%.1f", fDS, ftype.c_str(), fitMin, fitMax));

    // Load data from TChain with a cut string
    TChain *skimTree1 = new TChain("skimTree");
    if(fDS == 6) {
    skimTree1->Add("/mnt/mjdDisk1/Majorana/data/sandbox/latv4/lat/latSkimDS1_*.root" );
    skimTree1->Add("/mnt/mjdDisk1/Majorana/data/sandbox/latv4/lat/latSkimDS2_*.root" );
    skimTree1->Add("/mnt/mjdDisk1/Majorana/data/sandbox/latv4/lat/latSkimDS3_*.root" );
    skimTree1->Add("/mnt/mjdDisk1/Majorana/data/sandbox/latv4/lat/latSkimDS4_*.root" );
    skimTree1->Add("/mnt/mjdDisk1/Majorana/data/sandbox/latv4/lat/latSkimDS5_*.root" );
    }
    else if(fDS == 5){
      for(int i = 80; i < 113; i++) skimTree1->Add(Form("/mnt/mjdDisk1/Majorana/data/sandbox/latv4/lat/latSkimDS%d_%d_*.root", fDS, i));
    }
    else skimTree1->Add(Form("/mnt/mjdDisk1/Majorana/data/sandbox/latv4/lat/latSkimDS%d_*.root", fDS));
    cout << "Full Data entries: " << skimTree1->GetEntries() << endl;
    fitter1->LoadChainData(skimTree1, theCut0);

    // Construct PDF and do fit
    fitter1->ConstructPDF();
    fitter1->DoFit("Minuit");
    fitter1->DrawBasicShit(0.2, true, false, false);
    fitter1->GetFitResult()->Print("v");

    vector<double> valTrit1 = fitter1->GetVar("Tritium");
    double tritVal1 = valTrit1[0];
    double tritErr1 = valTrit1[1];

    vector<double> valGe1 = fitter1->GetVar("Ge68");
    double Ge68Val1 = valGe1[0];
    double Ge68Err1 = valGe1[1];

    vector<double> valFe1 = fitter1->GetVar("Fe55");
    double Fe55Val1 = valFe1[0];
    double Fe55Err1 = valFe1[1];

    cout << "Before and After cut comparison:" << endl;
    cout << "Tritium (Before): " << tritVal1 << " +/- " << tritErr1 << endl;
    cout << "Tritium (After): " << tritVal0 << " +/- " << tritErr0 << endl;
    cout << "Ge68 (Before): " << Ge68Val1 << " +/- " << Ge68Err1 << endl;
    cout << "Ge68 (After): " << Ge68Val0 << " +/- " << Ge68Err0 << endl;
    cout << "Fe55 (Before): " << Fe55Val1 << " +/- " << Fe55Err1 << endl;
    cout << "Fe55 (After): " << Fe55Val0 << " +/- " << Fe55Err0 << endl;
    cout << "Ratios: " << tritVal0/tritVal1 << " (Tritium) --- " << Ge68Val0/Ge68Val1 << " (Ge68) --- " << Fe55Val0/Fe55Val1 << " (Fe55)" << endl;
}

void RunBasicFit(int fDS, double fitMin, double fitMax, string ftype)
{
  // Basic fit, 1 round
	string inDir = "/Users/brianzhu/project/";
		string theCut0 = "";
    int bNat = 0;
    theCut0 += Form("trapENFCal>=%.2f&&trapENFCal<=%.2f", fitMin, fitMax); // Energy cut for fit range

    if(ftype == "Nat" || ftype == "nat"){
      theCut0 += "&&isNat"; // Set Enriched or Natural
      bNat = 1;
    }
    else if(ftype == "Enr" || ftype == "enr") {
      theCut0 += "&&isEnr";
      bNat = 0;
    }
		else if(ftype == "All" || ftype == "all")
		{
			theCut0 += "";
		}
    else theCut0 = ftype;
    GPXFitter *fitter0 = new GPXFitter(fDS, fitMin, fitMax);

    // Load data from TChain with a cut string
    TChain *skimTree0 = new TChain("skimTree");
    if(fDS == -1) {
    skimTree0->Add(Form("%s/cuts/corrfs_rn/corrfs_rn-DS1-*.root", inDir.c_str()) );
    skimTree0->Add(Form("%s/cuts/corrfs_rn/corrfs_rn-DS2-*.root", inDir.c_str()) );
    skimTree0->Add(Form("%s/cuts/corrfs_rn/corrfs_rn-DS3-*.root", inDir.c_str()) );
    skimTree0->Add(Form("%s/cuts/corrfs_rn/corrfs_rn-DS4-*.root", inDir.c_str()) );
    skimTree0->Add(Form("%s/cuts/corrfs_rn/corrfs_rn-DS5-*.root", inDir.c_str()) );
		fitter0->SetSavePrefix(Form("BasicFit_All_%s_%.1f_%.1f", ftype.c_str(), fitMin, fitMax));
		}
    else {
			skimTree0->Add(Form("%s/cuts/corrfs_rn/corrfs_rn-DS%d-*.root", inDir.c_str(), fDS) );
			fitter0->SetSavePrefix(Form("BasicFit_DS%d_%s_%.1f_%.1f", fDS, ftype.c_str(), fitMin, fitMax));
		}
    fitter0->LoadChainData(skimTree0, theCut0);

    // Construct PDF and do fit
    fitter0->ConstructPDF();
    fitter0->DoFit("Minuit");
    // fitter0->DoFitEff("Minuit");
    fitter0->DrawBasicShit(0.3, false, false, false);
    fitter0->GetFitResult()->Print("v");
    // fitter0->GetFitResultEff()->Print("v");
		// cout << "Extended Term: " << fitter0->GetWorkspace()->pdf("model")->extendedTerm(1) << endl;

    vector<double> valTrit0 = fitter0->GetVar("Tritium");
    double tritVal0 = valTrit0[0];
    double tritErr0 = valTrit0[1];

    vector<double> valGe0 = fitter0->GetVar("Ge68");
    double Ge68Val0 = valGe0[0];
    double Ge68Err0 = valGe0[1];

		vector<double> valZn0 = fitter0->GetVar("Zn65");
    double Zn65Val0 = valZn0[0];
    double Zn65Err0 = valZn0[1];

    vector<double> valFe0 = fitter0->GetVar("Fe55");
    double Fe55Val0 = valFe0[0];
    double Fe55Err0 = valFe0[1];

    cout << "Basic Fit Results DS" << fDS << " (" << ftype.c_str() << ")" << endl;
    cout << "Tritium: " << tritVal0 << " +/- " << tritErr0 << endl;
    cout << "Ge68: " << Ge68Val0 << " +/- " << Ge68Err0 << endl;
		cout << "Zn65: " << Zn65Val0 << " +/- " << Zn65Err0 << endl;
    cout << "Fe55: " << Fe55Val0 << " +/- " << Fe55Err0 << endl;
}
