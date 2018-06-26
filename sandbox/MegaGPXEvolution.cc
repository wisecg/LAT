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

void RunBasicFit(string fDS, double fitMin, double fitMax, string fMode);

int main(int argc, char** argv)
{
  	// gROOT->ProcessLine("gErrorIgnoreLevel = 3001;");
  	// gROOT->ProcessLine("RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);");

	if(argc <= 4)
	{
		cout << "Usage: " << argv[0] << " [DS (string)] [Fit Min] [Fit Max] [Cut Mode (string)]" << endl;
		return 0;
	}

	vector<string> dsList = {"0", "1", "2", "3", "4", "5A", "5B", "5C", "All", "LowBkg"};
	vector<string> modeList = {"All", "Nat", "Enr", "M1LowBkg", "M1All", "M2LowBkg"};

	string fDS = argv[1];
	float fitMin = atof(argv[2]);
	float fitMax = atof(argv[3]);
	string fMode = argv[4];

	bool bDS = std::find(dsList.begin(), dsList.end(), fDS) != dsList.end();
	bool bMode = std::find(modeList.begin(), modeList.end(), fMode) != modeList.end();
	if(bDS == 0)
	{
		cout << fDS << " is not an available dataset option!" << endl;
		cout << "Options are: 0, 1, 2, 3, 4, 5A, 5B, 5C, All, LowBkg" << endl;
		return 0;
	}
	if(bMode == 0)
	{
		cout << fMode << " is not an available mode option!" << endl;
		cout << "Options are: All, Nat, Enr, M1LowBkg, M1All, M2LowBkg" << endl;
		return 0;
	}

	gStyle->SetOptStat(0);
	RunBasicFit(fDS, fitMin, fitMax, fMode);

	return 0;
}

void RunBasicFit(string fDS, double fitMin, double fitMax, string fMode)
{
		// Calculated and saved from lat-expo.py
		// Reject C2P5D3 from good detector list in M2, it only contributes noise in DS5a, doesn't seem to exist in DS5b and is tiny in DS5c
		GPXFitter *fitter = new GPXFitter(fDS, fitMin, fitMax, fMode);

		std::map<std::string, std::vector<double>> expoFull;
		expoFull["0"] = {400.9074811876048, 157.92594939476592};
		expoFull["1"] = {631.539647097231, 29.523527363230723};
  	expoFull["2"] = {104.3238418726033, 5.177826450022776};
		expoFull["3"] = {334.9535551781866, 48.64815212750721};
		expoFull["4"] = {50.778998571351394, 23.149078616249657};
		expoFull["5A"] = {798.8446523114309, 337.4631592176636};
		expoFull["5B"] = {624.1459839987111, 226.37080226302683};
		expoFull["5C"] = {173.0972327703449, 63.354681716529285};
		expoFull["All"] = {3118.5913929874637, 891.613177148996};
		expoFull["LowBkg"] = {2717.683911799859, 733.6872277542301};

		std::map<std::string, std::vector<double>> expoM1;
		expoM1["0"] = {400.9074811876048, 157.92594939476592};
		expoM1["1"] = {631.539647097231, 29.523527363230723};
  	expoM1["2"] = {104.3238418726033, 5.177826450022776};
		expoM1["3"] = {334.9535551781866, 48.64815212750721};
		expoM1["4"] = {0., 0.};
		expoM1["5A"] = {643.9921013269255, 190.41003908459194};
		expoM1["5B"] = {481.69736138877846, 114.86835044020026};
		expoM1["5C"] = {137.6058203337168, 30.73016914149956};
		expoM1["All"] = {2735.0198083850464, 577.2840140018184};
		expoM1["LowBkg"] = {2334.1123271974416, 419.3580646070525};

		std::map<std::string, std::vector<double>> expoM2;
		expoM2["0"] = {0., 0.};
		expoM2["1"] = {0., 0.};
  	expoM2["2"] = {0., 0.};
		expoM2["3"] = {0., 0.};
		expoM2["4"] = {50.778998571351394, 23.149078616249657};
		expoM2["5A"] = {154.8525509845091, 147.05312013307218};
		expoM2["5B"] = {142.44862260993472, 111.50245182282691};
		expoM2["5C"] = {35.49141243662807, 32.62451257502965};
		expoM2["All"] = {383.5715846024233, 314.3291631471784};
		expoM2["LowBkg"] = {383.5715846024233, 314.3291631471784};

		string inDir = "/Users/brianzhu/project/LATv2/bkg/cut/final95";
		string theCut = "";
    theCut += Form("trapENFCal>=%.2f&&trapENFCal<=%.2f", fitMin, fitMax); // Energy cut for fit range

		// Set cut mode: Enr, Nat, All, or specific detector combo
		if(fMode == "Nat")
		{
      theCut += "&&isNat"; // Set Enriched or Natural
			fitter->SetExposureMap(expoFull);
		}
    else if(fMode == "Enr")
		{
      theCut += "&&isEnr";
			fitter->SetExposureMap(expoFull);
		}
		else if(fMode == "All")
		{
			theCut += "";
			fitter->SetExposureMap(expoFull);
		}
		else if(fMode == "M1LowBkg")
		{
			theCut += "&&(C==1&&P==1&&D==2)||(C==1&&P==1&&D==3)||(C==1&&P==1&&D==4)||(C==1&&P==2&&D==2)||(C==1&&P==2&&D==3)||(C==1&&P==3&&D==2)||(C==1&&P==3&&D==3)||(C==1&&P==3&&D==4)||(C==1&&P==5&&D==3)||(C==1&&P==6&&D==1)||(C==1&&P==6&&D==3)||(C==1&&P==6&&D==4)||(C==1&&P==7&&D==2)||(C==1&&P==7&&D==3)||(C==1&&P==7&&D==4)";
			fitter->SetExposureMap(expoM1);
		}
		else if(fMode == "M1All")
		{
			theCut += "&&(C==1&&P==1&&D==2)||(C==1&&P==1&&D==3)||(C==1&&P==1&&D==4)||(C==1&&P==2&&D==2)||(C==1&&P==2&&D==3)||(C==1&&P==3&&D==4)||(C==1&&P==5&&D==3)||(C==1&&P==6&&D==3)||(C==1&&P==7&&D==2)||(C==1&&P==7&&D==3)";
			fitter->SetExposureMap(expoM1);
		}
		else if(fMode == "M2LowBkg")
		{
			theCut += "&&(C==2&&P==1&&D==4)||(C==2&&P==3&&D==1)||(C==2&&P==3&&D==2)||(C==2&&P==6&&D==2)||(C==2&&P==7&&D==3)";
			fitter->SetExposureMap(expoM2);
		}
		else
		{
			cout << "Cut Mode not recognized, exiting" << endl;
			return;
		}

		// Load data from TChain with a cut string
    TChain *skimTree = new TChain("skimTree");
    if(fDS == "All") {
		skimTree->Add(Form("%s/final95_DS0*.root", inDir.c_str()) );
    skimTree->Add(Form("%s/final95_DS1*.root", inDir.c_str()) );
    skimTree->Add(Form("%s/final95_DS2*.root", inDir.c_str()) );
    skimTree->Add(Form("%s/final95_DS3*.root", inDir.c_str()) );
    skimTree->Add(Form("%s/final95_DS4*.root", inDir.c_str()) );
    skimTree->Add(Form("%s/final95_DS5*.root", inDir.c_str()) );
		fitter->SetSavePrefix(Form("BasicFit_%s_%s_%.1f_%.1f", fDS.c_str(), fMode.c_str(), fitMin, fitMax));
		}
		else if(fDS == "LowBkg")
		{
			skimTree->Add(Form("%s/final95_DS1*.root", inDir.c_str()) );
			skimTree->Add(Form("%s/final95_DS2*.root", inDir.c_str()) );
			skimTree->Add(Form("%s/final95_DS3*.root", inDir.c_str()) );
			skimTree->Add(Form("%s/final95_DS4*.root", inDir.c_str()) );
			skimTree->Add(Form("%s/final95_DS5*.root", inDir.c_str()) );
			fitter->SetSavePrefix(Form("BasicFit_%s_%s_%.1f_%.1f", fDS.c_str(), fMode.c_str(), fitMin, fitMax));
		}
		// Single Datasets
		else {
			skimTree->Add(Form("%s/final95_DS%s*.root", inDir.c_str(), fDS.c_str()) );
			fitter->SetSavePrefix(Form("BasicFit_DS%s_%s_%.1f_%.1f", fDS.c_str(), fMode.c_str(), fitMin, fitMax));
		}
    fitter->LoadChainData(skimTree, theCut);

    // Construct PDF and do fit
		bool bNoEff = false; // Turns on-off efficiency

		fitter->ConstructPDF(bNoEff);
		fitter->DoFit("Minuit");
		fitter->GetFitResult()->Print("v");
		fitter->DrawBasic(0.3, true, false, false);

		// vector<string> argTest = {"Tritium", "Ge68", "Zn65", "Fe55", "Mn54"};
		// auto LimitMap = fitter->ProfileNLL(argTest);

		// fitter->GetFitResult()->Print("v");
    // fitter->GetFitResultEff()->Print("v");
		// cout << "Extended Term: " << fitter->GetWorkspace()->pdf("model")->extendedTerm(1) << endl;

/*
    vector<double> valTrit0 = fitter->GetVar("Tritium");
    double tritVal0 = valTrit0[0];
    double tritErr0 = valTrit0[1];

    vector<double> valGe0 = fitter->GetVar("Ge68");
    double Ge68Val0 = valGe0[0];
    double Ge68Err0 = valGe0[1];

		vector<double> valZn0 = fitter->GetVar("Zn65");
    double Zn65Val0 = valZn0[0];
    double Zn65Err0 = valZn0[1];

    vector<double> valFe0 = fitter->GetVar("Fe55");
    double Fe55Val0 = valFe0[0];
    double Fe55Err0 = valFe0[1];

    cout << "Basic Fit Results DS" << fDS << " (" << fMode.c_str() << ")" << endl;
    cout << "Tritium: " << tritVal0 << " +/- " << tritErr0 << endl;
    cout << "Ge68: " << Ge68Val0 << " +/- " << Ge68Err0 << endl;
		cout << "Zn65: " << Zn65Val0 << " +/- " << Zn65Err0 << endl;
    cout << "Fe55: " << Fe55Val0 << " +/- " << Fe55Err0 << endl;
*/

}
