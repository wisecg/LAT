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
#include <iostream>
#include <fstream>
// #include <Python/Python.h>

using namespace RooFit;
using namespace std;

void RunBasicFit(string fDS, double fitMin, double fitMax, string fMode, string fCPD);

int main(int argc, char** argv)
{
  	// gROOT->ProcessLine("gErrorIgnoreLevel = 3001;");
  	// gROOT->ProcessLine("RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);");

	if(argc <= 4)
	{
		cout << "Usage: " << argv[0] << " [DS (string)] [Fit Min] [Fit Max] [Cut Mode (string)] [CPD (string, optional)]" << endl;
		return 0;
	}

	vector<string> dsList = {"0", "1", "2", "3", "4", "5A", "5B", "5C", "6", "All", "LowBkg"};
	vector<string> modeList = {"AllDet", "Nat", "Enr", "M1LowBkg", "M1All", "M2LowBkg", "M2LowCosmo", "LowBkg", "M1Enr", "M2Enr", "Det"};
	vector<string> detList = {"112", "113", "114", "122", "123", "132", "133", "134", "141", "142", "143", "144", "145", "151", "152", "153", "154", "161", "162", "163", "164", "171", "172", "173", "174", "222", "223", "231", "232", "241", "242", "244", "251", "254", "262", "273"};

	string fDS = argv[1];
	float fitMin = atof(argv[2]);
	float fitMax = atof(argv[3]);
	string fMode = argv[4];
	string fCPD = "";
	if(argc > 5)
	{
		fCPD = argv[5];
	}

	// Make sure the options are available
	bool bDS = std::find(dsList.begin(), dsList.end(), fDS) != dsList.end();
	bool bMode = std::find(modeList.begin(), modeList.end(), fMode) != modeList.end();
	bool bCPD = false;

	if(fMode == "Det" && argc == 6)
	{
		bCPD = std::find(detList.begin(), detList.end(), fCPD) != detList.end();
		if(bCPD == 0)
		{
			cout << fCPD << " is not an available detector!" << endl;
			cout << "Options are: 112, 113, 114, 122, 123, 132, 133, 134, 141, 142, 143, 144, 145, 151, 152, 153, 154, 161, 162, 163, 164, 171, 172, 173, 174, 222, 223, 231, 232, 241, 242, 244, 251, 254, 262, 273" << endl;
			return 0;
		}
	}
	if(bDS == 0)
	{
		cout << fDS << " is not an available dataset option!" << endl;
		cout << "Options are: 0, 1, 2, 3, 4, 5A, 5B, 5C, 6, All, LowBkg" << endl;
		return 0;
	}
	if(bMode == 0)
	{
		cout << fMode << " is not an available mode option!" << endl;
		cout << "Options are: AllDet, Nat, Enr, LowBkg, M1LowBkg, M1All, M2LowBkg, M2LowCosmo, Det" << endl;
		return 0;
	}

	gStyle->SetOptStat(0);
	RunBasicFit(fDS, fitMin, fitMax, fMode, fCPD);

	return 0;
}

void RunBasicFit(string fDS, double fitMin, double fitMax, string fMode, string fCPD)
{
		// Calculated and saved from lat-expo.py
		// Reject C2P5D3 from good detector list in M2, it only contributes noise in DS5a, doesn't seem to exist in DS5b and is tiny in DS5c
		GPXFitter *fitter = new GPXFitter(fDS, fitMin, fitMax, fMode);

		// Goes Dataset: {Enr, Nat}
		std::map<std::string, std::vector<double>> expoFull;
		expoFull["0"] = {400.9074811876048, 157.92594939476592};
		expoFull["1"] = {631.539647097231, 29.523527363230723};
  	expoFull["2"] = {104.3238418726033, 5.177826450022776};
		expoFull["3"] = {334.9535551781866, 48.64815212750721};
		expoFull["4"] = {50.778998571351394, 23.149078616249657};
		expoFull["5A"] = {798.8446523114309, 337.4631592176636};
		expoFull["5B"] = {624.1459839987111, 226.37080226302683};
		expoFull["5C"] = {173.0972327703449, 63.354681716529285};
		expoFull["6"] = {1144.5646365776347, 389.58469957417987};
		expoFull["All"] = {4288.611020012726, 1282.2269447556757};
		// Low Bkg is Only with selection of detectors for enriched
		expoFull["LowBkg"] = {3695.6489, 1282.2269447556757};

		std::map<std::string, std::vector<double>> expoM1;
		expoM1["0"] = {400.9074811876048, 157.92594939476592};
		expoM1["1"] = {631.539647097231, 29.523527363230723};
  	expoM1["2"] = {104.3238418726033, 5.177826450022776};
		expoM1["3"] = {334.9535551781866, 48.64815212750721};
		expoM1["4"] = {0., 0.};
		expoM1["5A"] = {643.9921013269255, 190.41003908459194};
		expoM1["5B"] = {481.69736138877846, 114.86835044020026};
		expoM1["5C"] = {137.6058203337168, 30.73016914149956};
		expoM1["6"] = {865.5945009994206, 205.83109592058972};
		expoM1["All"] = {3600.6143093844676, 783.1151099224081};
		expoM1["LowBkg"] = {3199.7068281968627, 625.1891605276422};
		// expoM1["LowBkg"] = {3090.6212, 0.}; // Real exposure

		std::map<std::string, std::vector<double>> expoM2;
		expoM2["0"] = {0., 0.};
		expoM2["1"] = {0., 0.};
  	expoM2["2"] = {0., 0.};
		expoM2["3"] = {0., 0.};
		expoM2["4"] = {50.778998571351394, 23.149078616249657};
		expoM2["5A"] = {154.8525509845091, 147.05312013307218};
		expoM2["5B"] = {142.44862260993472, 111.50245182282691};
		expoM2["5C"] = {35.49141243662807, 32.62451257502965};
		expoM2["6"] = {278.9701355782069, 183.75360365359083};
		expoM2["All"] = {662.5417201806301, 498.08276680076926};
		expoM2["LowBkg"] = {662.5417201806301, 498.08276680076926}; // One that works
		// expoM2["LowBkg"] = {605.0277, 0.}; // Real exposure

		string inDir = "/Users/brianzhu/project/LATv2/bkg/cut/final95";
		string theCut = "";
    theCut += Form("tOffset<4000&&trapENFCal>=%.2f&&trapENFCal<=%.2f", fitMin, fitMax); // Energy cut for fit range

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
		else if(fMode == "M1Enr")
		{
      theCut += "&&isEnr&&C==1";
			fitter->SetExposureMap(expoM1);
		}
		else if(fMode == "M2Enr")
		{
      theCut += "&&isEnr&&C==2";
			fitter->SetExposureMap(expoM2);
		}
		else if(fMode == "AllDet")
		{
			theCut += "";
			fitter->SetExposureMap(expoFull);
		}
		else if(fMode == "LowBkg")
		{
			theCut += "&&((C==1&&P==1&&D==2)||(C==1&&P==1&&D==3)||(C==1&&P==1&&D==4)||(C==1&&P==2&&D==2)||(C==1&&P==2&&D==3)||(C==1&&P==3&&D==2)||(C==1&&P==3&&D==3)||(C==1&&P==3&&D==4)||(C==1&&P==5&&D==3)||(C==1&&P==6&&D==1)||(C==1&&P==6&&D==3)||(C==1&&P==6&&D==4)||(C==1&&P==7&&D==2)||(C==1&&P==7&&D==3)||(C==1&&P==7&&D==4)&&(C==2&&P==1&&D==4)||(C==2&&P==3&&D==1)||(C==2&&P==3&&D==2)||(C==2&&P==6&&D==2)||(C==2&&P==7&&D==3))";
			fitter->SetExposureMap(expoFull);
		}
		else if(fMode == "M1LowBkg")
		{
			theCut += "&&((C==1&&P==1&&D==2)||(C==1&&P==1&&D==3)||(C==1&&P==1&&D==4)||(C==1&&P==2&&D==2)||(C==1&&P==2&&D==3)||(C==1&&P==3&&D==2)||(C==1&&P==3&&D==3)||(C==1&&P==3&&D==4)||(C==1&&P==5&&D==3)||(C==1&&P==6&&D==1)||(C==1&&P==6&&D==3)||(C==1&&P==6&&D==4)||(C==1&&P==7&&D==2)||(C==1&&P==7&&D==3)||(C==1&&P==7&&D==4))";
			fitter->SetExposureMap(expoM1);
		}
		else if(fMode == "M1All")
		{
			theCut += "&&((C==1&&P==1&&D==2)||(C==1&&P==1&&D==3)||(C==1&&P==1&&D==4)||(C==1&&P==2&&D==2)||(C==1&&P==2&&D==3)||(C==1&&P==3&&D==4)||(C==1&&P==5&&D==3)||(C==1&&P==6&&D==3)||(C==1&&P==7&&D==2)||(C==1&&P==7&&D==3))";
			fitter->SetExposureMap(expoM1);
		}
		else if(fMode == "M2LowBkg")
		{
			theCut += "&&((C==2&&P==1&&D==4)||(C==2&&P==3&&D==1)||(C==2&&P==3&&D==2)||(C==2&&P==6&&D==2)||(C==2&&P==7&&D==3))";
			fitter->SetExposureMap(expoM2);
		}
		// Rejecting 2 detectors with higher Ge68
		else if(fMode == "M2LowCosmo")
		{
			// theCut += "&&((C==2&&P==1&&D==4)||(C==2&&P==3&&D==1)||(C==2&&P==3&&D==2))";
			theCut += "&&((C==2&&P==3&&D==1)||(C==2&&P==3&&D==2))";
			fitter->SetExposureMap(expoM2);
		}
		else if(fMode == "Det")
		{
			theCut += Form("&&(C==%c&&P==%c&&D==%c)", fCPD.c_str()[0], fCPD.c_str()[1],fCPD.c_str()[2]);
			fitter->SetExposureMap(expoFull); // This doesn't matter because it won't be used
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
		skimTree->Add(Form("%s/final95_DS6*.root", inDir.c_str()) );
		}
		else if(fDS == "LowBkg")
		{
			skimTree->Add(Form("%s/final95_DS1*.root", inDir.c_str()) );
			skimTree->Add(Form("%s/final95_DS2*.root", inDir.c_str()) );
			skimTree->Add(Form("%s/final95_DS3*.root", inDir.c_str()) );
			skimTree->Add(Form("%s/final95_DS4*.root", inDir.c_str()) );
			skimTree->Add(Form("%s/final95_DS5*.root", inDir.c_str()) );
			skimTree->Add(Form("%s/final95_DS6*.root", inDir.c_str()) );
		}
		// Single Datasets
		else {
			skimTree->Add(Form("%s/final95_DS%s*.root", inDir.c_str(), fDS.c_str()) );
		}
		if(fMode == "Det")
		{
			fitter->SetSavePrefix(Form("C%cP%cD%c_%s_%s_%.1f_%.1f", fCPD.c_str()[0], fCPD.c_str()[1],fCPD.c_str()[2], fDS.c_str(), fMode.c_str(), fitMin, fitMax));
		}
		else
		{
			fitter->SetSavePrefix(Form("Cosmo_%s_%s_%.1f_%.1f", fDS.c_str(), fMode.c_str(), fitMin, fitMax));
		}
		fitter->LoadChainData(skimTree, theCut);

		// std::string effDir = "/Users/brianzhu/macros/code/LAT/data";
		// TFile *effFile = new TFile(Form("%s/lat-expo-efficiency_final95.root", effDir.c_str()));

    // Construct PDF and do fit
		bool bNoEff = false; // Turns on/off efficiency
		bool bWFMode = false; // Turns on/off WF collapse
		fitter->ConstructPDF(bNoEff, bWFMode, fCPD);
		fitter->DoFit("Minuit");
		fitter->DrawBasic(0.3, true, false, false);
		fitter->GetFitResult()->Print("v");

		vector<string> argTest = {"Tritium", "Ge68", "Ga68", "Zn65", "Fe55", "Mn54", "Pb210", "Bkg"};
		// auto LimitMap = fitter->ProfileNLL(argTest);

		ofstream outFile;
		outFile.open(Form("/Users/brianzhu/macros/code/MJDAnalysis/Wenqin/CosmoFits/C%cP%cD%c_%s_FitOutput.csv", fCPD.c_str()[0], fCPD.c_str()[1],fCPD.c_str()[2],fDS.c_str()));
		outFile << Form("C%cP%cD%c", fCPD.c_str()[0], fCPD.c_str()[1],fCPD.c_str()[2]) << endl;
		for(auto arg: argTest)
		{
			vector<double> valTemp = fitter->GetVar(arg);
			outFile << arg << "," << valTemp[0] << "," << valTemp[1] << endl;
		}
		outFile.close();

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
