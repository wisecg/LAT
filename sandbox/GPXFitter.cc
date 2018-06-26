#include "GPXFitter.hh"
#include "RooAddPdf.h"
#include "RooHistPdf.h"
#include "RooConstVar.h"
#include "RooDataSet.h"
#include "RooFitResult.h"

#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooArgList.h"
#include "RooArgSet.h"
#include "RooMinimizer.h"
#include "RooPlot.h"
#include "RooPolynomial.h"
#include "RooGaussian.h"
#include "RooRealVar.h"
#include "RooExtendPdf.h"
#include "RooWorkspace.h"
#include "RooDataHist.h"
#include "RooMCStudy.h"
#include "RooHist.h"
#include "RooEffProd.h"
#include "RooEfficiency.h"
#include "RooCategory.h"
#include "RooProdPdf.h"

#include "RooStats/ProfileLikelihoodCalculator.h"
#include "RooStats/LikelihoodInterval.h"
#include "RooStats/LikelihoodIntervalPlot.h"

#include "TAxis.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TChain.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TString.h"
#include "TLine.h"
#include "TTree.h"
#include "TStyle.h"
#include "TEntryList.h"
#include "TBranch.h"
#include "TPaveText.h"
#include "TText.h"
#include "TLegend.h"
#include <iostream>

using namespace RooFit;
using namespace std;

GPXFitter::GPXFitter() :
    fDS("LowBkg"),
    fFitMin(2.),
    fFitMax(50.),
    fMode("Enr"),
    fEnergy(nullptr),
    fRealData(nullptr),
    fMCData(nullptr),
    fMCStudy(nullptr),
    fModelPDF(nullptr),
    fNLL(nullptr),
    fProfileNLL(nullptr),
    fEffSpec(nullptr),
    fChiSquare(0),
    fMinimizer(nullptr),
    fFitResult(nullptr),
    fFitWorkspace(nullptr),
    fSavePrefix("FitResult")
{
  // Set default exposure as full exposure
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
  expoFull["LowBkg"] = {2717.683911799859, 733.6872277542301}; // Low Bkg is DS1-5
  fExposureMap = expoFull;

}

GPXFitter::~GPXFitter()
{
	delete fEnergy;
	fEnergy = nullptr;

  delete fRealData;
	fRealData = nullptr;

	delete fModelPDF;
	fModelPDF = nullptr;

  delete fMCData;
  fMCData = nullptr;

  delete fMCStudy;
  fMCStudy = nullptr;

  delete fEffSpec;
  fEffSpec = nullptr;

  delete fMinimizer;
  fMinimizer = nullptr;

  delete fNLL;
  fNLL = nullptr;

  delete fProfileNLL;
  fProfileNLL = nullptr;

	delete fFitResult;
	fFitResult = nullptr;

  delete fFitWorkspace;
  fFitWorkspace = nullptr;
}

// Constructs model PDF, use only after LoadData or else!
void GPXFitter::ConstructPDF(bool bNoEff)
{
  if (bNoEff) fSavePrefix += "NoEff";

  if(fRealData == nullptr)
	{
		std::cout << "Error: Data not loaded! Will Segfault!" << std::endl;
		return;
	}
    // RooWorkspace is necessary for the model and parameters to be persistent
    // this is necessary because we created a bunch of objects that aren't persistent here
    fFitWorkspace = new RooWorkspace("fFitWorkspace", "Fit Workspace");

    // Energy shift
    RooRealVar fDeltaE("fDeltaE", "fDeltaE", 0., -0.2, 0.2);
    fEnergyShift = new RooFormulaVar("fEnergyShift", "@0-@1", RooArgList(*fEnergy, fDeltaE));

    std::string tritDir = "/Users/brianzhu/macros/code/MJDAnalysis/Axion";
    TFile *tritFileBulk = new TFile(Form("%s/TritSpec_0.1keV.root", tritDir.c_str()));
    TH1D *tritSpec = dynamic_cast<TH1D*>(tritFileBulk->Get("tritHist"));

    TFile *axionFile = new TFile(Form("%s/axionHistos.root", tritDir.c_str()));
    TH1D *axionSpec = dynamic_cast<TH1D*>(axionFile->Get("hConv"));
    RooDataHist axionRooHist("axion", "Axion Histogram", *fEnergy, Import(*axionSpec));
    fEnergy->setRange(fFitMin, fFitMax);
    RooHistPdf axionPdf("axionPdf", "AxionPdf", *fEnergy, axionRooHist, 1);

    TFile *gausFile = new TFile(Form("%s/GausCosmoPDFs.root", tritDir.c_str()));
    TH1D *Ge68LSpec = dynamic_cast<TH1D*>(gausFile->Get("hGe68L"));
    TH1D *V49Spec = dynamic_cast<TH1D*>(gausFile->Get("hV49"));
    TH1D *Cr51Spec = dynamic_cast<TH1D*>(gausFile->Get("hCr51"));
    TH1D *Mn54Spec = dynamic_cast<TH1D*>(gausFile->Get("hMn54"));
    TH1D *Fe55Spec = dynamic_cast<TH1D*>(gausFile->Get("hFe55"));
    TH1D *Co57Spec = dynamic_cast<TH1D*>(gausFile->Get("hCo57"));
    TH1D *Zn65Spec = dynamic_cast<TH1D*>(gausFile->Get("hZn65"));
    TH1D *Ga68Spec = dynamic_cast<TH1D*>(gausFile->Get("hGa68"));
    TH1D *Ge68Spec = dynamic_cast<TH1D*>(gausFile->Get("hGe68"));
    TH1D *As73Spec = dynamic_cast<TH1D*>(gausFile->Get("hAs73"));
    TH1D *Pb210Spec = dynamic_cast<TH1D*>(gausFile->Get("hPb210"));

    TFile *pb210File = new TFile(Form("%s/Pb210PDFs_Full.root", tritDir.c_str()));
    TH1D *pb210Spec = dynamic_cast<TH1D*>(pb210File->Get("hPb210"));

    // Copy the pb210 (to get the same binning) for the flat background
    TH1D *bkgSpec = dynamic_cast<TH1D*>(pb210Spec->Clone("hBkg"));

    std::string effDir = "/Users/brianzhu/macros/code/LAT/data";
    TFile *effFile = new TFile(Form("%s/lat-expo-efficiency_final95.root", effDir.c_str()));

    if(fDS == "All" || fDS == "LowBkg")
    {
      if(fMode == "All")
      {
        fEffSpec = dynamic_cast<TH1D*>(effFile->Get("hDS1_Nat"));
        if(fDS == "All")
        {
          fEffSpec->Add(dynamic_cast<TH1D*>(effFile->Get("hDS0_Nat")));
          fEffSpec->Add(dynamic_cast<TH1D*>(effFile->Get("hDS0_Enr")));
        }
        fEffSpec->Add(dynamic_cast<TH1D*>(effFile->Get("hDS2_Nat")));
        fEffSpec->Add(dynamic_cast<TH1D*>(effFile->Get("hDS3_Nat")));
        fEffSpec->Add(dynamic_cast<TH1D*>(effFile->Get("hDS4_Nat")));
        fEffSpec->Add(dynamic_cast<TH1D*>(effFile->Get("hDS5A_Nat")));
        fEffSpec->Add(dynamic_cast<TH1D*>(effFile->Get("hDS5B_Nat")));
        fEffSpec->Add(dynamic_cast<TH1D*>(effFile->Get("hDS5C_Nat")));
        fEffSpec->Add(dynamic_cast<TH1D*>(effFile->Get("hDS1_Enr")));
        fEffSpec->Add(dynamic_cast<TH1D*>(effFile->Get("hDS2_Enr")));
        fEffSpec->Add(dynamic_cast<TH1D*>(effFile->Get("hDS3_Enr")));
        fEffSpec->Add(dynamic_cast<TH1D*>(effFile->Get("hDS4_Enr")));
        fEffSpec->Add(dynamic_cast<TH1D*>(effFile->Get("hDS5A_Enr")));
        fEffSpec->Add(dynamic_cast<TH1D*>(effFile->Get("hDS5B_Enr")));
        fEffSpec->Add(dynamic_cast<TH1D*>(effFile->Get("hDS5C_Enr")));
        fEffSpec->Scale(1./(fExposureMap[fDS][0] + fExposureMap[fDS][1]));
      }
      else if(fMode == "Nat" || fMode == "Enr")
      {
        fEffSpec = dynamic_cast<TH1D*>(effFile->Get(Form("hDS1_%s", fMode.c_str())));
        if(fDS == "All") fEffSpec->Add(dynamic_cast<TH1D*>(effFile->Get(Form("hDS0_%s", fMode.c_str()))));
        fEffSpec->Add(dynamic_cast<TH1D*>(effFile->Get(Form("hDS2_%s", fMode.c_str()))));
        fEffSpec->Add(dynamic_cast<TH1D*>(effFile->Get(Form("hDS3_%s", fMode.c_str()))));
        fEffSpec->Add(dynamic_cast<TH1D*>(effFile->Get(Form("hDS4_%s", fMode.c_str()))));
        fEffSpec->Add(dynamic_cast<TH1D*>(effFile->Get(Form("hDS5A_%s", fMode.c_str()))));
        fEffSpec->Add(dynamic_cast<TH1D*>(effFile->Get(Form("hDS5B_%s", fMode.c_str()))));
        fEffSpec->Add(dynamic_cast<TH1D*>(effFile->Get(Form("hDS5C_%s", fMode.c_str()))));
        if(fMode == "Enr") fEffSpec->Scale(1./(fExposureMap[fDS][0]));
        else if(fMode == "Nat") fEffSpec->Scale(1./(fExposureMap[fDS][1]));
      }
      else if (fMode == "M1All" || fMode == "M1LowBkg")
      {
        fEffSpec = dynamic_cast<TH1D*>(effFile->Get("hDS1_Enr_M1"));
        if(fMode == "M1LowBkg")fEffSpec->Add(dynamic_cast<TH1D*>(effFile->Get("hDS0_Enr_M1")));
        fEffSpec->Add(dynamic_cast<TH1D*>(effFile->Get("hDS2_Enr_M1")));
        fEffSpec->Add(dynamic_cast<TH1D*>(effFile->Get("hDS3_Enr_M1")));
        fEffSpec->Add(dynamic_cast<TH1D*>(effFile->Get("hDS5A_Enr_M1")));
        fEffSpec->Add(dynamic_cast<TH1D*>(effFile->Get("hDS5B_Enr_M1")));
        fEffSpec->Add(dynamic_cast<TH1D*>(effFile->Get("hDS5C_Enr_M1")));
        fEffSpec->Scale(1./(fExposureMap[fDS][0]));
      }
      else if (fMode == "M2LowBkg")
      {
        fEffSpec = dynamic_cast<TH1D*>(effFile->Get("hDS4_Enr_M2"));
        fEffSpec->Add(dynamic_cast<TH1D*>(effFile->Get("hDS5A_Enr_M2")));
        fEffSpec->Add(dynamic_cast<TH1D*>(effFile->Get("hDS5B_Enr_M2")));
        fEffSpec->Add(dynamic_cast<TH1D*>(effFile->Get("hDS5C_Enr_M2")));
        fEffSpec->Scale(1./(fExposureMap[fDS][0]));
      }
    }
    // Specific Dataset
    else
    {
      if (fMode == "All")
      {
        fEffSpec = dynamic_cast<TH1D*>(effFile->Get(Form("hDS%s_Enr", fDS.c_str())));
        fEffSpec = dynamic_cast<TH1D*>(effFile->Get(Form("hDS%s_Nat", fDS.c_str())));
        fEffSpec->Scale(1./(fExposureMap[fDS][0] + fExposureMap[fDS][1]));
      }
      else if (fMode == "Enr" || fMode == "Nat")
      {
        fEffSpec = dynamic_cast<TH1D*>(effFile->Get(Form("hDS%s_%s", fDS.c_str(), fMode.c_str() )));
        if(fMode == "Enr") fEffSpec->Scale(1./(fExposureMap[fDS][0]));
        else if(fMode == "Nat") fEffSpec->Scale(1./(fExposureMap[fDS][1]));
      }
    }

    // Pre-scale pdfs with energy dependence with Efficiency
    double xVal = 0;
    double effVal = 0;
    for(int i = 1; i < pb210Spec->GetNbinsX(); i++)
    {
      xVal = pb210Spec->GetBinCenter(i);
      effVal = fEffSpec->GetBinContent(fEffSpec->FindBin(xVal));

      // If this flag is activated, all PDFs are generated with no efficiency function!
      if (bNoEff) effVal = 1.;

      tritSpec->SetBinContent(i, tritSpec->GetBinContent(i)*effVal);
      pb210Spec->SetBinContent(i, pb210Spec->GetBinContent(i)*effVal);
      bkgSpec->SetBinContent(i, 1./((fFitMax-fFitMin)/0.1)*effVal);

      Ge68LSpec->SetBinContent(i, Ge68LSpec->GetBinContent(i)*effVal);
      V49Spec->SetBinContent(i, V49Spec->GetBinContent(i)*effVal);
      Cr51Spec->SetBinContent(i, Cr51Spec->GetBinContent(i)*effVal);
      Mn54Spec->SetBinContent(i, Mn54Spec->GetBinContent(i)*effVal);
      Fe55Spec->SetBinContent(i, Fe55Spec->GetBinContent(i)*effVal);
      Co57Spec->SetBinContent(i, Co57Spec->GetBinContent(i)*effVal);
      Zn65Spec->SetBinContent(i, Zn65Spec->GetBinContent(i)*effVal);
      Ga68Spec->SetBinContent(i, Ga68Spec->GetBinContent(i)*effVal);
      Ge68Spec->SetBinContent(i, Ge68Spec->GetBinContent(i)*effVal);
      As73Spec->SetBinContent(i, As73Spec->GetBinContent(i)*effVal);
      Pb210Spec->SetBinContent(i, Pb210Spec->GetBinContent(i)*effVal);
    }

    RooDataHist tritRooHist("tritBulk", "Tritium Histogram (Bulk)", *fEnergy, Import(*tritSpec));
    RooDataHist pb210RooHist("pb210", "Pb210 Histogram", *fEnergy, Import(*pb210Spec));
    RooDataHist bkgRooHist("bkg", "Background Histogram", *fEnergy, Import(*bkgSpec));

    RooDataHist Ge68LRooHist("Ge68L", "Ge68L Histogram", *fEnergy, Import(*Ge68LSpec));
    RooDataHist V49RooHist("V49", "V49 Histogram", *fEnergy, Import(*V49Spec));
    RooDataHist Cr51RooHist("Cr51", "Cr51 Histogram", *fEnergy, Import(*Cr51Spec));
    RooDataHist Mn54RooHist("Mn54", "Mn54 Histogram", *fEnergy, Import(*Mn54Spec));
    RooDataHist Fe55RooHist("Fe55", "Fe55 Histogram", *fEnergy, Import(*Fe55Spec));
    RooDataHist Co57RooHist("Co57", "Co57 Histogram", *fEnergy, Import(*Co57Spec));
    RooDataHist Zn65RooHist("Zn65", "Zn65 Histogram", *fEnergy, Import(*Zn65Spec));
    RooDataHist Ga68RooHist("Ga68", "Ga68 Histogram", *fEnergy, Import(*Ga68Spec));
    RooDataHist Ge68RooHist("Ge68", "Ge68 Histogram", *fEnergy, Import(*Ge68Spec));
    RooDataHist As73RooHist("As73", "As73 Histogram", *fEnergy, Import(*As73Spec));
    RooDataHist Pb210RooHist("Pb210", "Pb210 Histogram", *fEnergy, Import(*Pb210Spec));

    RooHistPdf tritPdf("tritPdf", "TritiumPdf", *fEnergy, tritRooHist, 1);
    RooHistPdf pb210Pdf("pb210Pdf", "Pb210Pdf", *fEnergy, pb210RooHist, 1);
    RooHistPdf bkgPdf("bkgPdf", "BkgPdf", *fEnergy, bkgRooHist, 1);

    RooHistPdf Ge68L_gauss("Ge68L_gauss", "Ge68L_gauss", *fEnergy, Ge68LRooHist, 1);
    RooHistPdf V49_gauss("V49_gauss", "V49_gauss", *fEnergy, V49RooHist, 1);
    RooHistPdf Cr51_gauss("Cr51_gauss", "Cr51_gauss", *fEnergy, Cr51RooHist, 1);
    RooHistPdf Mn54_gauss("Mn54_gauss", "Mn54_gauss", *fEnergy, Mn54RooHist, 1);
    RooHistPdf Fe55_gauss("Fe55_gauss", "Fe55_gauss", *fEnergy, Fe55RooHist, 1);
    RooHistPdf Co57_gauss("Co57_gauss", "Co57_gauss", *fEnergy, Co57RooHist, 1);
    RooHistPdf Zn65_gauss("Zn65_gauss", "Zn65_gauss", *fEnergy, Zn65RooHist, 1);
    RooHistPdf Ga68_gauss("Ga68_gauss", "Ga68_gauss", *fEnergy, Ga68RooHist, 1);
    RooHistPdf Ge68_gauss("Ge68_gauss", "Ge68_gauss", *fEnergy, Ge68RooHist, 1);
    RooHistPdf As73_gauss("As73_gauss", "As73_gauss", *fEnergy, As73RooHist, 1);
    RooHistPdf Pb210_gauss("Pb210_gauss", "Pb210_gauss", *fEnergy, Pb210RooHist, 1);

    // Normalization parameters
    // Make names pretty for plots
    RooRealVar num_trit("Tritium", "Tritium", 1000.0, 0.0, 50000.);
    RooRealVar num_axion("Axion", "Axion", 0.0, 0.0, 1000.);
    RooRealVar num_bkg("Bkg", "Bkg", 1500, 0.0, 10000.);
    RooRealVar num_V49("V49", "V49", 0.0, 0.0, 1000.);
    RooRealVar num_Cr51("Cr51", "Cr51", 0.0, 0.0, 1000.);
    RooRealVar num_Mn54("Mn54", "Mn54", 0.0, 0.0, 1000.);
    RooRealVar num_Fe55("Fe55", "Fe55", 3.0, 0.0, 1000.);
    RooRealVar num_Co57("Co57", "Co57", 0.0, 0.0, 1000.);
    RooRealVar num_Zn65("Zn65", "Zn65", 10.0, 0.0, 1000.);
    RooRealVar num_Ga68("Ga68", "Ga68", 1.0, 0.0, 1000.);
    RooRealVar num_Ge68("Ge68", "Ge68", 50.0, 0.0, 1000.);
    RooRealVar num_As73("As73", "As73", 0.0, 0.0, 1000.);
    RooRealVar num_Ge68L("Ge68L", "Ge68L", 0.0, 0.0, 10.);
    RooRealVar num_Pb210("Pb210", "Pb210", 300.0, 0.0, 5000.);

    // Extended PDF model -- use this to create an extended model
    RooExtendPdf tritPdfe("tritPdfe", "Extended Tritium", tritPdf, num_trit);
    RooExtendPdf axionPdfe("axionPdfe", "Extended axion", axionPdf, num_axion);
    RooExtendPdf BkgPolye("BkgPolye", "Extended BkgPoly", bkgPdf, num_bkg);
    RooExtendPdf V49_gausse("V49_gausse", "Extended V49_gauss", V49_gauss, num_V49);
    RooExtendPdf Cr51_gausse("Cr51_gausse", "Extended Cr51_gauss", Cr51_gauss, num_Cr51);
    RooExtendPdf Mn54_gausse("Mn54_gausse", "Extended Mn54_gauss", Mn54_gauss, num_Mn54);
    RooExtendPdf Fe55_gausse("Fe55_gausse", "Extended Fe55_gauss", Fe55_gauss, num_Fe55);
    RooExtendPdf Co57_gausse("Co57_gausse", "Extended Co57_gauss", Co57_gauss, num_Co57);
    RooExtendPdf Zn65_gausse("Zn65_gausse", "Extended Zn65_gauss", Zn65_gauss, num_Zn65);
    RooExtendPdf Ga68_gausse("Ga68_gausse", "Extended Ga68_gauss", Ga68_gauss, num_Ga68);
    RooExtendPdf Ge68_gausse("Ge68_gausse", "Extended Ge68_gauss", Ge68_gauss, num_Ge68);
    RooExtendPdf As73_gausse("As73_gausse", "Extended As73_gauss", As73_gauss, num_As73);
    RooExtendPdf Ge68L_gausse("Ge68L_gausse", "Extended Ge68L_gauss", Ge68L_gauss, num_Ge68L);
    // Using Gaussian
    RooExtendPdf Pb210_gausse("Pb210_gausse", "Extended Pb210_gauss", Pb210_gauss, num_Pb210);
    // Using MaGe PDF
    // RooExtendPdf Pb210_gausse("Pb210_gausse", "Extended Pb210_gauss", pb210Pdf, num_Pb210);

    RooArgList shapes(tritPdfe, BkgPolye, Mn54_gausse, Fe55_gausse, Zn65_gausse, Ge68_gausse);
    shapes.add(Ga68_gausse);
    if(fFitMin < 1.0) shapes.add(Ge68L_gausse);
    if(fFitMax > 48) shapes.add(Pb210_gausse);
    // Additional PDFs we don't use, but we can use
    // shapes.add(V49_gausse);
    // shapes.add(Cr51_gausse);
    // shapes.add(As73_gausse);

    RooAddPdf model("model", "total pdf", shapes);

    // Add model to workspace -- also adds all of the normalization constants
    // If this step isn't done, a lot of the later functions won`'t work!
    fFitWorkspace->import(RooArgSet(model));
    fModelPDF = fFitWorkspace->pdf("model");
    std::cout << "Generated PDFs" << std::endl;
}

void GPXFitter::DoFit(std::string Minimizer)
{
    // Create NLL (This is not a profile! When you draw it to one axis, it's just a projection!)
    fNLL = fFitWorkspace->pdf("model")->createNLL(*fRealData, Extended(), NumCPU(2));

    // Create minimizer, fit model to data and save result
    fMinimizer = new RooMinimizer(*fNLL);
    fMinimizer->setMinimizerType(Form("%s", Minimizer.c_str()));
    fMinimizer->setPrintLevel(-1);
    fMinimizer->setStrategy(2);
    fMinimizer->migrad();
    fMinimizer->hesse();
    fMinimizer->improve();
    fMinimizer->minos();

    fFitResult = fMinimizer->save();
    fFitWorkspace->import(*fFitResult);
}

void GPXFitter::DrawBasic(double binSize, bool drawLabels, bool drawResid, bool drawMatrix)
{

    TCanvas *cSpec = new TCanvas("cSpec", "cSpec", 1100, 800);
    // cSpec->SetLogy();
    RooPlot* frameFit = fEnergy->frame(Range(fFitMin, fFitMax), Bins((fFitMax - fFitMin)*1.0/binSize + 0.5));
    fRealData->plotOn(frameFit);
    fModelPDF->plotOn(frameFit, LineColor(kBlue));
    fModelPDF->plotOn(frameFit, Components("tritPdfe"), LineColor(kRed), LineStyle(kDashed));
    fModelPDF->plotOn(frameFit, Components("BkgPolye"), LineColor(kGreen+1), LineStyle(kDashed));
    fModelPDF->plotOn(frameFit, Components("Ge68_gausse"), LineColor(kMagenta), LineStyle(kDashed));
    fModelPDF->plotOn(frameFit, Components("Mn54_gausse"), LineColor(kMagenta+1), LineStyle(kDashed));
    fModelPDF->plotOn(frameFit, Components("Fe55_gausse"), LineColor(kYellow+1), LineStyle(kDashed));
    fModelPDF->plotOn(frameFit, Components("Zn65_gausse"), LineColor(kMagenta-1), LineStyle(kDashed));
    fModelPDF->plotOn(frameFit, Components("Pb210_gausse"), LineColor(kBlack), LineStyle(kDashed));
    frameFit->SetTitle("");

    // Get parameter values from first fit... these methods suck
    std::string tritDir = "/Users/brianzhu/macros/code/MJDAnalysis/Axion";

    TFile *tritFile = new TFile(Form("%s/TritSpec_0.1keV.root", tritDir.c_str()));
    TH1D *tritSpec = dynamic_cast<TH1D*>(tritFile->Get("tritHist"));
    tritSpec->Scale(1./tritSpec->Integral());

    double tritVal = dynamic_cast<RooRealVar*>(fFitResult->floatParsFinal().find("Tritium") )->getValV();
    double tritErr = dynamic_cast<RooRealVar*>(fFitResult->floatParsFinal().find("Tritium") )->getError();
    double tritValCorr = tritVal/(tritSpec->Integral(tritSpec->FindBin(fFitMin), tritSpec->FindBin(fFitMax)));
    double tritErrCorr = tritErr/(tritSpec->Integral(tritSpec->FindBin(fFitMin), tritSpec->FindBin(fFitMax)));
    double geVal = dynamic_cast<RooRealVar*>(fFitResult->floatParsFinal().find("Ge68"))->getValV();
    double geErr = dynamic_cast<RooRealVar*>(fFitResult->floatParsFinal().find("Ge68"))->getError();
    double gaVal = dynamic_cast<RooRealVar*>(fFitResult->floatParsFinal().find("Ga68"))->getValV();
    double gaErr = dynamic_cast<RooRealVar*>(fFitResult->floatParsFinal().find("Ga68"))->getError();
    double znVal = dynamic_cast<RooRealVar*>(fFitResult->floatParsFinal().find("Zn65"))->getValV();
    double znErr = dynamic_cast<RooRealVar*>(fFitResult->floatParsFinal().find("Zn65"))->getError();
    double feVal = dynamic_cast<RooRealVar*>(fFitResult->floatParsFinal().find("Fe55"))->getValV();
    double feErr = dynamic_cast<RooRealVar*>(fFitResult->floatParsFinal().find("Fe55"))->getError();
    double mnVal = dynamic_cast<RooRealVar*>(fFitResult->floatParsFinal().find("Mn54"))->getValV();
    double mnErr = dynamic_cast<RooRealVar*>(fFitResult->floatParsFinal().find("Mn54"))->getError();
    double pbVal = dynamic_cast<RooRealVar*>(fFitResult->floatParsFinal().find("Pb210"))->getValV();
    double pbErr = dynamic_cast<RooRealVar*>(fFitResult->floatParsFinal().find("Pb210"))->getError();
    double bkgVal = dynamic_cast<RooRealVar*>(fFitResult->floatParsFinal().find("Bkg"))->getValV();
    double bkgErr = dynamic_cast<RooRealVar*>(fFitResult->floatParsFinal().find("Bkg"))->getError();

    // Add chi-square to the plot - it's fairly meaningless as it's an unbinned fit but people will want it
    fChiSquare = frameFit->chiSquare(8);
    if(drawLabels)
    {
      TPaveText *leg = new TPaveText(0.50, 0.55, 0.88, .88, "NDC");
      leg->SetTextFont(133);
      leg->SetFillColor(0);
      leg->SetBorderSize(0);
      leg->SetTextSize(22);
      leg->AddText(Form("#chi^{2}/NDF = %.3f" ,fChiSquare ) );
      leg->AddText("Total Model");
      dynamic_cast<TText*>(leg->GetListOfLines()->Last())->SetTextColor(kBlue);
      leg->AddText(Form("Tritium (Uncorrected): %.3f #pm %.3f", tritVal, tritErr));
      dynamic_cast<TText*>(leg->GetListOfLines()->Last())->SetTextColor(kRed);
      leg->AddText(Form("Tritium (Corrected): %.3f #pm %.3f", tritValCorr, tritErrCorr));
      dynamic_cast<TText*>(leg->GetListOfLines()->Last())->SetTextColor(kRed);
      // leg->AddText(Form("Tritium (2-4 keV): %.3f #pm %.3f", tritValCorr*0.224487, tritErrCorr*0.224487));
      leg->AddText(Form("Ge68: %.3f #pm %.3f", geVal, geErr));
      dynamic_cast<TText*>(leg->GetListOfLines()->Last())->SetTextColor(kMagenta);
      leg->AddText(Form("Ga68: %.3f #pm %.3f", gaVal, gaErr));
      dynamic_cast<TText*>(leg->GetListOfLines()->Last())->SetTextColor(kBlue+1);
      leg->AddText(Form("Zn65: %.3f #pm %.3f", znVal, znErr));
      dynamic_cast<TText*>(leg->GetListOfLines()->Last())->SetTextColor(kMagenta-1);
      leg->AddText(Form("Fe55: %.3f #pm %.3f", feVal, feErr));
      dynamic_cast<TText*>(leg->GetListOfLines()->Last())->SetTextColor(kYellow+1);
      leg->AddText(Form("Mn54: %.3f #pm %.3f", mnVal, mnErr));
      dynamic_cast<TText*>(leg->GetListOfLines()->Last())->SetTextColor(kMagenta+1);
      leg->AddText(Form("Pb210: %.3f #pm %.3f", pbVal, pbErr));
      dynamic_cast<TText*>(leg->GetListOfLines()->Last())->SetTextColor(kBlack);
      leg->AddText(Form("Bkg: %.3f #pm %.3f", bkgVal, bkgErr));
      dynamic_cast<TText*>(leg->GetListOfLines()->Last())->SetTextColor(kGreen+1);
      frameFit->addObject(leg);
    }

    frameFit->Draw();
    cSpec->SaveAs(Form("./LATv2Result/%s_Spec.pdf", fSavePrefix.c_str()) );
    fFitWorkspace->import(*frameFit);

    if(drawMatrix)
    {
        TCanvas *cMatrix = new TCanvas("cMatrix", "cMatrix", 1100, 800);
        TH2D *fCorrMatrix = dynamic_cast<TH2D*>(fFitResult->correlationHist("Correlation Matrix"));
        fCorrMatrix->Draw("colz");
        cMatrix->SaveAs(Form("./LATv2Result/%s_CorrMatrix.pdf", fSavePrefix.c_str()) );
        fFitWorkspace->import(*fCorrMatrix);
    }

    if(drawResid)
    {
        TCanvas *cResidual = new TCanvas("cResidual", "cResidual", 1100, 800);
        // RooHist *hresid = frameFit->pullHist();
        RooHist *hresid = frameFit->residHist();
        RooPlot *frameResid = fEnergy->frame(Title("Fit Residuals"));
        frameResid->addPlotable(hresid, "P");
        frameResid->GetYaxis()->SetTitle("Residuals");
        frameResid->Draw();
        cResidual->SaveAs(Form("./LATv2Result/%s_Residual.pdf", fSavePrefix.c_str()) );
    }
}

void GPXFitter::DrawModels(double binSize)
{
  TCanvas *cModels = new TCanvas("cModels", "cModels", 1100, 800);
  RooPlot* frameFit = fEnergy->frame(Range(fFitMin, fFitMax), Bins((fFitMax - fFitMin)*1.0/binSize + 0.5));
  fFitWorkspace->pdf("model")->plotOn(frameFit, LineColor(kRed), LineStyle(kDashed));
  fModelPDF->plotOn(frameFit, LineColor(kBlue));
  frameFit->SetTitle("");
  frameFit->Draw();
  cModels->SaveAs(Form("./LATv2Result/%s_Models.pdf", fSavePrefix.c_str()) );
}

void GPXFitter::DrawContour(std::string argN1, std::string argN2)
{
    TCanvas *cContour = new TCanvas("cContour", "cContour", 1100, 800);
    RooPlot *frameContour = fMinimizer->contour( *fFitWorkspace->var(Form("%s", argN1.c_str())), *fFitWorkspace->var(Form("%s", argN2.c_str())), 1, 2, 3);
    frameContour->SetTitle(Form("Contour of %s vs %s", argN2.c_str(), argN1.c_str()) );

    // Get range for plot -- make stuff pretty
    double meanx = fFitWorkspace->var(Form("%s", argN1.c_str()))->getValV();
    double uperrx = fFitWorkspace->var(Form("%s", argN1.c_str()))->getErrorHi();
    double lowerrx = fFitWorkspace->var(Form("%s", argN1.c_str()))->getErrorLo(); // This value is negative!
    double meany = fFitWorkspace->var(Form("%s", argN2.c_str()))->getValV();
    double uperry = fFitWorkspace->var(Form("%s", argN2.c_str()))->getErrorHi();
    double lowerry = fFitWorkspace->var(Form("%s", argN2.c_str()))->getErrorLo(); // This value is negative!

    frameContour->GetXaxis()->SetRangeUser(meanx + 4*lowerrx, meanx + 4*uperrx);
    frameContour->GetYaxis()->SetRangeUser(meany + 4*lowerry, meany + 4*uperry);
    frameContour->Draw();

    // Save plots into workspace and pdf
    cContour->SaveAs(Form("./LATv2Result/%s_Contour_%svs%s.pdf", fSavePrefix.c_str(), argN2.c_str(), argN1.c_str()));
    fFitWorkspace->import(*frameContour);
}


// Use after constructing the model and minimization!
// This is a simple MC study only generating parameter information and pulls, the toy MC data can be saved
// I forget a lot of the crap I did here, a lot of it was for boring tests Jason and Reyco wanted
void GPXFitter::GenerateMCStudy(std::vector<std::string> argS, int nMC)
{
    // Right now I'm saving the fit output
    fMCStudy = new RooMCStudy(*fModelPDF, *fEnergy, Extended(), Silence(), FitOptions(Save()) );
    fMCStudy->generateAndFit(nMC);

    for(auto &argN: argS)
    {
        // Get parameter values from first fit... these methods suck but we have to use them
        double parVal = dynamic_cast<RooRealVar*>(fFitResult->floatParsFinal().find(Form("%s", argN.c_str())))->getValV();
        double parErr = dynamic_cast<RooRealVar*>(fFitResult->floatParsFinal().find(Form("%s", argN.c_str())))->getError();
        // Make test plots
        TCanvas *cMCStudy = new TCanvas("cMCStudy", "cMCStudy", 1100, 800);
        RooPlot *frame1 = fMCStudy->plotParam( *fFitWorkspace->var(Form("%s", argN.c_str())), Bins(50) );
        RooPlot *frame2 = fMCStudy->plotError( *fFitWorkspace->var(Form("%s", argN.c_str())), FrameRange(parErr-0.5*parErr, parErr+0.5*parErr), Bins(50));
        RooPlot *frame3 = fMCStudy->plotPull( *fFitWorkspace->var(Form("%s", argN.c_str())), FrameRange(-5, 5), Bins(50));
        RooPlot *frame4 = fMCStudy->plotNLL(Bins(50));

        // Add PaveTexts with values and such to make things pretty
        TPaveText *legParam = new TPaveText(0.60, 0.78, 0.89, 0.89, "NDC");
        legParam->SetTextFont(133);
        legParam->SetFillColor(0);
        legParam->SetBorderSize(1);
        legParam->SetTextSize(14);
        legParam->AddText(Form("Best Fit: %.3f #pm %.3f", parVal, parErr));
        frame1->addObject(legParam);

        // Workaround because fitting built into plotPull is terrible...
        // Get the histogram from the frame and then fit it myself
        RooHist *hist = frame3->getHist();
        hist->Fit("gaus", "ME");
        TF1 *gaus = hist->GetFunction("gaus");
        TPaveText *legpull = new TPaveText(0.60, 0.75, 0.89, 0.89, "NDC");
        legpull->SetTextFont(133);
        legpull->SetFillColor(0);
        legpull->SetBorderSize(1);
        legpull->SetTextSize(14);
        legpull->AddText(Form("Pull Mean: %.3f #pm %.3f", gaus->GetParameter(1), gaus->GetParError(1)) );
        legpull->AddText(Form("Pull Sigma: %.3f #pm %.3f", gaus->GetParameter(2), gaus->GetParError(2)) );
        frame3->addObject(legpull);

        TPaveText *legNLL = new TPaveText(0.60, 0.78, 0.89, 0.89, "NDC");
        legNLL->SetTextFont(133);
        legNLL->SetFillColor(0);
        legNLL->SetBorderSize(1);
        legNLL->SetTextSize(14);
        legNLL->AddText(Form("Best Fit NLL: %.3f", fFitResult->minNll()));
        frame4->addObject(legNLL);

        // Draw pretty lines
        TLine l1;
        l1.SetLineColor(kBlue);
        l1.SetLineWidth(2);
        l1.SetLineStyle(3);

        cMCStudy->Divide(2,2);
        cMCStudy->cd(1); gPad->SetLeftMargin(0.15); frame1->GetYaxis()->SetTitleOffset(1.4); frame1->Draw();
        // Draw a line at best fit position
        l1.DrawLine(parVal, frame1->GetMinimum(), parVal, frame1->GetMaximum());
        cMCStudy->cd(2); gPad->SetLeftMargin(0.15); frame2->GetYaxis()->SetTitleOffset(1.4); frame2->Draw();
        l1.DrawLine(parErr, frame2->GetMinimum(), parErr, frame2->GetMaximum());
        cMCStudy->cd(3); gPad->SetLeftMargin(0.15); frame3->GetYaxis()->SetTitleOffset(1.4); frame3->Draw();
        cMCStudy->cd(4); gPad->SetLeftMargin(0.15); frame4->GetYaxis()->SetTitleOffset(1.4); frame4->Draw();
        // Draw a line at minimum NLL position
        l1.DrawLine(fFitResult->minNll(), frame4->GetMinimum(), fFitResult->minNll(), frame4->GetMaximum());

        // Save MC Study to plot and workspace
        // fFitWorkspace->import(*fMCStudy);
        cMCStudy->SaveAs(Form("./LATv2Result/%s_%s_MCStudy.pdf", fSavePrefix.c_str(), argN.c_str()) );
    }
}

// Use after constructing the model and minimization!
// How useful is this when there's RooMCStudy?
void GPXFitter::GenerateToyMC(std::string fileName)
{
    TFile *fOut = new TFile( Form("./Data/%s_%s.root", fSavePrefix.c_str(), fileName.c_str()), "RECREATE" );
    fMCData = fModelPDF->generate( RooArgSet(*fEnergy), Name("Toy_dataset"), Extended());

    // Save data to workspace and write to a file
    RooWorkspace *workspace = new RooWorkspace("workspace", "workspace");
    workspace->import(*fMCData);
    fOut->cd();
    workspace->Write();
    fOut->Close();
}

// Returns best fit value as a vector with asymmetric errors {Value, ErrorHi, ErrorLo}
std::vector<double> GPXFitter::GetVar(std::string argN)
{
    double Val = dynamic_cast<RooRealVar*>(fFitResult->floatParsFinal().find( Form("%s", argN.c_str()) ))->getValV();
    double ErrHi = dynamic_cast<RooRealVar*>(fFitResult->floatParsFinal().find( Form("%s", argN.c_str()) ))->getErrorHi();
    double ErrLo = dynamic_cast<RooRealVar*>(fFitResult->floatParsFinal().find( Form("%s", argN.c_str()) ))->getErrorLo();

	std::vector<double> VarResult = {Val, ErrHi, ErrLo};
    return VarResult;
}

// Gets resolution, function and parameters from BDM PRL paper
// https://arxiv.org/abs/1612.00886
// In the future it should just be a convolution with all the other PDFs?
double GPXFitter::GetSigma(double energy)
{
    double p0, p1, p2;
    if(fDS=="0") {p0 = 0.147; p1=0.0173; p2=0.0003;}
    else if(fDS=="1") {p0 = 0.136; p1=0.0174; p2=0.00028;}
    else if(fDS=="2") {p0 = 0.143; p1=0.0172; p2=0.000284;}
    else if(fDS=="3") {p0 = 0.162; p1=0.0172; p2=0.000297;}
    else if(fDS=="4") {p0 = 0.218; p1=0.015; p2=0.00035;}
    // else if(fDS==5) {p0 = 0.2121; p1=0.01838; p2=0.00031137;}
    // else if(fDS==5) {p0 = 0.2592; p1=0.2057; p2=0.00030863;} // DS5b
    else if(fDS=="5") {p0 = 0.18148; p1=0.01690; p2=0.00031873;} // DS5b
    else {p0 = 0.2121; p1=0.01838; p2=0.00031137;} // Use DS5 numbers for any other DS
    double sig = std::sqrt(p0*p0 + p1*p1*energy + p2*p2*energy*energy );
	return sig;
}

// Assumes standard skim format -- converts stuff from vector<double> to scalar
// Also assumes trapENFCal is the energy parameter of choice
void GPXFitter::LoadChainData(TChain *skimTree, std::string theCut)
{
    fCutString = theCut;
    std::cout << Form("Found %lld entries before cuts", skimTree->GetEntries()) << std::endl;
    // First get TEntryList with TCut
    skimTree->Draw(">> elist", Form("%s", theCut.c_str() ), "entrylist goff");
    TEntryList *elist = dynamic_cast<TEntryList*>(gDirectory->Get("elist"));
    skimTree->SetEntryList(&*elist); // This works
    std::cout << Form("Using cut: %s", theCut.c_str() ) << std::endl;
    std::cout << Form("Found %lld entries passing cuts", elist->GetN()) << std::endl;

    // I found it easier to work like this rather than with a TTreeReader...
    std::vector<double> *ftrapENFCal = nullptr;
    std::vector<int> *fchannel = nullptr;
    skimTree->SetBranchAddress("trapENFCal", &ftrapENFCal);
    skimTree->SetBranchAddress("channel", &fchannel);

    // Create and fill a dummy TTree to load into the RooDataSet
    // I've only saved energy and channel so far... there's probably more useful parameters to keep around
    double trapENFCal = 0;
    int channel = 0;
    int treeNum = 0;
    TTree *dummyTree = new TTree("dummyTree", "Tree for RooDataSet");
    dummyTree->Branch("trapENFCal", &trapENFCal, "trapENFCal/D");
    dummyTree->Branch("channel", &channel, "channel/I");
    for(int i = 0; i < elist->GetN(); i++)
    {
        int treeEntry = elist->GetEntryAndTree(i, treeNum);
        skimTree->GetEntry( treeEntry + skimTree->GetTreeOffset()[treeNum] );
        if(i%5000==0) std::cout << "Processing event: " << i << std::endl;

        for(int j = 0; j < ftrapENFCal->size(); j++)
        {
            trapENFCal = ftrapENFCal->at(j);
            channel = fchannel->at(j);
            dummyTree->Fill();
        }
    }
    std::cout << "Filled Entries = " << dummyTree->GetEntries() << std::endl;

    // Can and perhaps should split the data up by channel in a more complicated fit
    fEnergy = new RooRealVar("trapENFCal", "trapENFCal", fFitMin, fFitMax, "keV");
    fRealData = new RooDataSet("data", "data", dummyTree, RooArgSet(*fEnergy));
}

// Implemented now in RooStats rather than RooFit
// Calculates profile likelihood and spits out limits
std::map<std::string, std::vector<double>> GPXFitter::ProfileNLL(std::vector<std::string> argS, double CL)
{
    std::map<std::string, std::vector<double>> LimitMap;
    for(auto &argN : argS)
    {
        cout << "Profiling " << argN << endl;
        // Draw Profile NLL and save as PDF
        TCanvas *cNLL = new TCanvas("cNLL", "cNLL", 900, 600);

        // Best fit value -- just using this to set range
        double parVal = dynamic_cast<RooRealVar*>(fFitResult->floatParsFinal().find(Form("%s", argN.c_str())))->getValV();

        RooStats::ProfileLikelihoodCalculator plc(*fRealData, *fModelPDF, RooArgSet(*fFitWorkspace->var(Form("%s", argN.c_str()))) );
        // Set 1 sigma interval
        plc.SetConfidenceLevel(CL);

        RooStats::LikelihoodInterval *interval = plc.GetInterval();
        double lowerLimit = interval->LowerLimit(*fFitWorkspace->var(Form("%s", argN.c_str())));
        double upperLimit = interval->UpperLimit(*fFitWorkspace->var(Form("%s", argN.c_str())));

        RooStats::LikelihoodIntervalPlot plot(interval);
        plot.SetRange(parVal - 1.5*(parVal - lowerLimit), parVal + 1.5*(upperLimit-parVal) );
        plot.Draw();
        cNLL->SaveAs(Form("./LATv2Result/%s_%sNLL.C", fSavePrefix.c_str(), argN.c_str()) );
        std::vector<double> Limits = {lowerLimit, upperLimit};
        LimitMap[argN.c_str()] = Limits;
    }

    return LimitMap;
}

void GPXFitter::SaveShit(std::string outfileName)
{
    TFile *fOut = new TFile( Form("./LATv2Result/%s_%s", fSavePrefix.c_str(), outfileName.c_str()), "RECREATE" );
    fOut->cd();
    fFitWorkspace->Write();
    fOut->Close();
}

void GPXFitter::SetFitRange(double fitMin, double fitMax)
{
    fFitMin = fitMin;
    fFitMax = fitMax;
}

// Must be done AFTER parameter is loaded
void GPXFitter::SetParameter(std::string arg, double val, bool fix)
{
    fFitWorkspace->var(Form("%s", arg.c_str()))->setVal(val);
    if(fix)fFitWorkspace->var(Form("%s", arg.c_str()))->setConstant(kTRUE);
    else fFitWorkspace->var(Form("%s", arg.c_str()))->setConstant(kFALSE);
}
