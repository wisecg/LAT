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
    fDS(1),
    fFitMin(2.),
    fFitMax(100.),
    fEnergy(nullptr),
    fRealData(nullptr),
    fMCData(nullptr),
    fMCStudy(nullptr),
    fModelPDF(nullptr),
    fNLL(nullptr),
    fNLLEff(nullptr),
    fProfileNLL(nullptr),
    fChiSquare(0),
    fMinimizer(nullptr),
    fMinimizerEff(nullptr),
    fFitResult(nullptr),
    fFitResultEff(nullptr),
    fFitWorkspace(nullptr),
    fModelPDFEff(nullptr),
    fModelPDFEffDraw(nullptr),
    fSavePrefix("FitResult")
{ }

GPXFitter::~GPXFitter()
{
	delete fEnergy;
	fEnergy = nullptr;

  delete fRealData;
	fRealData = nullptr;

	delete fModelPDF;
	fModelPDF = nullptr;

  delete fModelPDFEff;
	fModelPDFEff = nullptr;

  delete fMCData;
  fMCData = nullptr;

  delete fMCStudy;
  fMCStudy = nullptr;

  delete fMinimizer;
  fMinimizer = nullptr;

  delete fNLL;
  fNLL = nullptr;

  delete fProfileNLL;
  fProfileNLL = nullptr;

	delete fFitResult;
	fFitResult = nullptr;

  delete fFitResultEff;
	fFitResultEff = nullptr;

  delete fFitWorkspace;
  fFitWorkspace = nullptr;
}

// Constructs model PDF, use only after LoadData or else!
void GPXFitter::ConstructPDF(double enVal, bool bBDM)
{
	if(fRealData == nullptr)
	{
		std::cout << "Error: Data not loaded! Will Segfault!" << std::endl;
		return;
	}

    // RooWorkspace is necessary for the model and parameters to be persistent
    // this is necessary because we created a bunch of objects that aren't persistent here
    fFitWorkspace = new RooWorkspace("fFitWorkspace", "Fit Workspace");

    // Energy shift
    double fDeltaE = 0.00;

    // std::string tritDir = "/mnt/mjdDisk1/Majorana/users/psz/CUORE/MJDAnalysis/Wenqin/Data";
    std::string tritDir = "/Users/brianzhu/macros/code/MJDAnalysis/Axion";
    TFile *tritFileBulk = new TFile(Form("%s/TritSpec.root", tritDir.c_str()));
    TFile *tritFile = new TFile(Form("%s/TritSpecSim.root", tritDir.c_str()));
    TH1D *tritSpec = dynamic_cast<TH1D*>(tritFileBulk->Get("tritHist"));
    TH1D *tritSpecSurface = dynamic_cast<TH1D*>(tritFile->Get("tritHistSurface"));

    RooDataHist tritRooHist("tritBulk", "Tritium Histogram (Bulk)", *fEnergy, Import(*tritSpec));
    RooDataHist tritSurfaceRooHist("tritSurface", "Tritium Histogram (Transition Layer)", *fEnergy, Import(*tritSpecSurface));
    // Because Steve's histogram sucks
    // The range of the histogram is maxed out at 50 keV, so need to reset range after loading histogram
    fEnergy->setRange(fFitMin, fFitMax);
    RooHistPdf tritPdf("tritPdf", "TritiumPdf", *fEnergy, tritRooHist, 1);
    RooHistPdf tritSurfacePdf("tritSurfacePdf", "TritiumSurfacePdf", *fEnergy, tritSurfaceRooHist, 1);

    TFile *axionFile = new TFile(Form("%s/axionHistos.root", tritDir.c_str()));
    TH1D *axionSpec = dynamic_cast<TH1D*>(axionFile->Get("hConv"));
    RooDataHist axionRooHist("axion", "Axion Histogram", *fEnergy, Import(*axionSpec));
    fEnergy->setRange(fFitMin, fFitMax);
    RooHistPdf axionPdf("axionPdf", "AxionPdf", *fEnergy, axionRooHist, 1);

    // std::string effDir = "/mnt/mjdDisk1/Majorana/users/psz/CUORE/MJDAnalysis/Wenqin/Data";
    // std::string effDir = "/Users/brianzhu/macros/code/LAT/plots/spectra/PrelimSpectra";
    std::string effDir = "/Users/brianzhu/macros/code/LAT/plots/AThresh";
    TFile *effFile;
    TH1D *effSpec;
    if(fCutString.find("isNat") != std::string::npos) {
      if(fDS == -1){
        effFile = new TFile(Form("%s/Bkg_isNat_DS5.root", effDir.c_str()));
        effSpec = dynamic_cast<TH1D*>(effFile->Get("DS5_isNat_EffTot"));
      }
      else {
        effFile = new TFile(Form("%s/Bkg_isNat_DS%d.root", effDir.c_str(), fDS));
        effSpec = dynamic_cast<TH1D*>(effFile->Get(Form("DS%d_isNat_EffTot", fDS)));
      }
    }
    else if(fCutString.find("isEnr") != std::string::npos) {
      if(fDS == -1){
        effFile = new TFile(Form("%s/Bkg_isEnr_DS5.root", effDir.c_str()));
        effSpec = dynamic_cast<TH1D*>(effFile->Get("DS5_isEnr_EffTot"));
      }
      else {
        effFile = new TFile(Form("%s/Bkg_isEnr_DS%d.root", effDir.c_str(), fDS));
        effSpec = dynamic_cast<TH1D*>(effFile->Get(Form("DS%d_isEnr_EffTot", fDS)));
      }
    }
    else {
      effFile = new TFile(Form("%s/Bkg_isEnr_DS%d.root", effDir.c_str(), fDS));
      effSpec = dynamic_cast<TH1D*>(effFile->Get(Form("DS%d_isEnr_EffTot", fDS)));
    }

    // Empirical total efficiency that fits DS1 natural pretty well
    // RooFormulaVar effPdf("effPdf","0.5*(1+TMath::Erf((@fEnergy-2)/(TMath::Sqrt(2)*3) ))", *fEnergy);
    // TF1 *fEff1 = new TF1("fEff", "0.5*(1+TMath::Erf((x-[0])/(TMath::Sqrt(2)*[1]) ))", fFitMin, fFitMax);
    // fEff1->SetParameters(2, 3);
    // effSpec = dynamic_cast<TH1D*>(fEff1->CreateHistogram());

    // External pulsers from Clint
    TF1 *fEff1 = new TF1("fEff", "0.5*(1+TMath::Erf((x-[0])/(TMath::Sqrt(2)*[1]) ))", fFitMin, fFitMax);
    fEff1->SetParameters(0.36,1.26);

    effSpec->Multiply(fEff1);
    RooDataHist effRooHist("eff", "Efficiency Histogram", *fEnergy, Import(*effSpec));
    fEnergy->setRange(fFitMin, fFitMax);
    RooHistPdf effPdf("effPdf", "EffPdf", *fEnergy, effRooHist, 0);


    RooRealVar DeltaE("DeltaE", "DeltaE", 0., -2.0, 2.0);
    // Change this for polynomial background or linear
    // RooRealVar polySlope("polySlope", "Background Slope", 0.00002, -0.2, 0.2);
    // RooArgList polyList(polySlope);
    // RooPolynomial BkgPoly("Background", "Linear Background function", *fEnergy, polyList);
    RooPolynomial BkgPoly("Background", "Flat Background function", *fEnergy, RooArgList());

    RooRealVar Ge68L_mean("Ge68L_mean", "Ge68L_mean", 1.3 + fDeltaE);
    RooRealVar Ge68L_sigma("Ge68L_sigma", "Ge68L_sigma", GetSigma(1.3 + fDeltaE));
    RooGaussian Ge68L_gauss("Ge68L_gauss", "Ge68L Gaussian", *fEnergy, Ge68L_mean, Ge68L_sigma);

    RooRealVar V49_mean("V49_mean", "V49_mean", 4.97 + fDeltaE);
    RooRealVar V49_sigma("V49_sigma", "V49_sigma", GetSigma(4.97 + fDeltaE));
    RooGaussian V49_gauss("V49_gauss", "V49 Gaussian", *fEnergy, V49_mean, V49_sigma);

    RooRealVar Cr51_mean("Cr51_mean", "Cr51_mean", 5.46 + fDeltaE);
    RooRealVar Cr51_sigma("Cr51_sigma", "Cr51_sigma", GetSigma(5.46 + fDeltaE));
    RooGaussian Cr51_gauss("Cr51_gauss", "Cr51 Gaussian", *fEnergy, Cr51_mean, Cr51_sigma);

    RooRealVar Mn54_mean("Mn54_mean", "Mn54_mean", 5.99 + fDeltaE);
    RooRealVar Mn54_sigma("Mn54_sigma", "Mn54_sigma", GetSigma(5.99 + fDeltaE));
    RooGaussian Mn54_gauss("Mn54_gauss", "Mn54 Gaussian", *fEnergy, Mn54_mean, Mn54_sigma);

    RooRealVar Fe55_mean("Fe55_mean", "Fe55_mean", 6.54 + fDeltaE);
    RooRealVar Fe55_sigma("Fe55_sigma", "Fe55_sigma", GetSigma(6.54 + fDeltaE));
    RooGaussian Fe55_gauss("Fe55_gauss", "Fe55 Gaussian", *fEnergy, Fe55_mean, Fe55_sigma);

    RooRealVar Co57_mean("Co57_mean", "Co57_mean", 7.11 + fDeltaE);
    RooRealVar Co57_sigma("Co57_sigma", "Co57_sigma", GetSigma(7.11 + fDeltaE));
    RooGaussian Co57_gauss("Co57_gauss", "Co57 Gaussian", *fEnergy, Co57_mean, Co57_sigma);

    RooRealVar Zn65_mean("Zn65_mean", "Zn65_mean", 8.98 + fDeltaE);
    RooRealVar Zn65_sigma("Zn65_sigma", "Zn65_sigma", GetSigma(8.98 + fDeltaE));
    RooGaussian Zn65_gauss("Zn65_gauss", "Zn65 Gaussian", *fEnergy, Zn65_mean, Zn65_sigma);

    RooRealVar Ga68_mean("Ga68_mean", "Ga68_mean", 9.66 + fDeltaE);
    RooRealVar Ga68_sigma("Ga68_sigma", "Ga68_sigma", GetSigma(9.66 + fDeltaE));
    RooGaussian Ga68_gauss("Ga68_gauss", "Ga68 Gaussian", *fEnergy, Ga68_mean, Ga68_sigma);

    RooRealVar Ge68_mean("Ge68_mean", "Ge68_mean", 10.37 + fDeltaE);
    // RooRealVar Ge68_mean("Ge68_mean", "Ge68_mean", 10.37 + fDeltaE, 10.2, 10.45); // Floating mean
    RooRealVar Ge68_sigma("Ge68_sigma", "Ge68_sigma", GetSigma(10.37 + fDeltaE));
    RooGaussian Ge68_gauss("Ge68_gauss", "Ge68 Gaussian", *fEnergy, Ge68_mean, Ge68_sigma);

    RooRealVar As73_mean("As73_mean", "As73_mean", 11.3 + fDeltaE);
    RooRealVar As73_sigma("As73_sigma", "As73_sigma", GetSigma(11.3 + fDeltaE));
    RooGaussian As73_gauss("As73_gauss", "As73 Gaussian", *fEnergy, As73_mean, As73_sigma);

    RooRealVar Pb210_mean("Pb210_mean", "Pb210_mean", 46.54 + fDeltaE);
    RooRealVar Pb210_sigma("Pb210_sigma", "Pb210_sigma", GetSigma(46.54 + fDeltaE));
    RooGaussian Pb210_gauss("Pb210_gauss", "Pb210 Gaussian", *fEnergy, Pb210_mean, Pb210_sigma);

    // Normalization parameters
    // Make names pretty for plots
    RooRealVar num_trit("Tritium", "Tritium", 100.0, 0.0, 5000.);
    RooRealVar num_tritSurface("Surface Tritium", "Surface Tritium", 0.0, 0.0, 500.);
    RooRealVar num_axion("Axion", "Axion", 0.0, 0.0, 1000.);
    RooRealVar num_bkg("Bkg", "Bkg", 50.0, 0.0, 5000.);
    RooRealVar num_V49("V49", "V49", 0.0, 0.0, 500.);
    RooRealVar num_Cr51("Cr51", "Cr51", 0.0, 0.0, 500.);
    RooRealVar num_Mn54("Mn54", "Mn54", 0.0, 0.0, 500.);
    RooRealVar num_Fe55("Fe55", "Fe55", 3.0, 0.0, 500.);
    RooRealVar num_Co57("Co57", "Co57", 0.0, 0.0, 500.);
    RooRealVar num_Zn65("Zn65", "Zn65", 0.0, 0.0, 500.);
    RooRealVar num_Ga68("Ga68", "Ga68", 5.0, 0.0, 500.);
    RooRealVar num_Ge68("Ge68", "Ge68", 10.0, 0.0, 500.);
    RooRealVar num_As73("As73", "As73", 0.0, 0.0, 500.);
    RooRealVar num_Ge68L("Ge68L", "Ge68L", 0.0, 0.0, 500.);
    RooRealVar num_Pb210("Pb210", "Pb210", 5.0, 0.0, 500.);

    // RooRealVar num_Tot("Tot", "Tot", 5.0, 0.0, 100000.);

    // Extended PDF model -- use this to create an extended model
    RooExtendPdf tritPdfe("tritPdfe", "Extended Tritium", tritPdf, num_trit);
    RooExtendPdf tritSurfacePdfe("tritSurfacePdfe", "Extended Tritium Surface", tritSurfacePdf, num_tritSurface);
    RooExtendPdf axionPdfe("axionPdfe", "Extended axion", axionPdf, num_axion);
    RooExtendPdf BkgPolye("BkgPolye", "Extended BkgPoly", BkgPoly, num_bkg);
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
    RooExtendPdf Pb210_gausse("Pb210_gausse", "Extended Pb210_gauss", Pb210_gauss, num_Pb210);

    RooArgList shapes(tritPdfe, BkgPolye, Mn54_gausse, Fe55_gausse, Zn65_gausse, Ge68_gausse);
    shapes.add(Ga68_gausse);
    // shapes.add(tritSurfacePdfe);
    if(fFitMin <= 1.0) shapes.add(Ge68L_gausse);
    if(fFitMax > 48) shapes.add(Pb210_gausse);

    // Additional PDFs we don't use, but we can use
    // shapes.add(V49_gausse);
    // shapes.add(Cr51_gausse);
    // shapes.add(As73_gausse);

    RooAddPdf model("model", "total pdf", shapes);

    // Combine with efficiency -- Several methods here
    // RooProdPdf modelEff("modelEff", "model with efficiency", RooArgSet(model, effPdf)); // Can do with product also?
    // RooProdPdf modelEff("modelEff", "model with efficiency", model, effPdf);
    RooEffProd modelEff("modelEff", "model with efficiency", model, effPdf); // This looks pretty good
    // fModelPDFFinal = new RooEffProd("modelEffFinal", "model with efficiency 2", model, effPdf);
    // RooExtendPdf modelEffe("modelEffe", "Extended model with efficiency", modelEff, num_Tot);

    // Add model to workspace -- also adds all of the normalization constants
    // If this step isn't done, a lot of the later functions won`'t work!
    fFitWorkspace->import(RooArgSet(modelEff));
    fModelPDFEff = fFitWorkspace->pdf("modelEff");
    fModelPDF = fFitWorkspace->pdf("model");
    // fModelPDFEffFit = dynamic_cast<RooAbsPdf*>(modelEff2);
    // fModelPDFEffDraw = fFitWorkspace->pdf("modelEff");
}

void GPXFitter::DoFit(std::string Minimizer)
{
    // Create NLL (This is not a profile! When you draw it to one axis, it's just a projection!)
    fNLL = fFitWorkspace->pdf("model")->createNLL(*fRealData, Extended(), NumCPU(2));
    // fNLL = fModelPDF->createNLL(*fRealData, Extended(), NumCPU(4));

    // Create minimizer, fit model to data and save result
    fMinimizer = new RooMinimizer(*fNLL);
    fMinimizer->setMinimizerType(Form("%s", Minimizer.c_str()));
    fMinimizer->setPrintLevel(-1);
    fMinimizer->setStrategy(2);
    fMinimizer->migrad();
    fMinimizer->hesse();
    // Use all these, fk if I know they're any good
    fMinimizer->improve();
    fMinimizer->minos();

    fFitResult = fMinimizer->save();
    fFitWorkspace->import(*fFitResult);
}

void GPXFitter::DoFitEff(std::string Minimizer)
{
    // Create NLL (This is not a profile! When you draw it to one axis, it's just a projection!)
    // fNLLEff = fModelPDFEff->createNLL(*fRealData, Extended(), NumCPU(4));
    fNLLEff = fModelPDFEff->createNLL(*fRealData, NumCPU(2));

    // Create minimizer, fit model to data and save result
    fMinimizerEff = new RooMinimizer(*fNLLEff);
    fMinimizerEff->setMinimizerType(Form("%s", Minimizer.c_str()));
    fMinimizerEff->setPrintLevel(-1);
    fMinimizerEff->setStrategy(2);
    fMinimizerEff->migrad();
    fMinimizerEff->hesse();
    // Use all these, fk if I know they're any good
    fMinimizerEff->improve();
    fMinimizerEff->minos();

    fFitResultEff = fMinimizerEff->save();
    fFitWorkspace->import(*fFitResultEff);
}

void GPXFitter::DrawBasicShit(double binSize, bool drawLabels, bool drawResid, bool drawMatrix)
{

    // RooEffProd modelEffDraw("modelEffDraw", "model with efficiency", *fFitWorkspace->pdf("model"), *fFitWorkspace->pdf("effPdf"));

    TCanvas *cSpec = new TCanvas("cSpec", "cSpec", 1100, 800);
    RooPlot* frameFit = fEnergy->frame(Range(fFitMin, fFitMax), Bins((fFitMax - fFitMin)*1.0/binSize + 0.5));
    fRealData->plotOn(frameFit);
    // fModelPDF->plotOn(frameFit, Components("axionPdfe"), LineColor(kBlue), LineStyle(kDashed));
    fModelPDF->plotOn(frameFit, LineColor(kBlue));
    // modelEffDraw.plotOn(frameFit, LineColor(kRed));
    fModelPDFEff->plotOn(frameFit, LineColor(kRed));
    frameFit->SetTitle("");

    // Get parameter values from first fit... these methods suck
    // std::string tritDir = "/mnt/mjdDisk1/Majorana/users/psz/CUORE/MJDAnalysis/Wenqin/Data";
    std::string tritDir = "/Users/brianzhu/macros/code/MJDAnalysis/Axion";

    TFile *tritFile = new TFile(Form("%s/TritSpec.root", tritDir.c_str()));
    TH1D *tritSpec = dynamic_cast<TH1D*>(tritFile->Get("tritHist"));

    double tritVal = dynamic_cast<RooRealVar*>(fFitResult->floatParsFinal().find("Tritium") )->getValV();
    double tritErr = dynamic_cast<RooRealVar*>(fFitResult->floatParsFinal().find("Tritium") )->getError();
    double tritValCorr = tritVal/(tritSpec->Integral(tritSpec->FindBin(fFitMin), tritSpec->FindBin(fFitMax)));
    double tritErrCorr = tritErr/(tritSpec->Integral(tritSpec->FindBin(fFitMin), tritSpec->FindBin(fFitMax)));
    double geVal = dynamic_cast<RooRealVar*>(fFitResult->floatParsFinal().find("Ge68"))->getValV();
    double geErr = dynamic_cast<RooRealVar*>(fFitResult->floatParsFinal().find("Ge68"))->getError();
    double znVal = dynamic_cast<RooRealVar*>(fFitResult->floatParsFinal().find("Zn65"))->getValV();
    double znErr = dynamic_cast<RooRealVar*>(fFitResult->floatParsFinal().find("Zn65"))->getError();
    double feVal = dynamic_cast<RooRealVar*>(fFitResult->floatParsFinal().find("Fe55"))->getValV();
    double feErr = dynamic_cast<RooRealVar*>(fFitResult->floatParsFinal().find("Fe55"))->getError();
    double bkgVal = dynamic_cast<RooRealVar*>(fFitResult->floatParsFinal().find("Bkg"))->getValV();
    double bkgErr = dynamic_cast<RooRealVar*>(fFitResult->floatParsFinal().find("Bkg"))->getError();
    // double gelVal = dynamic_cast<RooRealVar*>(fFitResult->floatParsFinal().find("Ge68L"))->getValV();
    // double gelErr = dynamic_cast<RooRealVar*>(fFitResult->floatParsFinal().find("Ge68L"))->getError();

    // double tritVal = dynamic_cast<RooRealVar*>(fFitResultEff->floatParsFinal().find("Tritium") )->getValV();
    // double tritErr = dynamic_cast<RooRealVar*>(fFitResultEff->floatParsFinal().find("Tritium") )->getError();
    // double tritValCorr = tritVal/(tritSpec->Integral(tritSpec->FindBin(fFitMin), tritSpec->FindBin(fFitMax)));
    // double tritErrCorr = tritErr/(tritSpec->Integral(tritSpec->FindBin(fFitMin), tritSpec->FindBin(fFitMax)));
    // double geVal = dynamic_cast<RooRealVar*>(fFitResultEff->floatParsFinal().find("Ge68"))->getValV();
    // double geErr = dynamic_cast<RooRealVar*>(fFitResultEff->floatParsFinal().find("Ge68"))->getError();
    // double feVal = dynamic_cast<RooRealVar*>(fFitResultEff->floatParsFinal().find("Fe55"))->getValV();
    // double feErr = dynamic_cast<RooRealVar*>(fFitResultEff->floatParsFinal().find("Fe55"))->getError();
    // double bkgVal = dynamic_cast<RooRealVar*>(fFitResultEff->floatParsFinal().find("Bkg"))->getValV();
    // double bkgErr = dynamic_cast<RooRealVar*>(fFitResultEff->floatParsFinal().find("Bkg"))->getError();


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
      leg->AddText("Model without Efficiency function");
      dynamic_cast<TText*>(leg->GetListOfLines()->Last())->SetTextColor(kBlue);
      // leg->AddText("Model with Efficiency function");
      // dynamic_cast<TText*>(leg->GetListOfLines()->Last())->SetTextColor(kRed);
      leg->AddText(Form("Tritium (Uncorrected): %.3f #pm %.3f", tritVal, tritErr));
      leg->AddText(Form("Tritium (Corrected): %.3f #pm %.3f", tritValCorr, tritErrCorr));
      leg->AddText(Form("Tritium (2-4 keV): %.3f #pm %.3f", tritValCorr*0.224487, tritErrCorr*0.224487));
      leg->AddText(Form("Ge68K: %.3f #pm %.3f", geVal, geErr));
      // leg->AddText(Form("Ge68L: %.3f #pm %.3f", gelVal, gelErr));
      leg->AddText(Form("Zn65: %.3f #pm %.3f", znVal, znErr));
      leg->AddText(Form("Fe55: %.3f #pm %.3f", feVal, feErr));
      leg->AddText(Form("Bkg: %.3f #pm %.3f", bkgVal, bkgErr));
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
  fModelPDFEff->plotOn(frameFit, LineColor(kBlue));
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


    // Study done for Jason to convince him I know what I'm doing
    // double parVal0 = dynamic_cast<RooRealVar*>(fFitResult->floatParsFinal().find(Form("%s", argN.c_str())))->getValV();
    // double parErrHi0 = dynamic_cast<RooRealVar*>(fFitResult->floatParsFinal().find(Form("%s", argN.c_str())))->getErrorHi();
    // double parErrLo0 = dynamic_cast<RooRealVar*>(fFitResult->floatParsFinal().find(Form("%s", argN.c_str())))->getErrorLo();

    // TH1D *hMean = new TH1D("hMean", "Tritium Mean", 200, parVal0+2.5*parErrLo0, parVal0+2.5*parErrHi0);
    // TH1D *hErrLo = new TH1D("hErrLo", "Tritium Error Low", 200, parErrLo0+parErrLo0/4, parErrLo0+parErrHi0/4);
    // TH1D *hErrHi = new TH1D("hErrHi", "Tritium Error High", 200, parErrHi0+parErrLo0/4, parErrHi0+parErrHi0/4);

    // for(int i = 0; i < nMC; i++)
    // {
    //     const RooFitResult *fFitMCResult = fMCStudy->fitResult(i);
    //     double parVal2 = dynamic_cast<RooRealVar*>(fFitMCResult->floatParsFinal().find(Form("%s", argN.c_str())))->getValV();
    //     double parErrHi2 = dynamic_cast<RooRealVar*>(fFitMCResult->floatParsFinal().find(Form("%s", argN.c_str())))->getErrorHi();
    //     double parErrLo2 = dynamic_cast<RooRealVar*>(fFitMCResult->floatParsFinal().find(Form("%s", argN.c_str())))->getErrorLo();

    //     hMean->Fill(parVal2);
    //     hErrLo->Fill(parErrLo2);
    //     hErrHi->Fill(parErrHi2);
    // }

    // TCanvas *cM = new TCanvas("cM", "cM", 800, 600);
    // hMean->Draw();
    // l1.DrawLine(parVal0, 0, parVal0, hMean->GetBinContent(hMean->GetMaximumBin()) );
    // cM->SaveAs("./LATv2Result/MeanTest.pdf");

    // TCanvas *cLo = new TCanvas("cLo", "cLo", 800, 600);
    // hErrLo->Draw();
    // l1.DrawLine(parErrLo0, 0, parErrLo0, hErrLo->GetBinContent(hErrLo->GetMaximumBin()));
    // cLo->SaveAs("./LATv2Result/MeanTest_Low.pdf");

    // TCanvas *cHi = new TCanvas("cHi", "cHi", 800, 600);
    // hErrHi->Draw();
    // l1.DrawLine(parErrHi0, 0, parErrHi0, hErrHi->GetBinContent(hErrHi->GetMaximumBin()));
    // cHi->SaveAs("./LATv2Result/MeanTest_High.pdf");

    // double NLLmean = dynamic_cast<RooRealVar*>(fFitResult->floatParsFinal().find("NLL"))->getValV();
    // double NLLmean = dynamic_cast<RooRealVar*>(fFitResult->floatParsFinal().find("NLL"))->getError();

/*
    // Extract nLL variables from MCStudy -- maybe useful later
    RooDataSet MCFitData = fMCStudy->fitParDataSet();
    // RooDataSet *NLL = static_cast<RooDataSet*>(MCFitData.reduce(RooArgSet()));
    MCFitData.Print("v");
    const RooArgSet* row = MCFitData.get();
    row->Print("v");
    std::shared_ptr<TCanvas> cMCNLL( std::make_shared<TCanvas>("cMCNLL", "cMCNLL", 1100, 800) );
    RooRealVar *NLL = static_cast<RooRealVar*>(row->find("NLL"));
    RooFormulaVar sNLL("sNLL", "Shifted NLL", Form("NLL - %f", fFitResult->minNll()), *NLL);
    RooRealVar *SNLL = static_cast<RooRealVar*>(MCFitData.addColumn(sNLL));

    RooPlot* frameMC = SNLL->frame(Range((0.2*fFitResult->minNll()), -(0.2*fFitResult->minNll()) ));
    MCFitData.plotOn(frameMC);
    frameMC->Draw();
    cMCNLL->SaveAs("Test.pdf");
*/

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

// Returns best fit value as a vector with asymmetric errors {Value, ErrorHi, ErrorLo}
std::vector<double> GPXFitter::GetVarEff(std::string argN)
{
    double Val = dynamic_cast<RooRealVar*>(fFitResultEff->floatParsFinal().find( Form("%s", argN.c_str()) ))->getValV();
    double ErrHi = dynamic_cast<RooRealVar*>(fFitResultEff->floatParsFinal().find( Form("%s", argN.c_str()) ))->getErrorHi();
    double ErrLo = dynamic_cast<RooRealVar*>(fFitResultEff->floatParsFinal().find( Form("%s", argN.c_str()) ))->getErrorLo();

	std::vector<double> VarResult = {Val, ErrHi, ErrLo};
    return VarResult;
}


// Gets resolution, function and parameters from BDM PRL paper
// https://arxiv.org/abs/1612.00886
// In the future it should just be a convolution with all the other PDFs?
double GPXFitter::GetSigma(double energy)
{
    double p0, p1, p2;
    if(fDS==0) {p0 = 0.147; p1=0.0173; p2=0.0003;}
    else if(fDS==1) {p0 = 0.136; p1=0.0174; p2=0.00028;}
    else if(fDS==2) {p0 = 0.143; p1=0.0172; p2=0.000284;}
    else if(fDS==3) {p0 = 0.162; p1=0.0172; p2=0.000297;}
    else if(fDS==4) {p0 = 0.218; p1=0.015; p2=0.00035;}
    // else if(fDS==5) {p0 = 0.2121; p1=0.01838; p2=0.00031137;}
    // else if(fDS==5) {p0 = 0.2592; p1=0.2057; p2=0.00030863;} // DS5b
    else if(fDS==5) {p0 = 0.18148; p1=0.01690; p2=0.00031873;} // DS5b
    else {p0 = 0.2121; p1=0.01838; p2=0.00031137;} // Use DS5 numbers for any other DS
    double sig = std::sqrt(p0*p0 + p1*p1*energy + p2*p2*energy*energy );
	return sig;
}

// Assumes standard skim format -- converts stuff from vector<double> to scalar
// Also assumes trapENFCal is the energy parameter of choice
void GPXFitter::LoadChainData(TChain *skimTree, std::string theCut)
{
    fCutString = theCut;
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
        // Draw Profile NLL and save as PDF
        TCanvas *cNLL = new TCanvas("cNLL", "cNLL", 900, 600);

        // Best fit value -- just using this to set range
        double parVal = dynamic_cast<RooRealVar*>(fFitResultEff->floatParsFinal().find(Form("%s", argN.c_str())))->getValV();

        // RooStats::ProfileLikelihoodCalculator plc(*fRealData, *fModelPDF, RooArgSet(*fFitWorkspace->var(Form("%s", argN.c_str()))) );
        RooStats::ProfileLikelihoodCalculator plc(*fRealData, *fModelPDFEff, RooArgSet(*fFitWorkspace->var(Form("%s", argN.c_str()))) );
        // Set 1 sigma interval
        plc.SetConfidenceLevel(CL);

        RooStats::LikelihoodInterval *interval = plc.GetInterval();
        double lowerLimit = interval->LowerLimit(*fFitWorkspace->var(Form("%s", argN.c_str())));
        double upperLimit = interval->UpperLimit(*fFitWorkspace->var(Form("%s", argN.c_str())));

        RooStats::LikelihoodIntervalPlot plot(interval);
        plot.SetRange(parVal - 1.5*(parVal - lowerLimit), parVal + 1.5*(upperLimit-parVal) );
        plot.Draw();
        cNLL->SaveAs(Form("./LATv2Result/%s_%sNLL.pdf", fSavePrefix.c_str(), argN.c_str()) );
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
