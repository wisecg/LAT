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
#include <iostream>

using namespace std;
using namespace RooFit;

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
  fProfileNLL(nullptr),
  fChiSquare(0),
  fMinimizer(nullptr),
  fFitResult(nullptr),
  fFitWorkspace(nullptr),
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

  delete fFitWorkspace;
  fFitWorkspace = nullptr;
}

// Constructs model PDF, use only after LoadData or else!
void GPXFitter::ConstructPDF(double enVal, bool bBDM)
{
  if(fRealData == nullptr) {
    cout << "Error: Data not loaded! Will Segfault!" << endl;
    return;
  }

  // RooWorkspace is necessary for the model and parameters to be persistent
  // this is necessary because we created a bunch of objects that aren't persistent here
  fFitWorkspace = new RooWorkspace("fFitWorkspace", "Fit Workspace");

  // Energy shift
  double fDeltaE = 0.00;

  string tritDir = "/Users/wisecg/dev/spec-fit/data/";
  TFile *tritFile = new TFile(Form("%s/TritSpec.root", tritDir.c_str()));
  TH1D *tritSpec = dynamic_cast<TH1D*>(tritFile->Get("tritHist"));

  RooDataHist tritRooHist("trit", "Tritium Histogram", *fEnergy, Import(*tritSpec));
  // Because Steve's histogram sucks
  // The range of the histogram is maxed out at 50 keV, so need to reset range after loading histogram
  fEnergy->setRange(fFitMin, fFitMax);
  RooHistPdf tritPdf("tritPdf", "TritiumPdf", *fEnergy, tritRooHist, 2);

  TFile *axionFile = new TFile(Form("%s/axionHistos.root", tritDir.c_str()));
  TH1D *axionSpec = dynamic_cast<TH1D*>(axionFile->Get("hConv"));
  RooDataHist axionRooHist("axion", "Axion Histogram", *fEnergy, Import(*axionSpec));
  fEnergy->setRange(fFitMin, fFitMax);
  RooHistPdf axionPdf("axionPdf", "AxionPdf", *fEnergy, axionRooHist, 2);

  // Change this for polynomial background or linear
  // RooRealVar polySlope("polySlope", "Background Slope", 0.00002, -0.2, 0.2);
  // RooArgList polyList(polySlope);
  // RooPolynomial BkgPoly("Background", "Linear Background function", *fEnergy, polyList);
  RooPolynomial BkgPoly("Background", "Linear Background function", *fEnergy, RooArgList() );

  RooRealVar Ge68L_mean("Ge68L_mean", "Ge68L_mean", 1.2 + fDeltaE);
  RooRealVar Ge68L_sigma("Ge68L_sigma", "Ge68L_sigma", GetSigma(1.2 + fDeltaE));
  RooGaussian Ge68L_gauss("Ge68L_gauss", "Ge68L Gaussian", *fEnergy, Ge68L_mean, Ge68L_sigma);

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

  RooRealVar Ge68_mean("Ge68_mean", "Ge68_mean", 10.37 + fDeltaE);
  RooRealVar Ge68_sigma("Ge68_sigma", "Ge68_sigma", GetSigma(10.37 + fDeltaE));
  RooGaussian Ge68_gauss("Ge68_gauss", "Ge68 Gaussian", *fEnergy, Ge68_mean, Ge68_sigma);

  RooRealVar Pb210_mean("Pb210_mean", "Pb210_mean", 46.54 + fDeltaE);
  RooRealVar Pb210_sigma("Pb210_sigma", "Pb210_sigma", GetSigma(46.54 + fDeltaE));
  RooGaussian Pb210_gauss("Pb210_gauss", "Pb210 Gaussian", *fEnergy, Pb210_mean, Pb210_sigma);

  // Normalization parameters
  // Make names pretty for plots
  RooRealVar num_trit("Tritium", "Tritium", 6700.0, 0.0, 100000.);
  RooRealVar num_axion("Axion", "Axion", 0.0, 0.0, 10000.);
  RooRealVar num_bkg("Bkg", "Background", 50.0, 0.0, 100000.);
  RooRealVar num_Ge68L("Ge68L", "Ge68L", 5.0, 0.0, 50000.);
  RooRealVar num_Mn54("Mn54", "Mn54", 5.0, 0.0, 50000.);
  RooRealVar num_Fe55("Fe55", "Fe55", 5.0, 0.0, 50000.);
  RooRealVar num_Co57("Co57", "Co57", 5.0, 0.0, 50000.);
  RooRealVar num_Zn65("Zn65", "Zn65", 5.0, 0.0, 50000.);
  RooRealVar num_Ge68("Ge68", "Ge68", 180.0, 0.0, 50000.);
  RooRealVar num_Pb210("Pb210", "Pb210", 5.0, 0.0, 50000.);

  // Extended PDF model -- use this to create an extended model
  RooExtendPdf tritPdfe("tritPdfe", "Extended trit", tritPdf, num_trit);
  RooExtendPdf axionPdfe("axionPdfe", "Extended axion", axionPdf, num_axion);
  RooExtendPdf BkgPolye("BkgPolye", "Extended BkgPoly", BkgPoly, num_bkg);
  RooExtendPdf Ge68L_gausse("Ge68L_gausse", "Extended Ge68L_gauss", Ge68L_gauss, num_Ge68L);
  RooExtendPdf Mn54_gausse("Mn54_gausse", "Extended Mn54_gauss", Mn54_gauss, num_Mn54);
  RooExtendPdf Fe55_gausse("Fe55_gausse", "Extended Fe55_gauss", Fe55_gauss, num_Fe55);
  RooExtendPdf Co57_gausse("Co57_gausse", "Extended Co57_gauss", Co57_gauss, num_Co57);
  RooExtendPdf Zn65_gausse("Zn65_gausse", "Extended Zn65_gauss", Zn65_gauss, num_Zn65);
  RooExtendPdf Ge68_gausse("Ge68_gausse", "Extended Ge68_gauss", Ge68_gauss, num_Ge68);
  RooExtendPdf Pb210_gausse("Pb210_gausse", "Extended Pb210_gauss", Pb210_gauss, num_Pb210);

  // BDM PDF
  RooRealVar alpha("alpha", "Scale Correction", -0.0014, -0.01, -0.000001);
  RooRealVar E0("E0", "Offset Correction", -0.256, -0.4, -0.00001);
  RooRealVar envalRoo("envalRoo", "Real Energy", enVal); // Evaluated energy -- enVal changes
  RooArgList mAList(alpha, envalRoo, E0);
  RooFormulaVar peak_mean("mA", "envalRoo + alpha*(envalRoo - 95.0) + E0", mAList);
  RooRealVar peak_sigma("peak_sigma", "peak_sigma", GetSigma(enVal + fDeltaE));
  RooGaussian peakGaus("peak_gaus", "gaussian for DM signal", *fEnergy, peak_mean, peak_sigma);
  RooRealVar peakYieldInit("peakYieldInit", "yield signal peak", 0.1, 0.0, 10000.);

  RooArgList shapes(tritPdfe, axionPdfe, BkgPolye, Ge68L_gausse, Mn54_gausse, Fe55_gausse, Zn65_gausse, Ge68_gausse, Pb210_gausse);
  // RooArgList shapes(tritPdfe, BkgPolye, Mn54_gausse, Fe55_gausse, Zn65_gausse, Ge68_gausse);
  RooAddPdf model("model", "total pdf", shapes);

  // Add model to workspace -- also adds all of the normalization constants
  // If this step isn't done, a lot of the later functions won`'t work!
  fFitWorkspace->import(RooArgSet(model));
  fModelPDF = fFitWorkspace->pdf("model");
}

void GPXFitter::DoFit(string Minimizer)
{
  // Create NLL (This is not a profile! When you draw it to one axis, it's just a projection!)
  fNLL = fModelPDF->createNLL(*fRealData, Extended(), NumCPU(4));

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

void GPXFitter::DrawBasicShit(double binSize, bool drawResid, bool drawMatrix)
{
  TCanvas *cSpec = new TCanvas("cSpec", "cSpec", 1100, 800);
  RooPlot* frameFit = fEnergy->frame(Range(fFitMin, fFitMax), Bins((fFitMax - fFitMin)*1.0/binSize + 0.5));
  fRealData->plotOn(frameFit);

  // returns final value
  double axionScale = dynamic_cast<RooRealVar*>(fFitResult->floatParsFinal().find("Axion") )->getValV();

  RooAbsPdf* ax = fFitWorkspace->pdf("axionPdf");

  fModelPDF->plotOn(frameFit, Components("axionPdfe"), LineColor(kBlue), LineStyle(kDashed));

  fModelPDF->plotOn(frameFit, LineColor(kRed));
  frameFit->SetTitle("");

  // Get parameter values from first fit... these methods suck
  string tritDir = "/Users/wisecg/dev/spec-fit/data/";
  TFile *tritFile = new TFile(Form("%s/TritSpec.root", tritDir.c_str()));
  TH1D *tritSpec = dynamic_cast<TH1D*>(tritFile->Get("tritHist"));

  double tritVal = dynamic_cast<RooRealVar*>(fFitResult->floatParsFinal().find("Tritium") )->getValV();
  double tritErr = dynamic_cast<RooRealVar*>(fFitResult->floatParsFinal().find("Tritium") )->getError();
  double tritValCorr = tritVal/(tritSpec->Integral(tritSpec->FindBin(fFitMin), tritSpec->FindBin(fFitMax)));
  double tritErrCorr = tritErr/(tritSpec->Integral(tritSpec->FindBin(fFitMin), tritSpec->FindBin(fFitMax)));
  double geVal = dynamic_cast<RooRealVar*>(fFitResult->floatParsFinal().find("Ge68"))->getValV();
  double geErr = dynamic_cast<RooRealVar*>(fFitResult->floatParsFinal().find("Ge68"))->getError();
  double axVal = dynamic_cast<RooRealVar*>(fFitResult->floatParsFinal().find("Axion"))->getValV();
  double axErr = dynamic_cast<RooRealVar*>(fFitResult->floatParsFinal().find("Axion"))->getError();

  // Add chi-square to the plot - it's fairly meaningless as it's an unbinned fit but people will want it
  fChiSquare = frameFit->chiSquare(9);
  TPaveText *leg = new TPaveText(0.50, 0.75, 0.88, .88, "NDC");
  leg->SetTextFont(133);
  leg->SetFillColor(0);
  leg->SetBorderSize(1);
  leg->SetTextSize(22);
  leg->AddText(Form("#chi^{2}/NDF = %.3f" ,fChiSquare ) );
  leg->AddText(Form("Tritium (Uncorrected): %.3f #pm %.3f", tritVal, tritErr));
  leg->AddText(Form("Tritium (Corrected): %.3f #pm %.3f", tritValCorr, tritErrCorr));
  leg->AddText(Form("Ge68: %.3f #pm %.3f", geVal, geErr));
  leg->AddText(Form("Axions: %.3f #pm %.3f", axVal, axErr));
  frameFit->addObject(leg);
  frameFit->Draw();
  cSpec->SaveAs(Form("./plots/%s_Spec.pdf", fSavePrefix.c_str()) );
  fFitWorkspace->import(*frameFit);

  if(drawMatrix)
  {
    TCanvas *cMatrix = new TCanvas("cMatrix", "cMatrix", 1100, 800);
    TH2D *fCorrMatrix = dynamic_cast<TH2D*>(fFitResult->correlationHist("Correlation Matrix"));
    fCorrMatrix->Draw("colz");
    cMatrix->SaveAs(Form("./plots/%s_CorrMatrix.pdf", fSavePrefix.c_str()) );
    fFitWorkspace->import(*fCorrMatrix);
  }

  if(drawResid)
  {
    TCanvas *cResidual = new TCanvas("cResidual", "cResidual", 1100, 800);
    RooHist *hresid = frameFit->pullHist();
    RooPlot *frameResid = fEnergy->frame(Title("Normalized Fit Residuals"));
    frameResid->addPlotable(hresid, "P");
    frameResid->GetYaxis()->SetTitle("Normalized Residuals (#sigma)");
    frameResid->Draw();
    cResidual->SaveAs(Form("./plots/%s_Residual.pdf", fSavePrefix.c_str()) );
  }
}

void GPXFitter::DrawContour(string argN1, string argN2)
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
  cContour->SaveAs(Form("./plots/%s_Contour_%svs%s.pdf", fSavePrefix.c_str(), argN2.c_str(), argN1.c_str()));
fFitWorkspace->import(*frameContour);
}


// Use after constructing the model and minimization!
// This is a simple MC study only generating parameter information and pulls, the toy MC data can be saved
// I forget a lot of the crap I did here, a lot of it was for boring tests Jason and Reyco wanted
void GPXFitter::GenerateMCStudy(vector<string> argS, int nMC)
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
        cMCStudy->SaveAs(Form("./plots/%s_%s_MCStudy.pdf", fSavePrefix.c_str(), argN.c_str()) );
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
    // cM->SaveAs("./plots/MeanTest.pdf");

    // TCanvas *cLo = new TCanvas("cLo", "cLo", 800, 600);
    // hErrLo->Draw();
    // l1.DrawLine(parErrLo0, 0, parErrLo0, hErrLo->GetBinContent(hErrLo->GetMaximumBin()));
    // cLo->SaveAs("./plots/MeanTest_Low.pdf");

    // TCanvas *cHi = new TCanvas("cHi", "cHi", 800, 600);
    // hErrHi->Draw();
    // l1.DrawLine(parErrHi0, 0, parErrHi0, hErrHi->GetBinContent(hErrHi->GetMaximumBin()));
    // cHi->SaveAs("./plots/MeanTest_High.pdf");

    // double NLLmean = dynamic_cast<RooRealVar*>(fFitResult->floatParsFinal().find("NLL"))->getValV();
    // double NLLmean = dynamic_cast<RooRealVar*>(fFitResult->floatParsFinal().find("NLL"))->getError();

 /*
    // Extract nLL variables from MCStudy -- maybe useful later
    RooDataSet MCFitData = fMCStudy->fitParDataSet();
    // RooDataSet *NLL = static_cast<RooDataSet*>(MCFitData.reduce(RooArgSet()));
    MCFitData.Print("v");
    const RooArgSet* row = MCFitData.get();
    row->Print("v");
    shared_ptr<TCanvas> cMCNLL( make_shared<TCanvas>("cMCNLL", "cMCNLL", 1100, 800) );
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
void GPXFitter::GenerateToyMC(string fileName)
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

// Gets resolution, function and parameters from BDM PRL paper
// https://arxiv.org/abs/1612.00886
// In the future it should just be a convolution with all the other PDFs?
double GPXFitter::GetSigma(double energy)
{
    double p0, p1, p2;
    if(fDS==0) {p0 = 0.147; double p1=0.0173; double p2=0.0003;}
    else if(fDS==1) {p0 = 0.136; double p1=0.0174; double p2=0.00028;}
    else if(fDS==3) {p0 = 0.162; double p1=0.0172; double p2=0.000297;}
    else if(fDS==4) {p0 = 0.218; double p1=0.015; double p2=0.00035;}
    else if(fDS==5) {p0 = 0.2121; double p1=0.01838; double p2=0.00031137;}
    else if(fDS==6) {p0 = 0.2121; double p1=0.01838; double p2=0.00031137;}
    double sig = sqrt(p0*p0 + p1*p1*energy + p2*p2*energy*energy );

	return sig;
}

// Assumes standard skim format -- converts stuff from vector<double> to scalar
// Also assumes trapENFCal is the energy parameter of choice
void GPXFitter::LoadChainData(TChain *skimTree, string theCut)
{
  // First get TEntryList with TCut
  skimTree->Draw(">> elist", Form("%s", theCut.c_str()), "entrylist goff");
  TEntryList *elist = dynamic_cast<TEntryList*>(gDirectory->Get("elist"));
  skimTree->SetEntryList(&*elist); // This works
  cout << Form("Using cut: %s", theCut.c_str()) << endl;
  cout << Form("Found %lld entries passing cuts", elist->GetN()) << endl;

  // I found it easier to work like this rather than with a TTreeReader...
  vector<double> *ftrapENFCal = nullptr;
  vector<int> *fchannel = nullptr;
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
    if(i%5000==0) cout << "Processing event: " << i << endl;

    for(int j = 0; j < ftrapENFCal->size(); j++)
    {
      trapENFCal = ftrapENFCal->at(j);
      channel = fchannel->at(j);
      dummyTree->Fill();
    }
  }
  cout << "Dummy Tree filled entries: " << dummyTree->GetEntries() << endl;

  // Can and perhaps should split the data up by channel in a more complicated fit
  fEnergy = new RooRealVar("trapENFCal", "trapENFCal", fFitMin, fFitMax, "keV");
  // fEnergy = new RooRealVar("trapENFCal", "trapENFCal", 0, 250, "keV");
  fRealData = new RooDataSet("data", "data", dummyTree, RooArgSet(*fEnergy));
}

// Implemented now in RooStats rather than RooFit
// Calculates profile likelihood and spits out limits
map<string, vector<double>> GPXFitter::ProfileNLL(vector<string> argS, double CL)
{
  map<string, vector<double>> LimitMap;
  for(auto &argN : argS)
  {
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
    cNLL->SaveAs(Form("./plots/%s_%sNLL.pdf", fSavePrefix.c_str(), argN.c_str()) );
    vector<double> Limits = {lowerLimit, upperLimit};
    LimitMap[argN.c_str()] = Limits;
  }

return LimitMap;
}

void GPXFitter::SaveShit(string outfileName)
{
  TFile *fOut = new TFile( Form("./plots/%s_%s", fSavePrefix.c_str(), outfileName.c_str()), "RECREATE" );
  fOut->cd();
  fFitWorkspace->Write();
  fOut->Close();
}

void GPXFitter::SetFitRange(double fitMin, double fitMax)
{
  fFitMin = fitMin;
  fFitMax = fitMax;
}