// shiftFit: unbinned likelihood fitter.
// A weak-sauce version of GPXFitter by Brian.
//
// C. Wiseman, USC
// 8/22/17

#include <iostream>
#include <stdio.h>
#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "RooDataSet.h"
#include "RooFitResult.h"
#include "RooMinimizer.h"
#include "RooRealVar.h"
#include "RooAbsPdf.h"
#include "RooWorkspace.h"
#include "RooGaussian.h"
#include "RooExtendPdf.h"
#include "RooPolynomial.h"
#include "RooAddPdf.h"
#include "RooAbsReal.h"
#include "RooPlot.h"
using namespace std;
using namespace RooFit;

class ubFit
{
  public:
    ubFit(double eMin, double eMax, string mode, bool split, vector<double> aPks):
      fEMin(eMin),fEMax(eMax),fChiSquare(0.),fMode(mode),fSplit(split),fAxPks(aPks)
      // fEnergy(nullptr),fRealData(nullptr),fMinimizer(nullptr),fFitResult(nullptr),
      // fFitWorkspace(nullptr),fNLL(nullptr),fModelPDF(nullptr),fShapes(nullptr)
    {
      if (mode=="malbek") {
        fInputFile = "malbek_data.root";
        fTreeName = "malbek_wrt";
        fDrawVar = "energy_keV";
        fGePks = {  // Graham's thesis, pg 78
          4.97,       // 49V
          6.54,       // 55Fe
          1.10, 8.98, // 65Zn
          9.66,       // 68Ga
          1.30, 10.37 // 68Ge & 71Ge
        };
      }
      else if (mode=="mjd") {
        fInputFile = "mjd.root";
        fTreeName = "mjd_tree";
        fDrawVar = "trapENFCal";
        fGePks = {10.3}; // placeholder
      }
      fSavePrefix = Form("./plots/%s-%.0f-%.0f",mode.c_str(),eMin,eMax);
    }
    virtual ~ubFit() {};
    void LoadData();
    double GetSigma(double energy);
    void DeclarePDF();
    void RunFitter(string minType="Minuit", int nCPU=2);
    void DrawResults(double binSize=0.1, string outDir=".");
    RooFitResult *GetFitResult() { return fFitResult; }

  private:
    bool fSplit;
    double fEMin, fEMax, fChiSquare;
    vector<double> fAxPks, fGePks;
    string fSavePrefix, fInputFile, fMode, fTreeName, fDrawVar;
    RooRealVar* fEnergy;
    RooDataSet* fRealData;
    RooMinimizer* fMinimizer;
    RooFitResult* fFitResult;
    RooWorkspace* fFitWorkspace;
    RooAbsReal* fNLL;
    RooAbsPdf* fModelPDF;
    RooArgList* fShapes;
};

void ubFit::LoadData()
{
  TFile f(fInputFile.c_str());

  // Simple method - no tree merging
  if (!fSplit) {
    TTree* t = (TTree*)f.Get(fTreeName.c_str());
    fEnergy = new RooRealVar(fDrawVar.c_str(), fDrawVar.c_str(), fEMin, fEMax, "keV");
    fRealData = new RooDataSet("data","data", t, RooArgSet(*fEnergy));
  }

  // Frank method: duplicate and shift the input trees to line up the peaks.
  // Shift everything to line up w/ the last peak in the vector.
  else {
    TTree* tree0 = (TTree*)f.Get(fTreeName.c_str());
    cout << "Found " << tree0->GetEntries() << " entries.  Merging ...\n";
    double ene0, wt0, shift;
    tree0->SetBranchAddress("energy_keV",&ene0);
    tree0->SetBranchAddress("weight",&wt0);

    TList* tList = new TList();
    vector<TTree*> tVec(fAxPks.size());
    TFile *tmp = new TFile(Form("%s_splitData.root",fMode.c_str()),"RECREATE");
    for (size_t i = 0; i < fAxPks.size(); i++)
    {
      double ene;
      double shift = fAxPks[fAxPks.size()-1]-fAxPks[i];
      tVec[i] = new TTree(Form("t%zu",i),Form("t%zu",i));
      tVec[i]->Branch("energy_keV",&ene,"energy_keV/D");
      for (Long64_t j = 0; j < tree0->GetEntries(); j++) {
        tree0->GetEntry(j);
        ene = ene0*wt0;
        ene += shift;
        tVec[i]->Fill();
      }
      cout << Form("fAxPks: %i  tree %i - %llu entries  shift %.3f kev\n",(int)fAxPks.size(),(int)tVec.size(),tVec[i]->GetEntries(),shift);
      tList->Add(tVec[i]);
      tVec[i]->Write();
    }
    // remove(Form("%s_splitData.root",fMode.c_str())); // keep the file for checks
    TTree* t = TTree::MergeTrees(tList);
    cout << Form("Trees merged, with %lld entries total.\n",t->GetEntries());
    fEnergy = new RooRealVar(fDrawVar.c_str(), fDrawVar.c_str(), fEMin, fEMax, "keV");
    fRealData = new RooDataSet("data","data", t, RooArgSet(*fEnergy));
  }
}

double ubFit::GetSigma(double E)
{
  double sig, sig_e, expval, F;

  // MALBEK - pg. 92 of Graham's thesis. (he actually ended up letting sigma float)
  sig_e = 0.0698;
  F = 0.21;
  expval = 0.00296;
  sig = sqrt(sig_e * sig_e + expval * F * E);
  // cout << "sig is: " << sig << endl;
  return sig;

  // MJD - nominally Kris's PRL: https://arxiv.org/pdf/1612.00886.pdf
  // GPXFitter had some ds-specific values in it, maybe pinghan has a different function ...
}

void ubFit::DeclarePDF()
{
  // Polynomial or linear background
  RooRealVar polySlope("polySlope", "Slope", 0.00002, -0.2, 0.2);
  RooArgList polyList(polySlope);
  RooPolynomial pbkg("Background", "Poly Background", *fEnergy, polyList);

  RooPolynomial fbkg("Background", "Linear Background", *fEnergy, RooArgList() );
  RooRealVar num_bkg("bkg", "Background", 50.0, 0.0, 100000.);
  RooExtendPdf bkge("bkge", "Extended bkg", fbkg, num_bkg);

  // Gaussian peak, centered around last entry.  TODO: Multiple peaks
  RooRealVar pk_mean("pk_mean", "pk_mean", fAxPks[fAxPks.size()-1]);
  RooRealVar pk_sig("pk_sig", "pk_sig", GetSigma(fAxPks[fAxPks.size()-1]));
  RooGaussian pk_gaus("pk_gaus", "Peak Gaussian", *fEnergy, pk_mean, pk_sig);
  RooRealVar num_pk("pk", "pk", 5.0, 0.0, 50000.);
  RooExtendPdf pk_gausse("pk_gause", "Extended pk_gauss", pk_gaus, num_pk);

  if (fMode=="malbek") fShapes = new RooArgList(bkge, pk_gausse);
  else if (fMode=="mjd") fShapes = new RooArgList(bkge);

  // Declare the full model PDF and add everything to the workspace
  RooAddPdf model("model", "total pdf", *fShapes);
  fFitWorkspace = new RooWorkspace("fFitWorkspace", "Fit Workspace");
  fFitWorkspace->import(RooArgSet(model));
  fModelPDF = fFitWorkspace->pdf("model");
}

void ubFit::RunFitter(string minType, int nCPU)
{
  // Create NLL (This is not a profile! When you draw it to one axis, it's just a projection!)
  fNLL = fModelPDF->createNLL(*fRealData, Extended(), NumCPU(nCPU));

  // Create minimizer and run initial pass
  fMinimizer = new RooMinimizer(*fNLL);
  fMinimizer->setMinimizerType(minType.c_str());
  fMinimizer->setPrintLevel(-1);
  fMinimizer->setStrategy(2);
  fMinimizer->migrad();

  // do some refinement steps (might not need 'em all)
  fMinimizer->hesse();
  fMinimizer->improve();
  fMinimizer->minos();

  // save results
  fFitResult = fMinimizer->save();
  fFitWorkspace->import(*fFitResult);
}

void ubFit::DrawResults(double binSize, string outDir)
{
  TCanvas *c = new TCanvas("c", "c", 1100, 800);
  RooPlot* frameFit = fEnergy->frame(Range(fEMin, fEMax), Bins((fEMax - fEMin)*1.0/binSize + 0.5));
  fRealData->plotOn(frameFit);
  frameFit->SetTitle("");

  fModelPDF->plotOn(frameFit, LineColor(kRed));

  // double axionScale = dynamic_cast<RooRealVar*>(fFitResult->floatParsFinal().find("Axion") )->getValV();
  // RooAbsPdf* ax = fFitWorkspace->pdf("axionPdf");
  // fModelPDF->plotOn(frameFit, Components("axionPdfe"), LineColor(kBlue), LineStyle(kDashed));

  // these access methods fuuucking suck
  // double pkVal = dynamic_cast<RooRealVar*>(fFitResult->floatParsFinal().find("pk_gaus"))->getValV();
  // double pkErr = dynamic_cast<RooRealVar*>(fFitResult->floatParsFinal().find("pk_gaus"))->getError();
  // double bkgVal = dynamic_cast<RooRealVar*>(fFitResult->floatParsFinal().find("bkg"))->getValV();
  // double bkgErr = dynamic_cast<RooRealVar*>(fFitResult->floatParsFinal().find("bkg"))->getError();

  // Add chi-square to the plot - it's fairly meaningless as it's an unbinned fit but people will want it
  // fChiSquare = frameFit->chiSquare(9);
  // TPaveText *leg = new TPaveText(0.50, 0.75, 0.88, .88, "NDC");
  // leg->SetTextFont(133);
  // leg->SetFillColor(0);
  // leg->SetBorderSize(1);
  // leg->SetTextSize(22);
  // leg->AddText(Form("#chi^{2}/NDF = %.3f" ,fChiSquare ) );
  // frameFit->addObject(leg);

  frameFit->Draw();
  c->Print(Form("%s-test.pdf",fSavePrefix.c_str()));
}

// ===============================================================================================
int main(int argc, char** argv)
{
  if (argc < 3) {
    cout << "Usage: ./control [minE] [maxE]\n";
    return 0;
  }
  double eMin = atof(argv[1]), eMax = atof(argv[2]);
  gROOT->ProcessLine("RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);");
  gStyle->SetOptStat(0);
  gStyle->SetPalette(kRainBow);

  vector<double> axPeaks = {1.739,1.836,2.307,2.464};

  ubFit fitter(eMin, eMax, "malbek", true, axPeaks);
  fitter.LoadData();
  fitter.DeclarePDF();
  fitter.RunFitter();
  fitter.DrawResults();
  RooFitResult *fitResult = fitter.GetFitResult();
  fitResult->Print();
}
// ===============================================================================================
