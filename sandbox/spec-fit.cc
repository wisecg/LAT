// "General"-purpose MJD histogram fitter.
// Uses TMinuit.
// C. Wiseman USC/Majorana 2/24/17

#include <iostream>
#include <fstream>
#include <cmath>
#include "TFile.h"
#include "TROOT.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TRandom3.h"
#include "TMinuit.h"
#include "TMath.h"
#include "TGraph.h"
#include "TMatrixT.h"
#include "TLine.h"

using namespace std;

void generateSpec(double bins, double xlow, double xhi);
void fitSpec(double bins, double xlow, double xhi, double fitLo);
void getDetectorParams(double E, double &deltaE, double &sig);

char theCut[1000];
char standardCut[1000] = "trapENFCal < 50 && trapENFCal > 0.8 && channel!=596 && channel!=676 && channel!=644 && channel!=612 && gain==0 && trapETailMin < 0 && mH==1 && isGood && !isLNFill && !wfDCBits && !muVeto && !(time_s > 2192e3 && time_s < 2195e3) && !(time_s > 7370e3 && time_s <7371e3) && !(time_s > 7840e3 && time_s < 7860e3) && !(time_s > 8384e3 && time_s < 8387e3) && !(time_s > 8984e3 && time_s < 8985e3) && !(time_s > 9002e3 && time_s < 9005e3) && run != 13075 && run != 13093 && run != 13116 && run!=9648 && run!=10663 && run!=10745 && run!=11175 && run!=12445 && run!=12723 && run!=12735 && run!=12745 && run!=12746 && run!=12765 && run!=12766 && run!=12767 && run!=13004";
char ds1_toeCut[1000] = "(((kvorrT/trapENFCal) >= 1.2 && (kvorrT/trapENFCal) <= 2.1) || (channel==580 && ((kvorrT/trapENFCal) >= 0.9 && (kvorrT/trapENFCal) <= 2.1)) || (channel==664 && ((kvorrT/trapENFCal) >= 1.1 && (kvorrT/trapENFCal) <= 2.1)))";

int main()
{
  // gROOT->ProcessLine(".x ~/env/MJDClintPlotStyle.C");

  // Everything has to be binned the same way for the fit to work
  double binsPerkeV = 5;
  double xlow = 1, xhi = 50;
  double bins = (xhi - xlow)*binsPerkeV;
  // double fitLimit = 2;

  generateSpec(bins,xlow,xhi);
  // fitSpec(bins, xlow, xhi, fitLimit);
}

// Adjust gamma peaks in our model, using
// DS-0 low energy paper results.
void getDetectorParams(double E, double &deltaE, double &sig)
{
  // Energy scale correction, eq. 1
  // double alphaE = -0.0014, E_0 = -0.256;
  // deltaE = alphaE * (E - 95.) + E_0;
  // deltaE = -0.03; // brandon says fixed value for ds-1
  deltaE = -0.05; // brandon says fixed value for ds-1

  // Resolution widths, eq. 2
  double sig0 = 0.16, F=0.11, eps=0.00296;
  sig = sqrt(pow(sig0,2) + eps * F * (E + deltaE));
}

// Make a ROOT file with the spectrum model
void generateSpec(double bins, double xlow, double xhi)
{
  string outputFile = "./axionHistos.root";

  // ================= MJD DS-1 w/ T/E =================
  // TChain *cutSkim = new TChain("skimTree");
  // cutSkim->Add("./data/standardCut_DS1.root");
  // printf("Found %lli entries.\n",cutSkim->GetEntries());
  //
  // TH1D *hEnr = new TH1D("hEnr","hEnr",bins,xlow,xhi);
  // sprintf(theCut,"%s && %s && isEnr",standardCut,ds1_toeCut);
  // cutSkim->Project("hEnr","trapENFCal",theCut);
  //
  // TH1D *hNat = new TH1D("hNat","hNat",bins,xlow,xhi);
  // sprintf(theCut,"%s && !isEnr",standardCut);
  // cutSkim->Project("hNat","trapENFCal",theCut);

  // ================= Tritium =================
  vector<double> energies, tritVals;
  ifstream tritFile("./data/TritiumSpectrum.txt");
  double ener, trit;
  while(tritFile >> ener >> trit) {
    energies.push_back(ener);
    if (trit > 0) tritVals.push_back(trit);
    else tritVals.push_back(0.);
  }
  tritFile.close();
  TH1D *hTrit = new TH1D("hTrit","hTrit",bins,xlow,xhi);
  TAxis *xaxis = hTrit->GetXaxis();
  for (size_t i = 0; i < energies.size(); i++) {
    int binx = xaxis->FindBin(energies[i]);
    hTrit->SetBinContent(binx,tritVals[i]);
    // cout << binx << " " << tritVals[i] << endl;
  }
  hTrit->Scale(1/hTrit->Integral());

  // ================= Convolved Axion Spectrum =================
  vector<double> axEnergies, axVals;
  ifstream axFile("./data/RedondoSpectrum.txt");
  double axionFlux;
  while(axFile >> ener >> axionFlux) {
    axEnergies.push_back(ener);
    axVals.push_back(axionFlux);
  }
  axFile.close();
  vector<double> sigEne, ge76XS;
  ifstream xsFile("./data/xs.txt");
  double se, xs;
  while(xsFile >> se >> xs){
    sigEne.push_back(se);
    ge76XS.push_back(xs);
  }
  xsFile.close();
  TH1D *hAxion = new TH1D("hAxion","hAxion",bins,xlow,xhi);
  TH1D *hSig = new TH1D("hSig","hSig",bins,xlow,xhi);
  TH1D *hConv = new TH1D("hConv","hConv",bins,xlow,xhi);
  TH1D *hBrem = new TH1D("hBrem","hBrem",bins,xlow,xhi);
  TH1D *hBConv = new TH1D("hBConv","hBConv",bins,xlow,xhi);
  xaxis = hAxion->GetXaxis();
  for (size_t i = 0; i < energies.size(); i++)
  {
    int binx = xaxis->FindBin(energies[i]);
    hAxion->SetBinContent(binx,axVals[i]);

    // interpolate the axion cross section
    size_t idx = 0;
    for (size_t j = 0; j < sigEne.size(); j++){
      idx = j;
      if (sigEne[j] > energies[i]) break;
    }
    double sig_hi=ge76XS[idx], EH=sigEne[idx], E=energies[i], sig_lo=0, EL=0, sig=0;
    if (idx!=0) {
      sig_lo = ge76XS[idx-1];
      EL = sigEne[idx-1];
      sig = sig_lo + (sig_hi - sig_lo)*((E - EL)/(EH - EL));
    }
    else
      sig = sig_hi;
    // printf("%lu  %.1f  idx %lu  sigEne %.1f  XS %.1f  sig %.1f  E - EL %.2f  EH %.2f  EL %.2f\n" ,i,energies[i],idx,sigEne[idx],ge76XS[idx],sig,E-EL,EH,EL);

    // make brems curve: x^0.89 Exp[-0.7 x - 1.26 Sqrt[x]]
    double ene = energies[i];
    double brem = pow(ene,0.89) * exp(-0.7*ene - 1.26 * sqrt(ene));

    hSig->SetBinContent(binx,sig);
    hConv->SetBinContent(binx,axVals[i]*sig);
    hBrem->SetBinContent(binx,brem);
    hBConv->SetBinContent(binx,brem*sig);
  }
  // hAxion->Scale(1/hAxion->Integral("w"));
  // hConv->Scale(1/hConv->Integral("w"));
  // hSig->Scale(1/hSig->Integral("w"));
  // hBrem->Scale(1/hBrem->Integral("w"));
  // hBConv->Scale(1/hBConv->Integral("w"));

  // ================= noise curve =================
  TH1D *hNoise = new TH1D("hNoise","hNoise",bins,xlow,xhi);
  TAxis *ax2 = hNoise->GetXaxis();
  for (size_t i = 0; i < sigEne.size(); i++)
  {
    int binx = ax2->FindBin(sigEne[i]);

    double noise = 0;
    // if (sigEne[i] > 1. && sigEne[i] < 3.)
      // noise = 1/sigEne[i];
      noise = TMath::Gaus(sigEne[i], 2., 1.2); // x, mean, sigma

    hNoise->SetBinContent(binx,noise);
  }
  hNoise->Scale(1/hNoise->Integral("w"));

  // ================= Gamma peaks =================
  // cout << "Generating adjusted gamma peaks ...\n";
  vector<double> deltaEs(5);
  vector<double> widths(5);
  getDetectorParams(5.99,deltaEs[0],widths[0]);
  getDetectorParams(6.54,deltaEs[1],widths[1]);
  getDetectorParams(7.11,deltaEs[2],widths[2]);
  getDetectorParams(8.98,deltaEs[3],widths[3]);
  getDetectorParams(10.37,deltaEs[4],widths[4]);
  // cout << Form("For 5.99 keV, DeltaE %.2f  E+DeltaE %.2f  Width %.3f\n",deltaEs[0],5.99+deltaEs[0],widths[0])
  //      << Form("For 6.54 keV, DeltaE %.2f  E+DeltaE %.2f  Width %.3f\n",deltaEs[1],6.54+deltaEs[1],widths[1])
  //      << Form("For 7.11 keV, DeltaE %.2f  E+DeltaE %.2f  Width %.3f\n",deltaEs[2],7.11+deltaEs[2],widths[2])
  //      << Form("For 8.98 keV, DeltaE %.2f  E+DeltaE %.2f  Width %.3f\n",deltaEs[3],8.98+deltaEs[3],widths[3])
  //      << Form("For 10.37 keV, DeltaE %.2f  E+DeltaE %.2f  Width %.3f\n",deltaEs[4],10.37+deltaEs[4],widths[4]);
  TH1D *hFlat = new TH1D("hFlat","hFlat",bins,xlow,xhi);
  TH1D *pMn54 = new TH1D("pMn54","pMn54",bins,xlow,xhi);
  TH1D *pFe55 = new TH1D("pFe55","pFe55",bins,xlow,xhi);
  TH1D *pCo57 = new TH1D("pCo57","pCo57",bins,xlow,xhi);
  TH1D *pZn65 = new TH1D("pZn65","pZn65",bins,xlow,xhi);
  TH1D *pGe68 = new TH1D("pGe68","pGe68",bins,xlow,xhi);
  gRandom = new TRandom3();
  for (int i = 0; i < 2000; i++) {
    hFlat->SetBinContent(i,0.8);
    pMn54->Fill(gRandom->Gaus(5.99+deltaEs[0], widths[0]));
    pFe55->Fill(gRandom->Gaus(6.54+deltaEs[1], widths[1]));
    pCo57->Fill(gRandom->Gaus(7.11+deltaEs[2], widths[2]));
    pZn65->Fill(gRandom->Gaus(8.98+deltaEs[3], widths[3]));
    pGe68->Fill(gRandom->Gaus(10.37+deltaEs[4], widths[4]));
  }
  hFlat->Scale(1/hFlat->Integral("w"));
  pMn54->Scale(1/pMn54->Integral("w"));
  pFe55->Scale(1/pFe55->Integral("w"));
  pCo57->Scale(1/pCo57->Integral("w"));
  pZn65->Scale(1/pZn65->Integral("w"));
  pGe68->Scale(1/pGe68->Integral("w"));

  // ================= Write all output =================
  TFile *f = new TFile(outputFile.c_str(),"RECREATE");
  // hEnr->Write();
  // hNat->Write();
  hTrit->Write();
  hAxion->Write();
  hSig->Write();
  hConv->Write();
  hBrem->Write();
  hBConv->Write();
  hFlat->Write();
  pMn54->Write();
  pFe55->Write();
  pCo57->Write();
  pZn65->Write();
  pGe68->Write();
  hNoise->Write();
  f->Close();
}

// To use this fitter in any c++ program, copy this class
// and the function "globFCN" under it.  Make sure to
// add -lMinuit to the ROOT libraries in the makefile.
class histFitter : public TObject
{
  public:
  int fPars;
  double fChiSquare, fBins, fLo, fHi, fFitLo;
  vector<double> fParameters;
  vector<double> fParErrors;
  TH1D *hData;
  TH1D *hModelTot;
  vector<TH1D*> hModel;

  histFitter(int pars, double bins, double xl, double xh, double fitLo=0) {
    fPars = pars;
    fBins = bins;
    fLo = xl;
    fHi = xh;
    fFitLo = fitLo;
    fChiSquare = 0.;
    fParameters.resize(fPars);
    fParErrors.resize(fPars);
  };
  virtual ~histFitter() {
    delete hModelTot;
    delete hData;
    for (auto h : hModel) delete h;
  };
  void SetData(TH1D* hdat) {
    hModelTot = new TH1D("hModelTot", "sum of model components",fBins,fLo,fHi);
    hData = hdat;
  };
  void AddModelHist(TH1D *hComponent) {
    hModel.push_back(hComponent);
  };
  void SetParValue(int fIndex, double fValue) {
    fParameters[fIndex] = fValue;
  };
  void UpdateModel() {
    hModelTot->Reset();
    for (int im = 0; im<fPars; im++)
      hModelTot->Add(hModel[im],fParameters[im]);
  }
  double GetChiSquare()
  {
    double chiSquare = 0., datam1_i, modelm1_i;
    int fitFloor = (int)(fFitLo*(fBins/(fHi-fLo)));
    // for(int i = fitFloor; i < hData->GetNbinsX(); i++) {
    //   datam1_i = hData->GetBinContent(i)*hData->GetBinWidth(i);
    //   modelm1_i = hModelTot->GetBinContent(i)*hModelTot->GetBinWidth(i);
    //   if(modelm1_i != 0 && datam1_i != 0)
    //     chiSquare += 2 * (modelm1_i - datam1_i + datam1_i * TMath::Log(datam1_i/modelm1_i));
    //   if (datam1_i == 0) chiSquare += 2*modelm1_i;
    // }
    for(int i = fitFloor; i < 25; i++) {
      datam1_i = hData->GetBinContent(i)*hData->GetBinWidth(i);
      modelm1_i = hModelTot->GetBinContent(i)*hModelTot->GetBinWidth(i);
      if(modelm1_i != 0 && datam1_i != 0)
        chiSquare += 2 * (modelm1_i - datam1_i + datam1_i * TMath::Log(datam1_i/modelm1_i));
      if (datam1_i == 0) chiSquare += 2*modelm1_i;
    }
    for(int i = 30; i < 75; i++) {
      datam1_i = hData->GetBinContent(i)*hData->GetBinWidth(i);
      modelm1_i = hModelTot->GetBinContent(i)*hModelTot->GetBinWidth(i);
      if(modelm1_i != 0 && datam1_i != 0)
        chiSquare += 2 * (modelm1_i - datam1_i + datam1_i * TMath::Log(datam1_i/modelm1_i));
      if (datam1_i == 0) chiSquare += 2*modelm1_i;
    }
    for(int i = 80; i < hData->GetNbinsX(); i++) {
      datam1_i = hData->GetBinContent(i)*hData->GetBinWidth(i);
      modelm1_i = hModelTot->GetBinContent(i)*hModelTot->GetBinWidth(i);
      if(modelm1_i != 0 && datam1_i != 0)
        chiSquare += 2 * (modelm1_i - datam1_i + datam1_i * TMath::Log(datam1_i/modelm1_i));
      if (datam1_i == 0) chiSquare += 2*modelm1_i;
    }
    return chiSquare;
  };
  int GetNParam() { return fPars; }
};

// Minuit requires this global function to be able to call the
// fitter's method that calculates the quantity to minimize
void globFCN(int &n, double *grad, double &fval, double x[], int code)
{
  (void)n; (void)grad; (void)code;  // suppress unused parameter warnings
  histFitter* Obj = (histFitter*)gMinuit->GetObjectFit();
  for(int i = 0; i < Obj->GetNParam(); i++) Obj->SetParValue(i, x[i]);
  Obj->UpdateModel();
  fval = Obj->GetChiSquare();
}

// Fit the spectrum to data
void fitSpec(double bins, double xlow, double xhi, double fitLo)
{
  TFile *f = new TFile("./data/lowBGHistos.root");

  int npars = 8;
  histFitter *ftr = new histFitter(npars,bins,xlow,xhi,fitLo); // (# models, bins, xlo, xhi, fitLo)
  ftr->SetData((TH1D*)f->Get("hEnr"));
  ftr->AddModelHist((TH1D*)f->Get("hTrit"));
  ftr->AddModelHist((TH1D*)f->Get("hFlat"));
  ftr->AddModelHist((TH1D*)f->Get("pGe68"));
  ftr->AddModelHist((TH1D*)f->Get("pFe55"));
  ftr->AddModelHist((TH1D*)f->Get("pMn54"));
  ftr->AddModelHist((TH1D*)f->Get("pCo57"));
  ftr->AddModelHist((TH1D*)f->Get("pZn65"));
  ftr->AddModelHist((TH1D*)f->Get("hNoise"));
  // ftr->AddModelHist((TH1D*)f->Get("hConv"));

  TMinuit minuit(ftr->fPars);
  minuit.SetPrintLevel(0);           // -1 (Quiet), 0 (Normal), 1 (Verbose)
  minuit.Command("set strategy 2");  // 0, 1, or 2; 2 is the slowest but best
  minuit.SetFCN(globFCN);
  minuit.SetObjectFit(ftr);
  // Parameters: (#, name, init val, init err, lower lim, upper lim)
  minuit.DefineParameter(0, "trit", 0.1, 1., 0., 10000.);
  minuit.DefineParameter(1, "flat", 0.1, 1., 0., 300);
  minuit.DefineParameter(2, "ge68", 0.1, 1., 0., 50);
  minuit.DefineParameter(3, "fe55", 0.1, 0.1, 0, 50);
  minuit.DefineParameter(4, "mn54", 0.1, 1., 0., 50);
  minuit.DefineParameter(5, "co57", 0.1, 0.1, 0., 50);
  minuit.DefineParameter(6, "zn65", 0.1, 1., 0., 50);
  minuit.DefineParameter(7, "noise", 20., 1., 0., 300.);
  // minuit.DefineParameter(8, "axions", 1., 0.1, 0., 300);

  // ============ Do the minimization ============
  int status = minuit.Migrad();
  cout << "Fit Status: " << status << endl;

  // Access fit params and update model one final time
  for(int i = 0; i < ftr->fPars; i++)
    minuit.GetParameter(i, ftr->fParameters[i], ftr->fParErrors[i]);
  ftr->UpdateModel();
  double finalChiSquare = ftr->GetChiSquare();
  double NDF = bins - (fitLo * (bins/(xhi-xlow))) - npars; // (# bins 2-30 kev) - (# fit pars)

  // TMath::Prob(cs,ndf) gives you the probability that the chi squared your fit gives
  // exceeds the chisquared value by chance.  It should be close to 0.5, since the fit
  // should give a greater cs half the time from statistical fluctuations.
  cout << Form("Final Vals:\n   chiSq %.2f  NDF %.0f  chiSq/NDF %.2f  TMath::Prob(chiSq,NDF) %.2f\n", finalChiSquare,NDF,finalChiSquare/NDF,TMath::Prob(finalChiSquare,NDF));

  // ========== For convenience, pull histos out of the fitter ==========
  int nPar = ftr->fPars;
  vector<TH1D*> hModel = ftr->hModel;
  TH1D* hModelTot = ftr->hModelTot;

  // ============ Compute integrals ============
  double sumCts=0;
  vector<double> totalCts;
  for (int i = 0; i < nPar; i++) {
    hModel[i]->Scale(ftr->fParameters[i]);
    totalCts.push_back(hModel[i]->Integral(hModel[i]->GetXaxis()->FindBin(fitLo),bins));
    cout << "model " << i << ": " << totalCts[i] << endl;
  }
  for (auto c : totalCts) sumCts += c;
  cout << Form("model total (%.1f - 30 keV): %.1f\n",fitLo,sumCts);
  cout << Form("data total (%.1f - 30 keV): %.1f\n",fitLo,ftr->hData->Integral(ftr->hData->GetXaxis()->FindBin(fitLo),bins));

  // ============ Compute residuals ============
  vector<double> ener, residual;
  for (int i=0; i < bins; i++)
  {
    // for bin i: (data_i – model_i)/error_i
    double ene = ftr->hData->GetBinCenter(i);
    double data = ftr->hData->GetBinContent(i);
    double model = hModelTot->GetBinContent(i);
    double err = sqrt(data);
    if (data==0) err = 1.2;  // integral of the poisson probability for 1 sigma centered at 0
    double res = (data - model)/err;
    if (ene > fitLo) {
     ener.push_back(ene);
     residual.push_back(res);
    }
    // cout << Form("%i  ene %.2f  data %.1f  err %.2f  model %.2f  res %.2f\n",i,ene,data,err,model,res);

    // Set error on hData
    ftr->hData->SetBinError(i,err);
    // cout << "err: " << 1./err << endl;
  }

  // ========== Draw fit with residuals ==========
  TCanvas *c1 = new TCanvas("c1","c1",800,700);
  TPad *pad1 = new TPad("pad1","pad1",0,0.3,1,1);
  TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.3);
  pad1->SetBottomMargin(0.00001);
  pad1->SetBorderMode(0);
  // pad1->SetLogy();
  pad2->SetTopMargin(0.00001);
  pad2->SetBottomMargin(0.35);
  pad2->SetBorderMode(0);
  pad1->Draw();
  pad2->Draw();
  pad1->cd();

  ftr->hData->SetMinimum(0.1);
  // ftr->hData->SetMaximum(70);
  // ftr->hData->DrawCopy("hist");    // draw just the histo
  // hData->SetFillStyle(3018);
  ftr->hData->SetLineColor(kBlue);
  ftr->hData->Draw("e same");     // draw just the error bars
  ftr->hData->GetXaxis()->SetRangeUser(xlow,xhi);

  ftr->hData->GetXaxis()->SetTitle("Energy (keV)");
  ftr->hData->GetYaxis()->SetTitle(Form("Counts/%.2fkeV",(xhi-xlow)/bins));
  c1->Update();

  hModelTot->SetLineColorAlpha(kRed,0.7);
  hModelTot->SetLineWidth(2);
  hModelTot->Draw("hist same");

  hModel[0]->SetLineColor(kGreen+2);  // tritium
  hModel[0]->Draw("hist same");

  TH1D *hTritFlat = new TH1D("hTritFlat","hTritFlat",bins,xlow,xhi);
  hTritFlat = (TH1D*)hModel[0];
  hTritFlat->Add(hModel[1]); // add flat
  hTritFlat->Draw("hist same");

  // hModel[4]->SetLineColor(kCyan+2); // axion
  // hModel[4]->Draw("hist same");

  hModel[7]->SetLineColor(kOrange+4); // noise
  hModel[7]->Draw("hist same");

  TLegend* leg1 = new TLegend(0.45,0.6,0.87,0.92);
  leg1->AddEntry(ftr->hData,"DS-1 (w/ T/E)","l");
  leg1->AddEntry(hModelTot,Form("Model (c/n:%.2f p:%.1f)",finalChiSquare/NDF,TMath::Prob(finalChiSquare,NDF)),"l");
  leg1->AddEntry(hTritFlat,"tritium+flat","l");
  // leg1->AddEntry(hModel[4],"axion*sigma_pe","l");
  leg1->AddEntry(hModel[7],"noise","l");
  leg1->Draw("SAME");
  c1->Update();

  // Draw residuals
  pad2->cd();
  TGraph *gResid = new TGraph(ener.size(),&(ener[0]),&(residual[0]));
  gResid->GetXaxis()->SetTitle("Energy (keV)");
  gResid->GetXaxis()->SetTitleOffset(3);
  gResid->GetXaxis()->SetLimits(xlow,xhi);
  gResid->GetYaxis()->SetTitle("res");

  gResid->SetMarkerStyle(kFullDotMedium);
  gResid->SetMarkerColor(kBlue);
  gResid->SetLineColor(kRed);
  gResid->Draw("ALP");
  gResid->GetYaxis()->SetNdivisions(505);

  TLine *line = new TLine(xlow,0,xhi,0);
  line->SetLineWidth(2);
  line->SetLineColorAlpha(kBlack,0.5);
  line->Draw("same");

  TLegend *leg2 = new TLegend(0.68,0.38,0.87,0.52);
  leg2->AddEntry(gResid,"(data - model)/sqrt(N)","lp");
  leg2->Draw("same");

  c1->cd();
  c1->Update();
  c1->Print("./output/fitResiduals.pdf");

  // ========== Draw correlation matrix ==========
  TMatrixT<double> mCorrMatrix;
  mCorrMatrix.ResizeTo(nPar, nPar);
  minuit.mnemat(mCorrMatrix.GetMatrixArray(), nPar);
  for(int i = mCorrMatrix.GetRowLwb(); i <= mCorrMatrix.GetRowUpb(); i++)
    for(int j = mCorrMatrix.GetColLwb(); j <= mCorrMatrix.GetColUpb(); j++)
      mCorrMatrix(i,j) = mCorrMatrix(i,j)/(ftr->fParErrors[i]*ftr->fParErrors[j]);
  gStyle->SetPalette(55);
  TCanvas *c0 = new TCanvas("c0", "c0", 800, 600);
  c0->SetTopMargin(0.05);
  c0->SetBottomMargin(0.05);
  c0->SetLeftMargin(0.05);
  c0->SetRightMargin(0.15);
  c0->SetFillColor(0);
  c0->SetBorderMode(0);
  c0->SetBorderSize(0);
  c0->Draw();
  c0->cd();
  c0->SetGrid();
  c0->SetTicks();
  c0->SetFillStyle(4000);
  mCorrMatrix.Draw("colz");
  c0->Print("./output/corrMatrix.pdf");
}

