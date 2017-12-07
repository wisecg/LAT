// Run the evolved version of Wenqin

#include "GPXFitter.hh"
#include "TStyle.h"
#include "TChain.h"
#include "RooFitResult.h"
#include "RooWorkspace.h"
#include "RooAbsReal.h"
#include "RooArgList.h"
#include "RooAbsArg.h"
#include "RooMsgService.h"
#include "TROOT.h"

using namespace std;
using namespace RooFit;

int main(int argc, char** argv)
{
  gROOT->ProcessLine("RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);");

  if(argc <= 4) {
    cout << "Usage: " << argv[0] << " [DS] [Fit Min] [Fit Max] [Nat/Enr]" << endl;
    return 0;
  }
  int fDS = atoi(argv[1]);
  float fitMin = atof(argv[2]);
  float fitMax = atof(argv[3]);
  string ftype = argv[4];

  //// For drawing pretty shit
  gStyle->SetOptStat(0);
  gStyle->SetPalette(kRainBow);

  //// Set cuts here
  string theCut = "";

  if(ftype == "Nat") theCut += "isNat"; // Set Enriched or Natural
  else if(ftype == "Enr") theCut += "isEnr";



  if(fDS == 0) theCut += "&&!(run==6811&&(channel==600||channel==696)) && channel!=656";
  else if(fDS == 3) theCut += "&&channel!=592 && channel!=692";
  else if(fDS == 4) theCut += "&&!(run==60001692 && (channel==1144))&&channel!=1332";
  else if(fDS == 5) theCut += "&&channel!=1124";
  else if(fDS == 6) {
    theCut += " && ((run >= 2580 && run <= 6963 && !(run==6811 && (channel==600||channel==696)) && channel!=656)";
    theCut += " || (run >= 9422 && run <= 14502)";
    theCut += " || (run >= 16797 && run <= 17980 && channel!=592 && channel!=692)";
    theCut += " || (run >= 60000802 && run <= 60001888 && !(run==60001692 && (channel==1144)) && channel!=1332))";
    // theCut += " || (run >= 18623 && run <= 25671 && channel!=1124 ))";
  }

    // theCut += "&&!(run==6811&&(channel==600||channel==696)) && !(run==60001692 && (channel==1144))";

  theCut += Form("&& trapENFCal>=%.2f && trapENFCal<=%.2f", fitMin, fitMax); // Energy cut for fit range

  GPXFitter *fitter = new GPXFitter(fDS, fitMin, fitMax);
  // This is just a string for output files
  fitter->SetSavePrefix(Form("LowE_DS%d_%s_%.1f_%.1f", fDS, ftype.c_str(), fitMin, fitMax));

  // Load data from TChain with a cut string
  TChain *skimTree = new TChain("skimTree");
  if (fDS!=6)
    skimTree->Add(Form("~/project/latskim/latSkimDS%d_*.root", fDS) );
  else if (fDS==6) {
    skimTree->Add("~/project/latskim/latSkimDS0*");
    skimTree->Add("~/project/latskim/latSkimDS1*");
    skimTree->Add("~/project/latskim/latSkimDS3*");
    skimTree->Add("~/project/latskim/latSkimDS4*");
  }
  fitter->LoadChainData(skimTree, theCut);

  // Construct PDF and do fit
  fitter->ConstructPDF();
  fitter->DoFit();

  // This draws the spectrum as well as the covariance matrix and residuals if you want
  fitter->DrawBasicShit(0.2, true, true);

  // List of parameters we want to do more studies on
  // vector<string> argList = {"Axion"};
  //// Profile likelihood calculation
  // map<string, vector<double>> LimitMap = fitter->ProfileNLL(argList);

  // Generate Toy MC study -- input a list of parameters to look at from toy MC
  // This shit takes a long time
  // fitter->GenerateMCStudy(argList, 1000);

  // Print limits
  // for(auto &kv : LimitMap) cout << kv.first << " limit: " << kv.second[0] << " " << kv.second[1] << endl;

  RooFitResult *fitResult = fitter->GetFitResult();
  fitResult->Print();

  return 0;
}