// auto-thresh.cc
// Requires access to built data,
// output of process_mjd_cal (pass 1 gat)
// and the APDB (MkCookie must be run recently)
// B. Zhu, C. Wiseman
// v1. 2017/2/28
// v2. 2017/6/01 - changed to run over data subsets (as defined in DataSetInfo.hh)

#include <iostream>
#include <vector>
#include <string>
#include <map>
#include "TTree.h"
#include "TFile.h"
#include "TH1D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TEntryList.h"
#include "GATDataSet.hh"
#include "MGTWaveform.hh"
#include "MJDBUtilities.hh"
#include "MJDatabase.hh"
#include "MJCouchAccess.hh"
#include "MJDocument.hh"
#include "MJProvenance.hh"
#include "MJAnalysisDoc.hh"
#include "MJAnalysisParameters.hh"
#include "MJTRun.hh"
#include "DataSetInfo.hh"

using namespace std;
using namespace MJDB;

void FindThresholds(int dsNum, int subNum);
vector<double> ThreshTrapezoidalFilter( const vector<double>& anInput, double RampTime, double FlatTime, double DecayTime );

int main(int argc, char** argv)
{
  if (argc < 3) {
    cout << "Usage: ./auto-thresh [dsNum] [subNum] \n";
    return 1;
  }
  int dsNum = stoi(argv[1]);
  int subNum = stoi(argv[2]);
  FindThresholds(dsNum, subNum);
}

void FindThresholds(int dsNum, int subNum)
{
  string outputFile = Form("./threshDS%d_%d.root", dsNum, subNum);
  GATDataSet ds;
  LoadDataSet(ds, dsNum, subNum);

  // Set up output file
  TFile *fOutput = new TFile(outputFile.c_str(), "RECREATE");
  TTree *fThreshTree = new TTree("threshTree", "Tree with vector of channels and thresholds");

  int runMin = 0, runMax = 0;
  vector<int> channelList(0);
  vector<double> threshADC(0);
  vector<double> threshADCErr(0);
  vector<double> sigmaADC(0);
  vector<double> sigmaADCErr(0);
  vector<double> threshCal(0);
  vector<double> sigmaCal(0);
  vector<int> threshFitStatus(0);
  vector<int> sigmaFitStatus(0);
  vector<double> CalOffset(0);
  vector<double> CalScale(0);
  vector<int> numTrigger(0);
  vector<int> numNoise(0);

  fThreshTree->Branch("runMin", &runMin, "runMin/I");
  fThreshTree->Branch("runMax", &runMax, "runMax/I");
  fThreshTree->Branch("channelList", &channelList );
  fThreshTree->Branch("threshADC", &threshADC );
  fThreshTree->Branch("threshADCErr", &threshADCErr );
  fThreshTree->Branch("sigmaADC", &sigmaADC );
  fThreshTree->Branch("sigmaADCErr", &sigmaADCErr );
  fThreshTree->Branch("threshCal", &threshCal );
  fThreshTree->Branch("sigmaCal", &sigmaCal );
  fThreshTree->Branch("threshFitStatus", &threshFitStatus);
  fThreshTree->Branch("sigmaFitStatus", &sigmaFitStatus);
  fThreshTree->Branch("CalOffset", &CalOffset);
  fThreshTree->Branch("CalScale", &CalScale);
  fThreshTree->Branch("numTrigger", &numTrigger);
  fThreshTree->Branch("numNoise", &numNoise);

  MJTChannelSettings *chanSet= ds.GetChannelSettings();
  vector<uint32_t> en = chanSet->GetEnabledIDList();

  // Check if multisampling is on and get the waveform from built data
  TChain *builtChain = ds.GetBuiltChain(false);
  MJTRun *mjRun = 0;
  builtChain->SetBranchAddress("run",&mjRun);
  builtChain->GetEntry(0);
  bool isMS = mjRun->GetUseMultisampling();
  string wfBranchName = "fWaveforms";
  if (isMS) {
    cout << "Multisampling is active.  Getting fAuxWaveform ...\n";
    wfBranchName = "fAuxWaveforms";
  }
  TTreeReader bReader(builtChain);
  TTreeReaderValue<TClonesArray> wfBranch(bReader,wfBranchName.c_str());
  // Get channel and energy estimator
  TChain *gatChain = ds.GetGatifiedChain(false);
  TTreeReader gReader(gatChain);
  TTreeReaderValue<double> runIn(gReader, "run");
  TTreeReaderArray<double> wfChan(gReader,"channel");
  TTreeReaderArray<double> wfENF(gReader,"trapENM");
  int nEntries = gatChain->GetEntries();
  cout << "Scanning DS " << dsNum << " subNum " << subNum << " , " << nEntries << " entries.\n";

  // Initialize channel map and histogram vectors
  map<int,int> channelMap;
  vector<TH1D*> hTrigger;
  vector<TH1D*> hNoise;
  int nChannel = en.size();
  int nTrigger[nChannel];
  int nNoise[nChannel];
  for(int i = 0; i < nChannel; i++) {
    channelMap.insert({en[i], i});
    hTrigger.push_back(new TH1D(Form("hTrigger-%d", en[i]), Form("hTrigger-%d", en[i]), 1000, -30, 30));
    hNoise.push_back(new TH1D(Form("hNoise-%d", en[i]), Form("hNoise-%d", en[i]), 500, -30, 30));
    nTrigger[i] = 0;
    nNoise[i] = 0;
  }

  // Loop over entries
  int channel = 0;
  double NoiseSample = 0, TriggerSample = 0, trapENF = 0;
  vector<double> TrapFilter;
  vector<double> Waveform;
  cout << "Sampling from waveforms ...\n";
  // Set minimum run as run from the first entry
  gReader.SetEntry(0);
  runMin = *runIn;
  for(int i = 0; i < nEntries; i++)
  {
    bReader.SetEntry(i);
    gReader.SetEntry(i);
    int nWF = (*wfBranch).GetEntriesFast();
    int nCh = wfChan.GetSize();
    if (nWF != nCh) {
      cout << "Skipped entry " << i << endl;
      continue;
    }
    // Loop over waveforms
    for (int iWF = 0; iWF < nWF; iWF++)
    {
      Waveform.clear();
      TrapFilter.clear();
      shared_ptr<MGTWaveform> clone(dynamic_cast<MGTWaveform*>((*wfBranch).At(iWF)->Clone()));
      Waveform = clone->GetVectorData();
      TrapFilter = ThreshTrapezoidalFilter(Waveform, 400, 180, 0);

      NoiseSample = TrapFilter[1];   // 1st sample for noise
      TriggerSample = TrapFilter[9]; // 9th sample is the crossing
      channel = wfChan[iWF];
      trapENF = wfENF[iWF];

      // Low energy = flat signal => rough representation of threshold
      // Increased to 10 for DS4, higher noise? Early on pulsers weren't on
      if(trapENF > 0 && trapENF < 10)  {
        hTrigger[ channelMap[channel] ]->Fill(TriggerSample);
        nTrigger[ channelMap[channel] ]++;
      }

      // High energy = sharp rise => 1st sample good representation of noise
      if(trapENF > 50) {
        hNoise[ channelMap[channel] ]->Fill(NoiseSample);
        nNoise[ channelMap[channel] ]++;
      }
    }
    // if (i % 10000 == 0) cout << i << " entries saved so far.\n";
  }

  // Set maximum run as run from the last entry
  gReader.SetEntry(nEntries-1);
  runMax = *runIn;
  // Query database to get calibration parameters for enabled channels.
  // (The DB is only accessed for the first channel, then saves the run info in a buffer.)
  // Have you run MkCookie?

  MJAnalysisDoc findResult;
  EnergyCalibration mycalibration;
  mycalibration.SetPSource(kpsTrapENF);

  // Evaluate thresholds for this run
  cout << "Evaluating thresholds for each channel ...\n";
  TF1 *gaus1 = new TF1("gaus1", "gaus(0)", 0, 30); // Threshold value limited by 30 right now
  TF1 *gaus2 = new TF1("gaus2", "gaus(0)", -30, 30);
  for(auto i : en)
  {
    gaus1->SetParameters(0,0,0);
    gaus2->SetParameters(0,0,0);

    // Only try to fit if we have entries to fit, duh ...
    if (nTrigger[channelMap[i]] > 0 && (nNoise[channelMap[i]] > 0))
    {
      threshFitStatus.push_back(hTrigger[ channelMap[i] ]->Fit("gaus1", "qNR", "", 0.1, 10.0) );
      sigmaFitStatus.push_back(hNoise[ channelMap[i] ]->Fit("gaus2", "qNR+") );
    }
    else {
      threshFitStatus.push_back(999999);
      sigmaFitStatus.push_back(999999);
    }
    numTrigger.push_back(nTrigger[ channelMap[i] ]);
    numNoise.push_back(nNoise[ channelMap[i] ]);

    channelList.push_back(i);
    threshADC.push_back(gaus1->GetParameter(1));
    sigmaADC.push_back(gaus2->GetParameter(2));
    threshADCErr.push_back(gaus1->GetParError(1));
    sigmaADCErr.push_back(gaus2->GetParError(1));

    size_t Length = findResult.GetAnalysisParameter(runMin, i, mycalibration.GetPSource(), mycalibration.GetPType());

    if(Length == 0) {
      CalScale.push_back(0);
      CalOffset.push_back(0);
      threshCal.push_back(0);
      sigmaCal.push_back(0);
    }
    if(Length > 0) {
      MJAnalysisDoc temp = findResult[Length-1];
      mycalibration.GetDBValue(temp);
      double dScale = mycalibration.Scale.Value();
      double dOffset = mycalibration.Offset.Value();

      CalScale.push_back( dScale );
      CalOffset.push_back( dOffset );

      threshCal.push_back( gaus1->GetParameter(1)*dScale + dOffset);
      sigmaCal.push_back( gaus2->GetParameter(2)*dScale + dOffset);
    }
  }

  // If either fit failed, put the threshold at 99999 keV
  for (size_t i = 0; i < threshADC.size(); i++) {
    if (threshFitStatus[i] != 0 || sigmaFitStatus[i] != 0) {
      threshCal[i] = 99999;
      sigmaCal[i] = 99999;
    }
    // cout << Form("Ch %i  keV %.3f +/- %-8.3f   ADC %.2e +/- %-8.2e   Noise %.2e +/- %-8.2e  nADC %-5i  nNoise %-5i  Fit %i %i\n",
      // channelList[i], threshCal[i], sigmaCal[i], threshADC[i], threshADCErr[i], sigmaADC[i], sigmaADCErr[i],
      // nTrigger[i], nNoise[i], threshFitStatus[i], sigmaFitStatus[i]); //, CalScale[i], CalOffset[i]);
  }
  fThreshTree->Fill();
  fThreshTree->Write("",TObject::kOverwrite);

  // TF1 *fEff1 = new TF1("fEff", "0.5*(1+TMath::Erf((x-[0])/(TMath::Sqrt(2)*[1]) ))", gaus1->GetParameter(1), gaus2->GetParameter(2));
  // TF1 *fEff1 = new TF1("fEff", "0.5*(1+TMath::Erf((x-[0])/(TMath::Sqrt(2)*[1]) ))",1.3436,0.244441);
  // cout << "fits: " << gaus1->GetParameter(1) << "  " << gaus2->GetParameter(2) << endl;
  // fEff1->Draw();
  // for(auto i : en)
  // {
  //   gaus1->SetParameters(0,0,0);
  //   gaus2->SetParameters(0,0,0);
  //   hTrigger[ channelMap[i] ]->Fit("gaus1", "qNR", "", 0.1, 10.0);
  //   hNoise[ channelMap[i] ]->Fit("gaus2", "qNR+");
  // }

  // TCanvas *c1 = new TCanvas("c1","Bob Ross's Canvas",800,600);
  // TF1 *fEff1 = new TF1("fEff", "0.5*(1+TMath::Erf((x-[0])/(TMath::Sqrt(2)*[1]) ))",0,3);
  // fEff1->SetParameters(1.3436,0.244441);
  // fEff1->Draw();
  //
  // c1->Print("./plots/htrigger.pdf");

  fOutput->Close();
  cout << "Thresholds found.\n";

}

vector<double> ThreshTrapezoidalFilter( const vector<double>& anInput, double RampTime, double FlatTime, double DecayTime )
{
  double decayConstant = 0.0;
  if(DecayTime != 0) decayConstant = 1./(exp(1./DecayTime) - 1);
  double rampStep = RampTime;
  double flatStep = FlatTime;
  double baseline = 0; // No baseline for now
  double norm = rampStep;
  if(decayConstant != 0)norm *= decayConstant;

  vector<double> fVector;
  vector<double> anOutput;
  if(fVector.size() != anInput.size()) {
    fVector.resize(anInput.size());
    anOutput.resize(anInput.size());
  }

  fVector[0] = anInput[0] - baseline;
  anOutput[0] = (decayConstant+1.)*(anInput[0] - baseline);
  double scratch = 0.0;
  for(size_t i = 1; i < anInput.size(); i++)
  {
    // This is a little tricky with all the ternary operators, but it's faster
    // this way.  We must check the bounds.
    scratch = anInput[i]  - ((i>=rampStep) ? anInput[i-rampStep] : baseline)
      - ((i>=flatStep+rampStep) ? anInput[i-flatStep-rampStep] : baseline)
      + ((i>=flatStep+2*rampStep) ? anInput[i-flatStep-2*rampStep] : baseline);

    if(decayConstant != 0.0) {
        fVector[i] = fVector[i-1] + scratch;
        anOutput[i] = (anOutput[i-1] + fVector[i] + decayConstant*scratch);
    }
    else anOutput[i] = anOutput[i-1] + scratch;
  }

  // cout << "sz:" << anInput.size() << endl;

  for(size_t i = 2*rampStep+flatStep; i < anInput.size(); i++)
    anOutput[i-(2*rampStep+flatStep)] = anOutput[i];

  anOutput.resize(anOutput.size()-(2*rampStep+flatStep));

  // Rescale event by normalization factor
  for(size_t i = 1; i < anOutput.size(); i++)anOutput[i] = anOutput[i]/norm;

  return anOutput;
}
