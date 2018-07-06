// auto-thresh.cc
// Requires access to built data,
// output of process_mjd_cal (pass 1 gat)
// and the APDB (MkCookie must be run recently)
// B. Zhu, C. Wiseman
// v1. 2017/2/28
// v2. 2017/6/01 - changed to run over data subsets (as defined in DataSetInfo.hh)
// v3. 2018/4/06 - changed to run over specific runs.

#include <iostream>
#include <vector>
#include <string>
#include <map>
#include "TStyle.h"
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

void FindThresholds(int dsNum, int subNum, int runLo, int runHi, bool useDoubles, string outDir, bool bDebug);
vector<double> ThreshTrapezoidalFilter( const vector<double>& anInput, double RampTime, double FlatTime, double DecayTime );

int main(int argc, char** argv)
{
  if (argc < 3) {
    cout << "Usage: ./auto-thresh [ds] [sub]\n"
         << "   [-s [runLo] [runHi] : run limit mode (used to match TF boundaries)]\n"
         << "   [-d : use doubles in channel branches]\n"
         << "   [-o : specify output directory]\n"
         << "   [-v : verbose and debug mode, print out extra information and save histograms into 'plots' directory]\n";
    return 1;
  }
  int dsNum = stoi(argv[1]);
  int subNum = stoi(argv[2]);
  int runLo=-1, runHi=-1;
  bool useDoubles=0;
  string outDir = ".";
  bool bDebug = false;
  vector<string> opt(argv + 1, argv + argc);
  for (size_t i = 0; i < opt.size(); i++) {
    if (opt[i]=="-s") {
      runLo = stoi(opt[i+1]);
      runHi = stoi(opt[i+2]);
      cout << Form("Run limit mode active.  DS:%i Sub:%i runLo:%i runHi:%i\n",dsNum,subNum,runLo,runHi);
    }
    if (opt[i] == "-d") { useDoubles=1;    cout << "Using doubles in channel branches ...\n"; }
    if (opt[i] == "-o") { outDir=opt[i+1]; cout << "Writing to output directory: " << outDir << endl; }
    if (opt[i] == "-v") { bDebug = true; gStyle->SetOptStat(0); cout << "Verbose mode activated" << endl;}
  }

  // main routine
  FindThresholds(dsNum, subNum, runLo, runHi, useDoubles, outDir, bDebug);
}

void FindThresholds(int dsNum, int subNum, int runLo, int runHi, bool useDoubles, string outDir, bool bDebug)
{
  string outputFile = "";

  GATDataSet ds;
  if (runLo < 0 && runHi < 0) {
    LoadDataSet(ds, dsNum, subNum);
    outputFile = Form("%s/threshDS%d_%d.root", outDir.c_str(), dsNum, subNum);
  }
  else {
    GATDataSet tmp;
    LoadDataSet(tmp, dsNum, subNum);
    for (size_t i = 0; i < tmp.GetNRuns(); i++) {
      int run = tmp.GetRunNumber(i);
      if (run >= runLo && run <= runHi) {
        if(bDebug){cout << "Adding run " << run << endl;}
        ds.AddRunNumber(run);
      }
    }
    outputFile = Form("%s/threshDS%d_%d_%d_%d.root", outDir.c_str(), dsNum, subNum, runLo, runHi);
  }

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
  TTreeReaderArray<double> wfENF(gReader,"trapENM");

  // We have to get int/double right
  // you suck so much, doubles check.
  // TTreeReaderValue<double> runIn(gReader, "run");
  // TTreeReaderValue<vector<double>> wfChanIn(gReader, "channel");
  ROOT::Internal::TTreeReaderValueBase *wfChanIn = 0;
  ROOT::Internal::TTreeReaderValueBase *runIn = 0;
  if(!useDoubles) {
    wfChanIn = new TTreeReaderValue<vector<unsigned int>>(gReader, "channel");
    runIn = new TTreeReaderValue<int>(gReader,"run");
  }
  else {
    wfChanIn = new TTreeReaderValue<vector<double>>(gReader, "channel");
    runIn = new TTreeReaderValue<double>(gReader,"run");
  }

  int nEntries = gatChain->GetEntries();
  cout << "Scanning DS " << dsNum << " subNum " << subNum << " , " << nEntries << " entries.\n";

  // Initialize channel map and histogram vectors
  map<int,int> channelMap;
  vector<TH1D*> hTrigger;
  vector<TH1D*> hNoise;
  vector<TCanvas*> cDebug;
  vector<TCanvas*> cDebug2;
  TH1D *waveTrigger;
  TH1D *waveNoise;
  int nChannel = en.size();
  int nTrigger[nChannel];
  int nNoise[nChannel];
  for(int i = 0; i < nChannel; i++) {
    channelMap.insert({en[i], i});
    hTrigger.push_back(new TH1D(Form("hTrigger-%d", en[i]), Form("hTrigger-%d", en[i]), 1000, -30, 30));
    hNoise.push_back(new TH1D(Form("hNoise-%d", en[i]), Form("hNoise-%d", en[i]), 500, -30, 30));
    nTrigger[i] = 0;
    nNoise[i] = 0;
    if(bDebug)
    {
      // Only look at first
      cDebug.push_back(new TCanvas(Form("cDebug_ch%d", en[i]),"",800,600));
      cDebug2.push_back(new TCanvas(Form("cDebug2_ch%d", en[i]),"",800,600));
    }
  }

  // Loop over entries
  int channel = 0;
  double NoiseSample = 0, TriggerSample = 0, trapENF = 0;
  vector<double> TrapFilter;
  vector<double> Waveform;
  cout << "Sampling from waveforms ...\n";

  // Set minimum run as run from the first entry
  gReader.SetEntry(0);

  if (!useDoubles)
    runMin = **static_cast<TTreeReaderValue<int>*>(runIn);
  else
    runMin = (int)**static_cast<TTreeReaderValue<double>*>(runIn);


  for(int i = 0; i < nEntries; i++)
  {
    if(bDebug && i >= 5000000) break;
    bReader.SetEntry(i);
    gReader.SetEntry(i);
    int nWF = (*wfBranch).GetEntriesFast();

    if (fmod(100*(double)i/nEntries, 10.0) < 0.00001)
      cout << 100*(double)i/nEntries << " % done.\n";

    // Get channel list for event as vector<int>.
    // i hate you so much, doubles check
    vector<int> wfChan;
    if(useDoubles) {
      vector<double>& chList = **static_cast<TTreeReaderValue< vector<double> >* >(wfChanIn);
      for(double chan : chList) wfChan.push_back(int(chan));
    }
    else {
      vector<unsigned int>& chList = **static_cast<TTreeReaderValue< vector<unsigned int> >* >(wfChanIn);
      for(unsigned int chan : chList) wfChan.push_back(int(chan));
    }

    int nCh = wfChan.size();
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
        if(nTrigger[channelMap[channel]] < 100)
        {
          waveTrigger = new TH1D(Form("wtrig_ch%d_%d", channel, nTrigger[channelMap[channel]] ),"", 30, 0, 30);
          for(int bin = 1; bin < waveTrigger->GetNbinsX(); bin++)
          {
            waveTrigger->SetBinContent(bin, TrapFilter[bin+1]);
          }
          cDebug[channelMap[channel]]->cd();
          waveTrigger->SetLineColorAlpha(kBlue, 0.5);
          waveTrigger->Draw("SAME");
        }
      }

      // High energy = sharp rise => 1st sample good representation of noise
      if(trapENF > 50) {
        hNoise[ channelMap[channel] ]->Fill(NoiseSample);
        nNoise[ channelMap[channel] ]++;
        if(nNoise[channelMap[channel]] < 100)
        {
          waveNoise = new TH1D(Form("wnoise_ch%d_%d", channel, nNoise[channelMap[channel]] ),"", 30, 0, 30);
          for(int bin = 1; bin < waveNoise->GetNbinsX(); bin++)
          {
            waveNoise->SetBinContent(bin, TrapFilter[bin+1]);
          }
          cDebug[channelMap[channel]]->cd();
          waveNoise->SetLineColorAlpha(kRed, 0.5);
          waveNoise->Draw("SAME");
        }
      }
    }
    if(bDebug){ if (i % 100000 == 0) cout << i << " entries scanned so far.\n";}
  }

  // Set maximum run as run from the last entry
  gReader.SetEntry(nEntries-1);
  if (!useDoubles)
    runMax = **static_cast<TTreeReaderValue<int>*>(runIn);
  else
    runMax = (int)**static_cast<TTreeReaderValue<double>*>(runIn);

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

  if(bDebug)
  {
    for(auto i: channelMap)
    {
      if(i.first%2 == 0)
      {
        cDebug2[i.second]->cd();
        hNoise[i.second]->SetLineColor(kBlue);
        hTrigger[i.second]->SetLineColor(kRed);
        hNoise[i.second]->GetXaxis()->SetRangeUser(-5*sigmaADC[i.second], threshADC[i.second] + 5*sigmaADC[i.second]);
        hNoise[i.second]->DrawNormalized();
        hTrigger[i.second]->DrawNormalized("SAME");
        cDebug2[i.second]->SaveAs(Form("%s/plots/hTrigger_ch%d.pdf",outDir.c_str(), i.first));
        cDebug[i.second]->SaveAs(Form("%s/plots/hWaves_ch%d.pdf",outDir.c_str(), i.first));
      }
    }
  }
  // If either fit failed, put the threshold at 99999 keV
  for (size_t i = 0; i < threshADC.size(); i++) {
    if (threshFitStatus[i] != 0 || sigmaFitStatus[i] != 0) {
      threshCal[i] = 99999;
      sigmaCal[i] = 99999;
    }
    cout << Form("Ch %i  keV %.3f +/- %-8.3f   ADC %.2e +/- %-8.2e   Noise %.2e +/- %-8.2e  nADC %-5i  nNoise %-5i  Fit %i %i\n",
      channelList[i], threshCal[i], sigmaCal[i], threshADC[i], threshADCErr[i], sigmaADC[i], sigmaADCErr[i],
      nTrigger[i], nNoise[i], threshFitStatus[i], sigmaFitStatus[i]); //, CalScale[i], CalOffset[i]);
  }
  fThreshTree->Fill();
  fThreshTree->Write("",TObject::kOverwrite);

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
