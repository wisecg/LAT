// thresholds.cc
// Examines low-energy performance of MJD detectors.
// Requires MkCookie be run recently.
// B. Zhu (LANL) and C. Wiseman (USC)
// 1/27/2016

#include <iostream>
#include <vector>
#include <string>
#include <map>

#include "TTree.h"
#include "TFile.h"
#include "TROOT.h"
#include "TH1D.h"
#include "TStyle.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TBenchmark.h"
#include "TEntryList.h"
#include "TGraph.h"
#include "TLegend.h"

#include "GATDataSet.hh"
#include "DataSetInfo.hh"
#include "MGTWaveform.hh"
#include "MJDBUtilities.hh"
#include "MJDatabase.hh"
#include "MJCouchAccess.hh"
#include "MJDocument.hh"
#include "MJProvenance.hh"
#include "MJAnalysisDoc.hh"
#include "MJAnalysisParameters.hh"
#include "MJTRun.hh"

using namespace std;
using namespace katrin;
using namespace MJDB;

void RunDiagnostics();
void FindThresholds(vector<int> runs, int dsNum, int subNum=-1);
vector<double> ThreshTrapezoidalFilter( const vector<double>& anInput, double RampTime, double FlatTime, double DecayTime );
void ThresholdVsRun(string inputFile, string outputDir);
void FindTotalExposure();
void FindExposure(string inputFile, int dsNum);
void UpdateSkimFile();

int main(int argc, char** argv)
{
  if (argc < 2) {
    cout << "Usage: ./thresholds [-d [dsNum] [subNum] ]\n"
         << "                    [-l [run list file] ]\n"
         << "                    [-r [runNum] ]\n"
         << "                    [-t [thresh file] (threshold vs. run) ]\n"
         << "                    [-e [thresh file] [dsNum] (threshold vs exposure) ]\n"
         << "                    [-g (run diagnostics)]\n"
         << "                    [-u (update skim file)]\n";
    return 1;
  }
  bool d=0, l=0, r=0, t=0, e=0, g=0, u=0, f=0;
  string inputList, threshFile;
  string outputDir = "./plots/";
  int dsNum=0, subNum=0;
  vector<string> opt(argv + 1, argv + argc);
  for (size_t i = 0; i < opt.size(); i++) {
    if (opt[i] == "-d") { d=1; dsNum = stoi(opt[i+1]); subNum = stoi(opt[i+2]); }
    if (opt[i] == "-l") { l=1; inputList = opt[i+1]; }
    if (opt[i] == "-r") { r=1; subNum = stoi(opt[i+1]); }
    if (opt[i] == "-t") { t=1; threshFile = opt[i+1]; }
    if (opt[i] == "-e") { e=1; threshFile = opt[i+1]; dsNum = stoi(opt[i+2]); }
    if (opt[i] == "-g") { g=1; }
    if (opt[i] == "-u") { u=1; }
    if (opt[i] == "-f") { f=1; }
  }
  gROOT->ProcessLine( "gErrorIgnoreLevel = 2001;"); // suppress ROOT fit warning messages
  gROOT->ProcessLine(".x ~/env/MJDClintPlotStyle.C"); // load plot style
  // ================================================================

  // Scan thresholds for a subset of a DS
  if (d) {
    GATDataSet ds;
    vector<int> runs;
    LoadDataSet(ds,dsNum,subNum);
    for (size_t i = 0; i < ds.GetNRuns(); i++) runs.push_back(ds.GetRunNumber(i));
    // for (auto i : runs) cout << i << " "; cout << endl;
    FindThresholds(runs,dsNum,subNum);
  }
  // Scan thresholds for a list of runs (separate files)
  else if (l) {
    ifstream runFile(inputList);
    int rundummy;
    while (runFile >> rundummy) {
      vector<int> runs = {rundummy};
      FindThresholds(runs, rundummy);
    }
  }
  // Scan thresholds for a single run
  else if (r) {
    vector<int> runs = {subNum};
    FindThresholds(runs, subNum);
  }

  // Plot thresholds vs run number for each channel.
  else if (t) ThresholdVsRun(threshFile, outputDir);

  // Sanity Check: Estimate total DS exposure (no energy cuts).
  else if (f) FindTotalExposure();

  // Examine the exposure vs energy threshold value.
  else if (e) FindExposure(threshFile,dsNum);

  // Update a skim file with a threshold branch.
  else if (u) UpdateSkimFile();

  // Various diagnostics
  else if (g) RunDiagnostics();

  // ================================================================
}

// Option: -d [dsNum] [subNum], -l [file], -r [run]
void FindThresholds(vector<int> runs, int dsNum, int subNum)
{
  string outputFile;
  if (subNum != -1) outputFile = Form("./threshFiles/thresholdsDS%i_%i.root",dsNum,subNum);
  else outputFile = Form("./threshFiles/thresholds_run%i.root",dsNum);

  // Set up output file
  TFile *fOutput = new TFile(outputFile.c_str(), "RECREATE");
  TTree *fThreshTree = new TTree("threshTree", "Tree with vector of channels and thresholds");
  int fRun;
  double duration;
  vector<int> channelList;
  vector<double> threshADC;
  vector<double> threshADCErr;
  vector<double> sigmaADC;
  vector<double> sigmaADCErr;
  vector<double> threshCal;
  vector<double> sigmaCal;
  vector<int> threshFitStatus;
  vector<int> sigmaFitStatus;
  vector<double> CalOffset;
  vector<double> CalScale;
  vector<int> numTrigger;
  vector<int> numNoise;
  fThreshTree->Branch("run", &fRun);
  fThreshTree->Branch("duration",&duration);
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

  TBenchmark b;
  for (auto run : runs)
  {
    b.Start("runTimer");

    fRun = run;
    channelList.resize(0);
    threshADC.resize(0);
    threshADCErr.resize(0);
    sigmaADC.resize(0);
    sigmaADCErr.resize(0);
    threshCal.resize(0);
    sigmaCal.resize(0);
    threshFitStatus.resize(0);
    sigmaFitStatus.resize(0);
    CalOffset.resize(0);
    CalScale.resize(0);
    numTrigger.resize(0);
    numNoise.resize(0);

    // Access data for this run
    GATDataSet ds(run);
    duration = ds.GetRunTime()/1E9;
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
    TChain *gatChain = ds.GetGatifiedChain(false);
    TTreeReader gReader(gatChain);
    TTreeReaderArray<double> wfChan(gReader,"channel");
    TTreeReaderArray<double> wfENF(gReader,"trapENM");
    int nEntries = gatChain->GetEntries();
    cout << "Scanning run " << run << ", " << nEntries << " entries.\n";

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
      if (i % 10000 == 0) cout << i << " entries saved so far.\n";
    }

    // Query database to get calibration parameters for enabled channels.
    // (The DB is only accessed for the first channel, then saves the run info in a buffer.)
    // Have you run MkCookie?
    // Should you put in a check that the enabled channels are all found in the DB?

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

      threshFitStatus.push_back(hTrigger[ channelMap[i] ]->Fit("gaus1", "qNR", "", 0.1, 10.0) );
      sigmaFitStatus.push_back(hNoise[ channelMap[i] ]->Fit("gaus2", "qNR+") );
      numTrigger.push_back(nTrigger[ channelMap[i] ]);
      numNoise.push_back(nNoise[ channelMap[i] ]);

      channelList.push_back(i);
      threshADC.push_back(gaus1->GetParameter(1));
      sigmaADC.push_back(gaus2->GetParameter(2));
      threshADCErr.push_back(gaus1->GetParError(1));
      sigmaADCErr.push_back(gaus2->GetParError(1));

      size_t Length = findResult.GetAnalysisParameter(run, i, mycalibration.GetPSource(), mycalibration.GetPType());

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

    // If either fit failed, put the threshold at 9999 keV
    for (size_t i = 0; i < threshADC.size(); i++) {
      if (threshFitStatus[i] != 0 || sigmaFitStatus[i] != 0) {
        threshCal[i] = 9999;
        sigmaCal[i] = 9999;
      }
      cout << Form("Ch %i  keV %.3f +/- %-8.3f   ADC %.2e +/- %-8.2e   Noise %.2e +/- %-8.2e  nADC %-5i  nNoise %-5i  Fit %i %i\n",
        channelList[i], threshCal[i], sigmaCal[i], threshADC[i], threshADCErr[i], sigmaADC[i], sigmaADCErr[i],
        nTrigger[i], nNoise[i], threshFitStatus[i], sigmaFitStatus[i]); //, CalScale[i], CalOffset[i]);
    }
    fThreshTree->Fill();
    fThreshTree->Write("",TObject::kOverwrite);
    b.Show("runTimer");
    b.Reset();

  }
  fOutput->Close();
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

  for(size_t i = 2*rampStep+flatStep; i < anInput.size(); i++)
    anOutput[i-(2*rampStep+flatStep)] = anOutput[i];

  anOutput.resize(anOutput.size()-(2*rampStep+flatStep));

  // Rescale event by normalization factor
  for(size_t i = 1; i < anOutput.size(); i++)anOutput[i] = anOutput[i]/norm;

  return anOutput;
}

// Option: -t [file]
void ThresholdVsRun(string inputFile, string outputDir)
{
  // Load the input file
  TFile *inFile = new TFile(inputFile.c_str());
  TTree *threshTree = (TTree*)inFile->Get("threshTree");
	TTreeReader reader(threshTree);
	TTreeReaderValue<Int_t> run(reader, "run");
 	TTreeReaderArray<int> channelList(reader, "channelList");
 	TTreeReaderArray<double> threshADC(reader, "threshADC");
 	TTreeReaderArray<double> sigmaADC(reader, "sigmaADC");
 	TTreeReaderArray<double> threshADCErr(reader, "threshADCErr");
 	TTreeReaderArray<double> sigmaADCErr(reader, "sigmaADCErr");
 	TTreeReaderArray<double> threshCal(reader, "threshCal");
 	TTreeReaderArray<double> sigmaCal(reader, "sigmaCal");
 	TTreeReaderArray<int> threshFitStatus(reader, "threshFitStatus");
 	TTreeReaderArray<int> sigmaFitStatus(reader, "sigmaFitStatus");
 	TTreeReaderArray<double> CalOffset(reader, "CalOffset");
 	TTreeReaderArray<double> CalScale(reader, "CalScale");
 	TTreeReaderArray<int> numTrigger(reader, "numTrigger");
 	TTreeReaderArray<int> numNoise(reader, "numNoise");

  // Make a list of all unique channels in the input file,
  // and a map so we always write to the correct histogram index
  set<int> uniqueChans;
  map<int,int> chanMap;
  while (reader.Next()) {
    for (size_t i = 0; i < channelList.GetSize(); i++)
      uniqueChans.insert(channelList[i]);
  }
  vector<int> fullChanList(uniqueChans.begin(), uniqueChans.end());
  sort(fullChanList.begin(), fullChanList.end());

  cout << "Found " << fullChanList.size() << " unique channels:\n";
  for (size_t i = 0; i < fullChanList.size(); i++) {
    chanMap.insert( {fullChanList[i], i} );  // { key, value }
    cout << fullChanList[i] << " ";
  }
  cout << endl;

  // Create a histogram for every unique channel
  int nRuns = threshTree->GetEntries();
  vector<TH1D*> hThreshold;
  vector<TH1D*> hThresholdPS;
  vector<TH1D*> hThresholdNS;
  vector<TH1D*> hSigma;
	for(size_t i = 0; i < fullChanList.size(); i++) {
		hThreshold.push_back(new TH1D(Form("hThreshold-ch%d", fullChanList[i]), Form("hThreshold-ch%d", fullChanList[i]), nRuns, 0, nRuns));
		hThresholdPS.push_back(new TH1D(Form("hThresholdPS-ch%d", fullChanList[i]), Form("hThresholdPS-ch%d", fullChanList[i]), nRuns, 0, nRuns));
		hThresholdNS.push_back(new TH1D(Form("hThresholdNS-ch%d", fullChanList[i]), Form("hThresholdNS-ch%d", fullChanList[i]), nRuns, 0, nRuns));
		hSigma.push_back(new TH1D(Form("hSigma-ch%d", fullChanList[i]), Form("hSigma-ch%d", fullChanList[i]), nRuns, 0, nRuns));
	}

  // Reset the reader and loop over entries
  reader.SetTree(threshTree);
	for(int i = 0; i < nRuns; i++)
	{
		reader.SetEntry(i);
    int nChannels = channelList.GetSize();
  	cout << "Run: " << *run << "\t Channels: " << nChannels << endl;

    // Loop over channels
    for (size_t j=0; j < channelList.GetSize(); j++)
    {
      int k = chanMap[ channelList[j] ];  // histogram index
			if((threshCal[j] != threshCal[j]) || (sigmaCal[j] != sigmaCal[j])) {
				hThreshold[k]->SetBinContent(i+1, 0);
				hThresholdPS[k]->SetBinContent(i+1, 0);
				hThresholdNS[k]->SetBinContent(i+1, 0);
				hSigma[k]->SetBinContent(i+1, 0);
			}
			else {
				hThreshold[k]->SetBinContent(i+1, threshCal[j]);
				hThresholdPS[k]->SetBinContent(i+1, threshCal[j]+sigmaCal[j]);
				hThresholdNS[k]->SetBinContent(i+1, threshCal[j]-sigmaCal[j]);
				hSigma[k]->SetBinContent(i+1, sigmaCal[j]);;
			}
			if(i%25==0) {
				hThreshold[k]->GetXaxis()->SetBinLabel(i+1, Form("%d", *run) );
				hSigma[k]->GetXaxis()->SetBinLabel(i+1, Form("%d", *run));
			}
		}
	}
  // Plot thresholds vs. run
	TCanvas *c1 = new TCanvas("c1","Bob Ross's Canvas",1200,800);
	for(size_t i = 0; i < fullChanList.size(); i++) {
		hThreshold[i]->GetYaxis()->SetRangeUser(0, 7.5);
		hThreshold[i]->SetLineColor(kBlue);
		hThresholdPS[i]->SetLineColor(kRed);
		hThresholdPS[i]->SetLineStyle(2);
		hThresholdNS[i]->SetLineColor(kRed);
		hThresholdNS[i]->SetLineStyle(2);
		hThreshold[i]->Draw();
		hThresholdPS[i]->Draw("SAME");
		hThresholdNS[i]->Draw("SAME");
    c1->Print(Form("%s/%s.pdf",outputDir.c_str(),hThreshold[i]->GetTitle()));
	}
	// for(size_t i = 0; i < fullChanList.size(); i++) {
		// hSigma[i]->SetLineColor(kBlue);
		// hSigma[i]->Draw();
    // c1->Print(Form("%s/%s.pdf",outputDir.c_str(),hSigma[i]->GetTitle()));
	// }
  cout << "Printed up some pretty plots.\n";
}

// Option: -f
void FindTotalExposure()
{
  // Are you sure this is going to work?  How are you going to get the active detectors?
  // Make a list from the skim file?  Should you make the ActiveMasses map global?
  // What if the channel mapping changes during the dataset?

  int dsNum = 5;
  map<int,int> dsMap = {{0,76},{1,51},{3,24},{4,22},{5,80}};

  double totalLiveTime = 0;
  for (int i = 0; i <= dsMap[dsNum]; i++) {
    cout << "Loading DS-" << dsNum << " run sequence " << i << endl;
    GATDataSet ds;
    LoadDataSet(ds, dsNum, i);
    totalLiveTime += ds.GetRunTime()/1e9/86400;
  }
  cout << "Livetime: " << totalLiveTime << " days." << endl;
}

// Option: -e [file]
void FindExposure(string inputFile, int dsNum)
{
  // Active masses in kg, from Micah's document:
  // http://mjwiki.npl.washington.edu/pub/Majorana/AnalysisReports/ActiveMassCalcWithM1AndM2.pdf
  map<string,double> activeMasses = { {"C1P1D1",0.510}, {"C1P1D2",0.979}, {"C1P1D3",0.811}, {"C1P1D4",0.968}, {"C1P2D1",0.560}, {"C1P2D2",0.723}, {"C1P2D3",0.659}, {"C1P2D4",0.689}, {"C1P3D1",0.551}, {"C1P3D2",0.886}, {"C1P3D3",0.949}, {"C1P3D4",1.024}, {"C1P4D1",0.558}, {"C1P4D2",0.564}, {"C1P4D3",0.567}, {"C1P4D4",0.545}, {"C1P4D5",0.557}, {"C1P5D1",0.553}, {"C1P5D2",0.730}, {"C1P5D3",0.632}, {"C1P5D4",0.982}, {"C1P6D1",0.732}, {"C1P6D2",0.675}, {"C1P6D3",0.701}, {"C1P6D4",0.5722}, {"C1P7D1",0.561}, {"C1P7D2",0.710}, {"C1P7D3",0.5908}, {"C1P7D4",0.964}, {"C2P1D1",0.556}, {"C2P1D2",0.576}, {"C2P1D3",0.903}, {"C2P1D4",0.917}, {"C2P2D1",0.581}, {"C2P2D2",0.562}, {"C2P2D3",0.559}, {"C2P2D4",0.558}, {"C2P2D5",0.577}, {"C2P3D1",0.872}, {"C2P3D2",0.852}, {"C2P3D3",0.996}, {"C2P4D1",0.558}, {"C2P4D2",0.579}, {"C2P4D3",0.565}, {"C2P4D4",0.566}, {"C2P4D5",0.562}, {"C2P5D1",0.557}, {"C2P5D2",0.591}, {"C2P5D3",1.031}, {"C2P5D4",0.802}, {"C2P6D1",0.4622}, {"C2P6D2",0.775}, {"C2P6D3",0.821}, {"C2P6D4",0.778}, {"C2P7D1",0.566}, {"C2P7D2",0.968}, {"C2P7D3",0.562}, {"C2P7D4",0.567} };

  // vector<double> floors = {0.5,0.4};
  vector<double> floors = {100.0, 10.0, 7.5, 5.0, 4.0, 3.5, 3.0, 2.5, 2.0, 1.5, 1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1};
  vector<double> exposures(floors.size(),0);

  // Load threshold file
  TFile *f1 = new TFile(inputFile.c_str(),"UPDATE");  // open in update mode to save the exposure graph
  TTree *threshTree = (TTree*)f1->Get("threshTree");
  int run = 0;
  double duration = 0;
  vector<double> *threshKeV=0;
  threshTree->SetBranchAddress("run",&run);
  threshTree->SetBranchAddress("threshCal",&threshKeV);
  threshTree->SetBranchAddress("duration",&duration);

  for (size_t f = 0; f < floors.size(); f++)
  {
    string theCut = Form("threshCal < %.2f && threshCal > 0.1",floors[f]);

    // Assume the channel map is the same throughout this dataset.
    size_t n = threshTree->Draw("run",theCut.c_str(),"GOFF");
    if (n==0) continue;
    double *vRuns = threshTree->GetV1();
    GATDataSet ds(vRuns[0]);
    MJTChannelMap *map = ds.GetChannelMap();

    // Apply an entry list and start the loop
    string eListName = Form("elist_%.1f",floors[f]);
    threshTree->Draw(Form(">>%s",eListName.c_str()),theCut.c_str(), "entrylist");
    TEntryList *elist = (TEntryList*)gDirectory->Get(eListName.c_str());
    threshTree->SetEntryList(elist);
    for (size_t i = 0; i < (size_t)elist->GetN(); i++)
    {
      threshTree->GetEntry(i);
      string cut = Form("threshCal < %.2f && threshCal >= 0.1 && channelList %% 2 == 0", floors[f]);
      size_t n = threshTree->Draw("channelList:threshCal",cut.c_str(),"GOFF",1,i);
      if (n==0) continue;
      double* lChan = threshTree->GetV1();
      double* lThresh = threshTree->GetV2();
      vector<double> foundChans;
      double threshAvg = 0;
      for (size_t i = 0; i < n; i++) {
        string pos = map->GetDetectorPos(lChan[i]);
        foundChans.push_back(lChan[i]);
        exposures[f] += (duration/86400) * activeMasses[pos];
        threshAvg += lThresh[i];
      }
      threshAvg = threshAvg / (double)foundChans.size();

      // cout << Form("run %i  duration %.0f  exp (kg-d) %.2f  floor %.1f  avg %.2f  %lu chans: ",run,duration,exposures[f]/86400,floors[f],threshAvg,foundChans.size());
      // for (size_t j = 0; j < foundChans.size(); j++) cout << foundChans[j] << " ";
      // cout << endl;
    }
    cout << "Exposure for " << floors[f] << " keV floor: " << exposures[f]<< " kg * days.\n";
  }

  // Print final results
  string outputFile = Form("./plots/Exposure_DS%i.pdf",dsNum);
  TCanvas *c1 = new TCanvas("c1","Bob Ross's Canvas",800,600);
  c1->SetLogy(1);
  c1->SetLogx(1);
  TGraph *g1 = new TGraph(floors.size(),&(floors[0]),&(exposures[0]));
  g1->Write("",TObject::kOverwrite);  // save to thresholds file
  g1->SetMarkerStyle(kFullDotLarge);
  g1->SetMarkerColor(kRed);
  g1->SetLineColorAlpha(kBlue,0.5);
  g1->SetLineWidth(2);
  g1->GetXaxis()->SetTitle("Threshold (keV)");
  g1->GetYaxis()->SetTitle("Exposure (kg-days)");
  g1->Draw("ALP");
  c1->Print(outputFile.c_str());
  c1->SetLogy(0);
  c1->SetLogx(0);
  g1->Draw("ALP");
  c1->Print(TString::Format("./plots/Exposure_DS%i_lin.pdf",dsNum));

  f1->Close();
}

// Option: -g
void RunDiagnostics()
{
  // just a quick check to make sure we have all the branches we need
  // int dsNum=0, subNum=54;
  // GATDataSet dsList;
  // LoadDataSet(dsList,dsNum,subNum);
  // for (size_t i = 0; i < dsList.GetNRuns(); i++)
  // {
  //   int run = dsList.GetRunNumber(i);
  //
  //   GATDataSet ds(run);
  //   TChain *gatChain = ds.GetGatifiedChain(false);
  //   int nEntries = gatChain->GetEntries();
  //   cout << "Scanning run " << run << ", " << nEntries << " entries.\n";
  //
  //   static TString invalidBranch("trapENF");
  //   TBranch* br = (TBranch*)gatChain->GetListOfBranches()->FindObject(invalidBranch);
  //   if (!br) cout << "trapENF is dead\n";
  //
  //   TTreeReader gReader(gatChain);
  //   TTreeReaderArray<double> wfENF(gReader,"trapENF");   // why is this not trapENFCal?
  // }

  // Another quick check to see how many total runs are in a DS
  // Used this to compare to the final threshTree entry list
  // to make sure I didn't miss any runs.
  // int dsNum = 5;
  // map<int,int> dsMap = {{0,76},{1,51},{3,24},{4,22},{5,80}};
  // GATDataSet ds;
  // for (int i = 0; i <= dsMap[dsNum]; i++) LoadDataSet(ds,dsNum,i);
  // cout << "DS-" << dsNum << " runs: " << ds.GetNRuns() << endl;

  // Combine the Exposure TGraphs into one plot.
  TFile *f0 = new TFile("./final/thresholdsDS0.root");
  TFile *f1 = new TFile("./final/thresholdsDS1.root");
  TFile *f3 = new TFile("./final/thresholdsDS3.root");
  TFile *f4 = new TFile("./final/thresholdsDS4.root");
  TFile *f5 = new TFile("./final/thresholdsDS5.root");
  TGraph *g0 = (TGraph*)f0->Get("Graph");
  TGraph *g1 = (TGraph*)f1->Get("Graph");
  TGraph *g3 = (TGraph*)f3->Get("Graph");
  TGraph *g4 = (TGraph*)f4->Get("Graph");
  TGraph *g5 = (TGraph*)f5->Get("Graph");

  TCanvas *c1 = new TCanvas("c1","Bob Ross's Canvas",800,600);
  c1->SetLogx(1);
  g5->SetMarkerStyle(kFullDotLarge);
  g5->SetMarkerColor(kRed);
  g5->SetLineColorAlpha(kRed,0.5);
  g5->GetXaxis()->SetTitle("Threshold (keV)");
  g5->GetXaxis()->SetTitleOffset(1.1);
  g5->GetYaxis()->SetTitle("Exposure (kg-days)");
  g5->GetYaxis()->SetTitleOffset(1.2);
  g5->Draw("ALP");
  g0->SetMarkerStyle(kFullDotLarge);
  g0->SetMarkerColor(kBlue);
  g0->SetLineColorAlpha(kBlue,0.5);
  g0->Draw("SAME PLC");
  g1->SetMarkerStyle(kFullDotLarge);
  g1->SetMarkerColor(kGreen);
  g1->SetLineColorAlpha(kGreen,0.5);
  g1->Draw("SAME PLC");
  g3->SetMarkerStyle(kFullDotLarge);
  g3->SetMarkerColor(kMagenta);
  g3->SetLineColorAlpha(kMagenta,0.5);
  g3->Draw("SAME PLC");
  g4->SetMarkerStyle(kFullDotLarge);
  g4->SetMarkerColor(kOrange);
  g4->SetLineColorAlpha(kOrange,0.5);
  g4->Draw("SAME PLC");

  TLegend* leg1 = new TLegend(0.15,0.55,0.35,0.9);
	leg1->AddEntry(g0,"DS-0","l");
	leg1->AddEntry(g1,"DS-1","l");
	leg1->AddEntry(g3,"DS-3","l");
	leg1->AddEntry(g4,"DS-4","l");
	leg1->AddEntry(g5,"DS-5","l");
	leg1->Draw("SAME");
	c1->Update();

  c1->Print("./plots/CombinedExposure.pdf");
  c1->Print("./plots/CombinedExposure.png");
}

// Option: -u
void UpdateSkimFile()
{
    TFile *f2 = new TFile("./final/thresholdsDS3.root");
    TTree *threshTree = (TTree*)f2->Get("threshTree");
    TTreeReader threshReader(threshTree);
    TTreeReaderValue<int> runTh(threshReader, "run");
   	TTreeReaderArray<int> channelList(threshReader, "channelList");
   	TTreeReaderArray<double> threshCal(threshReader, "threshCal");
   	TTreeReaderArray<double> sigmaCal(threshReader, "sigmaCal");

    // Must use TFile, not TChain
    TFile *skimFile = new TFile("~/datasets/skim/skimDS3_0.root","UPDATE");
    TTree *skimTree = (TTree*)skimFile->Get("skimTree");

    int run=0;
    vector<int> *channel=0;
    vector<double> *thresh=0;
    vector<double> *threshSig=0;
    skimTree->SetBranchAddress("run",&run);
    skimTree->SetBranchAddress("channel",&channel);
    TBranch *thr = skimTree->Branch("threshKeV",&thresh);
    TBranch *sig = skimTree->Branch("threshSig",&threshSig);

    int prevRun = 0;
    map<int,int> threshMap;
    for (size_t i = 0; i < (size_t)skimTree->GetEntries(); i++)
    {
      skimTree->GetEntry(i);

      if (run!=prevRun)
      {
        skimTree->Write("",TObject::kOverwrite);

        int n = threshTree->Draw("Entry$",Form("run==%i",run),"GOFF");
        if (n==0) {
          cout << "Warning: No threshold data found for run " << run << ". Quitting ..." << endl;
          break;
        }
        double *lEntry = threshTree->GetV1();
        size_t thisEntry = (size_t)lEntry[0];
        threshReader.SetEntry(thisEntry);
        cout << "Found run " << run << endl;

        // Map channel to index -- threshCal and sigmaCal will have the same index
        threshMap.clear();
        for (size_t j = 0; j < channelList.GetSize(); j++) threshMap[ channelList[j] ] = j;
      }

      // Fill skim file threshold vectors
      thresh->resize(0);
      threshSig->resize(0);
      for (size_t j = 0; j < channel->size(); j++)
      {
        int chan = channel->at(j);
        double t = threshCal[ threshMap[chan] ];
        double s = sigmaCal[ threshMap[chan] ];
        thresh->push_back(t);
        threshSig->push_back(s);
      }
      thr->Fill();
      sig->Fill();

      // Save run for next entry
      prevRun=run;
    }
    cout << "skim Entries: " << skimTree->GetEntries() << "  threshBranch entries " << thr->GetEntries() << endl;

    skimTree->Write("",TObject::kOverwrite);
    skimFile->Close();
}
