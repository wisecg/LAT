// MJD Data Set livetime & exposure calculator.
// C. Wiseman, USC
// v.1 10/6/2016
// v.2 3/7/2017

/*
  Results (8 Mar 2016)

  DS-0:
  Total live time (days): 47.6217
  Reduction from veto (sec): 3683.46
  Reduced live time (days): 47.5791
  Active mass:  14.6000 kg   enriched 10.6900   natural 3.9100
  Exposure:  694.6549 kg-days total  natural 186.0343   enriched 508.6206

  DS-1:
  Total live time (days): 60.1582
  Reduction from veto (sec): 3435.08
  Reduced live time (days): 60.1184
  Active mass:  12.4300 kg   enriched 11.3100   natural 1.1200
  Exposure:  747.2720 kg-days total  natural 67.3326   enriched 679.9394

  DS-2:
  Total live time (days): 9.65801
  Reduction from veto (sec): 640.136
  Reduced live time (days): 9.6506
  Active mass:  13.0220 kg   enriched 11.9010   natural 1.1210
  Exposure:  125.6701 kg-days total  natural 10.8183   enriched 114.8518

  DS-3:
  Total live time (days): 29.9128
  Reduction from veto (sec): 1117.99
  Reduced live time (days): 29.8999
  Active mass:  15.4120 kg   enriched 12.6310   natural 2.7810
  Exposure:  460.8174 kg-days total  natural 83.1516   enriched 377.6657

  DS-4:
  Total live time (days): 23.6908
  Reduction from veto (sec): 9986.62
  Reduced live time (days): 23.5752
  Active mass:  9.4212 kg   enriched 5.4712   natural 3.9500
  Exposure:  222.1064 kg-days total  natural 93.1219   enriched 128.9845

  DS-5: (no veto reduction ... yet)
  Total live time (days): 82.5177
  Reduction from veto (sec): 0
  Reduced live time (days): 82.5177
  Active mass:  25.6530 kg   enriched 18.3420   natural 7.3110
  Exposure:  2116.8253 kg-days total  natural 603.2865   enriched 1513.5387

*/

#include <iostream>
#include <map>
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "GATDataSet.hh"
#include "DataSetInfo.hh"
#include "MJVetoEvent.hh"

using namespace std;

double vetoReduction(GATDataSet &ds, int dsNum);

int main(int argc, char** argv)
{
  // Get user args
  if (argc < 2) {
		cout << "Usage: ./ds_livetime [dataset number]\n";
		return 1;
	}
	int dsNum = stoi(argv[1]);

  // This must match the run sequences in DataSetInfo.hh
  map<int,int> dsMap = {{0,76},{1,51},{2,7},{3,24},{4,22},{5,80}};

  // Total active mass, taken from skim_mjd_data.cc
  double enrMass=0, natMass=0;
  if (dsNum == 0)      { enrMass = 10.69;  natMass = 3.91;  }
  else if (dsNum == 1) { enrMass = 11.31;  natMass = 1.12;  }
  else if (dsNum == 2) { enrMass = 11.901; natMass = 1.121; }
  else if (dsNum == 3) { enrMass = 12.631; natMass = 2.781; }
  else if (dsNum == 4) { enrMass = 5.4712; natMass = 3.950; }
  else if (dsNum == 5) { enrMass = 12.631; natMass = 3.352;    // M1
                         enrMass += 5.711; natMass += 3.959; } // M2

  // Calculate total live time and the reduction from muon veto
  cout << "DS-" << dsNum << ":\n";
  GATDataSet ds;
  double totalLiveTime = 0, reducedLiveTime=0, vetoDeadTime=0;
  for (int i = 0; i <= dsMap[dsNum]; i++) {
    LoadDataSet(ds, dsNum, i);
    cout << i << " " << ds.GetNRuns() << endl;
  }
  cout << "Calculating live time ...\n";
  totalLiveTime = ds.GetRunTime()/1e9/86400;
  if (dsNum != 5) vetoDeadTime = vetoReduction(ds, dsNum);
  reducedLiveTime = totalLiveTime - vetoDeadTime/86400;
  cout << "Total live time (days): " << totalLiveTime << endl;
  cout << "Reduction from veto (sec): " << vetoDeadTime << endl;
  cout << "Reduced live time (days): " << reducedLiveTime << endl;

  // Calculate exposures
  double enrExpo = reducedLiveTime * enrMass;
  double natExpo = reducedLiveTime * natMass;
  double totExpo = enrExpo + natExpo;
  cout << Form("Active mass:  %.4f kg   enriched %.4f   natural %.4f\n",enrMass+natMass, enrMass, natMass)
       << Form("Exposure:  %.4f kg-days total  natural %.4f   enriched %.4f\n", totExpo, natExpo, enrExpo);
}

double vetoReduction(GATDataSet &ds, int dsNum)
{
  TChain *vetoChain = NULL;

  // Make the muon list, exactly the same way as we do in skim_mjd_data
  if(vetoChain==NULL && dsNum!=4) {
    vetoChain = ds.GetVetoChain();
    cout << "Found " << vetoChain->GetEntries() << " veto entries.  Creating muon list ...\n";
  }
  vector<int> muRuns;
  vector<int> muTypes;
  vector<double> muRunTStarts;
  vector<double> muTimes;
  vector<double> muUncert;
  if (dsNum != 4)
  {
    TTreeReader vetoReader(vetoChain);
    TTreeReaderValue<MJVetoEvent> vetoEventIn(vetoReader,"vetoEvent");
    TTreeReaderValue<int> vetoRunIn(vetoReader,"run");
  	TTreeReaderValue<Long64_t> vetoStart(vetoReader,"start");
  	TTreeReaderValue<Long64_t> vetoStop(vetoReader,"stop");
  	TTreeReaderValue<double> xTime(vetoReader,"xTime");
    TTreeReaderValue<double> timeUncert(vetoReader,"timeUncert");
  	TTreeReaderArray<int> CoinType(vetoReader,"CoinType");	//[32]
    bool newRun=false;
  	int prevRun=0;
  	Long64_t prevStop=0;
  	while(vetoReader.Next())
  	{
      MJVetoEvent veto = *vetoEventIn;
      int run = *vetoRunIn;
  		if (run != prevRun) newRun=true;
  		else newRun = false;
  		int type = 0;
  		if (CoinType[0]) type=1;
  		if (CoinType[1]) type=2;	// overrides type 1 if both are true
  		if ((*vetoStart-prevStop) > 10 && newRun) type = 3;
      if (type > 0){
        muRuns.push_back(run);
        muRunTStarts.push_back(*vetoStart);
        muTypes.push_back(type);
        if (type!=3) muTimes.push_back(*xTime);
        else muTimes.push_back(*xTime); // time of the first veto entry in the run
        if (!veto.GetBadScaler()) muUncert.push_back(*timeUncert);
        else muUncert.push_back(8.0); // uncertainty for corrupted scalers
      }
  		prevStop = *vetoStop;  // end of entry, save the run and stop time
  		prevRun = run;
  	}
  }
  else if (dsNum==4) LoadDS4MuonList(muRuns,muRunTStarts,muTimes,muTypes,muUncert);
  size_t nMu = muTimes.size();
  if(nMu == 0) {
    cout << "couldn't load mu data" << endl;
    return 0;
  }
  cout << "Muon list has " << muRuns.size() << " entries.\n";

  double deadTimeTotal = 0;
  for (int i = 0; i < (int)muRuns.size(); i++)
  {
    // this matches what's in skim_mjd_data
    double deadTime = 1. + 2 * fabs(muUncert[i]); // 1 second window, uncertainties on both ends.
    if (dsNum==4) deadTime = 4. + 4. * fabs(muUncert[i]);  // larger ~10s window for ds-4
    deadTimeTotal += deadTime;
    // printf("%i  %i  %i  %.0f  %.3f +/- %.3f  dead %.3f\n",i,muRuns[i],muTypes[i],muRunTStarts[i],muTimes[i],muUncert[i],deadTimeTotal);
    // FIXME: ds-5 muon uncertainties are being reported as negative
  }
  // cout << "Total veto dead time: " << deadTimeTotal << " seconds from " << muRuns.size() << " muon candidate events.\n";

  return deadTimeTotal;
}