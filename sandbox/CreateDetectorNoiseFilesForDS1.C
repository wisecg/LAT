// Script originally from M. Buuck -- slightly changed cut to remove non-FT events
// Saves power spectra of each channel into a MGTWaveformFT
{
/* There are 4 runs with forced acquisition in DS1. We will use those runs to
 * build average Fourier Transforms for use with
 * MGDO/Transforms/MGWFAddNoiseFromFT.
 */

// Set up dataset and chain readers
// Only DS0, DS1, and DS4 have force trigger or delayed pulser runs
// GATDataSet ds(13071,13074);
GATDataSet ds(4201);

TChain* c=ds.GetChains();
vector<double>* channel=0;

// Micah uses trapEMinCal to distinguish some of the pulser events that snuck in
// I find it easier to use rawWFMax and rawWFMin
vector<double>* trapECal=0;
vector<double>* trapEMinCal=0;
std::vector<double> *rawWFMax = nullptr;
std::vector<double> *rawWFMin = nullptr;
MGTEvent* event=0;
c->SetBranchAddress("channel",&channel);
c->SetBranchAddress("event",&event);
c->SetBranchAddress("trapECal",&trapECal);
// c->SetBranchAddress("trapEMinCal",&trapEMinCal);
c->SetBranchAddress("rawWFMax",&rawWFMax);
c->SetBranchAddress("rawWFMin",&rawWFMin);
size_t nEntries=c->GetEntries();

// Need channel map and channel settings
MJTChannelMap* chmap=ds.GetChannelMap();
chmap->DumpChannelMap();
MJTChannelSettings* chset=ds.GetChannelSettings();

// DS1 and DS0 have different GRETINA card settings (actually only matters for run 4201)
vector<uint32_t> channelList=chset->GetEnabledIDList(); // DS1
// vector<uint32_t> channelList=chset->GetEnabledIDList("ORGretina4Model"); // DS0

vector<uint32_t> pulsers = chmap->GetPulserChanList();
//Remove pulsers from list of channels for these runs
for(auto pulser : pulsers) {
  channelList.erase(remove(channelList.begin(),channelList.end(),pulser));
}
//Remove P1D2LG and P7D3LG because they were not taking data
//channelList.erase(remove(channelList.begin(),channelList.end(),583));
//channelList.erase(remove(channelList.begin(),channelList.end(),595));
sort(channelList.begin(),channelList.end());

//We will store the FTs in a map indexed by channel
map<int,MGTWaveformFT> wfFTsByChan;
map<int,size_t> nGoodWfs;
MGTWaveformFT wfFT;
for(auto ch : channelList) {
  if(ch==595 || ch==583) continue; // Channels that weren't taking data
  //Set up the names of the FTs
  string wfFTName = "noiseWF_ch";
  wfFTName += to_string(ch);
  wfFTName += "_";
  wfFTName += chmap->GetDetectorName(ch);
  wfFTsByChan.emplace(piecewise_construct, make_tuple(ch), make_tuple());
  wfFTsByChan[ch].SetName(wfFTName.c_str());
  nGoodWfs[ch] = 0;
  cout << wfFTName.c_str() << endl;
}

//Do the FTs
MGWFFastFourierTransformFFTW fftw;
for(size_t iEntry=0; iEntry<nEntries; iEntry++) {
  c->GetEntry(iEntry);
  size_t nWfms=event->GetNWaveforms();
  if(iEntry%10000==0) cout << iEntry/double(nEntries)*100 << "\% done" << endl;
  for(size_t iWfm=0; iWfm<nWfms; iWfm++) {
    if(find(pulsers.begin(),pulsers.end(),channel->at(iWfm)) != pulsers.end()) continue;
    // if(abs(trapECal->at(iWfm))>5 || abs(trapEMinCal->at(iWfm))>5) continue; //Some pulser waveforms snuck into dataset
    if(rawWFMax->at(iWfm) - rawWFMin->at(iWfm) > 30) continue; //Some pulser waveforms snuck into dataset
    MGTWaveform* wf = event->GetWaveform(iWfm);
    wf->SetLength(2030); // Length of FT waveforms for DS0 is slightly longer
    fftw.PerformFFT(wf,&wfFT);
    wfFTsByChan[channel->at(iWfm)].MakeSimilarTo(wfFT,true);
    wfFTsByChan[channel->at(iWfm)].AddNormsSquared(wfFT);
    nGoodWfs[channel->at(iWfm)] += 1;
    // cout << nGoodWfs[channel->at(iWfm)] << endl;
  }
}

for(auto ch: pulsers)
{
  cout << ch << endl;
}

for(auto &ch: nGoodWfs)
{
  cout << ch.first << "\t" << ch.second << endl;
}

cout << "Writing to disk" << endl;
TFile file("/global/homes/b/bxyzhu/MJD/soft/noise_spectra_DS0.root","RECREATE");
for(auto ch : channelList) {
  if(ch==595 || ch==583) continue;
  // cout << ch << "\t" << nGoodWfs[ch] << endl;
  if(nGoodWfs[ch]==0) continue;
  wfFTsByChan[ch] /= sqrt(nGoodWfs[ch]);
  cout << "Divided channel " << ch << " by " << nGoodWfs[ch] << endl;
  wfFTsByChan[ch].Write();
}
}
