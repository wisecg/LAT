// Takes a sample file (with MGTWaveforms and parameters amplitude, rval, zval) and turns it into fake built file
// Then you can take the fake built file and gatify it... 
// Now we can make fake data to reprocess!
void GenerateMGTEvent(std::string InputFile)
{
	std::vector<MGTWaveform*> *fWFIn = nullptr;
	std::vector<double> *amplitude = nullptr;
	std::vector<double> *rval = nullptr;
	std::vector<double> *zval = nullptr;

	TChain *fIn = new TChain("skimTree");
	fIn->Add(InputFile.c_str());
	fIn->SetBranchAddress("MGTWaveforms", &fWFIn);
	fIn->SetBranchAddress("amplitude", &amplitude);
	fIn->SetBranchAddress("rval", &rval);
	fIn->SetBranchAddress("zval", &zval);

	TFile *fOut = new TFile(("MGTRigged" + InputFile).c_str(), "RECREATE");
	TTree *fTree = new TTree("MGTree", "MGTree");

	MGTEvent *fEvent = new MGTEvent();
    MJTRun *fRun = new MJTRun();
    
    double amp = 0, R = 0, Z = 0, S0 = 0, S1 = 0, S2 = 0, S3 = 0, S4 = 0, S5 = 0;

    fTree->Branch("event", &fEvent);
	fTree->Branch("run", &fRun);
	fTree->Branch("amplitude", &amp, "amplitude/D");
	fTree->Branch("rval", &R, "rval/D");
	fTree->Branch("zval", &Z, "zval/D");

    // Dummy run #
    fRun->SetRunNumber(9999);
    fRun->SetRunDescription("Simulation Run Scan");
    fRun->SetRunType( MGRunType::kData );
	fRun->SetStartTime( 0.);

    int nEntries = fIn->GetEntries();

	for(int i = 0; i < nEntries; i++)
	{
		if(i%1000 == 0) std::cout << "Processing Event: " << i << std::endl;
		
		fIn->GetEntry(i);

		fEvent->InitializeArrays("MJTGretina4DigitizerData"); // Sets name of digitizer class -- this may change...?
		// Fill Waveform
		MGTWaveform *wf = static_cast<MGTWaveform*>(fEvent->GetWaveforms()->ConstructedAt(0));
		wf->SetWFType(MGWaveform::kADC);
		wf->SetData(fWFIn->at(0)->GetVectorData());
		wf->SetSamplingFrequency(0.1); // This is the correct frequency! (Does it change with multi-sampled WFs?)
		
		// Set dummy digitizer values -- I don't know what the values mean!
		MGTVDigitizerData* dd = static_cast<MGTVDigitizerData*>(fEvent->GetDigitizerData()->ConstructedAt(0));
  		MJTGretina4DigitizerData* gdd = static_cast<MJTGretina4DigitizerData*>(dd);
  		gdd->SetEnergy(amplitude->at(0));
  		gdd->SetTimeStamp((double) i);
  		gdd->SetID(0, 0, 0);
  		gdd->SetIndex(i);
  		gdd->SetBoardSerialNumber(0.);
  		gdd->SetFifoIsHalfFull(0);
		gdd->SetFifoIsAlmostFull(0);

		fEvent->SetETotal(amplitude->at(0));
		fEvent->SetTime(0.);
		
		// Fill tree
		amp = amplitude->at(0);
		R = rval->at(0);
		Z = zval->at(0);
		fTree->Fill();

		// Must clear here
		fEvent->ClearEventData(true);
	}

	fTree->Print();

	fOut->cd();
	fTree->Write();
	fOut->Close();
}
