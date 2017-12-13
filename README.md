# LAT (Low-Energy Analysis Toolkit)

Written to perform a complete analysis of Majorana Demonstrator DS0-5 low-energy data.

The structure is mainly:
- job-panda (do file management, submit jobs, lots of misc things.)
- skim_mjd_data (produce low energy skim files w/ special options to reduce threshold)
- wave-skim (grab all waveforms for hits passing basic data cleaning cuts)
- lat, lat2, lat3 (perform secondary waveform processing)
- ds_livetime (calculate livetime of final analysis)
- specFit (unbinned likelihood fits to final spectra)

The bookeeping and many miscellaneous tasks are covered by 'job-panda.py', and most of the miscellaneous functions are stored in 'waveLibs.py'.  

Some general properties of the data are contained in 'DataSetInfo.py', while many more are stored in the calibration database file 'calDB.json'.  See waveLibs and LAT2 for examples of accessing and updating the database.

There are a lot of half-baked ideas, some useful, in the folder 'sandbox'.

## Authors

LAT was written by typing on keyboards.  Mainly by Clint Wiseman and Brian Zhu, with lots of input from the Majorana Collaboration.

## Data Production Procedures

Primary Data Production:
The primary data production workflow follows 6 steps, all batch job submission can and should be handled by `job-panda.py` on PDSF!
Pre-job submission checks
  -- If on PDSF, create SLURM submission script (eg: `./job-panda.py -makeSlurm`)
  -- Edit all input and output directories, make sure they're correct! (Primarily check directories in `slurm-job.sh` and `job-panda.py`)

1) Re-skim data with lower energy threshold (eg: `./job-panda.py -skim -ds [dsNum]`)
2) Grab waveforms for all skim data from built files and apply basic cut using wave-skim (eg: `./job-panda.py -wave -ds [dsNum]`)
3) Split skim files with waveforms into smaller chunks for faster processing (eg: `./job-panda.py -qsubSplit -ds [dsNum]`)
4) (Optional but recommended) Write basic cut to all split files (eg: `./job-panda.py -writeCut -ds [dsNum]`)
5) Run LAT for secondary waveform processing. LAT will do waveform fitting as well as wavelet packet transform so it takes a while! (eg: `./job-panda.py -lat -ds [dsNum]`)
6) Check all log files for errors! (eg: `./job-panda.py -checkLogs` and `./job-panda.py -checkLogs2`)

Cut Tuning:
After the secondary waveform processing is completed, `lat3.py` can be used to automatically tune a cut parameter. Cuts must be tuned on calibration data! Again, all workflow can and should be handled by `job-panda.py` on PDSF!
0) Create a directory ./LAT/plots/tuneCuts to store diagnostic plots
1) Tune cut parameter of choice, check the `tuneCuts` function in `job-panda.py` for all options (eg: `./job-panda.py -ds [dsNum] -tuneCuts '[argString]'`)
2) Check all diagnostic plots and logs to make sure cut values were set properly!
3) Apply cut parameters on data (eg: `./job-panda.py -lat3 [dsNum] [cutType]` or `./lat3.py -cut [dsNum] [cutType] [dataType]`)


## LAT Database Notes

The database `calDB.json` is a tinyDB whose items are nested dicts: `{"key":key, "vals":vals}`

### Cal Tables
Divides dataset into N indexes.
These tables give the run coverage for every `calIdx`.

    key: ds[DS]_calIdx
    vals: {[idx]:[cal lo, cal hi, cov lo, cov hi]}
Print one with `waveLibs.getDBCalTable(dsNum, verbose=True)`

### Cal Records
Calibration constants for each channel in each calIdx.
(Channel list comes from DataSetInfo.py.)

    key: ds[DS]_idx[n]
    vals: {[chan]:[trapENF, fitAmp, latAF, latAFC]}


### Cut Records

    key: [Name]_ds[i]_idx[j]_module[k]_[descriptor].
    vals: {[chan]:[1%, 5%, 90%, 95%, 99%]}

- Names are: "riseNoise", "fitSlo", "bcMax", "pol2", and "pol3"
- idx is the calIdx
- Descriptors: `Peak, Continuum, 50_90, 90_130, 130_170, 170_210`
- Continuum range is 5-50 keV
- Peak range is 236-240 keV
