# LAT (Low-Energy Analysis Toolkit)

Written to perform a complete analysis of Majorana Demonstrator DS0-6A low-energy data. Current analysis was performed with the LATv2 tag. 

## Dependencies: 
Base: ROOT, MGDO, and GAT.

Python Packages: Python (2.7 or 3.x), SciPy ecosystem (NumPy, SciPy, and Matplotlib), pyROOT, tinydb, and PyWavelets.

Optional Python packages: Pandas, Seaborn, pymc3, theano, PyTables

## Main Routines and data production workflow:
- `lat-jobs.py`: This is a script that wraps all the codes and submits all of the codes in order onto PDSF/CORI/etc. Also can run some diagnostics.
- `lat-checks.py`: This script runs basic file integrity checks after every stage in LAT up to the end of `lat.py`. This script verifies that no events are lost during the LAT production process.
- `skim_mjd_data.cc`: Standard MJD code that creates skim files, except we use the option to include additional parameters as well as lower the threshold
- `wave-skim.cc`: Matches skim data events with built data events to pull waveforms, applies basic data cleaning cuts (`"!(C==1&&isLNFill1) && !(C==2&&isLNFill2) && C!=0 && P!=0 && D!=0 && isGood && !muVeto && mH==1 && gain==0 && trapENFCal > 0.7”`). It also writes the applied cut into the ROOT file.
- `splitTree` and `writeCuts`: These are functions in `lat-jobs.py`. The function `splitTree` takes the skim files with waveforms and splits them into manageable chunks (~50 MB each). The function `writeCuts` writes the applied cut (in `wave-skim.cc`) into each split file; this is necessary as the ROOT automatically splitting does not write the cut into each split file. 
- `lat.py`: Evaluates additional PSA parameters, saves into skim files. 
- `lat-settings.py`: Scans the headers of all runs for changes in the HV settings and on-board trapezoid threshold settings. These settings determine which run ranges we need to evaluate thresholds.
- `auto-thresh.cc`: Standard MJD code that evaluates thresholds. This code is already run during auto-production per run, however, we run it again according to each subDS. This can be run in parallel with `lat.py`.
- `ds_livetime.cc`: Standard MJD code that evaluates exposure and livetime. However, we use an additional option to output granular data (run + channel + exposure of all data in the 0νββ analysis) into a ROOT file for additional analysis. This allows us to re-evaluate the exposure as necessary with any combination of run+channel selection. This can be run in parallel with `lat.py`.
- `lat2.py`: Tunes and applies PSA cuts (must be run after `lat.py`). Applies and saves threshold cut to skim files. Loads and saves all m2s238 calibration data, evaluates proper cut value for the slowness + high frequency cut for each detector. Generates new skim files with the PSA cuts applied.
- `lat3.py`: Calculates rate of each ch+subDS in order to perform outlier removal (burst cut). Makes skim files with addition of burst cut applied. Must be run after `lat2.py` and `ds_livetime.cc`.
- `lat-expo.py`: Evaluates the livetime and exposure of each ch+subDS with any combination of cuts. Evaluates efficiency and saves into files. Makes final skim files with `tOffset` cut applied. Makes movies of waveforms after all cuts are applied.
- `dsi.py`: Helper module that contains properties of the data sets, including background, calibration, special runs, and detector info
- `waveLibs.py`: Helper module that contains a variety of convenience functions (basic waveform processing, histogramming, various commonly used functions, simple filter, etc)
- `./data/runs*.json`: run lists for bkg runs (match `DataSetInfo.cc`), calibration, and special runs
- `spec-fit.py`: Final spectrum fits using RooFit

There are a lot of half-baked ideas, some useful, in `LAT/sandbox`.

### Authors
Clint Wiseman and Brian Zhu, with input from the Majorana Collaboration.

700 commits strong!
