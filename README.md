# LAT (Low-Energy Analysis Toolkit)

Written to perform a complete analysis of Majorana Demonstrator DS0-5 low-energy data.

## Main Routines:
- `skim_mjd_data`: produce low energy skim files w/ special options to reduce threshold
- `wave-skim`: grab all waveforms for hits passing basic data cleaning cuts
- `ds_livetime`: calculate livetime w/o low-e run and channel selection
- `lat-jobs`: control job submission
- `lat-settings`: save channel map and threshold/HV settings for all runs
- `lat`: waveform fitting and wavelet packet decomposition (PSA parameters)
- `lat2`: tune and apply PSA cuts
- `lat3`: tune and apply burst cut
- `lat-check`: data integrity checks for all files
- `lat-expo`: exposure and efficiency
- `spec-fit`: final spectrum fits
- `dsi.py`: properties of the data sets, including background, calibration, special runs, and detector info
- `./data/runs*.json`: run lists for bkg runs (match `DataSetInfo.cc`), calibration, and special runs
- `waveLibs.py`: convenience functions

There are a lot of half-baked ideas, some useful, in `LAT/sandbox`.

### Authors
Clint Wiseman and Brian Zhu, with input from the Majorana Collaboration.