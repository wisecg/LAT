# LAT (Low-Energy Analysis Toolkit)

Written to perform a complete analysis of Majorana Demonstrator DS0-5 low-energy data.

## Main Routines:
- `skim_mjd_data`: produce low energy skim files w/ special options to reduce threshold
- `wave-skim`: grab all waveforms for hits passing basic data cleaning cuts
- `ds_livetime`: calculate livetime w/o low-e run and channel selection
- `lat-jobs`: control job submission
- `lat`: waveform fitting
- `lat2`: tune and apply PSA cuts
- `lat3`: tune and apply burst cut
- `lat-expo`: exposure and efficiency
- `spec-fit`: final spectrum fits

Properties of the data sets, including background, calibration, special runs, and detector info is stored in `dsi.py`
Run lists come from `./data/runs*.json`
Convenience functions are in `waveLibs.py`
There are a lot of half-baked ideas, some useful, in `LAT/sandbox`.

### Authors
Clint Wiseman and Brian Zhu, with input from the Majorana Collaboration.