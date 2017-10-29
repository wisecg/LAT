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

### Authors

LAT was written by typing on keyboards.  Mainly by Clint Wiseman and Brian Zhu, with lots of input from the Majorana Collaboration.

### LAT Database Notes

The database `calDB.json` is a tinyDB whose items are nested dicts: `{"key":key, "vals":vals}`

Cal Tables (gives run coverages)

    key: ds[DS]_calIdx.
    vals: {[idx]:[cal lo, cal hi, cov lo, cov hi]}

- Print one with `waveLibs.getDBCalTable(dsNum, verbose=True)`

Cal Records (calib consts for each channel in each calIdx)
- Channel list comes from DataSetInfo.py

    key: ds[DS]_idx[n]
    vals: {[chan]:[trapENF, fitAmp, latAF, latAFC]}


Cut records

    key: [Name]_ds[i]_idx[j]_module[k]_[descriptor].
    vals: {[chan]:[1%, 5%, 90%, 95%, 99%]}

- Names are: "riseNoise", "fitSlo", "bcMax", "pol2", and "pol3"
- idx is the calIdx
- Descriptor is two numbers (energy range), "continuum", or "peak"
- descriptors: Peak, Continuum, 50_90, 90_130, 130_170, 170_210
- Continuum range is 5-50 kev
- Peak is 236-240

