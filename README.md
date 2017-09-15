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
