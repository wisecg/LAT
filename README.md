LAT (Low-Energy Analysis Toolkit) is intended to do a complete analysis of Majorana Demonstrator DS0-5 data.

The structure is mainly:
- skim_mjd_data (produce low energy skim files w/ special options to reduce threshold)
- wave-skim (grab all waveforms for hits passing basic data cleaning cuts)
- lat, lat2, lat3 (perform secondary waveform processing)
- ds_livetime (calculate livetime of final analysis)

The bookeeping and many miscellaneous tasks are covered by 'job-panda.py', and most of the miscellaneous functions are stored in 'waveLibs.py'.  

Some general properties of the data are contained in 'DataSetInfo.py', while many more are stored in the calibration database file 'calDB.json'.  See waveLibs and LAT2 for examples of accessing and updating the database.

LAT was written by Clint Wiseman and Brian Zhu, with input from the Majorana Collaboration.