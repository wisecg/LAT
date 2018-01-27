#!/bin/bash

echo "Job Start:"
date
echo "Node(s):  "$SLURM_JOB_NODELIST
echo "Job ID:  "$SLURM_JOB_ID

if [ -n "$SHIFTER_RUNTIME" ]; then
  echo "Shifter image active."
  echo "pwd: "`pwd`
  echo "gcc: "$CC
  echo "g++:"$CXX
  echo "Python:"`python --version`
  echo "ROOT:"`root-config --version`
fi

# This runs whatever commands job-panda.py passes to it.
echo "${@}"
echo "--------"
${@}

# ./skim_mjd_data -f 22513 -l -t 0.7 /global/projecta/projectdirs/majorana/users/wisecg/special/skim >& ./logs/specialSkim-DS5-22513.txt


echo "Job Complete:"
date
