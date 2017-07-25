#!/bin/bash
#$ -cwd
#$ -j y
#$ -o /global/homes/w/wisecg/skim-clean/logs/
#$ -P majorana
source /global/homes/w/wisecg/env/EnvBatch.sh # can also comment this out and run with qsub -V
cd /global/homes/w/wisecg/skim-clean

echo "Job Start:"
date
echo "Node:  "$HOSTNAME
echo "Job ID:  "$JOB_ID

# This runs whatever commands job-panda.py passes to it.
echo "${@}"
${@}

echo "Job Complete:"
date
