#!/bin/bash -l
#SBATCH -t 24:00:00  --ntasks=1
#SBATCH --mem 3400
#SBATCH --account=majorana
#SBATCH --workdir=/global/homes/w/wisecg/lat
#SBATCH --output=/global/homes/w/wisecg/lat/logs/slurm-%j.txt

echo "Job Start:"
date
echo "Node(s):  "$SLURM_JOB_NODELIST
echo "Job ID:  "$SLURM_JOB_ID

# This runs whatever commands job-panda.py passes to it.
echo "${@}"
${@}

echo "Job Complete:"
date