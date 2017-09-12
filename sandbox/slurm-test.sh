#!/bin/bash -l
#SBATCH -t 2:00:00  --ntasks=1
#SBATCH --mem 3400
#SBATCH --account=majorana
#SBATCH --workdir=/global/homes/w/wisecg/lat
#SBATCH --output=/global/homes/w/wisecg/lat/logs/latjob.o%j

#======START=====
echo "The current job ID is $SLURM_JOB_ID"
echo "Running on $SLURM_JOB_NUM_NODES nodes"
echo "Using $SLURM_NTASKS_PER_NODE tasks per node"
echo "A total of $SLURM_NTASKS tasks is used"
echo "Node list:"
sacct --format=JobID,NodeList%100 -j $SLURM_JOB_ID
# aprun -B ./a.out # same as aprun -n 160 ./a.out

echo $MJDDATADIR
echo $GATDIR
echo $ROOTSYS
which python
#=====END====