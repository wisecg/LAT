#!/bin/bash

# Micah's slurm/shifter script

#--- activate on Cori
#SBATCH --nodes=4 --ntasks=4
#SBATCH --partition=debug --time=00:30:00
#SBATCH --account=majorana
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mbuuck@uw.edu
#SBATCH --image=docker:buuck/mjsw:gat_psa

#--- activate on Cori KNL nodes
#-SBATCH --output=./logs/pulse1-knl.o%j
#-SBATCH --constraint=knl,quad,cache --job-name=pulse1-knl
#-SBATCH --cpus-per-task=272

#--- activate on Cori Haswell nodes
#SBATCH --output=./logs/pulse1-haswell.o%j
#SBATCH --constraint=haswell --job-name=pulse1-haswell
#SBATCH --cpus-per-task=64

#--- activate on PDSF-3
#-SBATCH --job-name=pbg-pdsf   --image=custom:pdsf-chos-sl64:v4

#tasks to be executed
job_sh=./mjdInShifter.sh

echo start-A
env|grep  SHIFTER_RUNTIME
ls -l  ${job_sh}
export OMP_NUM_THREADS=64
export OMP_PROC_BIND=true
export OMP_PLACES=threads

echo "enter cgroup ----"
srun -N 4 -n 4 -c 64 shifter  --volume=/global/project:/project /bin/bash ${job_sh}
echo "returned from cgroup ----"
echo end-A

# Jan's conditional stuff -- conditions are only allowed AFTER the SBATCH commands (sadly).
#!/bin/bash
#SBATCH -t 25:00   --ntasks=1 --account majorana

# Pick one of the following lines to toggle: chos or shifter or Cori
# (toggle  '#-SBATCH' vs. '#SBATCH'  )
#-SBATCH -J mjr-chos -p shared-chos
#-SBATCH -J mjr-shift -p shared --image=custom:pdsf-chos-sl64:v4
#SBATCH -J mjr-cori -p debug -N1 --image=custom:pdsf-chos-sl64:v2  -C haswell

#tasks to be executed
job_sh=mjrTask.sh

export RUN2=${1-16846}

echo "start-A "`hostname`"  RUN2="$RUN2
ls -l  ${job_sh}
if [[ $SLURM_JOB_PARTITION == *"-chos" ]]
then
  echo  run-in-chos
  CHOS=sl64 chos  ./${job_sh}
else
 echo  run-in-shifter
 shifter  --volume=/global/project:/project  /bin/bash ${job_sh}
fi
echo end-A


