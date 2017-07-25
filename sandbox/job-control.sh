#!/bin/bash
# Do some useful command line things.
# C. Wiseman, USC.
# DS map (from DataSetInfo.hh):
# map<int,int> dsMap = {{0,76},{1,51},{2,7},{3,24},{4,22},{5,112}};

function submitJobs()
{
  # run a test job
  # qsub -l h_vmem=2G qsub-job.sh 5 97
  # qsub -l h_vmem=2G qsub-job.sh 1 9503

  ##############################################################################
  # Main block.
  #
  # DID YOU CHECK qsub-job.sh???

  # Phase 1: submit skim_mjd_data jobs.
  # Put this in qsub-job.sh:
  # ./skim_mjd_data $1 $2 -t 0.8 -l -n ~/project/skim-v2/

  # Phase 2: submit wave-skim jobs.
  # Put this in qsub-job.sh:
  # ./wave-skim -r $1 $2

  # Phase 3: submit lat jobs.
  # Put this in qsub-job.sh:
  # ./lat.py -b -r $1 $2

  # for ((i=0; i<=76; i++)); do qsub -l h_vmem=2G qsub-job.sh 0 $i; done
  # for ((i=0; i<=51; i++)); do qsub -l h_vmem=2G qsub-job.sh 1 $i; done
  # for ((i=0; i<=7; i++)); do qsub -l h_vmem=2G qsub-job.sh 2 $i; done
  # for ((i=0; i<=24; i++)); do qsub -l h_vmem=2G qsub-job.sh 3 $i; done
  # for ((i=0; i<=22; i++)); do qsub -l h_vmem=2G qsub-job.sh 4 $i; done
  # for ((i=0; i<=112; i++)); do qsub -l h_vmem=2G qsub-job.sh 5 $i; done

  # Phase 4: Calibration data
  # ./lat.py -f 0 $i -b
  for((i=2931; i<=2940; i++)); do qsub -l h_vmem=2G qsub-job.sh 0 $i; done
  for((i=6854; i<=6863; i++)); do qsub -l h_vmem=2G qsub-job.sh 0 $i; done
  for((i=9497; i<=9503; i++)); do qsub -l h_vmem=2G qsub-job.sh 1 $i; done
  for((i=14149; i<=14155; i++)); do qsub -l h_vmem=2G qsub-job.sh 1 $i; done
  for((i=14568; i<=14574; i++)); do qsub -l h_vmem=2G qsub-job.sh 2 $i; done
  for((i=15789; i<=15794; i++)); do qsub -l h_vmem=2G qsub-job.sh 2 $i; done
  for((i=16911; i<=16920; i++)); do qsub -l h_vmem=2G qsub-job.sh 3 $i; done
  for((i=17950; i<=17959; i++)); do qsub -l h_vmem=2G qsub-job.sh 3 $i; done
  for((i=60001014; i<=60001023; i++)); do qsub -l h_vmem=2G qsub-job.sh 4 $i; done
  for((i=60001855; i<=60001864; i++)); do qsub -l h_vmem=2G qsub-job.sh 4 $i; done
  for((i=19055; i<=19064; i++)); do qsub -l h_vmem=2G qsub-job.sh 5 $i; done
  for((i=22513; i<=22523; i++)); do qsub -l h_vmem=2G qsub-job.sh 5 $i; done

  ##############################################################################
}

function mergeOutput()
{
  dsNum=2
  dsMax=7
  outFile="thresholdsDS${dsNum}.root"
  args=()
  for ((i=0; i<=${dsMax}; i++)); do
      args+=("thresholdsDS${dsNum}_$i.root")
  done
  echo "${args[@]}"
  thisDir=`pwd`
  cd data
  hadd $outFile ${args[@]}
  cd $thisDir
}

function findExposure()
{
  arr=(0 1 3 4 5)
  for i in "${arr[@]}"; do
    echo $i
    make -s && ./thresholds -e final/thresholdsDS$i.root $i
  done
}

# =================================
# MkCookie
submitJobs
# mergeOutput
# findExposure