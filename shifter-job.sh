#!/bin/bash
echo "Job Start:"
date
echo "Node(s):  "$SLURM_JOB_NODELIST
echo "Job ID:  "$SLURM_JOB_ID
echo inShifter:`env|grep  SHIFTER_RUNTIME`
echo $CHOS
echo `pwd`
echo $SHELL
echo $MJSWDIR
echo "homedir is:"$HOMEDIR

if [ "$SHIFTER_RUNTIME" ]; then
  echo "we're in shifter"
else
  echo "we're not in shifter"
fi

# maybe this source doesn't work in the shifter-job script ?
# thisDir=`pwd`
# source ~/env/EnvBatch.sh
# echo $ROOTSYS
# cd $ROOTSYS/bin
# source thisroot.sh
# cd $thisDir

# This runs whatever commands we pass to it.
# echo "${@}"
# ${@}

echo "made it back here."

echo $ROOTSYS

echo "attempting to run lat3:"
./lat3.py -cut 2 fs

echo "Job Complete:"
date
