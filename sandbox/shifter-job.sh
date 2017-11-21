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

echo "sourcing setup"
# source /global/homes/w/wisecg/env/EnvSlurm.sh
# source /global/homes/w/wisecg/env/EnvBatch.sh
# source /project/projectdirs/majorana/setupMajorana.sh
echo "rootsys:"$ROOTSYS
#
# # This runs whatever commands we pass to it.
# echo "${@}"
# ${@}
#
# echo "Job Complete:"
# date
