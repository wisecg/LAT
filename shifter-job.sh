#!/bin/bash
echo "Job Start:"
date
echo "Node(s): "$SLURM_JOB_NODELIST
echo "Job ID:  "$SLURM_JOB_ID
if [ -n "$SHIFTER_RUNTIME" ]; then
  echo "Shifter active"
  echo $CC
  echo $CXX
  which python
  which root
fi

# This runs whatever commands we pass to it.
# echo "${@}"
# ${@}

# echo "attempting to run lat3:"
# ./lat3.py -cut 2 fs

echo "Job Complete:"
date
