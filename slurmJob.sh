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

# This runs whatever commands lat-jobs.py passes to it.
echo "${@}"
${@}

echo "Job Complete:"
date
