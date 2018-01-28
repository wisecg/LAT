#!/bin/bash
# A port of "jobPump" by Jan Balewski, NERSC.  https://bitbucket.org/balewski/vm-eval-tool
# Designed to run many single-core jobs on a multi-core batch node.
# I only rewrote it so I could figure out what it was doing!
#
# USAGE: ./job-pump.sh [jobs.list] [program name, as seen by 'top']
# NOTE: lines in jobs.list can be ignored with the # character.
#
# Clint Wiseman, USC, Jan. 2018.

set -u ;         # exit if you try to use an uninitialized variable
set -e ;         # exit if any statement returns a non-true return value
set -o errexit ; # exit if any statement returns a non-true return value

# Pump Settings (takes 4 inputs from command line)
jobList=$1          # list of commands to execute
execName=$2         # program name as it appears in 'top'
maxJobLoad=$3       # equal to the nCores requested (#SBATCH -nXX)
peakWLoad=$4        # should be = (1.1 * nCores)
lowFreeRAM=5.0      # minimum free RAM allowed on the node
totJobLimit=9999999 # maximum number of jobs to attempt
jobDelay=10         # seconds between job launches (default: 2)
checkDelay=60       # seconds between checks (default: 20)


function countJobs {
  nRunA=`/usr/bin/pgrep $1 |wc -l `
  nRun="$(echo -e "${nRunA}" | tr -d '[:space:]')"
  # echo "Found $nRun $1 jobs running (limit: $maxJobLoad)"
}

function getFreeRAM {

  freeRAM=0
  if [ "`uname`" == "Darwin" ]; then
    freeRAM=`top -l 1 -s 0 | grep PhysMem | awk '{print $6/1000}'`

  elif [ "`uname`" = "Linux" ]; then
    freeRAM=`free -g | grep "cache:" |  awk  '{ print $4}'`

  fi
  # echo "Current free RAM (GB): $freeRAM (limit: $lowFreeRAM)"
}

function getLoad {
  # echo getLoad
  ln=`uptime | awk -F'[a-z]:' '{print $2}'`
  # echo ln=$ln
  w=`echo $ln |  cut -f1 -d\ `
  # echo w=$w
  wLoad=`echo $w | cut -f1 -d.`
  # echo "Current workload: $wLoad (limit: $peakWLoad)"
}

function printNodeHealth {
  echo "Checking health of node "`hostname -f`" on "`date`

  if [ "`uname`" == "Darwin" ]; then
    top -l 1 -s 0 | grep PhysMem
    top -l 1 | head -n 10
    top -l 1 | grep $execName | grep $USER | nl

  elif [ "`uname`" = "Linux" ]; then
    echo "free -g"
    free -g
    top bn1 | head -n 12
    top bn1 | grep $execName | grep $USER | nl
  fi
}

# ==========================================================================================

echo -e "Starting job pump ...
   Task List:           $jobList
   Task Name:           $execName
   Max Running Jobs:    $maxJobLoad
   Total Jobs Limit:    $totJobLimit
   Peak Load Limit:     $peakWLoad
   Free RAM Limit (GB): $lowFreeRAM
   New Job Delay (sec): $jobDelay
   Monitor Delay (sec): $checkDelay
   Running in:         "`pwd`"\n"

if [ ! -f $jobList ]; then
  echo "Task list '$jobList' not found!  Exiting..."
  exit
fi

# Begin loop over the job list (given as argument)
k=0
nQueue=0
myDelay=$jobDelay
while read line ; do

  if [[ $line == "#"* ]]; then
    echo "Skipping line:" $line
    continue
  fi
  if [ $nQueue -ge $totJobLimit ]; then
    echo "WARNING: total job limit ($nQueue) reached, skipping all remaining jobs in job list ..."
    break
  fi

  # Try to *not* submit new job - machine may be busy
  while true ; do
    sleep $myDelay
    countJobs $execName
    getLoad
    getFreeRAM

    echo "Check #$k. nQueue: $nQueue, just slept $myDelay seconds.  Date: "`date`
    k=$[ $k +1 ]

    # Every 10 checks, print the node health
    if [ $[ $k%10 ] -eq 5 ]; then
      echo " "
      printNodeHealth $k ;
      echo " "
    fi

    if [ $nRun -ge $maxJobLoad ]; then
      myDelay=$checkDelay
      echo -e "Already running max job load: $nRun ... sleeping $myDelay seconds.\n"
      continue
    fi
    if [ $wLoad  -ge $peakWLoad  ]; then
      myDelay=$checkDelay
      echo -e "Workload is too high: $wLoad ... sleeping $myDelay seconds.\n"
      continue
    fi
    if (( $(echo "$lowFreeRAM > $freeRAM" | bc -l) )); then
      myDelay=$checkDelay
      echo -e "Free RAM is too low: $freeRAM ... sleeping $myDelay seconds.\n"
      continue
  	fi
    break # no reason to slack any more, go ahead and start a new job.
  done

  # Submit a new job ('&' puts job into background)
  job=$line' &'
  echo "Starting job: "$job
  eval $job
  nQueue=$[ $nQueue + 1 ]

  # Compute a larger delay if nQueue is close to max workload or freeRAM
  myDelay=$checkDelay
  delay=$[  $peakWLoad - $wLoad ]
  if (( $(echo "$freeRAM < $delay" | bc -l) )); then
    delay=$freeRAM
  fi
  if (( $(echo "$delay > 0" | bc -l) )); then
    myDelay=$(echo $jobDelay + 60/$delay | bc)
  fi

  echo "Workload: $wLoad (lim $peakWLoad)  Free RAM: $freeRAM (lim $lowFreeRAM)  nRunning: $nRun (lim $maxJobLoad)  delay: $delay"
  echo "Now sleeping for $myDelay seconds ..."
  echo " "

done < $jobList

echo "Task list completed. $nQueue jobs submitted, now waiting for completion ..."
sleep 60  # for good measure
printNodeHealth $k

while [ true ] ; do
  countJobs $execName
  getLoad
  getFreeRAM
  k=$[ $k +1 ]
  if [ $[ $k%10 ] -eq 0 ]; then printNodeHealth $k ;  fi
  if [ $nRun -le 0 ]; then break; fi
  echo "Check #$k:  Finishing jobs ..."`date`"  sleeping $checkDelay seconds ..."
  sleep $checkDelay
  continue
done

echo "Finished job pump, date: "`date`
