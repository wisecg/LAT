#!/bin/bash
# Examine this to see how default sl64 environments are set:
# /global/project/projectdirs/majorana/setupMajorana.sh
# Source this file in .bashrc.ext
# Clint Wiseman, USC/Majorana

function sourceEnvMJD {
  umask "u=rwx,g=rx,o=r"

  export PATH=${HOMEDIR}:${PATH} # rmate
	export SWDIR=/global/project/projectdirs/majorana/software/sl64
  export MJSWDIR=${SWDIR}/mjsw/mjsw201712Prod
  export MJHOME=/global/project/projectdirs/majorana
  export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/usr/local/lib
  export MJDDATADIR=/global/project/projectdirs/majorana/data/mjd
  export MJSWDEV=$SWDIR/mjsw/mjsw201706Prod

  source ${SWDIR}/root/root-6.12.04/installp3/bin/thisroot.sh
  export CLHEP_BASE_DIR=${SWDIR}/CLHEP/2.4.0.1/CLHEP
	export CLHEP_INCLUDE_DIR=${CLHEP_BASE_DIR}/include
	export CLHEP_LIB_DIR=${CLHEP_BASE_DIR}/lib
	export CLHEP_LIB=CLHEP
	export PATH=$CLHEP_BASE_DIR/bin:${PATH}
	export LD_LIBRARY_PATH=${CLHEP_LIB_DIR}:${LD_LIBRARY_PATH}

  # export MGDODIR=${MJSWDIR}/MGDO
  # export GATDIR=${MJSWDIR}/GAT
  export MGDODIR=${HOMEDIR}/mgsw/MGDO
  export GATDIR=${HOMEDIR}/mgsw/GAT
  # export MGDODIR=${HOMEDIR}/mgsw/test/MGDO
  # export GATDIR=${HOMEDIR}/mgsw/test/GAT

  export TAMDIR=${MGDODIR}/tam
  export ORDIR=${MJSWDIR}/OrcaRoot
  export MJORDIR=${MJSWDIR}/MJOR
  export PATH=${MJSWDIR}/bin:${ORDIR}/Applications:${MJORDIR}:${GATDIR}/Apps:${GATDIR}/Scripts:${SIGGENDIR}:${MGDODIR}/bin:${PATH}
  export LD_LIBRARY_PATH=$MJSWDIR/lib:${ORDIR}/lib:${MGDODIR}/install/lib:${GATDIR}/lib:${LD_LIBRARY_PATH}:${MAGEDIR}/analysis
  export ROOT_INCLUDE_PATH=${CLHEP_INCLUDE_DIR}:${MGDODIR}/Base:${MGDODIR}/Gerda:${MGDODIR}/GerdaTransforms:${MGDODIR}/Majorana:${MGDODIR}/MJDB:${MGDODIR}/Root:${MGDODIR}/Tabree:${MGDODIR}/Tools:${MGDODIR}/Transforms:${TAMDIR}:${TAMDIR}/inc:${GATDIR}/BaseClasses:${GATDIR}/MGTEventProcessing:${GATDIR}/MGOutputMCRunProcessing:${GATDIR}/SiggenWrapper

  export LATDIR=${HOMEDIR}/lat
  export LATDATADIR=${HOMEDIR}/project
  export MAGERESULTS=/global/projecta/projectdirs/majorana/sim/MJDG41003Sims
  export PYTHONPATH=${LATDIR}:${PYTHONPATH}
}

function sourceEnvClint {
  export PS1="\[\033[38;5;14m\]\A\[$(tput sgr0)\]\[\033[38;5;15m\] \[$(tput sgr0)\]\[\033[38;5;10m\]\u@\h\[$(tput sgr0)\]\[\033[38;5;15m\]:\[$(tput bold)\]\[$(tput sgr0)\]\[\033[38;5;9m\]\W\[$(tput sgr0)\]\[$(tput sgr0)\]\[\033[38;5;15m\]\\$ \[$(tput sgr0)\]"
  alias rootmj="root -b -l"
  alias atom='rmate -p "$RMATE_PORT"'
  alias ls="ls --color"
  alias inter="salloc -t 10:00:00 -p shared"
  alias imj="shifter --image wisecg/mjsw:v2 bash"
}

function rpf {
  myarr=$(ps -u `whoami` | grep sshd | awk '{print $1}')
  kill $myarr
}

function loadModulesMJD {
  source /etc/profile.d/modules.sh
  export CC=/usr/common/usg/software/gcc/4.8.2/bin/gcc
  export CXX=/usr/common/usg/software/gcc/4.8.2/bin/g++
  module load python/3.4.3
  # module load python/3.6-anaconda-4.4
  module load gcc/4.8.2
  module load ImageMagick
  module load curl
  module load git
  module load cmake/3.5.2
  module load latex
  module load valgrind
  module load oprofile

  findrun () { find $MJDDATADIR -name "*$1*";}
  alias mjdiskspace="prjquota majorana | grep majorana | awk '{print \"prjquota: \" \$3-\$2\" GB\"}'; myquota -j /global/project/projectdirs/majorana | grep -e \"^/global\" | awk '{print \"myquota (/project): \" \$3-\$2 \" TB\"}'; myquota -j /global/projecta/projectdirs/majorana | grep -e \"^/global\" | awk '{print \"myquota (/projecta): \" \$3-\$2 \" TB\"}'"
  gatTagFromDecRev () { echo "obase=16; $1" | bc | xargs -I{} echo 0x{} | xargs printf %07x | xargs git --git-dir=$MJSWDEV/GAT/.git describe; }
  alias countInodesHere="find . -printf \"%h\n\" | cut -d/ -f-2 | sort | uniq -c | sort -rn"
}


# we're on pdsf (interactive node)
if [ "$NERSC_HOST" = "pdsf" ] && [[ $- == *i* ]] && [ -z "$SLURM_JOB_ID" ]; then
  if [ "${CHOS}" != "sl64" ]; then
    export CHOS=sl64
    chos
  fi
  sourceEnvClint
  loadModulesMJD
  sourceEnvMJD
fi

# we're on pdsf (slurm job)
if [ "$NERSC_HOST" = "pdsf" ] && [ -n "$SLURM_JOB_ID" ]; then
  sourceEnvClint
  sourceEnvMJD
fi

# do for cron jobs
if [ "$RUNCRON" = 1 ]; then
  source /etc/profile.d/modules.sh
  export CC=/usr/common/usg/software/gcc/4.8.2/bin/gcc
  export CXX=/usr/common/usg/software/gcc/4.8.2/bin/g++
  module load python/3.4.3
  module load gcc/4.8.2
  export CHOS=sl64
  chos
  sourceEnvMJD
fi

# we're on cori (interactive node)
if [ "$NERSC_HOST" = "cori" ] && [[ $- == *i* ]] && [ -z "$SLURM_JOB_ID" ]; then
  module load python/3.6-anaconda-4.4
  sourceEnvClint
  sourceEnvMJD
fi

# we're on cori (slurm job)
if [ "$NERSC_HOST" = "cori" ] && [ -n "$SLURM_JOBID" ]; then
  # module load python/3.5-anaconda
  module load python/3.6-anaconda-4.4
  sourceEnvClint
  sourceEnvMJD
fi

# we're on edison (interactive node)
if [ "$NERSC_HOST" = "edison" ] && [[ $- == *i* ]] && [ -z "$SLURM_JOB_ID" ]; then
  module load python/3.6-anaconda-4.4
  sourceEnvClint
  sourceEnvMJD
fi

# we're on edison (slurm job)
if [ "$NERSC_HOST" = "edison" ] && [ -n "$SLURM_JOBID" ]; then
  # module load python/3.5-anaconda
  module load python/3.6-anaconda-4.4
  sourceEnvClint
  sourceEnvMJD
fi

# do for only shifter jobs
# if [ -n "$SHIFTER_RUNTIME" ]; then
# fi

# do for pdsf slurm jobs
# if [ -n "$SLURM_JOB_ID" ]; then
# fi

# do for cori/edison slurm jobs
# if [ -n "$SLURM_JOB_ID" ]; then
# fi

# do for all interactive nodes
# if [[ $- == *i* ]]; then
# fi
