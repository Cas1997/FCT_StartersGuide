#!/usr/bin/env bash
now=$(date +%s)
seed=$((now + $1 - 1000000000)) # Get seed from unix time and add job number to it.

# Generate the .cfg file for pythia because the seed needs to be set independently
touch pythia8_conf.cfg
echo \#\#\# beams >> pythia8_conf.cfg
echo Beams:idA 2212			\# proton >> pythia8_conf.cfg
echo Beams:idB 2212 			\# proton >> pythia8_conf.cfg
echo Beams:eCM 14000. 		\# GeV >> pythia8_conf.cfg
echo  >> pythia8_conf.cfg
echo \#\#\# processes >> pythia8_conf.cfg
echo SoftQCD:nonDiffractive on >> pythia8_conf.cfg
echo  >> pythia8_conf.cfg
echo \#\#\# decays >> pythia8_conf.cfg
echo ParticleDecays:limitTau0 on >> pythia8_conf.cfg
echo ParticleDecays:tau0Max 10.	>> pythia8_conf.cfg
echo  >> pythia8_conf.cfg
echo \#\#\# seed >> pythia8_conf.cfg
echo Random:setSeed on >> pythia8_conf.cfg
echo Random:seed $seed >> pythia8_conf.cfg
echo \#\#\# >> pythia8_conf.cfg

# Run the o2sim

export ALIBUILD_WORK_DIR=<path to aliBuild work dir>
eval `/usr/local/bin/alienv shell-helper`
export PATH=${PATH}
source $o2EnvPath
set -x

MODULES="FCT TRK FT3 A3IP TF3 RCH"
SIGEVENTS=25000
NWORKERS=$2

export ALICE3_SIM_FIELD=ON
export ALICE3_MAGFIELD_MACRO=$magFieldPath

o2-sim --detectorList ALICE3 -m ${MODULES} -j ${NWORKERS} -n ${SIGEVENTS} -g external \
       --configKeyValues "GeneratorExternal.fileName=$lowPhotGenMPath; GeneratorExternal.funcName=GeneratorLowCleanFaster(); FCTBase.onlyChargedParticles=false; FCTBase.configFile=$fctConfig" \
       --configFile $o2simConfig
