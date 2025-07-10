#!/bin/bash
# sbatch.sh

#SBATCH --job-name="O2Sim_slurm"
#SBATCH --output="O2Sim_slurm.%j.out"
#SBATCH --error="O2Sim_slurm.%j.err"
#SBATCH --partition=main
#SBATCH --time=9:00:00
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=10000


export lowPhotGenHPath=<path to the low photon generator .h file>
export lowPhotGenMPath=<path to the low photon generator .macro file>
export o2EnvPath=<path to the file containing the o2 environment>
export o2simConfig=<path to o2sim_configuration.ini>
export magFieldPath=../field.C
export fctConfig=../FCT_SD.cfg

mkdir job$SLURM_ARRAY_TASK_ID
cd job$SLURM_ARRAY_TASK_ID
echo $SLURM_CPUS_PER_TASK
source ../simulation.sh $SLURM_ARRAY_TASK_ID $SLURM_CPUS_PER_TASK

mv ../O2Sim_slurm.${SLURM_JOB_ID}.out .
mv ../O2Sim_slurm.${SLURM_JOB_ID}.err .
