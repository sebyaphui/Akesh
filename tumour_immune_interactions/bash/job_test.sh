#!/bin/bash
#PBS -l select=1:ncpus=1:mem=1gb
#PBS -l walltime=00:01:00
#PBS -N hello_world
 
cd $PBS_O_WORKDIR
 
module load tools/prod
module load Python/3.11.2-GCCcore-12.2.0-bare
source sim_venv/bin/activate

python scripts/testing/hello_world.py > log.txt

mkdir $HOME/tumour_immune_interactions/job_data
cp * $HOME/tumour_immune_interactions/job_data