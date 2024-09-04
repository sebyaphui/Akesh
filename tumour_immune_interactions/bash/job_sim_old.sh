#!/bin/bash
#PBS -l select=1:ncpus=1:mem=1gb
#PBS -l walltime=01:00:00
#PBS -N simulation
 
cd $PBS_O_WORKDIR
mkdir $TMPDIR/sim
mkdir $TMPDIR/sim_data
mkdir $TMPDIR/outputs
cp sim/* $TMPDIR/sim
cp shared_sim_data/sim.pickle $TMPDIR/sim_data

module load tools/prod
module load Python/3.11.2-GCCcore-12.2.0-bare
source sim_venv/bin/activate

cd $TMPDIR
python sim/main.py -sf y -ow y -c Config37 > log.txt

mkdir $HOME/tumour_immune_interactions/job_data/$PBS_JOBID
cp log.txt $HOME/tumour_immune_interactions/job_data/$PBS_JOBID
cp outputs/* $HOME/tumour_immune_interactions/job_data/$PBS_JOBID/outputs -r
cp sim_data/* $HOME/tumour_immune_interactions/job_data/$PBS_JOBID/sim_data -r
