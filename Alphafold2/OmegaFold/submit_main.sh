#!/bin/bash
#
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -l mem_free=2G
#$ -l gpu_mem=5000
#$ -l scratch=2G
#$ -l h_rt=2:00:00

date
hostname

module load Sali cuda/11.5.0
conda activate env_omegafold

python main.py $1 $TMPDIR/omegafold_output

cp -r $TMPDIR/omegafold_output $2
