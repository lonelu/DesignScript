#!/bin/bash
#
#$ -S /bin/bash
#$ -o /wynton/home/degradolab/lonelu/DesignData/Chemodrugs/HB_RUC_mpnn/huggingface_sel/
#$ -e /wynton/home/degradolab/lonelu/DesignData/Chemodrugs/HB_RUC_mpnn/huggingface_sel/
#$ -cwd
#$ -j y
#$ -l mem_free=24G
#$ -l gpu_mem=24000
#$ -l scratch=24G
#$ -l h_rt=2:00:00
#$ -l compute_cap=61
#$ -q gpu.q
#$ -m bea
#$ -M lei.lu@ucsf.edu

date
hostname

module load Sali cuda/11.5.0
conda activate env_omegafold

python /wynton/home/degradolab/lonelu/GitHub_Design/OmegaFold/main.py /wynton/home/degradolab/lonelu/DesignData/Chemodrugs/HB_RUC_mpnn/huggingface_sel/huggingface_sel.fasta $TMPDIR/omegafold_output

cp -r $TMPDIR/omegafold_output /wynton/home/degradolab/lonelu/DesignData/Chemodrugs/HB_RUC_mpnn/huggingface_sel/

### qsub /wynton/home/degradolab/lonelu/GitHub_Design/OmegaFold/submit_main.sh
### qsub -q gpu.q /wynton/home/degradolab/lonelu/DesignData/Chemodrugs/submit_main_4.sh