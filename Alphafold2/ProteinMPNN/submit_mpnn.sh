#!/bin/bash
#
#$ -S /bin/bash
#$ -o /wynton/home/degradolab/lonelu/DesignData/Chemodrugs/HB_RUC2_MPNN/
#$ -e /wynton/home/degradolab/lonelu/DesignData/Chemodrugs/HB_RUC2_MPNN/
#$ -cwd
#$ -j y
#$ -q gpu.q
#$ -l mem_free=64G
#$ -l gpu_mem=64000
#$ -l scratch=32G
#$ -l h_rt=1:30:00
#$ -m bea
#$ -M lonelur@gmail.com

date
hostname

module load Sali cuda/11.5.0
conda activate SE3-nvidia


folder_with_pdbs="./inputs/"

output_dir="./outputs"
if [ ! -d $output_dir ]
then
    mkdir -p $output_dir
fi


path_for_parsed_chains=$output_dir"/parsed_pdbs.jsonl"
path_for_assigned_chains=$output_dir"/assigned_pdbs.jsonl"
path_for_fixed_positions=$output_dir"/fixed_pdbs.jsonl"
chains_to_design="A"
fixed_positions="47 127"


python ../helper_scripts/parse_multiple_chains.py --input_path=$folder_with_pdbs --output_path=$path_for_parsed_chains

python ../helper_scripts/assign_fixed_chains.py --input_path=$path_for_parsed_chains --output_path=$path_for_assigned_chains --chain_list "$chains_to_design"

python ../helper_scripts/make_fixed_positions_dict.py --input_path=$path_for_parsed_chains --output_path=$path_for_fixed_positions --chain_list "$chains_to_design" --position_list "$fixed_positions"

python ../protein_mpnn_run.py \
        --jsonl_path $path_for_parsed_chains \
        --chain_id_jsonl $path_for_assigned_chains \
        --fixed_positions_jsonl $path_for_fixed_positions \
        --out_folder $output_dir \
        --num_seq_per_target 100 \
        --sampling_temp "0.1" \
        --seed 37 \
        --batch_size 1
