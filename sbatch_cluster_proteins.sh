#!/bin/bash -l
#SBATCH --job-name=cluster_proteins
#SBATCH --mail-type=END         # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=gbirzu@stanford.edu # Where to send mail
#SBATCH --time=24:00:00 # Time limit hrs:min:sec (max 24 hrs on Sherlock normal queue)
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=256GB
#SBATCH --partition=bigmem
#SBATCH --output=slurm_cluster_proteins_%A.out

ml ncbi-blast+

main_dir='../results/single-cell/sscs_pangenome/'
#main_dir='../results/tests/pangenome_construction/large_cells/'
blastp_dir="${main_dir}_blastp_results/"
mcl_dir="${main_dir}_mcl_results/"
mkdir -p ${mcl_dir}

python3 pg_make_protein_identity_graph.py -I ${blastp_dir} -o ${mcl_dir}sscs_protein_identity_graph.abc -p sscs -n 8 --make_identity_graph
python3 pg_cluster_proteins.py -i ${mcl_dir}/sscs_protein_identity_graph.abc -A ../data/single-cell/filtered_annotations/sscs/ -O ${mcl_dir} -p sscs -n 8 --seq_type nucl
#python3 pg_align_seqs_and_build_trees.py -I ${output_dir}_mcl_results/ -O ${output_dir} -n 4

