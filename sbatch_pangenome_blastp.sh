#!/bin/bash -l
#SBATCH --job-name=run_blastp
#SBATCH --mail-type=END         # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=gbirzu@stanford.edu # Where to send mail
#SBATCH --time=24:00:00 # Time limit hrs:min:sec (max 24 hrs on Sherlock normal queue)
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16GB
#SBATCH --partition=hns
#SBATCH --output=slurm_run_blastp_%A.out
#SBATCH --array=1-10 	# Array range

ml ncbi-blast+

input_dir='../results/single-cell/sscs_pangenome/_blastp_results/'
input_files=($(find ${input_dir} -name '*.faa'))
in_database="${input_dir}sscs_protein_seqs.blastdb"

num_files=${#input_files[*]}
num_jobs=10
batch_size=$((${num_files} / ${num_jobs}))
#SLURM_ARRAY_TASK_ID=1
idx=$((${SLURM_ARRAY_TASK_ID} - 1))

i_start=$((${idx} * ${batch_size}))
if [ ${idx} -lt $((${num_jobs} - 1)) ]
then
    i_end=$(((${idx} + 1) * ${batch_size}))
else
    i_end=${num_files}
fi

for ((i=${i_start}; i<${i_end}; i++))
do
    in_file=${input_files[${i}]}
    sag_id=$(echo ${in_file} | sed 's/.*\///g' | sed 's/.faa$//g')
    python3 pg_make_protein_identity_graph.py -i ${in_file} -d ${in_database} -o "${input_dir}${sag_id}_blast_results.tab" -n 4
done

