#!/bin/bash -l
#SBATCH --job-name=align_orthogroups
#SBATCH --mail-type=END         # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=gbirzu@stanford.edu # Where to send mail
#SBATCH --time=6:00:00 # Time limit hrs:min:sec (max 24 hrs on Sherlock normal queue)
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16GB
#SBATCH --partition=hns
#SBATCH --output=slurm_align_orthogroups_%A.out
#SBATCH --array=1-10 	# Array range

ml ncbi-blast+

#input_dir='../results/tests/pangenome_construction/sscs_v1/filtered_orthogroups/'
#output_dir='../results/tests/pangenome_construction/sscs_v1/_aln_results/'
input_dir='../results/single-cell/sscs_pangenome/filtered_orthogroups/'
#input_dir='../results/single-cell/sscs_pangenome/fine_scale_og_clusters/'
output_dir='../results/single-cell/sscs_pangenome/_aln_results/'
input_files=($(find ${input_dir} -name '*.fna')) # Use unsorted order to spread out large alignments among jobs
#input_files=($(find ${input_dir} -name '*.fna' | sort)) # Use unsorted order to spread out large alignments among jobs
#input_files=($(find ${input_dir} -name 'YSG_1353.fna')) # Use unsorted order to spread out large alignments among jobs
mkdir -p ${output_dir}

num_files=${#input_files[*]}
num_jobs=10
#num_jobs=1
#SLURM_ARRAY_TASK_ID=1
batch_size=$((${num_files} / ${num_jobs}))
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
    num_seqs=$(cat ${in_file} | grep '^>' | wc -l)
    if [ ${num_seqs} -gt 1 ]
    then
        echo "Aligning ${in_file}"
        og_id=$(echo ${in_file} | sed 's/.*\///g' | sed 's/\.fna$//g')

        if [[ "${og_id}" =~ "rRNA" ]]
        then
            # Use nucleotide alignment for rRNA
            out_nucl_aln="${output_dir}${og_id}_aln.fna"
            mafft --thread 4 --quiet --reorder --auto ${in_file} > ${out_nucl_aln}
        else
            out_aa="${input_dir}${og_id}.faa"
            python3 pg_process_alignment_files.py -i ${in_file} -o ${out_aa} --translate_sequences

            out_aa_aln="${output_dir}${og_id}_aln.faa"
            mafft --thread 4 --quiet --reorder --auto ${out_aa} > ${out_aa_aln}

            out_nucl_aln="${output_dir}${og_id}_aln.fna"
            python3 pg_process_alignment_files.py -I ${input_dir} -i ${out_aa_aln} -o ${out_nucl_aln} --back_translate
        fi

        out_tree="${output_dir}${og_id}_aln.nwk"
        FastTree -nj -noml -nt ${out_nucl_aln} > ${out_tree}

        if [ $? -eq 0 ]
        then
            # Clean up
            rm -f ${out_aa}
            rm -f ${out_aa_aln}
        fi
    else
        echo "${in_file} has only ${num_seqs} sequences. No alignment is necessary."
    fi

    echo
done

