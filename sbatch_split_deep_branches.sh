#!/bin/bash -l
#SBATCH --job-name=split_deep_branches
#SBATCH --mail-type=END         # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=gbirzu@stanford.edu # Where to send mail
#SBATCH --time=12:00:00 # Time limit hrs:min:sec (max 24 hrs on Sherlock normal queue)
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16GB
#SBATCH --partition=hns
#SBATCH --output=slurm_split_deep_branches_%A.out
#SBATCH --array=1-10 	# Array range

data_dir='../results/single-cell/sscs_pangenome/'
seqs_dir="${data_dir}filtered_orthogroups/"
aln_dir="${data_dir}_aln_results/"
tree_files=($(find ${aln_dir} -name '*.nwk'))
#data_dir='../results/tests/pangenome_construction/sscs_v1/'
#tree_files=($(find ${aln_dir} -name '*.nwk' | sort))
#tree_files=("${aln_dir}YSG_1376_aln.nwk")
#tree_files=("${aln_dir}YSG_1353_aln.nwk")

num_files=${#tree_files[*]}
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
    in_file=${tree_files[${i}]}
    og_id=$(echo ${in_file} | sed 's/.*\///g' | sed 's/_aln\.nwk$//g')
    echo "Analyzing ${og_id}..."
    out_file="${seqs_dir}${og_id}_subclusters.txt"
    updates_file="${seqs_dir}${og_id}_og_updates.dat"
    if [ -f "${out_file}" ]
    then
        rm -f ${out_file}
    fi

    if [[ "${og_id}" =~ "rRNA" ]]
    then
        branch_cutoff=0.05
    else
        branch_cutoff=0.3
    fi
    python3 pg_split_deep_branches.py -D ${data_dir} -i ${in_file} -f ${out_file} -u ${updates_file} -b ${branch_cutoff} -s nucl -e fna -p sscs

    subcluster_list=($(cat ${out_file}))
    if [ ${#subcluster_list[*]} -gt 0 ]
    then
        has_subclusters=1
    else
        has_subclusters=0
    fi

    i_subcluster=0
    while [ ${#subcluster_list[*]} -gt 0 ]
    do
        subcluster_id=${subcluster_list[${i_subcluster}]}
        f_subcluster=${seqs_dir}${subcluster_id}.fna
        num_seqs=$(cat ${f_subcluster} | grep '^>' | wc -l)

        if [ "${num_seqs}" -gt 1 ]
        then
            sog_id=$(echo ${f_subcluster} | sed 's/.*\///g' | sed 's/\.fna$//g')
            out_nucl_aln="${aln_dir}${sog_id}_aln.fna"

            if [[ "${og_id}" =~ "rRNA" ]]
            then
                # Use nucleotide alignment for RNA
                mafft --thread 4 --quiet --reorder --auto ${f_subcluster} > ${out_nucl_aln}
            else
                # Use codone-aware alignment for proteins
                out_aa="${seqs_dir}${sog_id}.faa"
                python3 pg_process_alignment_files.py -i ${f_subcluster} -o ${out_aa} --translate_sequences
                out_aa_aln="${aln_dir}${sog_id}_aln.faa"
                mafft --thread 4 --quiet --reorder --auto ${out_aa} > ${out_aa_aln}
                python3 pg_process_alignment_files.py -I ${seqs_dir} -i ${out_aa_aln} -o ${out_nucl_aln} --back_translate
            fi

            out_tree="${aln_dir}${sog_id}_aln.nwk"
            FastTree -nj -noml -nt ${out_nucl_aln} > ${out_tree}

            tail -n +2 ${out_file} > temp.txt
            mv temp.txt ${out_file}
            python3 pg_split_deep_branches.py -D ${data_dir} -i ${out_tree} -f ${out_file} -u ${updates_file} -s nucl -e fna -p sscs -v

            # Clean up
            rm -f ${out_aa}
            rm -f ${out_aa_aln}
        else
            tail -n +2 ${out_file} > temp.txt
            mv temp.txt ${out_file}
        fi

        subcluster_list=($(cat ${out_file}))
    done

    # Clean up
    rm -f ${out_file}
    if [ ${has_subclusters} -eq 0 ]
    then
        rm -f ${updates_file}
    fi

done

