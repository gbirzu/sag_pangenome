#!/bin/bash
set -e
set -u

usage_message="Usage: split_deep_branches.sh -i INPUT_TREES -D DATA_DIR -S SEQS_DIR -A ALN_DIR [-t NUM_THREADS]"
num_threads=1

# Read input parameters
# TODO: Add scripts dir
while getopts "i:D:S:A:t:" OPTION
do
  case $OPTION in
    i)
      input_trees=${OPTARG}
      ;;
    D)
      data_dir=${OPTARG}
      ;;
    S)
      seqs_dir=${OPTARG}
      ;;
    A)
      aln_dir=${OPTARG}
      ;;
    t)
      num_threads=${OPTARG}
      ;;
    :)
      echo "Error: -${OPTARG} requires an argument."
      echo $usage_message 1>&2
      exit
      ;;
    ?)
      echo $usage_message 1>&2
      exit
      ;;
  esac
done

# Check input parameters are valid
if [ $OPTIND -lt 9 ]
then
    echo "Error: Some parameters were not specified; -i -D -S -A are all required."
    echo $usage_message 1>&2
    exit
fi


while read in_file
do
    echo $in_file
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

    #python3 pg_split_deep_branches.py -D ${data_dir} -i ${in_file} -f ${out_file} -u ${updates_file} -b ${branch_cutoff} -s nucl -e fna -p sscs

done <$input_trees

