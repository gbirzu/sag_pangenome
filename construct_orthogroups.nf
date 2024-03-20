params.assemblies = "$projectDir/data/*.gff"
params.data_dir = "$projectDir/data"
params.scripts_dir = "$projectDir/scripts"
params.out_dir = "$projectDir/results"
params.mcl_inflation = 1.5

process makeProteinSeqFiles {
    input:
    tuple val(sample_id), path(assembly)

    output:
    tuple val(sample_id), path("${sample_id}_genes.faa")

    script:
    """
    python3 ${params.scripts_dir}/export_gff_fasta_records.py -i $assembly -f genes -s prot -o ${sample_id}_genes.faa
    """

}

process makeBlastDatabase {
    input:
    path(genes_files)

    output:
    path "genes_blastdb*"

    script:
    """
    cat $genes_files > genes.faa
    makeblastdb -dbtype prot -out genes_blastdb -in genes.faa 
    """
}

process blastProteins {
    cpus 4

    input:
    tuple(val(sample_id), path(proteins))
    val(db)

    output:
    path "${sample_id}_blast_results.tsv"

    script:
    """
    blastp -db $db -num_threads $task.cpus -word_size 3 -evalue 1E-3 -outfmt '6 std qlen slen' -qcov_hsp_perc 75 -out ${sample_id}_blast_results.tsv -query $proteins
    """
}

process mergeBlastResults {
    input:
    path blast_results

    output:
    path "merged_blast_results.tsv"

    script:
    """
    cat ${blast_results} > merged_blast_results.tsv
    """
}

process clusterProteins {
    cpus 8

    input:
    path blast_results

    output:
    tuple path("mcl_cluster_seqs"), path("raw_orthogroup_table.tsv")

    script:
    """
    python3 ${params.scripts_dir}/make_bitscore_graph.py -i $blast_results -o bitscore_graph.abc
    mcl bitscore_graph.abc --abc -I ${params.mcl_inflation} -te ${task.cpus} -o mcl_clusters.tsv
    python3 ${params.scripts_dir}/construct_orthogroup_table.py -A "${params.data_dir}/" -i mcl_clusters.tsv -o raw_orthogroup_table.tsv
    mkdir mcl_cluster_seqs
    python3 ${params.scripts_dir}/extract_orthogroup_sequences.py -A "${params.data_dir}/" -O mcl_cluster_seqs/ -g raw_orthogroup_table.tsv -s nucl
    """
}

process filterClusterLengths {
    input:
    tuple path(mcl_dir), path(orthogroup_table)

    output:
    path "filtered_orthogroups"

    script:
    """
    mkdir -p filtered_orthogroups
    python3 ${params.scripts_dir}/filter_sequence_fragments.py -I $mcl_dir/ -O filtered_orthogroups/ -g ${orthogroup_table} -o filtered_orthogroup_table.tsv

    mkdir -p $projectDir/results/filtered_orthogroups
    #cp filtered_cluster_seqs/*.fna $projectDir/results/filtered_orthogroups
    """
}

process getSequenceFiles {
    input:
    path seqs_dir

    output:
    path "seq_files.txt"

    script:
    """
    find -L "\$PWD/${seqs_dir}" -name '*.fna' > seq_files.txt
    """
}

process batchAlignSeqsAndConstructTrees {
    cpus 2

    input:
    path in_files

    output:
    path "_aln_results"

    script:
    """
    mkdir -p _aln_results
    for f_in in ${in_files.join(' ')}
    do
        num_seqs=\$(grep '^>' \$f_in | wc -l)
        if [ \$num_seqs -gt 1 ]
        then
            og_id=\$(echo \$f_in | awk -F'/' '{print \$NF}' | awk -F'.' '{print \$1}')
            out_aa="\${og_id}.faa"
            python3 ${params.scripts_dir}/process_alignments.py -i \$f_in -o \$out_aa

            out_aa_aln="\${og_id}_aln.faa"
            mafft --thread ${task.cpus} --quiet --reorder --auto \${out_aa} > \${out_aa_aln}

            out_aln="_aln_results/\${og_id}_aln.fna"
            python3 ${params.scripts_dir}/process_alignments.py -i \${out_aa_aln} -n \${f_in} -o \${out_aln} -d backward

            out_tree="_aln_results/\${og_id}_tree.nwk"
            FastTree -nj -noml -nt \${out_aln} > \${out_tree}

        fi
    done

    # Clean up
    if [ \$? -eq 0 ]
    then
        rm -f *.faa
    fi
    """
}

process batchSplitDeepBranches {
    cpus 2

    input:
    path tree_files

    script:
    """
    seqs_dir=${params.out_dir}/filtered_orthogroups/
    for f_tree in ${tree_files.join(' ')}
    do
        og_id=\$(echo \$f_tree | awk -F'/' '{print \$NF}' | awk -F'_' '{print \$1"_"\$2}')
        out_file="\${og_id}_subclusters.txt"
        updates_file="\${og_id}_og_updates.dat"
        python3 ${params.scripts_dir}/pg_split_deep_branches.py -S \${seqs_dir} -i \${f_tree} -f \${out_file} -u \${updates_file} -b 0.3 -s nucl -e fna -p sscs

        subcluster_list=(\$(cat \${out_file}))
        if [ \${#subcluster_list[*]} -gt 0 ]
        then
            has_subclusters=1
        else
            has_subclusters=0
        fi

        i_subcluster=0
        while [ \${#subcluster_list[*]} -gt 0 ]
        do
            subcluster_id=\${subcluster_list[\${i_subcluster}]}
            f_subcluster=\${seqs_dir}\${subcluster_id}.fna
            num_seqs=\$(cat \${f_subcluster} | grep '^>' | wc -l)

            if [ "\${num_seqs}" -gt 1 ]
            then
                out_aa="\${subcluster_id}.faa"
                python3 ${params.scripts_dir}/process_alignments.py -i \$f_subcluster -o \$out_aa
                out_aa_aln="\${subcluster_id}_aln.faa"
                mafft --thread ${task.cpus} --quiet --reorder --auto \${out_aa} > \${out_aa_aln}
                out_aln="\${subcluster_id}_aln.fna"
                python3 ${params.scripts_dir}/process_alignments.py -i \${out_aa_aln} -n \${f_subcluster} -o \${out_aln} -d backward
                out_tree="\${subcluster_id}_tree.nwk"
                FastTree -nj -noml -nt \${out_aln} > \${out_tree}

                tail -n +2 \${out_file} > \${og_id}_temp.txt
                mv \${og_id}_temp.txt \${out_file}
                python3 ${params.scripts_dir}/pg_split_deep_branches.py -S \${seqs_dir} -i \${out_tree} -f \${out_file} -u \${updates_file} -b 0.3 -s nucl -e fna -p sscs

                # Clean up
                rm -f \${out_aa}
                rm -f \${out_aa_aln}
            else
                tail -n +2 \${out_file} > \${og_id}_temp.txt
                mv \${og_id}_temp.txt \${out_file}
            fi

            subcluster_list=(\$(cat \${out_file}))
        done

        # Clean up
        rm -f \${out_file}
        if [ \${has_subclusters} -eq 0 ]
        then
            rm -f \${updates_file}
        fi

    done
    """
}

process splitDeepBranches {
    cpus 2

    input:
    path _aln_results
    path filtered_orthogroups

    script:
    """
    tree_files=(\$(find -L _aln_results -name '*.nwk'))
    for f_tree in \${tree_files[*]}
    do
        echo \$f_tree
        og_id=\$(echo \$f_tree | awk -F'/' '{print \$NF}' | awk -F'_' '{print \$1"_"\$2}')
        out_file="\${og_id}_subclusters.txt"
        updates_file="\${og_id}_og_updates.dat"
        python3 ${params.scripts_dir}/pg_split_deep_branches.py -S filtered_orthogroups/ -i \${f_tree} -f \${out_file} -u \${updates_file} -b 0.3 -s nucl -e fna -p sscs

        subcluster_list=(\$(cat \${out_file}))
        if [ \${#subcluster_list[*]} -gt 0 ]
        then
            has_subclusters=1
        else
            has_subclusters=0
        fi

        i_subcluster=0
        while [ \${#subcluster_list[*]} -gt 0 ]
        do
            subcluster_id=\${subcluster_list[\${i_subcluster}]}
            f_subcluster=filtered_orthogroups/\${subcluster_id}.fna
            num_seqs=\$(cat \${f_subcluster} | grep '^>' | wc -l)

            if [ "\${num_seqs}" -gt 1 ]
            then
                out_aa="\${subcluster_id}.faa"
                python3 ${params.scripts_dir}/process_alignments.py -i \$f_subcluster -o \$out_aa
                
                out_aa_aln="\${subcluster_id}_aln.faa"
                mafft --thread ${task.cpus} --quiet --reorder --auto \${out_aa} > \${out_aa_aln}

                out_aln="\${subcluster_id}_aln.fna"
                python3 ${params.scripts_dir}/process_alignments.py -i \${out_aa_aln} -n \${f_subcluster} -o \${out_aln} -d backward

                out_tree="\${subcluster_id}_tree.nwk"
                FastTree -nj -noml -nt \${out_aln} > \${out_tree}

                tail -n +2 \${out_file} > \${og_id}_temp.txt
                mv \${og_id}_temp.txt \${out_file}
                python3 ${params.scripts_dir}/pg_split_deep_branches.py -S filtered_orthogroups/ -i \${out_tree} -f \${out_file} -u \${updates_file} -b 0.3 -s nucl -e fna -p sscs

                # Clean up
                rm -f \${out_aa}
                rm -f \${out_aa_aln}
            else
                tail -n +2 \${out_file} > \${og_id}_temp.txt
                mv \${og_id}_temp.txt \${out_file}
            fi

            subcluster_list=(\$(cat \${out_file}))
        done

        # Clean up
        rm -f \${out_file}
        if [ \${has_subclusters} -eq 0 ]
        then
            rm -f \${updates_file}
        fi
    done

    # Copy any new alignments and trees
    nwk_files=\$(find ./ -name '*.nwk' | wc -l)
    if [ \${nwk_files} -gt 0 ]
    then
        cp *.nwk _aln_results
    fi
        
    aln_files=\$(find ./ -name '*_aln.fna' | wc -l)
    if [ \${aln_files} -gt 0 ]
    then
        cp *_aln.fna _aln_results
    fi
    """
}


workflow {
    Channel
        .fromPath(params.assemblies)
        .map { tuple(it.baseName.split('.gff')[0], it) }
        .set { assemblies_ch }

    seqs_ch = makeProteinSeqFiles(assemblies_ch)
    blastdb_ch = makeBlastDatabase(seqs_ch.map { it[1] }.collect())
    db_path = blastdb_ch
        .flatten()
        .map { "${it.parent}/${it.baseName}" }
        .unique()
        .first() // make value channel to BLAST all sequence files against database
    blast_ch = blastProteins(seqs_ch, db_path)
    merged_blast_ch = mergeBlastResults(blast_ch.collect())
    mcl_ch = clusterProteins(merged_blast_ch)
    seq_clusters_ch = filterClusterLengths(mcl_ch)

    seq_files_ch = getSequenceFiles(seq_clusters_ch)
    seq_files_ch
        .splitCsv()
        .map { row -> file(row[0]) }
        .buffer(size: 10, skip: 1000, remainder: true)
        .set { batched_files_ch }
    aln_ch = batchAlignSeqsAndConstructTrees(batched_files_ch)
    aln_ch.view()
    splitDeepBranches(aln_ch, seq_clusters_ch)

    //Channel
    //    .fromPath("${params.out_dir}/_aln_results/*.nwk")
    //    .set { tree_files_ch }
    //batchSplitDeepBranches(tree_files_ch.buffer(size: 10, remainder: true))
}
