params.assemblies = "$projectDir/data/*.gff"
params.scripts_dir = "$projectDir/scripts"
params.out_dir = "$projectDir/results"

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
    //tuple val(sample_id), path(proteins)
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
    blast_ch.view()
}
