import argparse
import subprocess
import glob
import pickle
import utils
import pandas as pd
import seq_processing_utils as seq_utils
import pangenome_utils as pg_utils
import time
from pangenome_utils import PangenomeMap
from Bio import SeqIO
from Bio import AlignIO
from Bio.SeqRecord import SeqRecord

def test_all():
    start_time = time.time()
    output_dir = '../results/tests/pangenome_construction/n4/'
    #pangenome_map = PangenomeMap('../data/single-cell/jgi_gff/')
    #pangenome_map.export_protein_seqs(f'{output_dir}closely_related_sag_proteins.faa')
    pangenome_map = PangenomeMap('../data/tests/pangenome_construction/n4/')
    #pangenome_map.export_protein_seqs(f'{output_dir}n4/sag_protein_seqs.faa')
    pangenome_map.export_rrna_seqs(f'{output_dir}filtered_orthogroups/')

    #subprocess.call(['mkdir', '-p', f'{output_dir}_blastp_results/'])
    #pg_utils.calculate_pairwise_identities(f'{output_dir}closely_related_sag_proteins.faa', output_dir=f'{output_dir}_blastp_results/', num_threads=4, ext='.faa')
    #run_time = utils.timeit(start_time)
    #print(run_time)

def merge_blast_results(blast_results_files, f_out):
    with open(f_out, 'w') as out_handle:
        for f_results in blast_results_files:
            with open(f_results, 'r') as in_handle:
                for line in in_handle.readlines():
                    out_handle.write(line)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input_file', help='FASTA file with protein sequences.')
    parser.add_argument('-I', '--input_dir', help='Directory with SAG input files.')
    parser.add_argument('-d', '--database', default=None, help='Path to BLASTP reference database.')
    parser.add_argument('-b', '--blastp', default='blastp', help='Path to BLASTP binary.')
    parser.add_argument('-n', '--num_threads', type=int, default=1, help='Number of CPU threads.')
    parser.add_argument('-O', '--output_dir', help='Output directory. Note, by default all results will be stored in output_dir/_blastp_results/.')
    parser.add_argument('-o', '--output_file', help='Output file (for single query file only).')
    parser.add_argument('-p', '--prefix', default='sscs', help='Output files prefix.')
    parser.add_argument('--make_sequence_files', action='store_true', help='Create reference database with all protein sequences.')
    parser.add_argument('--make_identity_graph', action='store_true', help='Merge BLASTP results and create .abc graph.')
    parser.add_argument('-t', '--test_all', action='store_true', help='Run tests.')
    args = parser.parse_args()

    if args.test_all == True:
        test_all()

    elif args.make_sequence_files == True:
        pangenome_map = PangenomeMap(args.input_dir)
        pangenome_map.export_protein_seqs(f'{args.output_dir}{args.prefix}_protein_seqs.faa')
        subprocess.call(['mkdir', '-p', f'{args.output_dir}filtered_orthogroups/'])
        pangenome_map.export_rrna_seqs(f'{args.output_dir}filtered_orthogroups/')

        subprocess.call(['mkdir', '-p', f'{args.output_dir}_blastp_results/'])
        pangenome_map.export_sag_protein_seqs(f'{args.output_dir}_blastp_results/')

        # Make BLASTP database
        if args.database:
            f_blastdb = args.database
        else:
            f_blastdb = f'{args.output_dir}_blastp_results/{args.prefix}_protein_seqs.blastdb'
        makeblastdb_out = subprocess.run(['makeblastdb', '-dbtype', 'prot', '-out', f_blastdb, '-in', f'{args.output_dir}{args.prefix}_protein_seqs.faa'], stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.STDOUT)

    elif args.make_identity_graph == True:
        blast_results_files = glob.glob(f'{args.input_dir}*_blast_results.tab')
        merge_blast_results(blast_results_files, f'{args.input_dir}sag_protein_blast_results.tab')
        pg_utils.write_protein_graph(f'{args.input_dir}sag_protein_blast_results.tab', args.output_file)

    else:
        if args.database:
            f_blastdb = args.database
        else:
            f_blastdb = f'{args.output_dir}_blastp_results/sag_protein_seqs.blastdb'

        blastp_out = subprocess.run([args.blastp, '-db', f_blastdb, '-num_threads', f'{args.num_threads}', '-word_size', f'3', '-evalue', f'1E-3', '-outfmt', '6 std qlen slen', '-qcov_hsp_perc', f'75', '-out', args.output_file, '-query', args.input_file], stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.STDOUT)

