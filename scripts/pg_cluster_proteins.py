import argparse
import subprocess
import glob
import pickle
import utils
import numpy as np
import pandas as pd
import pangenome_utils as pg_utils
import time
from pangenome_utils import PangenomeMap
from Bio import SeqIO
from Bio import AlignIO
from Bio.SeqRecord import SeqRecord

def test_all():
    print('Running pg_cluster_proteins.py tests...')
    data_dir = '../results/tests/pangenome_construction/large_cells/_mcl_results/'
    annotations_dir = '../data/tests/pangenome_construction/large_cells/'

    cluster_index = 1
    cluster_dict, cluster_idx = pg_utils.process_mcl_clusters(f'{data_dir}/sag_protein_clusters.tsv', idx=cluster_index)
    pangenome_map = PangenomeMap(annotations_dir)
    og_table = make_orthogroup_table(cluster_dict, pangenome_map)

    print(og_table)
    print(og_table.loc[og_table['seqs_per_cell'] > 1, :])

def make_orthogroup_table(cluster_dict, pangenome_map):
    og_ids = list(cluster_dict.keys())
    sag_ids = sorted([sag_id.strip() for sag_id in list(pangenome_map.cell_contigs.keys())])
    gene_id_map = pangenome_map.get_gene_id_map()
    og_table = pd.DataFrame(index=og_ids, columns=['num_seqs', 'num_cells', 'seqs_per_cell', 'avg_length'] + sag_ids)

    for og_id in cluster_dict:
        gene_ids = cluster_dict[og_id]
        og_table.loc[og_id, 'num_seqs'] = len(gene_ids)
        og_table.loc[og_id, 'avg_length'] = pg_utils.calculate_mean_gene_length(gene_ids)

        sag_grouped_ids = pg_utils.group_gene_ids(gene_ids, gene_id_map)
        og_table.loc[og_id, 'num_cells'] = len(sag_grouped_ids)

        for sag_id in sag_grouped_ids:
            og_table.loc[og_id, sag_id] = ';'.join(sag_grouped_ids[sag_id])

    og_table['seqs_per_cell'] = og_table['num_seqs'] / og_table['num_cells']

    return og_table


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input_graph', help='Input graph in .abc format.')
    parser.add_argument('-s', '--sequences', help='Sequences file used to generate graph.')
    parser.add_argument('-A', '--annotations_dir', help='Directory with SAG annotations in .gff format.')
    parser.add_argument('-m', '--mcl', default='mcl', help='Path to mcl binary.')
    parser.add_argument('-n', '--num_threads', type=int, default=1, help='Number of CPU threads.')
    parser.add_argument('-O', '--output_dir', help='Output directory.')
    parser.add_argument('-p', '--prefix', default='sscs', help='Output files prefix.')
    parser.add_argument('--seq_type', default='nucl', help='Output sequence type ["prot" | "nucl"].')
    parser.add_argument('-t', '--test_all', action='store_true', help='Run tests.')
    args = parser.parse_args()

    if args.test_all == True:
        test_all()
    else:
        print(f'Clustering {args.input_graph} protein graph...')
        start_time = time.time()

        subprocess.call(['mkdir', '-p', f'{args.output_dir}'])
        pg_utils.run_mcl(args.input_graph, f'{args.output_dir}{args.prefix}_protein_clusters.tsv', num_threads=args.num_threads)
        cluster_index = 1
        cluster_dict, cluster_idx = pg_utils.process_mcl_clusters(f'{args.output_dir}{args.prefix}_protein_clusters.tsv', idx=cluster_index)
        pangenome_map = PangenomeMap(args.annotations_dir)

        # Save OG table
        og_table = make_orthogroup_table(cluster_dict, pangenome_map)
        og_table.to_csv(f'{args.output_dir}{args.prefix}_orthogroup_presence.tsv', sep='\t')

        # Write protein cluster FASTA files
        if args.seq_type == 'nucl':
            ext = 'fna'
        else:
            protein_seqs = utils.read_fasta(args.sequences)
            ext = 'faa'

        for cluster_id in cluster_dict:
            cluster_genes = cluster_dict[cluster_id]
            seq_records = []
            for seq_id in cluster_genes:
                if args.seq_type == 'nucl':
                    record = pangenome_map.extract_gene_record(seq_id)
                else:
                    record = protein_seqs[seq_id]
                seq_records.append(record)
            SeqIO.write(seq_records, f'{args.output_dir}{cluster_id}.{ext}', 'fasta')

        run_time = utils.timeit(start_time)
        print(f'Done! Total pg_cluster_proteins.py runtime: {run_time}')
