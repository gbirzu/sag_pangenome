import argparse
import subprocess
import os
import glob
import pickle
import utils
import numpy as np
import pandas as pd
import seq_processing_utils as seq_utils
import time
import ete3
import pangenome_utils as pg_utils
from pangenome_utils import PangenomeMap
from Bio import Phylo
from Bio import SeqIO
from Bio import AlignIO
from Bio.SeqRecord import SeqRecord

def read_trees_and_split_subclusters(args):
    subcluster_dict = {}
    tree_files = sorted(glob.glob(f'{args.data_dir}_aln_results/*.nwk'))
    for f_tree in tree_files:
        subclusters = find_subclusters(f_tree, args.branch_cutoff)
        if subclusters:
            cluster_id = f_tree.split('/')[-1].replace('_aln.nwk', '')
            subcluster_dict[cluster_id] = subclusters
    return subcluster_dict

def find_subclusters(f_tree, branch_cutoff):
    tree = ete3.PhyloTree(f_tree, format=0)
    cluster_id = f_tree.split('/')[-1].replace('_aln.nwk', '')
    tree_subclusters = split_deep_branches(tree, branch_cutoff)
    if tree_subclusters is not None:
        subclusters = label_subclusters(tree_subclusters, cluster_id)
        return subclusters
    else:
        return None

def split_deep_branches(tree, bc):
    '''
    Splits off subtrees with branches longer than bc into
    separate clusters and returs list of subclusters.
    '''

    # Get candidate subcluster leaves
    subtree_leaves = {}
    node_dict = tree.get_cached_content()
    for node in node_dict:
        if node.dist > bc:
            # Add subtree to cluster
            node_leaves = set([leaf.name for leaf in node_dict[node]])
            subtree_leaves[node] = node_leaves

    if len(subtree_leaves) > 0:
        tree_leaves = set([leaf.name for leaf in tree.get_leaves()])

        # Assign leaves to subclusters
        subclusters = {}
        for node in node_dict:
            if node in subtree_leaves:
                subcluster_leaves = tree_leaves.intersection(subtree_leaves[node])
                if len(subcluster_leaves) > 0:
                    subclusters[node] = list(subcluster_leaves)
                    tree_leaves = tree_leaves - subcluster_leaves

        # Add remaining leaves to 'root'
        if len(tree_leaves) > 0:
            subclusters['root'] = list(tree_leaves)

    else:
        subclusters = None

    return subclusters

def label_subclusters(subcluster_dict, cluster_id, labeling_scheme='alphabetic'):
    labeled_subclusters = {}
    counter = 0
    for node in subcluster_dict:
        if labeling_scheme == 'alphabetic':
            subcluster_id = f'{cluster_id}{utils.get_alphabetic_index(counter)}'
        else:
            subcluster_id = f'{cluster_id}-{counter:03d}'
        labeled_subclusters[subcluster_id] = subcluster_dict[node]
        counter += 1
    return labeled_subclusters

def export_subcluster_sequences(split_clusters, seqs_dir, pangenome_map=None, seq_type='prot', cleanup=False):
    if seq_type == 'prot':
        ext = 'faa'
    else:
        ext = 'fna'

    for cluster_id in split_clusters:
        if seqs_dir:
            cluster_records = utils.read_fasta(f'{seqs_dir}{cluster_id}.{ext}')
        for subcluster_id in split_clusters[cluster_id]:
            gene_ids = split_clusters[cluster_id][subcluster_id]

            if pangenome_map:
                subcluster_records = [pangenome_map.extract_gene_record(gene_id, type=seq_type) for gene_id in gene_ids]
            else:
                subcluster_records = [cluster_records[gene_id] for gene_id in gene_ids]
            SeqIO.write(subcluster_records, f'{seqs_dir}{subcluster_id}.{ext}', 'fasta')

    if cleanup == True:
        for cluster_id in split_clusters:
            subprocess.call(['rm', '-f', f'{seqs_dir}{cluster_id}.{ext}'])

def write_subcluster_seqs(subcluster_gene_ids, subcluster_id, cluster_id, seqs_dir, seq_type='prot'):
    if seq_type == 'prot':
        ext = 'faa'
    else:
        ext = 'fna'
    cluster_records, seq_ids = utils.read_fasta(f'{seqs_dir}{cluster_id}.{ext}', return_seq_order=True)
    ordered_subcluster_ids = [gene_id for gene_id in seq_ids if gene_id in subcluster_gene_ids]
    subcluster_records = []
    for gene_id in ordered_subcluster_ids:
        subcluster_records.append(cluster_records[gene_id])
    SeqIO.write(subcluster_records, f'{seqs_dir}{subcluster_id}.{ext}', 'fasta')

def update_orthogroup_table(subcluster_dict, args, pangenome_map=None):
    og_table = pd.read_csv(f'{args.data_dir}filtered_orthogroups/{args.prefix}_filtered_orthogroup_presence.tsv', sep='\t', index_col=0)
    if pangenome_map:
        gene_id_map = pangenome_map.get_gene_id_map()

    for cluster_id in subcluster_dict:
        subclusters = subcluster_dict[cluster_id]
        for og_id in subclusters:
            gene_ids = subclusters[og_id]
            og_table.loc[og_id, 'num_seqs'] = len(gene_ids)
            og_table.loc[og_id, 'avg_length'] = pg_utils.calculate_mean_gene_length(gene_ids)

            if pangenome_map:
                sag_grouped_ids = pg_utils.group_gene_ids(gene_ids, gene_id_map)
            else:
                sag_grouped_ids = pg_utils.group_og_table_genes(gene_ids, og_table)
            og_table.loc[og_id, 'num_cells'] = len(sag_grouped_ids)

            for sag_id in sag_grouped_ids:
                og_table.loc[og_id, sag_id] = ';'.join(sag_grouped_ids[sag_id])
        og_table['seqs_per_cell'] = og_table['num_seqs'] / og_table['num_cells']

    # Remove old clusters
    og_index = list(og_table.index)
    old_cluster_ids = [cluster_id for cluster_id in list(subcluster_dict.keys()) if cluster_id in og_index]
    og_table = og_table.drop(index=old_cluster_ids)
    print(og_table)
    og_table.to_csv(f'{args.data_dir}_mcl_results/{args.prefix}_orthogroup_presence.tsv', sep='\t')

def save_orthogroup_update(subcluster_dict, args):
    og_table = pd.read_csv(f'{args.data_dir}filtered_orthogroups/{args.prefix}_filtered_orthogroup_presence.tsv', sep='\t', index_col=0)
    updates_dict = {}
    for cluster_id in subcluster_dict:
        subclusters = subcluster_dict[cluster_id]
        updates_dict[cluster_id] = []
        for og_id in subclusters:
            table_row = pg_utils.make_subcluster_og_table_row(og_id, subclusters[og_id], og_table)
            updates_dict[cluster_id].append(table_row)
    with open(f'{args.data_dir}filtered_orthogroups/{args.updates_file}', 'wb') as out_handle:
        pickle.dump(updates_dict, out_handle)

def append_to_updates_file(f_updates, subcluster_id, gene_ids):
    if os.path.exists(f_updates):
        cluster_updates = pickle.load(open(f_updates, 'rb'))
        cluster_updates[subcluster_id] = gene_ids
    else:
        cluster_updates = {subcluster_id:gene_ids}
    pickle.dump(cluster_updates, open(f_updates, 'wb'))
    return cluster_updates


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-D', '--data_dir', help='Main directory containing pangenome output.')
    parser.add_argument('-S', '--seqs_dir', help='Directory containing orthogroup sequences.')
    parser.add_argument('-b', '--branch_cutoff', default=0.3, type=float, help='Cutoff length for long branches.')
    parser.add_argument('-e', '--ext', default='fna', help='Sequence files extension.')
    parser.add_argument('-f', '--subclusters_file', default=None, help='Text file where subcluster IDs are saved.')
    parser.add_argument('-i', '--input_file', default=None, help='Input FASTA file with sequences.')
    parser.add_argument('-p', '--prefix', default='sscs', help='Output files prefix.')
    parser.add_argument('-s', '--seq_type', default='nucl', help='Output sequence type ["prot" | "nucl"].')
    parser.add_argument('-u', '--updates_file', default=None, help='File where gene IDs of subclusters are saved.')
    parser.add_argument('-v', '--verbose', action='store_true', help='Run in verbose mode.')
    args = parser.parse_args()

    if args.input_file:
        #cluster_id = args.input_file.split('/')[-1].replace('_aln.nwk', '')
        cluster_id = args.input_file.split('/')[-1].replace('_tree.nwk', '')
        tree = ete3.PhyloTree(args.input_file, format=0)
        tree_subclusters = split_deep_branches(tree, args.branch_cutoff)

        with open(args.subclusters_file, 'a') as out_handle:
            if tree_subclusters is not None:
                aln_subclusters = label_subclusters(tree_subclusters, cluster_id)
                for sog_id in aln_subclusters:
                    gene_ids = aln_subclusters[sog_id]
                    write_subcluster_seqs(gene_ids, sog_id, cluster_id, f'{args.seqs_dir}', seq_type=args.seq_type)
                    out_handle.write(f'{sog_id}\n')

                    # Add single sequence clusters to updates file
                    if len(gene_ids) == 1 and args.updates_file is not None:
                        append_to_updates_file(args.updates_file, sog_id, gene_ids)

            else:
                gene_ids = []
                for node in tree.get_leaves():
                    gene_ids.append(node.name)
                if args.updates_file is not None:
                    cluster_updates = append_to_updates_file(args.updates_file, cluster_id, gene_ids)

                    if args.verbose:
                        print(cluster_updates)


