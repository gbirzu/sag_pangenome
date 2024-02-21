import argparse
import pandas as pd
import pangenome_utils as pg_utils

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
    parser.add_argument('-A', '--annotations_dir', help='Directory with SAG annotations in .gff format.')
    parser.add_argument('-i', '--input_cluster_table', help='Input graph in .abc format.')
    parser.add_argument('-o', '--output_file', help='Output file.')
    args = parser.parse_args()

    cluster_dict, cluster_idx = pg_utils.process_mcl_clusters(args.input_cluster_table, idx=1)
    pangenome_map = pg_utils.PangenomeMap(args.annotations_dir)
    og_table = make_orthogroup_table(cluster_dict, pangenome_map)
    og_table.to_csv(args.output_file, sep='\t')
