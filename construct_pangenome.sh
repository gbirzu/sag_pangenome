#!/bin/bash

#jout=$(sbatch sbatch_filter_cells.sh)
#jid0=$(echo "${jout}" | awk '{print $NF}')

#jout=$(sbatch --dependency=afterok:${jid0} sbatch_make_pangenome_sequence_files.sh)
#jout=$(sbatch sbatch_make_pangenome_sequence_files.sh)
#jid1=$(echo "${jout}" | awk '{print $NF}')

#jout=$(sbatch --dependency=afterok:${jid1} sbatch_pangenome_blastp.sh)
#jid2=$(echo "${jout}" | awk '{print $NF}')

#jout=$(sbatch --dependency=afterok:${jid2} sbatch_cluster_proteins.sh)
#jid3=$(echo "${jout}" | awk '{print $NF}')

#jout=$(sbatch --dependency=afterok:${jid3} sbatch_filter_fragments.sh)
#jout=$(sbatch sbatch_filter_fragments.sh)
#jid3a=$(echo "${jout}" | awk '{print $NF}')

#jout=$(sbatch --dependency=afterok:${jid3a} sbatch_align_orthogroups_and_make_trees.sh)
#jid4=$(echo "${jout}" | awk '{print $NF}')

#jout=$(sbatch --dependency=afterok:${jid4} sbatch_split_deep_branches.sh)
#jid5=$(echo "${jout}" | awk '{print $NF}')

#jout=$(sbatch --dependency=afterok:${jid5} sbatch_update_orthogroups.sh)
#jid6=$(echo "${jout}" | awk '{print $NF}')

#jout=$(sbatch --dependency=afterok:${jid6} sbatch_calculate_pairwise_divergences.sh)
jout=$(sbatch sbatch_calculate_pairwise_divergences.sh)
jid7=$(echo "${jout}" | awk '{print $NF}')

jout=$(sbatch --dependency=afterok:${jid7} sbatch_find_fine_scale_orthogroups.sh)
jid8=$(echo "${jout}" | awk '{print $NF}')

jout=$(sbatch --dependency=afterok:${jid8} sbatch_map_ogs_to_references.sh)
jid9=$(echo "${jout}" | awk '{print $NF}')

