from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment

#######################################
# Biopython wrappers
#######################################

def write_seqs_dict(seq_dict, output_file, format='fasta'):
    records = list(seq_dict.values())
    SeqIO.write(records, output_file, format)

def write_seqs(seq_records, output_file, format='fasta'):
    SeqIO.write(seq_records, output_file, format)

def make_Seq(seq_str, **kwargs):
    '''
    Wrapper for Biopython Seq constructor
    '''
    return Seq(seq_str, **kwargs)

def make_SeqRecord(seq, **kwargs):
    '''
    Wrapper for Biopython SeqRecord constructor
    '''
    return SeqRecord(seq, **kwargs)


def read_alignment(input_file, format='fasta'):
    return AlignIO.read(open(input_file, 'r'), format)

def write_alignment(aln, output_file, format='fasta'):
    AlignIO.write(aln, output_file, format)

#######################################
# Alignment manipulation
#######################################

def back_translate(aa_aln, original_seqs, replace_ambiguous_chars=False):
    nucl_aln = []
    for rec in aa_aln:
        num_gaps = 0
        aln_seq = []
        nucl_seq = original_seqs[rec.id]
        for i in range(len(rec.seq)):
            if rec[i] == '-':
                aln_seq.append('---')
                num_gaps += 1
            else:
                idx = 3 * (i - num_gaps)
                site_codon_str = str(nucl_seq[idx:(idx + 3)].seq)
                if replace_ambiguous_chars == True and 'N' in site_codon_str:
                    aln_seq.append('---')
                else:
                    aln_seq.append(site_codon_str)
        nucl_aln.append(SeqRecord(Seq(''.join(aln_seq)), id=rec.id, description=rec.description))
    return MultipleSeqAlignment(nucl_aln)



#######################################
# String manipulation
#######################################

def split_gene_id(gene_id):
    id_comps = gene_id.split('_')
    return '_'.join(id_comps[:-2]), (int(id_comps[-2]), int(id_comps[-1]))

if __name__ == '__main__':
    print('No tests to run...')

