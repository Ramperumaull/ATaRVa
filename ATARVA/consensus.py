import pyabpoa as pa
import statistics as stats
from collections import Counter
from itertools import chain

def seq_organiser(seqs):
    seq_len_dict = {}
    for each_seq in seqs:
        lent = len(each_seq)
        if lent in seq_len_dict:
            seq_len_dict[lent].append(each_seq)
        else:
            seq_len_dict[lent] = [each_seq]
    sorted_mode_seqs = sorted(seq_len_dict.values(), key=lambda x : len(x), reverse=True)
    flat_list = list(chain.from_iterable(sorted_mode_seqs))
    del seq_len_dict, sorted_mode_seqs
    return flat_list

def consensus_seq_poa(seqs):
    if len(seqs)<7:
        cons_algrm='MF'
    else:
        cons_algrm='HB'

    median_len = len(stats.median_high(seqs))
    if median_len > 10000:
        fseqs = [seq for seq in seqs if len(seq) == median_len]
        return ''.join(Counter(bases).most_common(1)[0][0] for bases in zip(*fseqs))

    # sorted_seqs = sorted(seqs, key=lambda x : len(x), reverse=True)
    sorted_seqs = seq_organiser(seqs)
    abpoa = pa.msa_aligner(cons_algrm=cons_algrm)
    result = abpoa.msa(sorted_seqs, out_cons=True, out_msa=False)
    return result.cons_seq[0]