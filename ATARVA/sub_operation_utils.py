import hdbscan
import numpy as np
import warnings
from ATARVA.consensus import *
from ATARVA.decomp_utils import motif_decomposition


def dbscan(data, hap_reads):
    data = np.array(data).reshape(-1, 1)
    min_samples = max(10, round(0.2*len(data))) # min 20% of the data or 10 reads
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=FutureWarning)
        clusterer = hdbscan.HDBSCAN(min_cluster_size=min_samples)
        cluster_labels = clusterer.fit_predict(data)
    unique_labels = set(cluster_labels)
    
    if len(unique_labels)==1: # cluster case = (0), (-1)
        return [False,None,None] # proceed with Kmeans
        
    elif (len(unique_labels)==2) and (-1 in unique_labels): # cluster case = (0,-1)
        return [False,None,None] # proceed with Kmeans
        
    elif len(unique_labels)>=2: # cluster case = (0,1), (0,1,-1), (0,1,2)
        main_label = unique_labels-{-1}

        main_clusters = {}
        
        for label in main_label:
            c_label = [i for i, x in enumerate(cluster_labels) if x == label]
            alen = [data[i][0] for i in c_label]
            if len(c_label) in main_clusters:
                main_clusters[len(c_label)+1] = [c_label, alen]
            else:
                main_clusters[len(c_label)] = [c_label, alen]
            
        top2_clus_idx = [v for _,v in sorted(main_clusters.items(), reverse=True)[:2]] # getting top 2 cluster with more support

        new_haplotypes = [[hap_reads[idx] for idx in top2_clus_idx[0][0]], [hap_reads[idx] for idx in top2_clus_idx[1][0]]] # getting respective read ids

        new_alen = [top2_clus_idx[0][1], top2_clus_idx[1][1]]

        if set(new_alen[0])==set(new_alen[1]):
            return [False,None,None]
        
        return [True, new_haplotypes, new_alen]
    
def mm_tag_extract(pos_qual, meth_start, meth_end, read_sequence, meth_cutoff, frwd_strand):
    read_meth_range = []
    last_index = len(read_sequence)-1
    if (meth_start!=None) and (meth_end!=None):
        for each_pos in pos_qual:
            if (each_pos[1]/255) < meth_cutoff:
                continue
            meth_pos = each_pos[0]
            meth_chunk_start = meth_pos if frwd_strand else meth_pos-1 # to check the meth context, start index
            meth_chunk_end = meth_pos+2 if frwd_strand else meth_pos+1 # to check the meth context, end index
            if meth_start <= meth_pos <= meth_end:
                if (meth_pos+1 <= last_index) and (read_sequence[meth_chunk_start : meth_chunk_end]=='CG'):
                    read_meth_range.append(each_pos)
    return read_meth_range
            
def methylation_calc(hap_reads, global_loci_variations, locus_key):
    meth_reads = 0
    hap_meth = 0
    locus_read_meth = global_loci_variations[locus_key]['read_meth']
    for read_id in hap_reads:
        if locus_read_meth[read_id]:
            meth_reads += 1
            hap_meth += locus_read_meth[read_id]
    if meth_reads > 0:
        return [round(hap_meth/meth_reads, 2), meth_reads]
    else:
        return [None, None]
    
def confidence_interval(data):
    data = np.array(data)
    ci = np.percentile(data, [2.5, 97.5])
    return [round(ci[0]), round(ci[1])]

def alt_sequence(read_seqs, hap_reads, amplicon, motif_size):
    seqs = [seq for seq in [read_seqs[read_id][0] for read_id in hap_reads] if seq!='']
    if len(seqs)>0:
        ALT = consensus_seq_poa(seqs)
        allele_length = len(ALT)
    else:
        ALT = '<DEL>'
        allele_length = 0

    decomp_seq = ''
    repeativity = True
    if amplicon and allele_length and (motif_size<=10):
        decomp_seq, nonrep_percent = motif_decomposition(ALT, motif_size)
        # nonrep_percent = non_repeat_length(decomp_seq)/len(ALT)
        if nonrep_percent > 0.30: # if more than 30% of the sequence is non-repeat, repeativity = False
            repeativity = False
    return [ALT, allele_length, decomp_seq, repeativity]