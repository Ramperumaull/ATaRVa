import hdbscan
import numpy as np
import warnings
import statistics
from ATARVA.consensus import *
from ATARVA.decomp_utils import motif_decomposition

methviz_tag = False
def set_methviz_tag(value):
    global methviz_tag
    methviz_tag = value

# base64 encodings
encode64_dict = {0: 'A', 1: 'B', 2: 'C', 3: 'D', 4: 'E', 5: 'F', 6: 'G', 7: 'H', 8: 'I', 9: 'J', 10: 'K', 11: 'L', 12: 'M', 13: 'N', 14: 'O', 15: 'P', 16: 'Q',17: 'R', 18: 'S', 19: 'T', 20: 'U', 21: 'V', 22: 'W',
                23: 'X', 24: 'Y', 25: 'Z', 26: 'a', 27: 'b', 28: 'c', 29: 'd', 30: 'e', 31: 'f', 32: 'g', 33: 'h', 34: 'i', 35: 'j', 36: 'k', 37: 'l', 38: 'm', 39: 'n', 40: 'o', 41: 'p', 42: 'q', 43: 'r', 44: 's',
                45: 't', 46: 'u', 47: 'v', 48: 'w', 49: 'x', 50: 'y', 51: 'z', 52: '0', 53: '1', 54: '2', 55: '3', 56: '4', 57: '5', 58: '6', 59: '7', 60: '8', 61: '9',62: '+', 63: '/', 64: '/'}

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
    # upper_bound = meth_cutoff
    # lower_bound = 1 - meth_cutoff
    read_meth_range = []
    last_index = len(read_sequence)-1
    if (meth_start!=None) and (meth_end!=None):
        for each_pos in pos_qual:
            # current_prob = each_pos[1]/255
            # if lower_bound < current_prob < upper_bound: # skip the MM if it is within the lower and upper bound; only store the extremies
            #     continue
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
    encrypted_meth = None
    locus_read_meth = global_loci_variations[locus_key]['read_meth']
    matrix = []
    max_meth = 0
    for read_id in hap_reads:
        if locus_read_meth[read_id] is not None:
            meth_reads += 1
            hap_meth += locus_read_meth[read_id][0]
            meth_probs = locus_read_meth[read_id][1]
            if len(meth_probs) > max_meth: # getting max length to pad the matrix
                max_meth = len(meth_probs)
            matrix.append(meth_probs)
            
    if meth_reads > 0:
        if methviz_tag:
            encrypted_meth = methylation_encoding(matrix, max_meth)
        return [round(hap_meth/meth_reads, 2), meth_reads, encrypted_meth]
    else:
        return [None, None, None]
    
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

def methylation_encoding(matrix, max_meth):
    encryted_meth = ''
    for each_row in matrix:
        if len(each_row) < max_meth:
            diff = max_meth - len(each_row)
            each_row.extend([-2]*diff)
    for col in zip(*matrix):
        col_array = np.array(col)
        mode = statistics.mode(col_array)
        if mode == -2:
            pass # skipping the positions where there is error call in few reads
        elif mode == -1:
            encryted_meth += '-' #adding - for ambiguous calls
        else:
            col_array = col_array[col_array != -1]
            col_mean = round(np.mean(col_array), 2) * 100
            col_mean/=1.5625
            encryted_meth += encode64_dict[col_mean]
    return encryted_meth