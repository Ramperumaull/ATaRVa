import numpy as np
from scipy.spatial.distance import squareform
import stringzilla as sz
from scipy.cluster.hierarchy import linkage, dendrogram
from ATARVA.sub_operation_utils import *

def build_distance_matrix(tot_seqs):
    N = len(tot_seqs)
    distance_matrix = []
    for i in range(N):
        for j in range(N):
            if i >= j:
                continue
            a = tot_seqs[i]; b = tot_seqs[j]
            dist = sz.edit_distance(a, b)
            max_len = max(len(a), len(b))
            score = dist/max_len
            distance_matrix.append(score)
    return np.array(distance_matrix)

def correlation_clustering(read_seqs, read_indices, motif_size, global_loci_variations, locus_key):

    seq_list = []
    for read_id in read_indices:
        seq_list.append(read_seqs[read_id][0])

    seq_len = [len(i) for i in seq_list]
    distance_matrix = build_distance_matrix(seq_list)
    linkage_matrix = linkage(distance_matrix, method='complete')

    # # 3rd col of Linkage has the distance values and 85th percentile is calculated to cut the dendrogram to get tight unique clusters
    # dist_percentile = np.percentile(linkage_matrix[:, 2], 85)

    # Finiding the largest jump in the merge distance
    linkage_dist = np.sort(linkage_matrix[:,2]) # sorting the distances
    jumps = np.diff(linkage_dist) # calculating the difference between consecutive distances
    jump_index = np.argmax(jumps) # finding the max jump index
    if jump_index < len(linkage_dist)-1: # ensuring that the jump index is within bounds
        jump_index += 1 # moving to the next index to get the distance after the jump, which will be the cutoff for clustering
    cutoff = linkage_dist[jump_index]

    # Generating dendrogram details without plotting
    den_detail = dendrogram(linkage_matrix, count_sort='descending', color_threshold=cutoff, no_plot=True)
    # Grouping reads based on cluster colors
    color = den_detail['leaves_color_list']
    leaves = den_detail['leaves']
    cluster_dict = {}
    for idx, each_color in enumerate(color):
        if each_color == 'C0': continue # skipping the outlier cluster
        if each_color not in cluster_dict:
            cluster_dict[each_color] = [read_indices[leaves[idx]]] # getting the read_id from leave_index of dendrogram
        else:
            cluster_dict[each_color].append(read_indices[leaves[idx]])

    if len(cluster_dict)==0: # no valid clusters found
        return [False, 6, {}]
    
    haplotypes = list(cluster_dict.values())

    alt_seq_lens = set()
    alt_seqs = set()

    genotype_dict = {}
    for hap_reads in haplotypes:

        if len(hap_reads)<3: continue

        ALT, allele_length, decomp_seq, repeativity = alt_sequence(read_seqs, hap_reads, True, motif_size) # true for amplicon, to check the repetitiveness in the sequence


        # skip if the cluster has less reads and is not a repeat
        if not repeativity:
            continue

        if allele_length in alt_seq_lens:
            if ALT not in alt_seqs:
                alt_seqs.add(ALT)
            else:
                continue
        else:
            alt_seq_lens.add(allele_length)
            alt_seqs.add(ALT) 

        ci = confidence_interval([len(read_seqs[read_id][0]) for read_id in hap_reads])
        meth_info = methylation_calc(hap_reads, global_loci_variations, locus_key, ALT)
        if allele_length not in genotype_dict:
            genotype_dict[allele_length] = (ALT, ci, decomp_seq, len(hap_reads), meth_info)
        else:
            genotype_dict[str(allele_length)] = (ALT, ci, decomp_seq, len(hap_reads), meth_info)

    del read_seqs, alt_seqs, alt_seq_lens, haplotypes, cluster_dict, distance_matrix, linkage_matrix, den_detail
    return [True, 10, genotype_dict]