from ATARVA.snp_utils import haplocluster_reads
from ATARVA.vcf_writer import *
from ATARVA.consensus import *
import numpy as np
from scipy.stats import nbinom
from sklearn.cluster import KMeans
import warnings
from threadpoolctl import threadpool_limits
from ATARVA.decomp_utils import motif_decomposition


def confidence_interval(data):
    data = np.array(data)
    mean = np.mean(data)
    var = np.var(data)
    if var == 0:
        return [int(mean), int(mean)]
    n = (mean**2)/((var**2) - mean) # number of success
    p = mean/(var**2) # prob of single success
    ci = nbinom.interval(0.95, n, p, loc=0)

    if np.any(np.isnan(ci)):
        return [int(mean), int(mean)]
    else:
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
    
def length_genotyper(hallele_counter, global_loci_info, global_loci_variations, locus_key, read_indices, contig, locus_start, locus_end, ref, out, male, log_bool, decomp, read_seqs, amplicon):

    read_indices = sorted(read_indices)
    locus_read_allele = global_loci_variations[locus_key]['read_allele']
    unique_alen = list(hallele_counter.keys())
    motif_size = int(float(global_loci_info[locus_key][4])) # <= 10 # boolean for motif-decomp check
    
    # if not amplicon:
    #     alen_with_1read = [item[0] for item in hallele_counter.items() if item[1]==1] # allele with 1 read contribution
    #     # if more than 10% of the reads support, empty the alen_with_1read list
    #     # if (len(alen_with_1read) / len(read_indices)) >= 0.15:
    #     #     alen_with_1read = []
    # else:
    #     alen_with_1read = []
    alen_with_1read = [item[0] for item in hallele_counter.items() if item[1]==1] # allele with 1 read contribution
    alen_with_gread = set(unique_alen) - set(alen_with_1read) # allele with more than 1 read contribution
    main_read_id = []
    alen_data = []
    
    for id in read_indices:
        if locus_read_allele[id][0] in alen_with_1read: # checking if the '1 read - allele' is nearby any of other 'good read - allele'
            num = locus_read_allele[id][0]
            for i in set(unique_alen): #for i in alen_with_gread:
                if i == num: continue
                window = round(0.1*i)
                if (i-window) <= num <= (i+window): # '1 read - allele' is considered if other allele are within 10% on either of the side
                    alen_data.append(num)
                    main_read_id.append(id)
                    break
        else:
            alen_data.append(locus_read_allele[id][0])
            main_read_id.append(id)

    if len(alen_data) < 3:
        return [False, 6]

    data = np.array(alen_data)
    data = data.reshape(-1, 1)
    with threadpool_limits(limits=1):
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=UserWarning)
            kmeans = KMeans(n_clusters=2, init='k-means++', n_init=5, random_state=0).fit(data)
    cluster_labels = kmeans.labels_  
    c1 = [i for i, x in enumerate(cluster_labels) if x == 0]
    c2 = [i for i, x in enumerate(cluster_labels) if x == 1]

    alen_c1 = [alen_data[i] for i in c1]
    alen_c2 = [alen_data[i] for i in c2]

    haplotypes = ([main_read_id[idx] for idx in c1], [main_read_id[idx] for idx in c2])
    cutoff = 0.15*len(alen_data) # 15%

    br = False
    if c1 and c2:
        def process_conditions(alen_x, alen_y):
            nonlocal br, cutoff
            max_val = max(alen_y)
            slide = max(max_val*0.1, 10)
            min_bound = min(alen_y)-slide
            max_bound = max_val+slide
            for min_al in alen_x:
                if min_bound <= min_al <= max_bound:
                    br = True
                    break

            if not br:
                cutoff = max(0.05, len(alen_x) / len(alen_data)) * len(alen_data) # min 5 % of total reads should be in the cluster
                if amplicon:
                    cutoff = max(5, cutoff) # min 5 reads should be in the cluster if amplicon

        if len(c1) < cutoff and len(c2) >= cutoff:
            process_conditions(alen_c1, alen_c2)
                               
        elif len(c2) < cutoff and len(c1) >= cutoff:
            process_conditions(alen_c2, alen_c1)


    if male:
        cluster_len = [len(c1), len(c2)]
        cidx = cluster_len.index(max( cluster_len ))
        if cluster_len[cidx]>=cutoff:
            mac = haplotypes[cidx]
            mal = alen_c1 if cidx==0 else alen_c2 # major allele length cluster
            lower,upper = confidence_interval(mal)
            allele_range = f'{lower}-{upper}'
            ALT, allele_length, decomp_seq, repeativity = alt_sequence(read_seqs, mac, amplicon, motif_size)
            if repeativity:
                vcf_homozygous_writer(ref, contig, locus_key, global_loci_info, allele_length, global_loci_variations, len(mac), out, ALT, log_bool, 'kmeans', decomp, hallele_counter, True, allele_range, decomp_seq)
            else:
                return [False, 6]
    
    elif (c1!=[] and len(c1)>=cutoff) and (c2!=[] and len(c2)>=cutoff):
        phased_read = ['.','.']
        chosen_snpQ = '.'
        snp_num = '.'        

        genotypes = []
        allele_count = {}
        ALT_seqs = []
        repeativity_list = []
        decomp_seq_list = []
        for hap_reads in haplotypes:
            ALT, allele_length, decomp_seq, repeativity = alt_sequence(read_seqs, hap_reads, amplicon, motif_size)
            repeativity_list.append(repeativity)
            decomp_seq_list.append(decomp_seq)
            ALT_seqs.append(ALT)
            genotypes.append(allele_length)
            if allele_length not in allele_count:
                allele_count[allele_length] = len(hap_reads)
            else:
                allele_count[str(allele_length)] = len(hap_reads)

        lower1,upper1 = confidence_interval(alen_c1)
        lower2,upper2 = confidence_interval(alen_c2)
        allele_range = f'{lower1}-{upper1},{lower2}-{upper2}'

        if all(repeativity_list):
            vcf_heterozygous_writer(contig, genotypes, locus_start, global_loci_variations, locus_end, allele_count, len(read_indices), global_loci_info, ref, out, chosen_snpQ, phased_read, snp_num, ALT_seqs, log_bool, 'kmeans', decomp, hallele_counter, allele_range, decomp_seq_list)
        elif any(repeativity_list):
            if repeativity_list[0]:
                allele_range = f'{lower1}-{upper1},{lower1}-{upper1}'
                vcf_homozygous_writer(ref, contig, locus_key, global_loci_info, genotypes[0], global_loci_variations, len(haplotypes[0]), out, ALT_seqs[0], log_bool, 'kmeans', decomp, hallele_counter, False, allele_range, decomp_seq_list[0])
            else:
                allele_range = f'{lower2}-{upper2},{lower2}-{upper2}'
                vcf_homozygous_writer(ref, contig, locus_key, global_loci_info, genotypes[1], global_loci_variations, len(haplotypes[1]), out, ALT_seqs[1], log_bool, 'kmeans', decomp, hallele_counter, False, allele_range, decomp_seq_list[1])
        else:
            return [False, 6]

    elif c1!=[] and len(c1)>=cutoff:
        lower1,upper1 = confidence_interval(alen_c1)
        allele_range = f'{lower1}-{upper1},{lower1}-{upper1}'
        ALT, allele_length, decomp_seq, repeativity = alt_sequence(read_seqs, haplotypes[0], amplicon, motif_size)
        if repeativity:
            vcf_homozygous_writer(ref, contig, locus_key, global_loci_info, allele_length, global_loci_variations, len(haplotypes[0]), out, ALT, log_bool, 'kmeans', decomp, hallele_counter, False, allele_range, decomp_seq)
        else:
            return [False, 6]

    elif c2!=[] and len(c2)>=cutoff:
        lower2,upper2 = confidence_interval(alen_c2)
        allele_range = f'{lower2}-{upper2},{lower2}-{upper2}'
        ALT, allele_length, decomp_seq, repeativity = alt_sequence(read_seqs, haplotypes[1] , amplicon, motif_size)
        if repeativity:
            vcf_homozygous_writer(ref, contig, locus_key, global_loci_info, allele_length, global_loci_variations, len(haplotypes[1]), out, ALT, log_bool, 'kmeans', decomp, hallele_counter, False, allele_range, decomp_seq)
        else:
            return [False, 6]
        
    else:
        return [False, 6] # write allele distribution with only one read supporting to it in vcf
    
    return [True, 10]


def analyse_genotype(contig, locus_key, global_loci_info,
                     global_loci_variations, global_read_variations, global_snp_positions, hallele_counter,
                     ref, out, sorted_global_snp_list, snpQ, snpC, snpD, snpR, phasingR, maxR, max_limit, male, log_bool, decomp, amplicon):
            
    locus_start = int(global_loci_info[locus_key][1])
    locus_end = int(global_loci_info[locus_key][2])
    motif_size = int(float(global_loci_info[locus_key][4]))

    state = False

    if max_limit == 0:
        read_indices = global_loci_variations[locus_key]['reads']
    else:
        read_indices = global_loci_variations[locus_key]['reads'][:maxR]

    read_seqs = global_loci_variations[locus_key]['read_sequence']

    if male or amplicon: 
        state, skip_point = length_genotyper(hallele_counter, global_loci_info, global_loci_variations, locus_key, read_indices, contig, locus_start, locus_end, ref, out, male, log_bool, decomp, read_seqs, amplicon)
        return [state, skip_point]


    snp_positions = set()
    for rindex in read_indices:
        snp_positions |= (global_read_variations[rindex]['snps'])

    snp_positions = sorted(list(filter(lambda x: (x in global_snp_positions) and (global_snp_positions[x]['cov'] >= 3) and
                                                    (locus_start - snpD < x < locus_end + snpD),
                            snp_positions)))


    snp_allelereads = {}
    read_indices = set(read_indices)
    non_ref_snp_cov = {}
    for pos in snp_positions:
        c_point=0
        coverage = set()
        non_ref_nucs = [nucleotides for nucleotides in global_snp_positions[pos] if nucleotides not in ['cov', 'Qval', 'r']]
        for each_nuc in non_ref_nucs:
            reads_of_nuc = global_snp_positions[pos][each_nuc].intersection(read_indices)
            if len(reads_of_nuc) == 0: continue
            coverage.add(len(reads_of_nuc))

            if (sum([global_snp_positions[pos]['Qval'][read_idx] for read_idx in reads_of_nuc])/len(reads_of_nuc)) <= 13:
                c_point=1
                break
        if (len(coverage)==0) or (c_point==1): continue
        else: non_ref_snp_cov[pos] = max(coverage)
            
        snp_allelereads[pos] = { 'cov': 0, 'reads': set(), 'alleles': {}, 'Qval': {} }
        for nuc in global_snp_positions[pos]:
            if (nuc == 'cov') or (nuc == 'Qval'): continue
            snp_allelereads[pos]['alleles'][nuc] = global_snp_positions[pos][nuc].intersection(read_indices)
            snp_allelereads[pos]['cov'] += len(snp_allelereads[pos]['alleles'][nuc])
            if nuc!='r':
                snp_allelereads[pos]['Qval'].update(dict([(read_idx,global_snp_positions[pos]['Qval'][read_idx]) for read_idx in snp_allelereads[pos]['alleles'][nuc]]))

    del_positions = list(filter(lambda x: snp_allelereads[x]['cov'] < 5, snp_allelereads.keys()))

    for pos in del_positions:
        del snp_allelereads[pos]


    ordered_snp_on_cov = sorted(snp_allelereads.keys(), key = lambda item : non_ref_snp_cov[item], reverse = True)


    haplotypes, min_snp, skip_point, chosen_snpQ, phased_read, snp_num = haplocluster_reads(snp_allelereads, ordered_snp_on_cov, read_indices, snpQ, snpC, snpR, phasingR) # SNP ifo and supporting reads for specific locus are given to the phasing function

    if haplotypes == (): # if the loci has no significant snps
        state, skip_point = length_genotyper(hallele_counter, global_loci_info, global_loci_variations, locus_key, read_indices, contig, locus_start, locus_end, ref, out, male, log_bool, decomp, read_seqs, False)
        del read_seqs
        return [state, skip_point]
    
    if min_snp != -1:
        min_idx = sorted_global_snp_list.index(min_snp)
        del sorted_global_snp_list[:min_idx]
        del_snps = set()
        for pos in global_snp_positions:
            if pos < min_snp: del_snps.add(pos)
        for pos in del_snps:
            del global_snp_positions[pos]


    genotypes = []
    allele_count = {}
    ALT_seqs = []
    alen_list = []
    for hap_reads in haplotypes:
        ALT, allele_length,_,_ = alt_sequence(read_seqs, hap_reads, False, motif_size)
        alen_list.append([len(read_seqs[read_id][0]) for read_id in hap_reads])
        ALT_seqs.append(ALT)
        genotypes.append(allele_length)
        if allele_length not in allele_count:
            allele_count[allele_length] = len(hap_reads)
        else:
            allele_count[str(allele_length)] = len(hap_reads)

    del read_seqs
    lower1,upper1 = confidence_interval(alen_list[0])
    lower2,upper2 = confidence_interval(alen_list[1])
    allele_range = f'{lower1}-{upper1},{lower2}-{upper2}'
    vcf_heterozygous_writer(contig, genotypes, locus_start, global_loci_variations, locus_end, allele_count, len(read_indices), global_loci_info, ref, out, chosen_snpQ, phased_read, snp_num, ALT_seqs, log_bool, 'SNP', decomp, hallele_counter, allele_range, [None])
    state = True
    return [state, skip_point]
    