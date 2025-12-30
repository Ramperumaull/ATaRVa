import sys
import pysam
# from ATARVA.consensus import consensus_seq_poa
from ATARVA.decomp_utils import motif_decomposition

info_opt_tag = 'ID'
def set_info_opt_tag(tag):
    global info_opt_tag
    info_opt_tag = tag

def vcf_writer(out, bam, bam_name):

    vcf_header = pysam.VariantHeader()

    # command
    vcf_header.add_line(f"##command=ATaRVa_0.5.0 {' '.join(sys.argv)}")

    for contig in bam.header['SQ']:
        vcf_header.contigs.add(contig['SN'], length=contig['LN'])
    #sample_name
    vcf_header.add_sample(bam_name)
    # FILTER
    vcf_header.filters.add('LESS_READS', number=None, type=None, description="Read depth below threshold")
    # INFO
    vcf_header.info.add("AC", number='A', type="Integer", description="Number of alternate alleles in called genotypes")
    vcf_header.info.add("AN", number=1, type="Integer", description="Number of alleles in called genotypes")
    vcf_header.info.add("MOTIF", number=1, type="String", description="Repeat motif")
    vcf_header.info.add("START", number=1, type="Integer", description="Start position of the repeat region in 0-based coordinate system")
    vcf_header.info.add("END", number=1, type="Integer", description="End position of the repeat region")
    vcf_header.info.add("CT", number=1, type="String", description="Cluster type")
    vcf_header.info.add("EAC", number=1, type="String", description="Each Allele Count")
    # FORMAT
    vcf_header.formats.add("GT", number=1, type="String", description="Genotype")
    vcf_header.formats.add("AL", number=2, type="Integer", description="Allele length in base pairs")
    vcf_header.formats.add("AR", number='.', type="String", description="Allele length range")
    vcf_header.formats.add("SD", number='.', type="Integer", description="Number of reads supporting for the alleles")
    vcf_header.formats.add("PC", number=2, type="Integer", description="Number of reads in the phased cluster for each allele")
    vcf_header.formats.add("DP", number=1, type="Integer", description="Number of the supporting reads for the repeat locus")
    vcf_header.formats.add("SN", number='.', type="Integer", description="Number of SNPs used for phasing")
    vcf_header.formats.add("SQ", number='.', type="Float", description="Phred-scale qualities of the SNPs used for phasing")
    vcf_header.formats.add("MM", number='.', type="Float", description="Mean methylation level for each allele")
    vcf_header.formats.add("MR", number='.', type="Integer", description="Number of reads providing methylation info for each allele")
    vcf_header.formats.add("DS", number='A', type="String", description="Motif decomposed sequence")

    out.write(str(vcf_header))

def vcf_homozygous_writer(ref, contig, locus_key, global_loci_info, homozygous_allele, reads_len, DP, out, ALT_read, log_bool, tag, decomp, hallele_counter, haploid_state, allele_range, decomp_seq, meth_info):

    locus_start = int(global_loci_info[locus_key][1])
    locus_end = int(global_loci_info[locus_key][2])
    
    if len(global_loci_info[locus_key]) > 5:
        optional_tag = f';{info_opt_tag}={global_loci_info[locus_key][5]}'
    else:
        optional_tag = ''

    # if type(reads_len) == list:
    #     reads_len = len(reads_len) #removable

    meth_prob = meth_info[0] #methylation probability
    meth_reads = str(meth_info[1]) if meth_info[1] is not None else '.' #number of methylated reads
    meth_prob = [str(meth_prob) if meth_prob is not None else '.']*2 # for homozygous, make it two same values to keep the format consistent
    
    ref_allele_length = locus_end - locus_start
    # DP = len(global_loci_variations[locus_key]['reads'])

    AC = 0; AN = 2; GT = '0/0'; ALT = '.'; alt_state = False
    MM = ','.join(meth_prob)
    if homozygous_allele != ref_allele_length:
        AC = 2
        GT = '1/1'

        ALT = ALT_read
        if ALT[0]!='<': alt_state = True


    if log_bool:
        eac = sorted(hallele_counter.items(), key = lambda x: x[1], reverse=True)
        INFO = 'AC=' + str(AC) + ';AN=' + str(AN) + ';MOTIF=' + str(global_loci_info[locus_key][3]) + ';START=' + str(locus_start) + ';END=' + str(locus_end) + optional_tag + ';CT=' + tag + ';EAC=' + str(eac)
    else:
        INFO = 'AC=' + str(AC) + ';AN=' + str(AN) + ';MOTIF=' + str(global_loci_info[locus_key][3]) + ';START=' + str(locus_start) + ';END=' + str(locus_end) + optional_tag

    if decomp:
        motif_size = int(float(global_loci_info[locus_key][4]))

        FORMAT = 'GT:AL:AR:SD:DP:SN:SQ:MM:MR:DS'
        if alt_state & (motif_size<=10):
            if decomp_seq:
                deseq = decomp_seq
            else:
                deseq,_ = motif_decomposition(ALT, motif_size)
        else:
            deseq = '.'

        if not haploid_state:
            SAMPLE = str(GT) + ':' + str(homozygous_allele) + ',' + str(homozygous_allele) + ':' + allele_range + ':' + str(reads_len) + ':' + str(DP) + ':.:.' + ':' + MM + ':' + meth_reads + ':' + deseq
        else:
            SAMPLE = GT[0] + ':' + str(homozygous_allele) + ':' + allele_range + ':' + str(reads_len) + ':' + str(DP) + ':.:.' + ':' + meth_prob[0] + ':' + meth_reads + ':' + deseq
    else:
        # FORMAT = 'GT:AL:SD:PC:DP:SN:SQ'
        FORMAT = 'GT:AL:AR:SD:DP:SN:SQ:MM:MR'
        # SAMPLE = str(GT) + ':' + str(homozygous_allele) + ',' + str(homozygous_allele) + ':' + str(reads_len) + ':.:' + str(DP) + ':.:.'
        if not haploid_state:
            SAMPLE = str(GT) + ':' + str(homozygous_allele) + ',' + str(homozygous_allele) + ':' + allele_range + ':' + str(reads_len) + ':' + str(DP) + ':.:.' + ':' + MM + ':' + meth_reads
        else:
            SAMPLE = GT[0] + ':' + str(homozygous_allele) + ':' + allele_range + ':' + str(reads_len) + ':' + str(DP) + ':.:.' + ':' + meth_prob[0] + ':' + meth_reads

    
    print(*[contig, locus_start+1, '.',  ref.fetch(contig, locus_start, locus_end), ALT , 0, 'PASS', INFO, FORMAT, SAMPLE], file=out, sep='\t')
    del ALT_read
    del global_loci_info[locus_key]


def vcf_heterozygous_writer(contig, genotypes, locus_start, locus_end, allele_count, DP, global_loci_info, ref, out, chosen_snpQ, phased_read, snp_num, ALT_reads, log_bool, tag, decomp, hallele_counter, allele_range, decomp_seq, meth_info):

    locus_key = f'{contig}:{locus_start}-{locus_end}'

    if len(global_loci_info[locus_key]) > 5:
        optional_tag = f';{info_opt_tag}={global_loci_info[locus_key][5]}'
    else:
        optional_tag = ''

    final_allele = set(genotypes)
    heterozygous_allele = ''
    AC = 'AC'
    AN = 2
    GT = 'GT'
    SD = 'SD'
    PC = 'PC'
    ALT = '.'
    alt_seqs = []

    ref_allele_length = locus_end - locus_start

    meth_prob = []
    meth_reads = []
    for each_meth in meth_info:
        meth_prob.append(str(each_meth[0]) if each_meth[0] is not None else '.') #methylation probability
        meth_reads.append(str(each_meth[1]) if each_meth[1] is not None else '.') #number of methylated reads

    # meth_prob = [str(i) for i in meth_prob]
    # meth_reads = [str(i) for i in meth_reads]

    if len(final_allele) == 1:

        if ref_allele_length == tuple(final_allele)[0]:
            AC = 0
            GT = '0|0'
            heterozygous_allele+=str(ref_allele_length)+','+str(ref_allele_length)
            SD = str(allele_count[ref_allele_length])+','+str(allele_count[str(ref_allele_length)])
            alt_seqs.append('')
        else:
            AC = 2; GT = '1|1'
            heterozygous_allele+=str(tuple(final_allele)[0])+','+str(tuple(final_allele)[0])
            SD = str(allele_count[tuple(final_allele)[0]])+','+str(allele_count[str(tuple(final_allele)[0])])

            ALT = ALT_reads[0]
            if ALT[0]!='<': alt_seqs.append(ALT)
            else: alt_seqs.append('')
        PC = str(phased_read[0])+','+str(phased_read[1])
        MM = ','.join(meth_prob)
        MR = ','.join(meth_reads)
    else:

        if len(set((ref_allele_length,)) & final_allele) == 1:
            AC = 1
            GT = '0|1'
            heterozygous_allele+=str(ref_allele_length)+','+str(tuple(final_allele-{ref_allele_length})[0])
            SD = str(allele_count[ref_allele_length])+','+str(allele_count[tuple(final_allele-{ref_allele_length})[0]])
            if genotypes.index(ref_allele_length) == 0:
                PC = str(phased_read[0])+','+str(phased_read[1])

                alt_seqs.append(None) # dummy added for ref, to keep the length of alt_seqs as 2
                ALT = ALT_reads[1]
                if ALT[0]!='<': alt_seqs.append(ALT)
                else: alt_seqs.append('')
                MM = ','.join(meth_prob)
                MR = ','.join(meth_reads)
            else:
                PC = str(phased_read[1])+','+str(phased_read[0])

                ALT = ALT_reads[0]
                if ALT[0]!='<': alt_seqs.append(ALT)
                else: alt_seqs.append('')
                alt_seqs.append(None) # dummy added for ref, to keep the length of alt_seqs as 2
                allele_range = ','.join(allele_range.split(',')[::-1]) # reverse the allele range to keep the order consistent with GT
                MM = ','.join(meth_prob[::-1]) # reverse the meth_prob to keep the order consistent with GT
                MR = ','.join(meth_reads[::-1])
        else:
            AC = '1,1'
            GT = '1|2'
            heterozygous_allele+=str(genotypes[0])+','+str(genotypes[1])
            SD = str(allele_count[genotypes[0]])+','+str(allele_count[genotypes[1]])
            PC = str(phased_read[0])+','+str(phased_read[1])

            ALT1 = ALT_reads[0]
            if ALT1[0]!='<': alt_seqs.append(ALT1)
            else: alt_seqs.append('')
                
            ALT2 = ALT_reads[1]
            if ALT2[0]!='<': alt_seqs.append(ALT2)
            else: alt_seqs.append('')

            ALT = ALT1 + ',' + ALT2
            MM = ','.join(meth_prob)
            MR = ','.join(meth_reads)


    if PC == '.,.': PC = '.' # due to length genotyper
    if log_bool:
        eac = sorted(hallele_counter.items(), key = lambda x: x[1], reverse=True)
        INFO = 'AC='+str(AC)+';AN='+str(AN)+';MOTIF=' + str(global_loci_info[locus_key][3]) + ';START=' + str(locus_start) + ';END='+str(locus_end) + optional_tag + ';CT=' + tag + ';EAC=' + str(eac)
    else:
        INFO = 'AC='+str(AC)+';AN='+str(AN)+';MOTIF=' + str(global_loci_info[locus_key][3]) + ';START=' + str(locus_start) + ';END='+str(locus_end) + optional_tag

    if decomp:
        motif_size = int(float(global_loci_info[locus_key][4]))
        # FORMAT = 'GT:AL:SD:PC:DP:SN:SQ:DS'
        FORMAT = 'GT:AL:AR:SD:DP:SN:SQ:MM:MR:DS'
        if motif_size>10:
            deseq = ','.join(['.']*len(alt_seqs))
        else:
            ds = []
            for index,iseq in enumerate(alt_seqs):
                if iseq:
                    if all(decomp_seq):
                        ds.append(decomp_seq[index])
                    else:
                        i_deseq,_ = motif_decomposition(iseq, motif_size)
                        ds.append(i_deseq)
                elif iseq=='':
                    ds.append('.')
            deseq = ','.join(ds)
        # SAMPLE = str(GT)+':'+heterozygous_allele+':' + SD + ':' + PC + ':' + str(DP) + ':' + str(snp_num) + ':' + chosen_snpQ + ':' + deseq
        SAMPLE = str(GT)+':'+heterozygous_allele+':' + allele_range + ':' + SD + ':' + str(DP) + ':' + str(snp_num) + ':' + chosen_snpQ + ':' + MM + ':' + MR + ':' + deseq
    else: 
        # FORMAT = 'GT:AL:SD:PC:DP:SN:SQ'
        FORMAT = 'GT:AL:AR:SD:DP:SN:SQ:MM:MR'
        # SAMPLE = str(GT)+':'+heterozygous_allele+':' + SD + ':' + PC + ':' + str(DP) + ':' + str(snp_num) + ':' + chosen_snpQ
        SAMPLE = str(GT)+':'+heterozygous_allele+':' + allele_range + ':' + SD + ':' + str(DP) + ':' + str(snp_num) + ':' + chosen_snpQ + ':' + MM + ':' + MR

    del ALT_reads
    del alt_seqs

    print(*[contig, locus_start+1, '.',  ref.fetch(contig, locus_start, locus_end), ALT, 0, 'PASS', INFO, FORMAT, SAMPLE], file=out, sep='\t')
    del global_loci_info[locus_key]

def vcf_fail_writer(contig, locus_key, global_loci_info, ref, out, DP, skip_point):

    locus_start = int(global_loci_info[locus_key][1])
    locus_end = int(global_loci_info[locus_key][2])

    if len(global_loci_info[locus_key]) > 5:
        optional_tag = f';{info_opt_tag}={global_loci_info[locus_key][5]}'
    else:
        optional_tag = ''

    if skip_point == 0:
        FILTER = 'LESS_READS'
          
    locus_key = f'{contig}:{locus_start}-{locus_end}'

    INFO = 'AC=0;AN=0;MOTIF=' + str(global_loci_info[locus_key][3]) + ';START=' + str(locus_start) + ';END=' + str(locus_end) + optional_tag
    FORMAT = 'GT:AL:AR:SD:DP:SN:SQ:MM:MR'
    SAMPLE = '.:.:.:.:.:.:.:.:.'

    print(*[contig, locus_start+1, '.',  ref.fetch(contig, locus_start, locus_end), '.', 0, FILTER, INFO, FORMAT, SAMPLE], file=out, sep='\t')
    del global_loci_info[locus_key]

def vcf_multizygous_writer(contig, genotype_dict, locus_start, locus_end, DP, global_loci_info, ref, out, log_bool, decomp, hallele_counter):

    locus_key = f'{contig}:{locus_start}-{locus_end}'

    tag = "correlation_clustering"

    if len(global_loci_info[locus_key]) > 5:
        optional_tag = f';{info_opt_tag}={global_loci_info[locus_key][5]}'
    else:
        optional_tag = ''

    GT_dict = {}
    gt_idx = 0
    ref_allele_length = locus_end - locus_start
    ref_seq = ref.fetch(contig, locus_start, locus_end)
    for each_genotype in genotype_dict:
        current_gt = genotype_dict[each_genotype]
        if int(each_genotype) == ref_allele_length:
            if ref_seq == current_gt[0]:
                GT_dict[0] = (current_gt[0], str(each_genotype), current_gt[3], f'{current_gt[1][0]}-{current_gt[1][1]}', current_gt[4][0], current_gt[4][1], current_gt[2])
            else:
                gt_idx += 1
                GT_dict[gt_idx] = (current_gt[0], str(each_genotype), current_gt[3], f'{current_gt[1][0]}-{current_gt[1][1]}', current_gt[4][0], current_gt[4][1], current_gt[2])
        else:
            gt_idx += 1
            GT_dict[gt_idx] = (current_gt[0], str(each_genotype), current_gt[3], f'{current_gt[1][0]}-{current_gt[1][1]}', current_gt[4][0], current_gt[4][1], current_gt[2])
    del genotype_dict
    
    GT = []
    if gt_idx> 0:
        AN = (gt_idx + 1) if (0 in GT_dict) else gt_idx
    else:
        AN = 2
    # AN = (gt_idx + 1) if gt_idx>0 else 2
    AC = ','.join(['1']*gt_idx) if gt_idx>0 else 0
    ALT = []
    AL = []
    AR = []
    SD = []
    MM = []
    MR = []
    deseq = []
    for gt_key in sorted(GT_dict.keys()):
        GT.append(str(gt_key))
        if gt_key != 0:
            ALT.append(GT_dict[gt_key][0])
            deseq.append(GT_dict[gt_key][6] if GT_dict[gt_key][6] else '.')
        AL.append(GT_dict[gt_key][1])
        SD.append(str(GT_dict[gt_key][2]))
        AR.append(GT_dict[gt_key][3])
        MM.append(str(GT_dict[gt_key][4]) if GT_dict[gt_key][4] is not None else '.')
        MR.append(str(GT_dict[gt_key][5]) if GT_dict[gt_key][5] is not None else '.')
    del GT_dict

    GT = '/'.join(GT)
    ALT = ','.join(ALT) if ALT else '.'
    AL = ','.join(AL)
    AR = ','.join(AR)
    SD = ','.join(SD)
    MM = ','.join(MM)
    MR = ','.join(MR)
        
    if log_bool:
        eac = sorted(hallele_counter.items(), key = lambda x: x[1], reverse=True)
        INFO = 'AC='+str(AC)+';AN='+str(AN)+';MOTIF=' + str(global_loci_info[locus_key][3]) + ';START=' + str(locus_start) + ';END='+str(locus_end) + optional_tag + ';CT=' + tag + ';EAC=' + str(eac)
    else:
        INFO = 'AC='+str(AC)+';AN='+str(AN)+';MOTIF=' + str(global_loci_info[locus_key][3]) + ';START=' + str(locus_start) + ';END='+str(locus_end) + optional_tag

    if decomp:
        FORMAT = 'GT:AL:AR:SD:DP:SN:SQ:MM:MR:DS'
        deseq = ','.join(deseq) if deseq else '.'
        SAMPLE = GT + ':' + AL + ':' + AR + ':' + SD + ':' + str(DP) + ':.:.:' + MM + ':' + MR + ':' + deseq
    else:
        FORMAT = 'GT:AL:AR:SD:DP:SN:SQ:MM:MR'
        SAMPLE = GT + ':' + AL + ':' + AR + ':' + SD + ':' + str(DP) + ':.:.:' + MM + ':' + MR

    print(*[contig, locus_start+1, '.',  ref_seq, ALT, 0, 'PASS', INFO, FORMAT, SAMPLE], file=out, sep='\t')

    del GT, ALT, AL, AR, SD, MM, MR, deseq, global_loci_info[locus_key]
