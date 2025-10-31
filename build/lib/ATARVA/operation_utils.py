from ATARVA.soft_clip_utils import detect_flank


def convert_eqx_read(chrom, read_start, cigar_tuples, read_seq, ref):
    new_read_seq = ''
    qpos = 0
    for cigar in cigar_tuples:
        if cigar[0] == 4:
            new_read_seq += read_seq[qpos:qpos+cigar[1]]
            qpos += cigar[1]
        elif cigar[0] == 1:
            new_read_seq += read_seq[qpos:qpos+cigar[1]]
            qpos += cigar[1]
        elif cigar[0] == 2:
            read_start += cigar[1]
        elif cigar[0] == 0:
            ref_end = read_start + cigar[1]
            inter_ref = ref.fetch(chrom, read_start, ref_end)
            current_segment = read_seq[qpos:qpos+cigar[1]]
            read_segment = ''
            for idx,base in enumerate(current_segment):
                if base!='=':
                    read_segment += base
                else:
                    read_segment += inter_ref[idx]
            new_read_seq += read_segment
            read_start = ref_end
            qpos += cigar[1]
        elif cigar[0] == 7:
            ref_end = read_start + cigar[1]
            inter_ref = ref.fetch(chrom, read_start, ref_end)
            new_read_seq += inter_ref
            read_start = ref_end
            qpos += cigar[1]
        elif cigar[0] == 8:
            new_read_seq += read_seq[qpos:qpos+cigar[1]]
            read_start += cigar[1]
            qpos += cigar[1]
            
    return new_read_seq

def update_homopolymer_coords(ref_seq, locus_start, homopoly_positions):
    """
    Record all the homopolymer stretches of at least 3 bases within the repeat coordinates
    """
    prev_N = ref_seq[0]; start = -1
    for i, n in enumerate(ref_seq[1:]):
        if n == prev_N:
            if start == -1: start = i
        else:
            if start != -1 and (i+1)-start >= 4:
                for l,c in enumerate(range(locus_start+start, locus_start+i+1)):
                    homopoly_positions[c] = (i-start+1)-l
            start = -1
        prev_N = n

    if start != -1 and (i+1)-start >= 4:
        for l,c in enumerate(range(locus_start+start, locus_start+i+1)):
            # for each position in the homopolymer stretch we record the length of the 
            # homopolymer nucleotides on the right
            homopoly_positions[c] = (i-start+1)-l


def match_jump(rpos, repeat_index, loci_keys, loci_coords, tracked, locus_qpos_range, qpos, match_len, loci_flank_qpos_range, flank_track, left_flank, right_flank, amp_right_flank_list, amp_left_flank_list, out_insertion_qpos_ranges_right, out_insertion_qpos_ranges_left, right_ins_rpos, left_ins_rpos, amplicon_variables):
    """
    Return the number of repeat indices to jump when scanning through a match segment
    """

    if amp_left_flank_list:
        chrom, ref, query_sequence, flank_length, qpos_start, qpos_end = amplicon_variables

    previous_rpos = rpos - match_len
    r = 0 
    for r,coord in enumerate(loci_coords[repeat_index:]):
        coord_start, coord_end = coord
        
        if rpos < coord_start: break
        
        if previous_rpos > coord_end: continue
            
        locus_key = loci_keys[r+repeat_index]
        if not tracked[r+repeat_index]:

            if coord_start <= rpos:
                
                locus_qpos_range[r+repeat_index][0] = qpos - (rpos - coord_start)
                if amp_left_flank_list:
                    lstart, lend = ref_repeat(locus_key)
                    Record_left_out_ins(amp_left_flank_list, r, repeat_index, chrom, ref, query_sequence, flank_length, lstart, lend, locus_qpos_range, loci_flank_qpos_range, out_insertion_qpos_ranges_left, left_ins_rpos, flank_track, qpos_start, qpos_end)

            if coord_end <= rpos:
                
                locus_qpos_range[r+repeat_index][1] = qpos - (rpos - coord_end)
                if amp_left_flank_list:
                    lstart, lend = ref_repeat(locus_key)
                    Record_right_out_ins(amp_right_flank_list, r, repeat_index, chrom, ref, query_sequence, flank_length, lstart, lend, locus_qpos_range, loci_flank_qpos_range, out_insertion_qpos_ranges_right, right_ins_rpos, flank_track, qpos_start, qpos_end)

            tracked[r+repeat_index] = True 

            # for storing repeat qpos ranges
            if (not flank_track[r+repeat_index][0]) and (coord_start+left_flank[r+repeat_index] <= rpos):
                loci_flank_qpos_range[r+repeat_index][0] = qpos - (rpos - coord_start)+left_flank[r+repeat_index]
                flank_track[r+repeat_index][0] = True
            if (not flank_track[r+repeat_index][1]) and (coord_end-right_flank[r+repeat_index] <= rpos):
                loci_flank_qpos_range[r+repeat_index][1] = qpos - (rpos - coord_end)-right_flank[r+repeat_index]
                if rpos > coord_end-right_flank[r+repeat_index]: flank_track[r+repeat_index][1] = True
                
                

        elif coord_end <= rpos:
            
            locus_qpos_range[r+repeat_index][1] = qpos - (rpos -coord_end)
            if amp_left_flank_list:
                lstart, lend = ref_repeat(locus_key)
                Record_right_out_ins(amp_right_flank_list, r, repeat_index, chrom, ref, query_sequence, flank_length, lstart, lend, locus_qpos_range, loci_flank_qpos_range, out_insertion_qpos_ranges_right, right_ins_rpos, flank_track, qpos_start, qpos_end)

        # for storing repeat qpos ranges
        if not flank_track[r+repeat_index][0]:
            if coord_start+left_flank[r+repeat_index] <= rpos:
                loci_flank_qpos_range[r+repeat_index][0] = qpos - (rpos -coord_start)+left_flank[r+repeat_index]
                flank_track[r+repeat_index][0] = True 
            if coord_end-right_flank[r+repeat_index] <= rpos:
                loci_flank_qpos_range[r+repeat_index][1] = qpos - (rpos - coord_end)-right_flank[r+repeat_index]
                if rpos > coord_end-right_flank[r+repeat_index]: flank_track[r+repeat_index][1] = True
            
        elif (not flank_track[r+repeat_index][1]) and (coord_end-right_flank[r+repeat_index] <= rpos):
            loci_flank_qpos_range[r+repeat_index][1] = qpos - (rpos -coord_end)-right_flank[r+repeat_index]
            if rpos > coord_end-right_flank[r+repeat_index]: flank_track[r+repeat_index][1] = True

    
    jump = 0    # jump beyond the repeat where all positions are tracked
    if loci_coords[repeat_index + r - 1][1] < rpos:
        for f in loci_coords[repeat_index:]:
            if f[1] < rpos: jump += 1
            else: break
    
    return jump


def deletion_jump(deletion_length, rpos, repeat_index, loci_keys, tracked, loci_coords, homopoly_positions, read_loci_variations, locus_qpos_range, qpos, loci_flank_qpos_range, flank_track, left_flank, right_flank, amp_right_flank_list, amp_left_flank_list, out_insertion_qpos_ranges_right, out_insertion_qpos_ranges_left, right_ins_rpos, left_ins_rpos, amplicon_variables):
    """
    Return the number of repeat indices to jump when scanning through a deletion segment.
    The function tracks specifically if the deletion is segment has complete repeats in them
    or segments of the repeat is deleted.
    """

    if amp_left_flank_list:
        chrom, ref, query_sequence, flank_length, qpos_start, qpos_end = amplicon_variables

    # rpos - corresponds to the position in the reference after tracking the deletion
    r = 0   # required to be initialised outside the loop
    for r, coord in enumerate(loci_coords[repeat_index:]):
        coord_start, coord_end = coord
        # if rpos is before the start of the repeat; repeat is unaffected
        if rpos < coord_start: break

        # actual position in the reference where the deletion is occurring
        del_pos = rpos - deletion_length
        if del_pos > coord_end: continue

        locus_key = loci_keys[r+repeat_index]
        if not tracked[r+repeat_index]:
            # if the locus is not tracked
            # deletion is encountered beyond
            if coord_start <= rpos:    
                locus_qpos_range[r+repeat_index][0] = qpos        
                tracked[r+repeat_index] = True    # set tracked as true
                if amp_left_flank_list:
                    lstart, lend = ref_repeat(locus_key)
                    Record_left_out_ins(amp_left_flank_list, r, repeat_index, chrom, ref, query_sequence, flank_length, lstart, lend, locus_qpos_range, loci_flank_qpos_range, out_insertion_qpos_ranges_left, left_ins_rpos, flank_track, qpos_start, qpos_end)

            if coord_end < rpos:
                
                locus_qpos_range[r+repeat_index][1] = qpos
                if amp_left_flank_list:
                    lstart, lend = ref_repeat(locus_key)
                    Record_right_out_ins(amp_right_flank_list, r, repeat_index, chrom, ref, query_sequence, flank_length, lstart, lend, locus_qpos_range, loci_flank_qpos_range, out_insertion_qpos_ranges_right, right_ins_rpos, flank_track, qpos_start, qpos_end)

            # for storing repeat qpos ranges
            if (not flank_track[r+repeat_index][0]) and (coord_start+left_flank[r+repeat_index] <= rpos):
                loci_flank_qpos_range[r+repeat_index][0] = qpos
                flank_track[r+repeat_index][0] = True
            if (not flank_track[r+repeat_index][1]) and (coord_end-right_flank[r+repeat_index] < rpos):
                loci_flank_qpos_range[r+repeat_index][1] = qpos
                flank_track[r+repeat_index][1] = True

        elif coord_end < rpos:
            
            locus_qpos_range[r+repeat_index][1] = qpos
            if amp_left_flank_list:
                lstart, lend = ref_repeat(locus_key)
                Record_right_out_ins(amp_right_flank_list, r, repeat_index, chrom, ref, query_sequence, flank_length, lstart, lend, locus_qpos_range, loci_flank_qpos_range, out_insertion_qpos_ranges_right, right_ins_rpos, flank_track, qpos_start, qpos_end)

        # for storing repeat qpos ranges
        if not flank_track[r+repeat_index][0]:
            if coord_start+left_flank[r+repeat_index] <= rpos:
                loci_flank_qpos_range[r+repeat_index][0] = qpos
                flank_track[r+repeat_index][0] = True 
            if coord_end-right_flank[r+repeat_index] < rpos:
                loci_flank_qpos_range[r+repeat_index][1] = qpos
                flank_track[r+repeat_index][1] = True
        elif (not flank_track[r+repeat_index][1]) and (coord_end-right_flank[r+repeat_index] <= rpos):
            loci_flank_qpos_range[r+repeat_index][1] = qpos
            if rpos > coord_end-right_flank[r+repeat_index]: flank_track[r+repeat_index][1] = True

        # updating the allele with the deletion considered
        # read_loci_variations[locus_key][rpos] = f'D|{deletion_length}'
        
        # del_len = min(coord[1], rpos) - max(coord[0], del_pos)
        del_len = min(coord_end-right_flank[r+repeat_index], rpos) - max(coord_start+left_flank[r+repeat_index], del_pos)
        if (rpos >= coord_start+left_flank[r+repeat_index]) and (del_pos <= coord_end-right_flank[r+repeat_index]): # introduced to include length only if it comes inside repeat region
            if del_pos not in homopoly_positions:
                read_loci_variations[locus_key]['alen'] -= del_len
                read_loci_variations[locus_key]['halen'] -= del_len
            else:
                if del_len <= homopoly_positions[del_pos]:
                    # if the deletion is only limited to the homopolymer positions
                    read_loci_variations[locus_key]['halen'] -= del_len
                else:
                    read_loci_variations[locus_key]['alen'] -= del_len
                    read_loci_variations[locus_key]['halen'] -= del_len


    jump = 0    # jump beyond the repeat where all positions are tracked
    if loci_coords[repeat_index + r - 1][1] < rpos:
        for f in loci_coords[repeat_index:]:
            if f[1] < rpos: jump += 1
            else: break

    return jump


def insertion_jump(insertion_length, insert, rpos, repeat_index, loci_keys, tracked, loci_coords, homopoly_positions, read_loci_variations, locus_qpos_range, qpos, loci_flank_qpos_range, flank_track, left_flank, right_flank, out_insertion_qpos_ranges_left, out_insertion_qpos_ranges_right, left_ins_rpos, right_ins_rpos, amp_right_flank_list, amp_left_flank_list, amplicon_variables):
    """
    Return the number of repeat indices to jump when scanning through a insertion segment.
    The function tracks specifically if the deletion is segment has complete repeats in them
    or segments of the repeat is deleted.
    """
    if amp_left_flank_list:
        chrom, ref, query_sequence, flank_length, qpos_start, qpos_end = amplicon_variables

    r = 0   # required to be initialised outside the loop
    for r, coord in enumerate(loci_coords[repeat_index:]):
        coord_start, coord_end = coord
        # if rpos is before the start of the repeat; repeat is unaffected
        if rpos < coord_start: break

        # if the insertion is happening beyond, the repeat in unaffected
        if rpos > coord_end: continue

        locus_key = loci_keys[r+repeat_index]
        if not tracked[r+repeat_index]:
            # if the locus is not tracked
            # deletion is encountered beyond
            if coord_start <= rpos:
                locus_qpos_range[r+repeat_index][0] = qpos-insertion_length
                tracked[r+repeat_index] = True    # set tracked as true
                if amp_left_flank_list:
                    lstart, lend = ref_repeat(locus_key)
                    Record_left_out_ins(amp_left_flank_list, r, repeat_index, chrom, ref, query_sequence, flank_length, lstart, lend, locus_qpos_range, loci_flank_qpos_range, out_insertion_qpos_ranges_left, left_ins_rpos, flank_track, qpos_start, qpos_end)

            if coord_end == rpos:
               
                locus_qpos_range[r+repeat_index][1] = qpos
                if amp_left_flank_list:
                    lstart, lend = ref_repeat(locus_key)
                    Record_right_out_ins(amp_right_flank_list, r, repeat_index, chrom, ref, query_sequence, flank_length, lstart, lend, locus_qpos_range, loci_flank_qpos_range, out_insertion_qpos_ranges_right, right_ins_rpos, flank_track, qpos_start, qpos_end)

                # here jump can be done

            # for storing repeat qpos ranges
            if (not flank_track[r+repeat_index][0]) and (coord_start+left_flank[r+repeat_index]-1 <= rpos):
                loci_flank_qpos_range[r+repeat_index][0] = qpos-insertion_length
                flank_track[r+repeat_index][0] = True
            if (not flank_track[r+repeat_index][1]) and (coord_end-right_flank[r+repeat_index] <= rpos):
                loci_flank_qpos_range[r+repeat_index][1] = qpos
                if rpos > coord_end-right_flank[r+repeat_index]: flank_track[r+repeat_index][1] = True


        elif coord_end == rpos:
            
            locus_qpos_range[r+repeat_index][1] = qpos
            if amp_left_flank_list:
                lstart, lend = ref_repeat(locus_key)
                Record_right_out_ins(amp_right_flank_list, r, repeat_index, chrom, ref, query_sequence, flank_length, lstart, lend, locus_qpos_range, loci_flank_qpos_range, out_insertion_qpos_ranges_right, right_ins_rpos, flank_track, qpos_start, qpos_end)

        # for storing repeat qpos ranges
        if not flank_track[r+repeat_index][0]:
            if coord_start+left_flank[r+repeat_index] <= rpos:
                loci_flank_qpos_range[r+repeat_index][0] = qpos-insertion_length
                flank_track[r+repeat_index][0] = True 
            if coord_end-right_flank[r+repeat_index] <= rpos:
                loci_flank_qpos_range[r+repeat_index][1] = qpos
                if rpos > coord_end-right_flank[r+repeat_index]: flank_track[r+repeat_index][1] = True
        elif (not flank_track[r+repeat_index][1]) and (coord_end-right_flank[r+repeat_index] <= rpos):
            loci_flank_qpos_range[r+repeat_index][1] = qpos
            if rpos > coord_end-right_flank[r+repeat_index]: flank_track[r+repeat_index][1] = True

        # read_loci_variations[locus_key][rpos] = f'I|{insertion_length}'
        if coord_start+left_flank[r+repeat_index] <= rpos <= coord_end-right_flank[r+repeat_index]: # introduced to include length only if it comes inside repeat region
            if rpos not in homopoly_positions:
                read_loci_variations[locus_key]['alen'] += insertion_length
                read_loci_variations[locus_key]['halen'] += insertion_length
            else:
                if len(set(insert)) == 1:
                    # only if the insertion is a homopolymer; consider it as homopolymer insertion
                    read_loci_variations[locus_key]['halen'] += insertion_length
                else:
                    read_loci_variations[locus_key]['alen'] += insertion_length
                    read_loci_variations[locus_key]['halen'] += insertion_length

        if coord_start <= rpos <= coord_start+left_flank[r+repeat_index]-1: # -1 is included so ins near the start pos is not taken into account as it is already added
            try:
                out_insertion_qpos_ranges_left[r+repeat_index].append((qpos-insertion_length, qpos))
                left_ins_rpos[r+repeat_index].append(rpos)
            except AttributeError:
                pass
        elif coord_end-right_flank[r+repeat_index]+1 <= rpos <= coord_end: # +1 is included so ins near the end pos is not taken into account as it is already added
            try:
                out_insertion_qpos_ranges_right[r+repeat_index].append((qpos-insertion_length, qpos))
                right_ins_rpos[r+repeat_index].append(rpos)
            except AttributeError:
                pass
    jump = 0    # jump beyond the repeat where all positions are tracked
    if loci_coords[repeat_index + r - 1][1] < rpos:
        for f in loci_coords[repeat_index:]:
            if f[1] < rpos: jump += 1
            else: break

    return jump

# def Record_right_out_ins(amp_right_flank_list, r, repeat_index, locus_qpos_range, out_insertion_qpos_ranges_right, qpos, right_ins_rpos, coord_end, seq_len):
#     needed_len = amp_right_flank_list[r+repeat_index]
#     if needed_len>0:
#         locus_qpos_range[r+repeat_index][1] += needed_len
#         # soft_ins_len = needed_len
#         if qpos < seq_len:
#             out_insertion_qpos_ranges_right[r+repeat_index].append((qpos, qpos + needed_len))
#             right_ins_rpos[r+repeat_index].append(coord_end+1)

# def Record_left_out_ins(amp_left_flank_list, r, repeat_index, locus_qpos_range, out_insertion_qpos_ranges_left, left_ins_rpos, coord_start):
#     needed_len = amp_left_flank_list[r+repeat_index]
#     if needed_len>0:
#         mod_flank = locus_qpos_range[r+repeat_index][0] - needed_len
#         if mod_flank>0:
#             locus_qpos_range[r+repeat_index][0] = mod_flank
#             soft_ins_len = needed_len
#         else:
#             locus_qpos_range[r+repeat_index][0] = 0
#             soft_ins_len = needed_len + mod_flank

#         current_qpos = locus_qpos_range[r+repeat_index][0]
#         if current_qpos != (current_qpos + soft_ins_len):
#             out_insertion_qpos_ranges_left[r+repeat_index].append((current_qpos, current_qpos + soft_ins_len))
#             left_ins_rpos[r+repeat_index].append(coord_start - 1)

def ref_repeat(locus_key):
    lstart = int(locus_key[locus_key.index(':')+1 : locus_key.index('-')])
    lend = int(locus_key[locus_key.index('-')+1:])
    return lstart, lend

def Record_left_out_ins(amp_left_flank_list, r, repeat_index, chrom, ref, query_sequence, flank_length, coord_start, coord_end, locus_qpos_range, loci_flank_qpos_range, out_insertion_qpos_ranges_left, left_ins_rpos, flank_track, qpos_start, qpos_end):
    needed_len = amp_left_flank_list[r+repeat_index]
    if needed_len>0:
        updated_seq_start, status = detect_flank(chrom, query_sequence, ref, flank_length, coord_start, coord_end, qpos_start, qpos_end, True)
        if status:
            locus_qpos_range[r+repeat_index][0] = updated_seq_start # setting the read flank start
            flank_track[r+repeat_index][0] = True 
            loci_flank_qpos_range[r+repeat_index][0] = updated_seq_start # setting the repeat start in read same as its flank start
            out_insertion_qpos_ranges_left[r+repeat_index] = () # restricting the left insertion ranges as empty
            left_ins_rpos[r+repeat_index] = () # restricting the left insertion ref pos as empty
        else:
            out_insertion_qpos_ranges_left[r+repeat_index] = (None,) # tagging the read as None to ignore this read for genotyping, as it's softclip did not have sufficient ref flank

def Record_right_out_ins(amp_right_flank_list, r, repeat_index, chrom, ref, query_sequence, flank_length, coord_start, coord_end, locus_qpos_range, loci_flank_qpos_range, out_insertion_qpos_ranges_right, right_ins_rpos, flank_track, qpos_start, qpos_end):
    needed_len = amp_right_flank_list[r+repeat_index]
    if needed_len>0:
        updated_seq_end, status = detect_flank(chrom, query_sequence, ref, flank_length, coord_start, coord_end, qpos_start, qpos_end, False)
        if status:
            locus_qpos_range[r+repeat_index][1] = updated_seq_end # setting the read flank end
            flank_track[r+repeat_index][1] = True
            loci_flank_qpos_range[r+repeat_index][1] = updated_seq_end # setting the repeat end in read same as its flank end
            out_insertion_qpos_ranges_right[r+repeat_index] = () # restricting the right insertion ranges as empty
            right_ins_rpos[r+repeat_index] = () # restricting the right insertion ref pos as empty
        else:
            out_insertion_qpos_ranges_right[r+repeat_index] = (None,) # tagging the read as None to ignore this read for genotyping, as it's softclip did not have sufficient ref flank