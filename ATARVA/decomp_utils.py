from collections import Counter
from bitarray import bitarray

def convert_to_bitset(seq):
    lbit = {'A': '0', 'C': '0', 'G': '1', 'T': '1', 'N': '1'}
    rbit = {'A': '0', 'C': '1', 'G': '0', 'T': '1', 'N': '1'}
    
    lbitseq = bitarray()
    rbitseq = bitarray()
    
    for s in seq:
        lbitseq.extend(lbit.get(s, '1'))
        rbitseq.extend(rbit.get(s, '1'))
    
    return lbitseq, rbitseq

def shift_and_match(seq):
    shift_list = []
    # best_shift = motif_length
    max_matches = 0

    shift_values = set(range(1, 11))

    for shift in sorted(shift_values):
        if shift < 1:
            continue

        lbitseq, rbitseq = convert_to_bitset(seq)

        lmatch = ~(lbitseq ^ (lbitseq >> shift))
        rmatch = ~(rbitseq ^ (rbitseq >> shift))
        match = lmatch & rmatch
        shift_list.append(match)

    return shift_list

def kmp_search_non_overlapping(text, pattern):
    def compute_lps(pattern):
        lps = [0] * len(pattern)
        length = 0
        i = 1
        while i < len(pattern):
            if pattern[i] == pattern[length]:
                length += 1
                lps[i] = length
                i += 1
            else:
                if length != 0:
                    length = lps[length - 1]
                else:
                    lps[i] = 0
                    i += 1
        return lps

    lps = compute_lps(pattern)
    result = []
    i = 0
    j = 0  

    while i < len(text):
        if pattern[j] == text[i]:
            i += 1
            j += 1

        if j == len(pattern):  
            result.append(i - j)
            j = 0  
        elif i < len(text) and pattern[j] != text[i]:
            if j != 0:
                j = lps[j - 1]
            else:
                i += 1

    return result

def get_most_frequent_motif(sequence, motif_size, motif):
    cyc_motif = get_canonical_motif(motif) if motif else None
    if len(sequence) < motif_size:  
        return sequence, None  

    repeating_units = [sequence[i:i + motif_size] for i in range(len(sequence) - motif_size + 1)]
    motif_counts = Counter(repeating_units)
    
    
    if not motif_counts:  
        return sequence, None  

    most_common = []
    for each_motif,count in motif_counts.most_common():
        if len(most_common) >=2: break
        if get_canonical_motif(each_motif) == cyc_motif:
            if motif not in most_common: most_common.append(motif)
        else:
            most_common.append(each_motif)

    primary_motif = most_common[0] if most_common else sequence
    secondary_motif = most_common[1] if len(most_common) > 1 else None

    return primary_motif, secondary_motif
    
def get_canonical_motif(motif):
    return min(motif[i:] + motif[:i] for i in range(len(motif)))

def max_match(shift_list, gap_regions):
    gap_wise_shift = []
    
    for each_gap in gap_regions:
        max_matches = 0
        best_shift = -1
        for idx, each_shift in enumerate(shift_list):
            if each_shift==0: continue
            each_shift = each_shift[each_gap[0]: each_gap[1]]
            count = 0
            count += each_shift.count() if each_shift.count() > (idx+1) else 0
    
            if count > max_matches:
                max_matches = count
                best_shift = idx+1
        gap_wise_shift.append(best_shift)
    return gap_wise_shift

    
slide_threshold = {1:1, 2:2, 3:3, 4:4, 5:5, 6:6, 7:7, 8:8}
def window_scan(shift_list, motif_size, sequence, sequential_decomp, sequential_part, overall_boundary, current_gap):
    shift_seq = shift_list[motif_size-1]#[current_gap[0]:current_gap[1]]

    slide_size = slide_threshold.get(motif_size, 8)
    
    i=current_gap[0]
    start = i; end = current_gap[1]
    start_track = start
    initial = True
    while i<(current_gap[1]-slide_size + 1):
        words = shift_seq[i:i+slide_size]

        if (words.count()/len(words)) >= 0.9:
            if initial:
                calc_start = i-motif_size if i-motif_size >= start_track else start_track
                if calc_start>-1:
                    for b in overall_boundary:
                        if not (b[0] <= calc_start < b[1]):
                            pass
                        else:
                            start = b[1]
                            break
                    else:
                        start = calc_start
                else:
                    start = i #i-motif_size if i-motif_size > -1 else i
                initial = False
            end = i+motif_size
            i+=motif_size

            continue
        # i+=slide_size
        else:
            i+=1
            if not initial:
                calc_end = end+motif_size
                for b in overall_boundary:
                    if (start< calc_end <= b[0]) or (calc_end > start >= b[1]):
                        pass
                    else:
                        end = b[0]
                        break
                else:
                    end = calc_end    
                if start!=end:
                    start_track = decomposer([start, end], sequence, motif_size, sequential_decomp, sequential_part, shift_list, overall_boundary)
                initial = True
    if not initial:
        calc_end = end+motif_size
        for b in overall_boundary:
            if (start< calc_end <= b[0]) or (calc_end > start >= b[1]):
                pass
            else:
                end = b[0]
                break
        else:
            end = calc_end 
        if start!=end:
            start_track = decomposer([start, end], sequence, motif_size, sequential_decomp, sequential_part, shift_list, overall_boundary)


def shift_decomp(seq, motif_size, motif, boundary, state):
    decomposed_parts = []
    count = 1  
    primary_motif, secondary_motif = get_most_frequent_motif(seq, motif_size, motif)
    positions = kmp_search_non_overlapping(seq, primary_motif)
    
    if not positions:
        return seq, boundary


    seq = seq[positions[0] : positions[-1]+motif_size]
    
    if state:
        b1 = boundary[0]
        boundary = [b1+positions[0], b1+positions[-1]+motif_size]

        
    for i in range(1, len(positions)):
        if positions[i] == positions[i - 1] + len(primary_motif):
            count += 1
        else:
            decomposed_parts.append(f"({primary_motif}){count}")
            interspersed = seq[positions[i - 1] + len(primary_motif):positions[i]]
            if interspersed:
                primary_motif, secondary_motif = get_most_frequent_motif(seq, motif_size, '')
                secondary_decomp, boundary = shift_decomp(interspersed, motif_size, '', boundary, False)
                if secondary_decomp:
                    decomposed_parts.append(secondary_decomp)

            count = 1
    decomposed_parts.append(f"({primary_motif}){count}")
    last_motif_end = positions[-1] + len(primary_motif)
    leftover_sequence = seq[last_motif_end:]
    if leftover_sequence:
        decomposed_parts.append(leftover_sequence)
    return "-".join(decomposed_parts), boundary

def decomposer(tuple_bound, sequence, best_shift, sequential_decomp, sequential_part, shift_list, overall_boundary):
    processed_seq, tuple_bound = shift_decomp(sequence[tuple_bound[0]:tuple_bound[1]], best_shift, '', tuple_bound, True)
    b1 = tuple_bound[0]; b2 = tuple_bound[1]
    sequential_decomp[b1] = processed_seq
    sequential_part[b1:b2] = 0
    for id in range(len(shift_list)):
        if shift_list[id] == 0: continue
        shift_list[id][b1:b2] = 0
    overall_boundary.append(tuple_bound)
    return b2

def gap_boundaries(overall_boundary, seq_len, chunk_state):
    overall_boundary = sorted(overall_boundary)
    if chunk_state:
        chunk_boundary = overall_boundary[0][0]
        chunk_boundary_end = overall_boundary[-1][1]
    else:
        chunk_boundary = 0
        chunk_boundary_end = seq_len
        
    last_start = overall_boundary[-1][0]
    gap_regions = []
    for id,occupied_reg in enumerate(overall_boundary):
    
        if id == 0:
            if occupied_reg[0] != chunk_boundary:
                gap_regions.append([chunk_boundary, occupied_reg[0]])
            elif len(overall_boundary) == 1:
                if occupied_reg[1] <= (seq_len-1):
                    gap_regions.append([occupied_reg[1],seq_len])
            current_end = occupied_reg[1]
    
        else:
            if occupied_reg[0] != current_end:
                gap_regions.append([current_end, occupied_reg[0]])
            current_end = occupied_reg[1]

            
    if current_end != chunk_boundary_end:
        gap_regions.append([current_end, chunk_boundary_end])
        
    return gap_regions

def motif_decomposition(sequence):

    shift_list = shift_and_match(sequence)
    overall_boundary = []
    sequential_decomp = {}
    sequential_part = bitarray([1]*len(sequence))
    seq_len = len(sequence)
    gap_regions = [[0, seq_len]]
    rounds = 0
    while any(sequential_part):
        if rounds >= 1:
            shift_list[0][:] = 0
        gap_wise_shift = max_match(shift_list, gap_regions)

        for index, best_shift in enumerate(gap_wise_shift):
            
            current_gap = gap_regions[index]
            
            if best_shift!=-1:
                previous_ovbound = overall_boundary.copy()
                window_scan(shift_list, best_shift, sequence, sequential_decomp, sequential_part, overall_boundary, current_gap)
                if previous_ovbound == overall_boundary:
                    c1 = current_gap[0]; c2 = current_gap[1]
                    covered_chunks = [i for i in overall_boundary if (i[0]>=c1 and i[1]<=c2) or ((i[1]==c1) or (i[0]==c2))]
                    if covered_chunks == []:
                        sequential_part[c1 : c2] = 0
                        sequential_decomp[c1] = sequence[c1 : c2]
                        overall_boundary.append([c1,c2])
                        continue

                    gap_in_chunks = gap_boundaries(covered_chunks, c2, True)
                    if not gap_in_chunks:
                        if (c1 == 0) or (c2 == seq_len):
                            gap_in_chunks = [current_gap]
                    
                    for each_gap in gap_in_chunks:
                        sequential_part[each_gap[0] : each_gap[1]] = 0
                        for shift_id in range(len(shift_list)):
                            shift_list[shift_id][each_gap[0] : each_gap[1]] = 0
                        sequential_decomp[each_gap[0]] = sequence[each_gap[0] : each_gap[1]]

                        overall_boundary.append([each_gap[0] , each_gap[1]])

                
            else:
                if sequence[current_gap[0] : current_gap[1]]:
                    sequential_decomp[current_gap[0]] = sequence[current_gap[0] : current_gap[1]] 
                sequential_part[current_gap[0] : current_gap[1]] = 0
                overall_boundary.append([current_gap[0] , current_gap[1]])


        gap_regions = gap_boundaries(overall_boundary, seq_len, False)

        rounds += 1
        if rounds >= 20 :
            break
        
    fseq = [v for k,v in sorted(sequential_decomp.items(), key = lambda x : x[0])]

    return "-".join(fseq)
