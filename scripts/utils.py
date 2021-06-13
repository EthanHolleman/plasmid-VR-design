import numpy as np
import pandas as pd


nuc_cache = {}


def nuc_count_calculator(skew, content, seq_len, closest=True):
    '''Calculate the number of G and C nucleotides to include in a DNA
    sequence given a sequence length, GC skew and content. 

    Args:
        gc_skew (float): Level of GC skew, between 0 and 1.
        gc_content (float): Level of GC content between 0 and 1.
        seq_len (int): Length of sequence in nucleotides.
        closest (bool, optional): If no exact result returns the closest. Defaults to True.
    '''
    closest_skew, closest_content = 0, 0
    best_nuc_combo = (0, 0)
    if seq_len not in nuc_cache:
        nuc_combos = []
        for k in range(2, seq_len):
            nuc_combos += [(i, k-i) for i in range(0, k+1)]
        nuc_cache[seq_len] = nuc_combos

    
    for nuc_combo in nuc_cache[seq_len]:
        
        cur_content, cur_skew = calculate_content(nuc_combo, seq_len), calculate_skew(nuc_combo)
        
        cur_content_diff = abs(cur_content-content)
        cur_skew_diff = abs(cur_skew-content)
        closest_content_diff = abs(closest_content-content)
        closest_skew_diff = abs(closest_skew-skew)

        if cur_content_diff < closest_content_diff and cur_skew_diff < closest_skew_diff:
            
            closest_skew = cur_skew
            closest_content = cur_content
            print(cur_skew, cur_content, nuc_combo, skew, content)
            best_nuc_combo = nuc_combo

    
    return best_nuc_combo


def calculate_skew(nuc_combo):
    return (nuc_combo[1] - nuc_combo[0]) / sum(nuc_combo)


def calculate_content(nuc_combo, seq_len):
    return sum(nuc_combo) / seq_len
        
    
def range_is_occupied(occupied_coords, start, end):
    sites = np.arange(start, end, dtype=int)
    occupied_coords = occupied_coords.astype(int)
    if any(np.take(occupied_coords, sites)) != 0:
        return True  # is occupied
    else:
        return False


def random_range_of_length_n(length_seq, range_length):
    start = int(np.random.choice(np.arange(0, length_seq-range_length), 1)[0])
    end = start + range_length
    return start, end


def find_available_random_range(occupied_coords, range_length):
    if longest_unoccupied_gap(occupied_coords) >= range_length:
        while True:
            start, end = random_range_of_length_n(
                len(occupied_coords), range_length
                )
            if range_is_occupied(occupied_coords, start, end):
                continue
            else:
                return int(start), int(end)
    else:
        # no possible ranges
        return False


def longest_unoccupied_gap(occupied_coords):
    # helper function to find longest gap (run of false values) in a boolean
    # array. Used to determine is there is still space in a list for another
    # cluster
    return len(max(''.join([str(i) for i in occupied_coords]).split('1'), key=lambda s: len(s)))


def get_int_half_length(length):
    half_len = int(length / 2)
    other_half = length - half_len
    return half_len, other_half


def read_variable_region_config_file(file_path):
    table = pd.read_table(file_path)
    table = table.to_dict(orient='records')
    for row in table:
        for key, val in row.items():
            if val == 'NA' or pd.isna(val):
                row[key] = None

    return table


