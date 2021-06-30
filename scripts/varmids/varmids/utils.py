import numpy as np
import pandas as pd
import csv
from Bio.Restriction import *
from Bio.Seq import Seq

SEED = 12311997  # turn this into a parameter in snakemake somewhere in future

RAND_GEN = np.random.default_rng(SEED)

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
    best_nuc_combo = (1, 1)

    number_nucs = int(seq_len * content)

    nuc_combos = [(i, number_nucs-i) for i in range(1, number_nucs+1)]
    skew_content_dict = {
        tuple(nuc_combo): (calculate_content(nuc_combo, seq_len), calculate_skew(nuc_combo))
        for nuc_combo in nuc_combos
    }
    best_nuc_combo, cur_distance = None, float('inf')

    for nuc_combo in skew_content_dict:
        combo_content, combo_skew = skew_content_dict[nuc_combo]
        distance = abs(combo_content - content) + abs(combo_skew - skew)
        if distance == 0:
            best_nuc_combo = nuc_combo
            break
        elif distance < cur_distance:
            best_nuc_combo, cur_distance = nuc_combo, distance

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
    start = int(RAND_GEN.choice(np.arange(0, length_seq-range_length), 1)[0])
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


def read_variable_region_dataframe(dataframe):
    # convert pandas input to clean dictionary
    table = dataframe.to_dict(orient='records')
    for row in table:
        for key, val in row.items():
            if val == 'NA' or pd.isna(val):
                row[key] = None

    return table


def read_variable_region_config_file(file_path):
    table = pd.read_table(file_path)
    table = table.to_dict(orient='records')
    for row in table:
        for key, val in row.items():
            if val == 'NA' or pd.isna(val):
                row[key] = None

    return table


def write_sequence_list_to_output_files(seq_list, fasta_path, tsv_path):
    seq_rows = []
    with open(fasta_path, 'w') as fasta_handle:
        for seq in seq_list:
            fasta_handle.write(seq.as_fasta_entry())
            seq_rows.append(seq.to_dict())

    pd.DataFrame(seq_rows).to_csv(tsv_path, sep='\t', index=False, na_rep='NA')

    return fasta_path, tsv_path


def fasta_seq_printer(seq):
    n = 80
    return '\n'.join([seq[i:i+n] for i in range(0, len(seq), n)])


def count_nucleotides_in_seq(seq):
    nucs = ['A', 'T', 'G', 'C']
    counts = {each_nuc: seq.count(each_nuc) for each_nuc in nucs}
    return counts


class RestrictionSiteChecker():

    all_enzymes_dict = {str(e): e for e in AllEnzymes}

    def __init__(self, variable_region, cutters):
        self.variable_region = variable_region
        self.cutters = cutters

    @property
    def cutters(self):
        return self._cutters

    @cutters.setter
    def cutters(self, new_cutters):
        print(new_cutters)
        approved_cutters = []
        for each_new_cutter in new_cutters:
            if isinstance(each_new_cutter, str):
                if each_new_cutter in RestrictionSiteChecker.all_enzymes_dict:
                    approved_cutters.append(
                        RestrictionSiteChecker.all_enzymes_dict[each_new_cutter]
                    )
                else:
                    raise TypeError('Enzyme not present in collection')
            else:
                raise TypeError(
                    'Restriction enzymes must be specified by their name as a string')
        self._cutters = RestrictionBatch(approved_cutters)

    def _check_for_cutter(self, sequence):
        a = Analysis(self.cutters, Seq(sequence), linear=False).full()
        for cutter, cuts in a.items():
            if len(cuts) > 0:
                return False
        return True

    def generate_RE_free_sequence(self):
        seq = self.variable_region.generate_sequence()
        while self._check_for_cutter(seq.nuc_seq) == False:
            seq = self.variable_region.generate_sequence()
        return seq
