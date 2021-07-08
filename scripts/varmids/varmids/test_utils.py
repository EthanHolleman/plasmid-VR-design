import pytest
import numpy as np
import os
import random

from Bio.Restriction import *
from Bio.Seq import Seq

from varmids.utils import *
from varmids.test_variable_region_maker import all_VRs
from varmids.variable_region_maker import Sequence



# @pytest.fixture
# def known_skew_contents():
#     return {
        
#     }

# def test_nuc_count_calculator(known_skew_contents):
#     for params, vals in known_skew_content.items():
#         nuc_combo = nuc_count_calculator(*params)
#         assert nuc_combo == vals

@pytest.fixture
def seq_list(all_VRs):
    return [vr.generate_sequence() for vr in all_VRs]


@pytest.fixture
def fasta_path():
    cur_dir = os.path.dirname(os.path.realpath(__file__))
    return os.path.join(cur_dir, 'test_files/test.fasta')


@pytest.fixture
def tsv_path():
    cur_dir = os.path.dirname(os.path.realpath(__file__))
    return os.path.join(cur_dir, 'test_files/test.tsv')


@pytest.fixture
def variable_region_file():
    cur_dir = os.path.dirname(os.path.realpath(__file__))
    return os.path.join(cur_dir, 'test_files/initiation_plamids.tsv')

@pytest.fixture
def random_cutters():
    n_cutters = 5
    all_enzymes = list(AllEnzymes)
    return [
        str(all_enzymes[random.randint(0, len(all_enzymes))])
        for e in range(n_cutters)
        ]


@pytest.fixture
def rsc(all_VRs, random_cutters):
    for vr in all_VRs:
        return RestrictionSiteChecker(vr, random_cutters)


@pytest.fixture
def known_skews():
    return {
        (20, 20): 0,
        (200, 200): 0,
        (150, 250): -0.25,
        (250, 150): 0.25,
        (1500, 2500): -0.25,
        (2500, 1500): 0.25,
        (0, 100): -1,
        (100, 0): 1
    }


@pytest.fixture
def known_contents():
    kc = {}
    for _ in range(0, 100):
        count_a, count_b = np.random.randint(0, 1000, 2)
        length = np.random.randint(count_a + count_b, count_a + count_b + 1000)
        content = (count_a + count_b) / length
        kc[((count_a, count_b), length)] = content
    return kc


@pytest.fixture
def known_occupied():
    ko = {}
    
    for _ in range(0, 100):
        range_len = np.random.randint(1, 11)
        seq_len = np.random.randint(15, 500)
        sites = np.random.choice([0, 1], seq_len, p=[0.75, 0.25])
        s = np.random.randint(0, seq_len-range_len)
        e = s + range_len
        if sum(sites[s:e]) == 0:
            answer = True
        else:
            answer = False
        ko[answer] = (sites, s, e)
    
    return ko

    
@pytest.fixture
def known_range_lengths():
    rl = {}
    for _ in range(0, 100):
        seq_len = np.random.randint(101, 1000)
        range_len = np.random.randint(1, 100)
        rl[(seq_len, range_len)] = range_len
    return rl


def test_calculate_skew(known_skews):
    for params, val in known_skews.items():
        calculated_skew = calculate_skew(params)
        assert calculated_skew == val, f'{params}: {val} != {calculated_skew}'


def test_calculate_content(known_contents):
    for params, val in known_contents.items():
        calculated_content = calculate_content(*params)
        assert calculated_content == val


def test_random_range_of_length_n(known_range_lengths):
    for params, val in known_range_lengths.items():
        s, e = random_range_of_length_n(*params)
        assert e <= params[0]
        assert e > s
        assert e - s == val


# def test_range_is_occupied(known_occupied):
#     for val, params in known_occupied.items():
#         is_occupied = range_is_occupied(*params)
#         assert is_occupied == val

def test_get_int_half_length():
    for _ in range(0, 100):
        rand_int = np.random.randint(0, 10000)
        half = int(rand_int / 2)
        half_2 = rand_int - half
        calculated_halves = get_int_half_length(rand_int)
        assert calculated_halves == (half, half_2)
        assert sum(calculated_halves) - rand_int <= 1


def test_read_variable_region_config_file(variable_region_file):
    table = read_variable_region_config_file(variable_region_file)
    assert isinstance(table, list)
    assert(len(table)) > 0
    
    lengths = []

    for item in table:
        assert isinstance(item, dict)
        lengths.append(len(item))
        for key, val in item.items():
            assert not val == 'NA'
    
    assert len(set(lengths)) == 1


def test_write_sequence_list_to_output_files(seq_list, fasta_path, tsv_path):
    write_fasta, write_tsv =  write_sequence_list_to_output_files(
                                seq_list, fasta_path, tsv_path
                            )
    
    assert write_fasta == fasta_path
    assert write_tsv == tsv_path

    assert os.path.isfile(write_fasta)
    assert os.path.isfile(write_tsv)

    assert len(open(write_fasta).readlines()) > 0
    assert len(open(write_tsv).readlines()) > 0

    os.remove(write_fasta)
    os.remove(write_tsv)


def test_init_restrictionSiteChecker(rsc):
    assert isinstance(rsc, RestrictionSiteChecker)
    assert isinstance(rsc.cutters, RestrictionBatch)


def test_generate_RE_free_sequence(rsc):
    seq_id = random.randint(1, 1000)
    sequence = rsc.generate_RE_free_sequence(seq_id)
    assert sequence.id_num == seq_id
    assert isinstance(sequence, Sequence)
    assert isinstance(rsc.cutters, RestrictionBatch)
    ana_full = Analysis(rsc.cutters, Seq(sequence.nuc_seq)).full()
    for cutter, cuts in ana_full.items():
        assert len(cuts) == 0
    
    
        







