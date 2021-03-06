import os
import pytest

from varmids.variable_region_maker import *
from varmids.utils import *
import random

import os

cur_dir = os.path.dirname(os.path.realpath(__file__))

TEST_PLASMIDS_PATH = os.path.join(cur_dir, 'test_files/initiation_plamids.tsv')
TEST_PLASMIDS_DICT = read_variable_region_config_file(TEST_PLASMIDS_PATH)


@pytest.fixture
def all_VRs():
    vrs = []
    for p in TEST_PLASMIDS_DICT:
        vr = VariableRegion.init_from_csv_row(p)
        vrs.append(vr)
    return vrs


@pytest.fixture
def test_VR():
    p = TEST_PLASMIDS_DICT[0]
    return VariableRegion.init_from_csv_row(p)


@pytest.fixture
def test_VR_with_at():
    p = TEST_PLASMIDS_DICT[1]
    # second variable region also has AT content and skew
    # set
    return VariableRegion.init_from_csv_row(p)


@pytest.fixture
def test_VR_cluster():
    p = TEST_PLASMIDS_DICT[2]
    return VariableRegion.init_from_csv_row(p)


def test_VR_creation_from_file(all_VRs):
    for each_vr in all_VRs:
        assert isinstance(each_vr, VariableRegion)


def test_calculate_nucleotide_counts(all_VRs):
    for each_vr in all_VRs:
        each_vr.calculate_nucleotide_counts()
        assert sum(each_vr.gc_count) > 1
        assert sum(each_vr.at_count) > 1
        assert sum(each_vr.nuc_dict.values()) == each_vr.length


def test_generate_sequence(all_VRs):
    for each_vr in all_VRs:
        seq = each_vr.generate_sequence()
        assert isinstance(seq, Sequence)
        assert len(seq.nuc_seq) == each_vr.length
        assert isinstance(seq.nuc_seq, str)
        assert len(set(seq.nuc_seq)) <= 4  # no more than 4 diff nucleotides


def test_nuc_seq_nucleotide_count(all_VRs):
    for each_vr in all_VRs:
        seq = each_vr.generate_sequence()
        for nucleotide, count in each_vr.nuc_dict.items():
            assert seq.nuc_seq.count(nucleotide) - count <= 2
            # tolerance of 1 nuc diff to account for rounding


def test_seq_name(all_VRs):
    seq_id_num = random.randint(0, 1000)
    for each_vr in all_VRs:
        assert isinstance(each_vr._sequence_name(seq_id_num), str)
        assert each_vr._sequence_name(seq_id_num)


def test_all_gc_content(all_VRs):
    for each_vr in all_VRs:
        if each_vr.gc_content:
            seq = each_vr.generate_sequence()
            assert abs(each_vr.gc_content - seq.gc_content) <= 0.06


def test_all_at_content(all_VRs):
    for each_vr in all_VRs:
        if each_vr.at_content:
            seq = each_vr.generate_sequence()
            assert abs(each_vr.at_content - seq.at_content) <= 0.06


def test_all_gc_skew(all_VRs):
    for each_vr in all_VRs:
        if each_vr.gc_skew:
            seq = each_vr.generate_sequence()
            assert abs(each_vr.gc_skew - seq.gc_skew) <= 0.06


def test_count_vs_nuc_dict(all_VRs):
    for each_vr in all_VRs:
        seq = each_vr.generate_sequence()
        for each_nuc in seq.nuc_counts:
            assert seq.count_dict[each_nuc] - seq.nuc_counts[each_nuc] <= 2


def test_all_at_skew(all_VRs):
    for each_vr in all_VRs:
        if each_vr.at_skew:
            seq = each_vr.generate_sequence()
            assert abs(each_vr.at_skew - seq.at_skew) <= 0.05, each_vr.name


# def test_reverse_complement(all_VRs):
#     for each_vr in all_VRs:
#         seq = each_vr.generate_sequence()
#         rc = seq.reverse_complement()
#         assert isinstance(rc, Sequence)
#         assert len(seq) == len(rc)
#         # gc and at content should not change
#         assert seq.gc_content == rc.gc_content
#         assert seq.at_content == rc.at_content

#         # gc and at skew should be same but opposite
#         assert seq.gc_skew == -rc.gc_skew
#         assert seq.at_skew == -rc.at_skew
