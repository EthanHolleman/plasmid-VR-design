import os
import pytest

from variable_region_maker import *
from utils import *

import os

cur_dir = os.path.dirname(os.path.realpath(__file__))

TEST_PLASMIDS_PATH = os.path.join(cur_dir, 'test_files/initiation_plamids.csv')
TEST_PLASMIDS_DICT = read_variable_region_config_file(TEST_PLASMIDS_PATH)


@pytest.fixture
def all_VRs():
    vrs = []
    for p in TEST_PLASMIDS_DICT:
        vr = VairableRegion.init_from_csv_row(p)
        vrs.append(vr)
    return vrs


@pytest.fixture
def test_VR():
    p = TEST_PLASMIDS_DICT[0]
    return VairableRegion.init_from_csv_row(p)


@pytest.fixture
def test_VR_with_at():
    p = TEST_PLASMIDS_DICT[1] 
    # second variable region also has AT content and skew
    # set
    return VairableRegion.init_from_csv_row(p)


@pytest.fixture
def test_VR_cluster():
    p = TEST_PLASMIDS_DICT[2]
    return VairableRegion.init_from_csv_row(p)


def test_VR_creation_from_file(all_VRs):
    for each_vr in all_VRs:
        assert isinstance(each_vr, VairableRegion)


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
            assert seq.nuc_seq.count(nucleotide) - count <= 1
            # tolerance of 1 nuc diff to account for rounding


def test_seq_name(all_VRs):
    for each_vr in all_VRs:
        assert isinstance(each_vr._sequence_name(), str)
        assert each_vr._sequence_name()


def test_all_gc_content(all_VRs):
    for each_vr in all_VRs:
        if each_vr.gc_content:
            seq = each_vr.generate_sequence()
            calculated_gc_count = (seq.nuc_seq.count('G'), seq.nuc_seq.count('C'))
            calculated_gc_content = calculate_content(calculated_gc_count, len(seq))
            assert abs(calculated_gc_content - each_vr.gc_content) < 0.05


def test_all_at_content(all_VRs):
    for each_vr in all_VRs:
        if each_vr.at_content:
            seq = each_vr.generate_sequence()
            calculated_at_count = (seq.nuc_seq.count('A'), seq.nuc_seq.count('T'))
            calculated_at_content = calculate_content(calculated_at_count, len(seq))
            assert abs(calculated_at_content - each_vr.at_content) < 0.05


def test_all_gc_skew(all_VRs):
    for each_vr in all_VRs:
        if each_vr.gc_skew:
            seq = each_vr.generate_sequence()
            calculated_gc_count = (seq.nuc_seq.count('G'), seq.nuc_seq.count('C'))
            calculated_gc_skew = calculate_skew(calculated_gc_count)
            assert abs(calculated_gc_skew - each_vr.gc_skew) < 0.05


def test_all_at_skew(all_VRs):
    for each_vr in all_VRs:
        if each_vr.at_skew:
            seq = each_vr.generate_sequence()
            calculated_at_count = (seq.nuc_seq.count('A'), seq.nuc_seq.count('T'))
            calculated_at_skew = calculate_skew(calculated_at_count)
            assert abs(calculated_at_skew - each_vr.at_skew) < 0.05




        
    