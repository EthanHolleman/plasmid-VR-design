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
        vrs.append(VairableRegion.init_from_csv_row(p))
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


def test_VR_creation_from_file(test_VR, test_VR_with_at, test_VR_cluster):
    assert isinstance(test_VR, VairableRegion)
    assert isinstance(test_VR_with_at, VairableRegion)
    assert isinstance(test_VR_cluster, VairableRegion)


def test_calculate_nucleotide_counts(test_VR):
    test_VR.calculate_nucleotide_counts()
    assert sum(test_VR.gc_count) > 1 
    assert sum(test_VR.at_count) > 1
    assert sum(test_VR.nuc_dict.values()) == len(test_VR)


def test_calc_nuc_counts_cluster(test_VR_cluster):
    test_VR_cluster.calculate_nucleotide_counts()
    assert sum(test_VR_cluster.gc_count) > 1 
    assert sum(test_VR_cluster.at_count) > 1
    assert sum(test_VR_cluster.nuc_dict.values()) == len(test_VR_cluster)


def test_generate_sequence_cluster(test_VR_cluster):
    seq = test_VR_cluster.generate_sequence()
    assert isinstance(seq, Sequence)
    assert len(seq.nuc_seq) == len(test_VR_cluster)
    assert isinstance(seq.nuc_seq, str)


def test_calc_nuc_counts_at(test_VR_with_at):
    test_VR_with_at.calculate_nucleotide_counts()
    assert sum(test_VR_with_at.gc_count) > 1 
    assert sum(test_VR_with_at.at_count) > 1
    assert sum(test_VR_with_at.nuc_dict.values()) == len(test_VR_with_at)


def test_generate_sequence(test_VR):
    seq = test_VR.generate_sequence()
    assert isinstance(seq, Sequence)
    assert len(seq.nuc_seq) == len(test_VR)
    assert isinstance(seq.nuc_seq, str)


def test_generate_sequence_at(test_VR_with_at):
    seq = test_VR_with_at.generate_sequence()
    assert isinstance(seq, Sequence)
    assert len(seq.nuc_seq) == len(test_VR_with_at)
    assert isinstance(seq.nuc_seq, str) 


def test_nuc_seq_nucleotide_count(test_VR):
    seq = test_VR.generate_sequence()
    for nucleotide, count in test_VR.nuc_dict.items():
        assert seq.nuc_seq.count(nucleotide) == count


def test_seq_name(test_VR):
    assert isinstance(test_VR._sequence_name(), str)
    assert test_VR._sequence_name()


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




        
    