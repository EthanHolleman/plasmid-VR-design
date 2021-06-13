import os
import pytest

from variable_region_maker import *
from utils import *

import os

cur_dir = os.path.dirname(os.path.realpath(__file__))

TEST_PLASMIDS_PATH = os.path.join(cur_dir, 'test_files/initiation_plamids.csv')
TEST_PLASMIDS_DICT = read_variable_region_config_file(TEST_PLASMIDS_PATH)


@pytest.fixture
def test_VR():
    p = TEST_PLASMIDS_DICT.pop()
    return VairableRegion.init_from_csv_row(p)


def test_VR_creation_from_file(test_VR):
    assert isinstance(test_VR, VairableRegion)


def test_calculate_nucleotide_counts(test_VR):
    test_VR.calculate_nucleotide_counts()
    assert sum(test_VR.gc_count) > 1 
    assert sum(test_VR.at_count) > 1
    assert 0 == 1