import pytest
import os
import sys
from pathlib import Path
from pydna.genbankrecord import GenbankRecord

test_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(test_dir)

from assemble_inserts import *

VR_DIR = Path(test_dir).joinpath('test_files/variable_regions')

@pytest.fixture
def variable_region_tsv_files():
    return list(Path(VR_DIR).iterdir())

@pytest.fixture
def downstream_filepaths():
    return [Path(test_dir).joinpath('test_files/insert/3_prime_homology_arm.gb')]

@pytest.fixture
def upstream_filepaths():
    return [Path(test_dir).joinpath('test_files/insert/3_prime_homology_arm.gb')]


@pytest.fixture
def vr_dicts(variable_region_tsv_files):
    dicts = []
    for each_file in variable_region_tsv_files:
        dicts.append(read_variable_region_tsv(each_file))
    return dicts


def test_read_variable_region_tsv(variable_region_tsv_files):
    for each_file in variable_region_tsv_files:
        vr_dict = read_variable_region_tsv(each_file)
        assert isinstance(vr_dict, dict)
        assert len(vr_dict) > 0
        assert 'name' in vr_dict
        assert 'GC_skew' in vr_dict
        assert 'GC_content' in vr_dict


def test_variable_region_to_labeled_record(vr_dicts):
    id_val = 'test_id'
    for each_dict in vr_dicts:
        assert isinstance(each_dict, dict)
        record = variable_region_to_labeled_record(each_dict, id_val)
        assert record.locus == each_dict['name']
        assert len(record.seq) > 0
        assert len(record.features) == 1


def test_read_records(upstream_filepaths, downstream_filepaths):
    for filepath_list in (upstream_filepaths, downstream_filepaths):
        records = read_records(filepath_list)
        assert len(records) == len(filepath_list)
        for each_record in records:
            assert isinstance(each_record, GenbankRecord)



# def test_assemble_insert(upstream_features, variable_region_records, downstream_records):
#     for each_vr in variable_region_records:
#         insert = assemble_inserts(*upstream_features, each_vr, *downstream_records)
#         assert isinstance(insert, Genbankrecord)



