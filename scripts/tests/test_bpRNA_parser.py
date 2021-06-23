import pytest
import os
import sys
from pathlib import Path

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from bpRNA_parser import *

cur_dir = os.path.dirname(os.path.realpath(__file__))

@pytest.fixture
def st_path():
    return os.path.join(cur_dir, '../test_files/init-1.st')

@pytest.fixture
def sf_instance(st_path):
    return StructureFile(st_path)

@pytest.fixture
def test_output_tsv(st_path):
    return str(Path(st_path).with_suffix('.test.tsv'))


def test_init_StructureFile(st_path):
    sf = StructureFile(st_path)
    assert isinstance(sf, StructureFile)
    assert sf.dot_bracket != None
    assert sf.struct_str != None
    assert sf.seq != None
    assert sf.header != None


def test_struct_descrip_parsing(sf_instance):
    for sd in sf_instance.struct_descrips:
        assert isinstance(sd, StructDescrip)
        assert sd.code in StructDescrip.structure_codes


# should be three from looking at init-1.st
def test_number_hairpins(sf_instance):
    assert sf_instance.num_hairpins == 3


# should be 0.79
def test_prop_unpaired(sf_instance):
    assert sf_instance.prop_unpaired == 0.79

# 20 ribos in hairpin / length 200
def test_prop_hairpin(sf_instance):
    assert sf_instance.prop_hairpin == 0.1


# the struct_span (length) of all hairpin structures should be parsable
def test_hairpin_spans(sf_instance):
    for each_sd in sf_instance.struct_descrips:
        if each_sd.code == 'H':  # is a hairpin
            assert each_sd.struct_range
            assert each_sd.struct_range[0] < each_sd.struct_range[1]

def test_write_to_csv(sf_instance, test_output_tsv):
    sf_instance.to_tsv(test_output_tsv)
    assert os.path.isfile(test_output_tsv)
    os.remove(test_output_tsv)


