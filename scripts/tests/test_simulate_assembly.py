import os
import sys
from pathlib import Path
import random
import numpy as np

import pytest
from pydna.genbankrecord import GenbankRecord
from pydna.readers import read
from pydna.dseqrecord import Dseqrecord

test_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(test_dir)

from simulate_assembly import *


TEST_FILES = Path(test_dir).joinpath('test_files')
INSERT_DIR = TEST_FILES.joinpath('assembled_inserts')
T7_INIT_BACKBONE = TEST_FILES.joinpath('backbones/pFC9.gb')
T7_TERM_BACKBONE = TEST_FILES.joinpath('backbones/pFC8.gb')


def read_genbank_record(filepath):
    return GenbankRecord(
        read(filepath)
    )

@pytest.fixture
def dummy_output():
    return TEST_FILES.joinpath('dummy.gb')

@pytest.fixture
def dummy_output_dir():
    return TEST_FILES.joinpath('dummy')

@pytest.fixture
def assembled_inserts():
    inserts = []
    for each_ai in INSERT_DIR.iterdir():
        inserts.append(
            read_genbank_record(each_ai)
        )
    return inserts

@pytest.fixture
def assembled_insert_paths():
      
    return list(INSERT_DIR.iterdir())


@pytest.fixture
def strong_initiator():
    # strong initiator seq not known unitl experiments with initiators are
    # done. Seq doesn't matter as long as starts with HindIII site and ends
    # with an EcoRI site so generate it it here.
    si = Dseqrecord(
        ''.join(
            [
                HindIII.site,
                ''.join(list(np.random.choice(
                    ['A', 'T', 'C', 'G'], 200, replace=True))),
                EcoRI.site
            ]
        )
    )
    return si


@pytest.fixture
def T7_init_backbone():
    return read_genbank_record(T7_INIT_BACKBONE)


@pytest.fixture
def T7_term_backbone():
    return read_genbank_record(T7_TERM_BACKBONE)


@pytest.fixture
def random_name():
    return f'RAND_NAME_{random.randint(0, 1000)}'


@pytest.fixture
def series_instances(T7_init_backbone, T7_term_backbone, strong_initiator,
                     assembled_inserts):
    random_insert = assembled_inserts[random.randint(
        0, len(assembled_inserts)-1)]
    return [
        T7InitiationSeries('T7-init', T7_init_backbone, random_insert),
        T7TerminationSeries(
            'T7-term', T7_term_backbone, random_insert, strong_initiator
        )
    ]


def test_series_init(random_name, T7_init_backbone, assembled_inserts):
    for each_insert in assembled_inserts:
        s = Series(random_name, T7_init_backbone, each_insert)
        assert isinstance(s, Series)


def test_large_fragment(series_instances):
    for each_series in series_instances:
        assert isinstance(each_series, Series)
        lf = each_series.large_fragment
        assert isinstance(lf, Dseqrecord)


def test_assembly(series_instances):
    for each_series in series_instances:
        assembly = each_series.assemble_circular()
        assert isinstance(assembly, Dseqrecord)
        assert len(assembly.seq) < len(each_series.backbone) * 2


# def test_write_assembly(series_instances, dummy_output):
#     for each_series in series_instances:
#         each_series.write_assembly(
#             dummy_output, 'v1.0'
#         )
#         assert dummy_output.is_file()
#         contents = open(str(dummy_output)).read()
#         assert len(contents) > 10  # should have stuff
    
#     os.remove(str(dummy_output))
# this gets done below many times


def test_write_assembly_dir_t7_initiation(assembled_insert_paths,
                                          dummy_output_dir):

    output_files = T7InitiationSeries.write_assemblies_from_insert_list(
        assembled_insert_paths, str(T7_INIT_BACKBONE), dummy_output_dir
    )
    for each_file in output_files:
        assert each_file.is_file()
        os.remove(str(each_file))


def test_write_assemble_dir_t7_termination(assembled_insert_paths,
                                           dummy_output_dir, strong_initiator):
    
    output_files = T7TerminationSeries.write_assemblies_from_insert_list(
        assembled_insert_paths, str(T7_TERM_BACKBONE), strong_initiator,
        dummy_output_dir
        )
    for each_file in output_files:
        assert each_file.is_file()
