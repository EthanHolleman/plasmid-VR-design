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
TAC_BACKBONE = TEST_FILES.joinpath('backbones/pFC53tacT1T2.gb')
INITIATOR = TEST_FILES.joinpath('insert/placeholder_initiator.gb')


def read_genbank_record(filepath):
    return GenbankRecord(
        read(str(filepath))
    )

@pytest.fixture
def dummy_output():
    return TEST_FILES.joinpath('dummy.gb')

@pytest.fixture
def dummy_output_dir():
    p = TEST_FILES.joinpath('temp')
    p.mkdir(exist_ok=True, parents=True)
    return p

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
    return read_genbank_record(INITIATOR)


@pytest.fixture
def T7_init_backbone():
    return read_genbank_record(T7_INIT_BACKBONE)


@pytest.fixture
def T7_term_backbone():
    return read_genbank_record(T7_TERM_BACKBONE)

@pytest.fixture
def tac_backbone():
    return TAC_BACKBONE


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


@pytest.fixture
def t7_init_series(random_name, T7_init_backbone, assembled_inserts):
    constructs = []
    for i, each_insert in enumerate(assembled_inserts):
        name = f'test_construct_{i}'
        constructs.append(
            T7InitiationSeries(name, T7_init_backbone, each_insert)
        )
    
    return constructs


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


def test_write_assembly_dir_t7_initiation(assembled_insert_paths,
                                          dummy_output_dir):

    output_files = T7InitiationSeries.write_assemblies_from_insert_list(
        assembled_insert_paths, str(T7_INIT_BACKBONE), dummy_output_dir
    )
    for each_file in output_files:
        assert each_file.is_file()
        record = read_genbank_record(str(each_file))
        os.remove(str(each_file))


def test_write_assemble_dir_t7_termination(assembled_insert_paths,
                                           dummy_output_dir, strong_initiator):
    
    output_files = T7TerminationSeries.write_assemblies_from_insert_list(
        assembled_insert_paths, str(T7_TERM_BACKBONE), strong_initiator,
        dummy_output_dir
        )
    for each_file in output_files:
        assert each_file.is_file()
        record = read_genbank_record(str(each_file))
        # make sure strong init is included
        assert record.seq.find(str(strong_initiator.seq)) != -1, f'{record.seq.find(strong_initiator.seq)}'
        os.remove(str(each_file))


def test_tac_initiation_series_init(t7_init_series, tac_backbone, dummy_output_dir):
    output_file = Path(dummy_output_dir).joinpath('test_tac_assembly.gb')
    for each_construct in t7_init_series:

        each_construct.write_assembly(str(output_file), 'test_tac_assembly', 1, 'test')

        assert os.path.isfile(str(output_file))

        tac_init_series = TacInitiationSeries(
            'test_tac', str(output_file), str(tac_backbone)
            )
        assert isinstance(tac_init_series, TacInitiationSeries)
        assert tac_init_series.parent_construct
        assert isinstance(tac_init_series.parent_construct, GenbankRecord)


@pytest.fixture
def t7_initiation_series_assemblies(assembled_insert_paths, dummy_output_dir):
    assemblies = []

    t7_dir = Path(dummy_output_dir).joinpath('t7_init_assemblies')
    t7_dir.mkdir(exist_ok=True, parents=True)
    output_files = T7InitiationSeries.write_assemblies_from_insert_list(
        assembled_insert_paths, str(T7_INIT_BACKBONE), t7_dir
    )
    return output_files
  

@pytest.fixture
def t7_termination_series_assemblies(assembled_insert_paths,
                                           dummy_output_dir, strong_initiator):
    t7_dir = Path(dummy_output_dir).joinpath('t7_term_assemblies')
    t7_dir.mkdir(exist_ok=True, parents=True)
    output_files = T7TerminationSeries.write_assemblies_from_insert_list(
        assembled_insert_paths, str(T7_TERM_BACKBONE), strong_initiator,
        t7_dir
        )
    return output_files


# def test_tac_init_write_assemblies(t7_initiation_series_assemblies, 
#                                    tac_backbone, dummy_output_dir):

#     out_dir = Path(dummy_output_dir).joinpath('Tac_init_constructs')
#     out_dir.mkdir(exist_ok=True, parents=True)
#     primer_path = out_dir.joinpath('primers.fa')
#     output_files = TacInitiationSeries.write_assemblies_from_parent_constructs(
#         t7_initiation_series_assemblies, str(tac_backbone), out_dir
#     )
#     for each_file in output_files:
#         record = GenbankRecord(read(each_file))
#         assert len(record) > 0
#         #assert len(record.features) == 13
#         feature_types = [each_feature.type for each_feature in record.features]
#         assert 'promoter' in feature_types
#         assert 'terminator' in feature_types
    
#     TacInitiationSeries.write_primers(primer_path)

# def test_tac_term_write_assemblies(t7_termination_series_assemblies, 
#                                     tac_backbone, dummy_output_dir):
#     out_dir = Path(dummy_output_dir).joinpath('Tac_term_constructs')
#     out_dir.mkdir(exist_ok=True, parents=True)
#     primer_path = out_dir.joinpath('primers.fa')
#     output_files = TacTerminationSeries.write_assemblies_from_parent_constructs(
#         t7_termination_series_assemblies, tac_backbone, out_dir
#     )
#     for each_file in output_files:
#         record = GenbankRecord(read(each_file))
#         assert len(record) > 0
#         assert len(record.features) == 12
#         feature_types = [each_feature.type for each_feature in record.features]
#         assert 'promoter' in feature_types
#         assert 'terminator' in feature_types
    
#     TacTerminationSeries.write_primers(primer_path)



