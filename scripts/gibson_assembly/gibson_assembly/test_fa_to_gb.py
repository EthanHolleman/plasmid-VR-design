import pytest
from pathlib import Path
import pandas as pd
import os
import sys

CUR_DIR = os.path.dirname(os.path.abspath(__file__))

sys.path.append(CUR_DIR)

from gibson_assembly.fa_to_gb import *


@pytest.fixture
def genbank_record():
    return str(Path(CUR_DIR).joinpath('test_files/test_gb.gb'))


@pytest.fixture
def fasta_record():
    return str(Path(CUR_DIR).joinpath('test_files/fasta.fa'))


@pytest.fixture
def record_data():
    return {
        'label': 'testRecord',
        'locus': 'testLocus',
        'definition': 'testDef',
        'misc': 'This locus can fit 2 cats'
    }


@pytest.fixture
def vr_def_row():
    df = pd.read_csv(
        'scripts/gibson_assembly/gibson_assembly/test_files/initiation_plamids_short.RC.tsv',
        sep='\t')
    for i, row in df.iterrows():
        return row


@pytest.fixture
def vr_genbank_file(fasta_record, vr_def_row, genbank_record):
    return convert_variable_region_fasta_to_genbank(
        fasta_record, vr_def_row, genbank_record, 
    )


def test_convert_variable_region_fasta_to_genbank(vr_genbank_file):
    assert os.path.isfile(vr_genbank_file)
    content = open(vr_genbank_file).read()
    assert len(content) > 1


@pytest.fixture
def ftgb(fasta_record):
    return fastaToGenbank(fasta_record)


@pytest.fixture
def ftgb_data(fasta_record, record_data):
    return fastaToGenbank(fasta_record, record_data)


@pytest.fixture
def feature_label():
    return 'TEST_LABEL 1'

@pytest.fixture
def ftgb_label(ftgb_data, feature_label):
    ftgb_data.add_label_feature(feature_label)
    return ftgb_data


def test_fastaToGenBank_init(ftgb):
    assert ftgb
    assert isinstance(ftgb, fastaToGenbank)
    assert ftgb.path


def test_fastaToGenBank_init_data(ftgb_data, record_data):
    assert ftgb_data
    assert isinstance(ftgb_data, fastaToGenbank)
    assert ftgb_data.path
    assert ftgb_data.data == record_data


def test_fastaToGenBank_init_data_properties(ftgb_data, record_data):
    
    assert ftgb_data.label == record_data['label']
    assert ftgb_data.locus == record_data['locus']
    assert ftgb_data.definition == record_data['definition']


def test_write_record(ftgb, genbank_record):
    ftgb.write_record(genbank_record)
    assert os.path.isfile(genbank_record)
    os.remove(genbank_record)


def test_write_record_data(ftgb_data, genbank_record):
    ftgb_data.write_record(genbank_record)
    assert os.path.isfile(genbank_record)
    os.remove(genbank_record)


def test_add_label_feature(ftgb_label, genbank_record, feature_label):
    assert isinstance(ftgb_label, fastaToGenbank)
    assert len(ftgb_label.record.features) > 0
    ftgb_label.write_record(genbank_record)
    assert os.path.isfile(genbank_record)
    gb_text = open(genbank_record).read()
    assert feature_label in gb_text  # dirty way to make sure feature was written





    




