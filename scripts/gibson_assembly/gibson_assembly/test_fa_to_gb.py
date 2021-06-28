import pytest
from pathlib import Path
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





    




