import pytest
from pathlib import Path
import os
import sys

CUR_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

sys.path.append(CUR_DIR)

from test_fa_to_gb import *

@pytest.fixture
def genbank_record():
    return str(Path(CUR_DIR).join('test_files/genbank.gb')

@pytest.fixture
def fasta_record():
    return str(Path(CUR_DIR).join('test_files/fasta.fa')

@pytest.fixture
def record_data():
    return {
        'label': 'testRecord',
        'locus': 'testLocus',
        'definition': 'testDef',
        'misc': 'This locus can fit 2 cats'
    }

@pytest.fixture
def fTgb(fasta_record):
    return fastaToGenbank(fasta_record)


@pytest.fixture
def fTgb_data(fasta_record, record_data):
    return fastaToGenbank(fasta_record, record_data)


def test_fastaToGenBank_init(fTgb):
    assert fTgb
    assert isinstance(fTgb, fastaToGenbank)
    assert fTgb.path


def test_fastaToGenBank_init_data(fTgb_data, record_data):
    assert fTgb_data
    assert isinstance(fTgb_data, fastaToGenbank)
    assert fTgb.path
    assert fTgb_data.data == record_data


def test_fastaToGenBank_init_data_properties(fTgb_data, record_data):
    
    assert fTgb_data.label == record_data['label']
    assert fTgb_data.locus == record_data['locus']
    assert fTgb_data.definition == record_data['definition']


def test_write_record(ftgb_data, genbank_record):

    ftgb_data.write_record(genbank_record)
    assert os.path.isfile(genbank_record)
    
    




