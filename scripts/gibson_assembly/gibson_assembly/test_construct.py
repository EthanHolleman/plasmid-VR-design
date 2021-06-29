import pytest
import os
import pandas as pd
from pathlib import Path
import shutil

from gibson_assembly.construct import *
from gibson_assembly.test_fa_to_gb import vr_genbank_file, fasta_record, vr_def_row
from gibson_assembly.fa_to_gb import *
from gibson_assembly import INSERT_KEYWORDS
from pydna.genbankfile import GenbankFile
from pydna.assembly import Assembly
from Bio import Restriction


@pytest.fixture
def genbank_record():
    return 'scripts/gibson_assembly/gibson_assembly/test_files/test_construct.gb'


@pytest.fixture
def construct_yaml():
    return 'scripts/gibson_assembly/gibson_assembly/test_files/test_constructs.yml'

@pytest.fixture
def assembly_dir():
    return 'scripts/gibson_assembly/gibson_assembly/test_files/test_assembly'

@pytest.fixture
def constructs(construct_yaml):
    return Construct.init_from_yaml(construct_yaml).values()


@pytest.fixture
def vr_constructs(constructs, vr_genbank_file):
    vr_constructs = []
    for construct in constructs:
        vr_constructs.append(
            construct.specify_variable_region(vr_genbank_file))
    return vr_constructs


def test_init_from_yaml(constructs):
    for each_construct in constructs:
        assert isinstance(each_construct, Construct)
        assert each_construct.name
        assert isinstance(each_construct.backbone, Backbone)
        for each_insert in each_construct.inserts:
            assert isinstance(
                each_insert, GenbankFile) or each_insert in INSERT_KEYWORDS


def test_init_backbone(constructs):
    for construct in constructs:
        backbone = construct.backbone
        assert isinstance(backbone.genbank, GenbankFile)
        assert os.path.isfile(backbone.filepath)


def test_closest_downstream_unique_RS(constructs):
    for construct in constructs:
        backbone = construct.backbone
        cutter = backbone.get_closest_downstream_unique_RS(
            construct.downstream_of
        )
        assert isinstance(cutter, tuple)
        enzyme, cut_details, downstream_of_feat = cutter

        # make sure cut is orientated relative to sense of cut downstream of
        # feature
        if downstream_of_feat.strand == -1:
            assert cut_details['cut_site'] < int(
                downstream_of_feat.location.start)
        else:
            assert cut_details['cut_site'] > int(
                downstream_of_feat.location.start)


def test_specify_variable_region(vr_genbank_file, constructs):
    for construct in constructs:
        vr_construct = construct.specify_variable_region(vr_genbank_file)
        assert vr_construct.inserts != construct.inserts
        if INSERT_KEYWORDS[0] in construct.inserts:
            assert INSERT_KEYWORDS[0] not in vr_construct.inserts
            insert_index = construct.inserts.index(INSERT_KEYWORDS[0])
            assert isinstance(vr_construct.inserts[insert_index], GenbankFile)


# def test_insert_fragments(vr_genbank_file, constructs):
#     for construct in constructs:
#         vr_construct = construct.specify_variable_region(vr_genbank_file)
#         backbone = vr_construct.backbone

#         assembly = backbone.insert_fragments(
#             vr_construct.inserts, vr_construct.downstream_of)
#         assert isinstance(assembly['assembly'], Assembly)
#         assert isinstance(assembly, dict)
#         assert len(assembly['fragments']) == 2 + len(vr_construct.inserts)


# def test_write_assembly(vr_constructs, assembly_dir):

#     if os.path.isdir(assembly_dir):
#         shutil.rmtree(str(assembly_dir))
#     os.mkdir(assembly_dir)
#     for each_construct in vr_constructs:
#         each_construct.write_assembly(str(assembly_dir))
#         files = []
#         for each_file in Path(assembly_dir).iterdir():
#             files.append(each_file)
        
#         assert len(files) > 1, files
    


