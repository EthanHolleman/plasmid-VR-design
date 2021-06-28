import pytest
import os

from gibson_assembly.construct import *
from pydna.genbankfile import GenbankFile
from Bio import Restriction


@pytest.fixture
def construct_yaml():
    return 'scripts/gibson_assembly/gibson_assembly/test_files/test_constructs.yml'

@pytest.fixture
def constructs(construct_yaml):
    return Construct.init_from_yaml(construct_yaml)


def test_init_from_yaml(construct_yaml):
    constructs = Construct.init_from_yaml(construct_yaml)
    assert constructs
    for each_construct in constructs:
        assert isinstance(each_construct, Construct)
        assert each_construct.name
        assert isinstance(each_construct.backbone, Backbone)


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
            assert cut_details['cut_site'] < int(downstream_of_feat.location.start)
        else:
            assert cut_details['cut_site'] > int(downstream_of_feat.location.start)