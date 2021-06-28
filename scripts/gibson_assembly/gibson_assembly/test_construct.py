import pytest

from gibson_assembly.construct import *
from pydna.genbankfile import GenbankFile


@pytest.fixture
def construct_yaml():
    return 'scripts/gibson_assembly/gibson_assembly/test_files/test_constructs.yml'


def test_init_from_yaml(construct_yaml):
    constructs = Construct.init_from_yaml(construct_yaml)
    assert constructs
    for each_construct in constructs:
        assert isinstance(each_construct, Construct)
        assert each_construct.name
        assert isinstance(each_construct.backbone, GenbankFile)



    
