import pytest
import os
import pandas as pd
from pathlib import Path
import shutil

from gibson_assembly.test_fa_to_gb import vr_genbank_file, fasta_record, vr_def_row
from gibson_assembly.test_construct import genbank_record, construct_yaml, constructs, vr_constructs

from gibson_assembly.assembly import GibsonAssembler
from gibson_assembly.construct import Construct, Backbone

@pytest.fixture
def gib_assembly(vr_constructs):
    assemblers = []
    for construct in vr_constructs:
        assemblers.append(GibsonAssembler(construct))
    return assemblers

def test_gibson_assembler_init(gib_assembly):
    for ga in gib_assembly:
        assert isinstance(ga, GibsonAssembler)
        assert isinstance(ga.construct, Construct)
        assert isinstance(ga.construct.backbone, Backbone)


def test_literal_inserts(gib_assembly):
    for ga in gib_assembly:
        for i, lit_insert in enumerate(ga.literal_inserts):
            if ga.downstream_of_feature().strand == -1:
                assert lit_insert == ga.construct.inserts[i].reverse_complement()
            else:
                assert lit_insert == ga.construct.inserts[i]


def test_create_fragment_list(gib_assembly):
    for ga in gib_assembly:
        ga._create_fragment_list(limit=20)
        assert ga.fragment_list
        assert len(ga.fragment_list) == len(ga.inserts) + 2
        # add 2 because circular so include linear backbone at start and end
        # of fragment list