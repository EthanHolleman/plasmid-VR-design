from pydna.assembly import Assembly
from pydna.dseqrecord import Dseqrecord
from pydna.design import primer_design
from pydna.design import assembly_fragments
import time

from pathlib import Path


class GibsonAssembler():

    def __init__(self, construct):
        self.construct = construct
        self.fragment_list = None
        self.assembly = None

    
    def downstream_of_feature(self):
        d = self.construct.backbone.feature_labels_dict

        return d[
            self.construct.downstream_of]

    @property
    def literal_inserts(self):
        # get orientation of downstream_of feature
        inserts = []
        downstream_feature = self.downstream_of_feature()
        if downstream_feature.strand == -1:  # need to take reverse complement of inserts
            for insert in self.construct.inserts:
                inserts.append(insert.reverse_complement())

        return inserts

    def _convert_to_dseqrecord(self, a):
        # helper function to convert objects to Dseqrecord if not one already

        if not isinstance(a, Dseqrecord):
            return Dseqrecord(a)
        return a

    def _create_fragment_list(self, limit):
        # layout for circular assembly
        linear_backbone = self.construct.backbone.linearize(self.construct.downstream_of)
        assembly_layout = [linear_backbone] + \
            self.literal_inserts + [linear_backbone]
        fragment_list = []
        last_added_index = None  # keeps track last index added to frag list
        for i in range(len(assembly_layout)-1):
            # identify existing homology between adjacent fragments
            
            a, b = self._convert_to_dseqrecord(
                assembly_layout[i]), self._convert_to_dseqrecord(assembly_layout[i+1])
            if self._identify_existing_homology(a, b, limit):
                additional_frags = [a, b]
            else:
                # no homology so need to create primers
                additional_frags = assembly_fragments((a, primer_design(b)))

            # check if current index has already been added
            if last_added_index == i:
                fragment_list.append(additional_frags[-1])
            else:
                fragment_list += additional_frags

            last_added_index = i+1

        self.fragment_list = fragment_list
        return self.fragment_list

    def create_assembly(self, limit=15):
        fragment_list = self._create_fragment_list(limit=limit)

        self.assembly = Assembly(fragment_list)
        return self.assembly

    def _write_primers(self, output_path):
        amplicons = []
        for frag in self.fragment_list:
            if isinstance(frag, Amplicon):
                amplicons.append(frag)

        primers = [(y.forward_primer, y.reverse_primer) for y in amplicons]

        with open(str(output_path), 'w') as handle:
            for forward, reverse in primers:
                handle.write(forward.format('fasta'))
                handle.write(reverse.format('fasta'))
        return output_path

    def _write_ascii_fig(self, output_path):
        with open(str(output_path), 'w') as handle:
            output_path.write(self.assembly[0].figure())

    def _write_ascii_fig_detailed(self, output_path):
        with open(str(output_path), 'w') as handle:
            output_path.write(self.assembly[0]).detailed_figure()

    def _write_assembly_gb(self, output_path):
        self.assembly[0].write(str(output_path))

    def _write_literal_inserts(self, output_dir):
        for insert in self.literal_inserts:
            insert_path = Path(output_dir).with_name(Path(insert.path).name)
            insert.write(str(insert_path))

    def write_assembly(self, output_dir, limit=15):
        if not self.assembly:
            self.create_assembly(limit=limit)
        output_dir = Path(output_dir)

        primer_path = str(output_dir.joinpath('gibson_primers.fasta'))
        assembly_ascii_detailed = str(
            output_dir.joinpath('detailed_figure.txt'))
        assembly_ascii = str(output_dir.joinpath('figure.txt'))
        expect_assembly = str(output_dir.joinpath('expected_construct.gb'))
        literal_insert_dir = str(output_dir.joinpath('literal_inserts'))

        self._write_ascii_fig(assembly_ascii)
        self._write_ascii_fig_detailed(assembly_ascii_detailed)
        self._write_primers(primer_path)
        self._write_assembly_gb(expect_assembly)
        self._write_literal_inserts(literal_insert_dir)

        return output_dir

    def _identify_existing_homology(self, a, b, limit):
        test_assembly = Assembly((a, b))
        linear_assembly = test_assembly.assemble_linear()

        if linear_assembly:
            return True
        else:
            return False

    # this object should take over all assembly functions
    # also should implement reverse complmentation of inserted
    # regions if required as that is currently not implemented
    # additionally should search for regions of existing homology
    # before converting seqrecords to amplicon instances.
