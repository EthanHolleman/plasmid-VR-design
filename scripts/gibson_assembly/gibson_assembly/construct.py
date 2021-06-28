import os
from pandas.core.reshape.tile import cut
import yaml
import copy
from pathlib import Path

from gibson_assembly import INSERT_KEYWORDS

from pydna.genbankfile import GenbankFile
from pydna.readers import read
from pydna.design import primer_design
from pydna.design import assembly_fragments
from pydna.assembly import Assembly
from Bio.Restriction import Analysis, RestrictionBatch
from pydna.amplicon import Amplicon
from Bio import SeqIO


class Construct():

    @classmethod
    def init_from_yaml(cls, yml_path):
        with open(yml_path, 'r') as stream:
            try:
                yaml_dict = yaml.safe_load(stream)
            except yaml.YAMLError as exc:
                raise exc

        constructs = {}
        for construct, attributes in yaml_dict.items():
            constructs[construct] = cls(
                name=construct,
                backbone=attributes['backbone'],
                downstream_of=attributes['insert_downstream_of'],
                inserts=attributes['contents']
            )
        return constructs

    def __init__(self, name, backbone, downstream_of, inserts):
        self.name = name
        self.backbone = backbone
        self.downstream_of = downstream_of
        self.inserts = inserts

    @property
    def backbone(self):
        return self._backbone

    @backbone.setter
    def backbone(self, new_bone):
        assert os.path.isfile(new_bone)
        self._backbone = Backbone(new_bone)

    @property
    def downstream_of(self):
        return self._downstream_of

    @downstream_of.setter
    def downstream_of(self, new_label):
        self._downstream_of = new_label

    @property
    def inserts(self):
        return self._inserts

    @inserts.setter
    def inserts(self, new_inserts):
        insert_instances = []
        for each_insert in new_inserts:
            if not isinstance(each_insert, GenbankFile) and os.path.isfile(each_insert):
                # attempt to read filepaths as genbank entries, these represent
                # constant regions. Otherwises should be a keyword represeting
                # type of insert region.
                print('EACH INSERT')
                print(each_insert)
                each_insert = read(each_insert)
            else:
                assert each_insert in INSERT_KEYWORDS

            insert_instances.append(
                each_insert
            )
        self._inserts = insert_instances

    def write_assembly(self, output_dir):
        assembly = self.backbone.insert_fragments(
            self.inserts, self.downstream_of)

        loci = '_'.join(set([self.backbone.genbank.locus] +
                        [insert.locus for insert in self.inserts]))

        assembly['candidate'].locus = loci

        write_assembly_dir(output_dir, assembly)
        return output_dir

    def specify_variable_region(self, vr_genbank_path):
        assert os.path.isfile(vr_genbank_path)
        vr_genbank = read(vr_genbank_path)

        if INSERT_KEYWORDS[0] in self.inserts:  # variable region keyword
            var_index = self.inserts.index(INSERT_KEYWORDS[0])
            var_construct = copy.deepcopy(self)
            var_construct.inserts[var_index] = vr_genbank
            # deepcopy to avoid retaining reference to nested
            # objects

            return var_construct
        else:
            raise TypeError(
                'This construct does not contain a variable region!')


class Backbone():

    def __init__(self, filepath):
        self.filepath = filepath
        self.genbank = read(filepath)

    @property
    def feature_labels_dict(self):
        # labels of all features
        features_by_labels = {}
        for feature in self.genbank.features:
            label = feature.qualifiers['label']
            if label:
                assert len(label) == 1
                label = label.pop()
                features_by_labels[label] = feature
        return features_by_labels

    @property
    def feature_labels(self):
        feat_labs = set([])
        for feature in self.genbank.features:
            label = feature.qualifiers['label']
            if label:
                assert len(label) == 1
                label = label.pop()
                feat_labs.add(label)
        return feat_labs

    def get_closest_downstream_unique_RS(self, feature_label):
        seqFeature = self.feature_labels_dict[feature_label]
        rb = RestrictionBatch([], ['C'])
        cut_analysis = Analysis(rb, self.genbank.seq, linear=False)
        unique_cutters = {
            enzyme: cut_loc for enzyme, cut_loc in cut_analysis.full().items()
            if len(cut_loc) == 1
        }
        def unique_cutters_dist_func(
            loc): return loc[0] - int(seqFeature.location.start)
        distances = {enzyme: unique_cutters_dist_func(
            loc) for enzyme, loc in unique_cutters.items()}
        distances = {enzyme: dist for enzyme,
                     dist in distances.items() if dist * seqFeature.strand > 0}

        # get enzyme with RS min distance from seqFeature
        best_enzyme = min(distances, key=lambda e: abs(distances.get(e)))
        cut_details = {
            'distance': distances[best_enzyme],
            'cut_site': unique_cutters[best_enzyme].pop()
        }
        return best_enzyme, cut_details, seqFeature

    def insert_fragments(self, inserts, insert_downstream_of):
        cutter = self.get_closest_downstream_unique_RS(insert_downstream_of)
        enzyme = cutter[0]
        linear = self.genbank.linearize(enzyme)
        amplicons = self._inserts_to_amplicons(inserts)
        frag_list = assembly_fragments(
            [linear] + amplicons + [linear]
        )
        primers = [(y.forward_primer, y.reverse_primer) for y in amplicons]
        assembly_final = Assembly(frag_list[:-1])
        candidate = assembly_final.assemble_circular()[0]
        return {
            'assembly': assembly_final,
            'fragments': frag_list,
            'primers': primers,
            'candidate': candidate
        }

    def _inserts_to_amplicons(self, inserts):
        return [primer_design(each_insert) for each_insert in inserts]


def write_primer_list(primers, output_path):
    with open(str(output_path), 'w') as handle:
        for pair in primers:
            handle.write(pair[0].format("fasta"))
            handle.write(pair[1].format("fasta"))
    return output_path


def write_assembly_dir(output_dir, assembly_dict):
    output_dir = Path(output_dir)
    primes_path = str(output_dir.joinpath('gibson_primers.fasta'))
    expect_construct = str(output_dir.joinpath('expected_construct.gb'))
    write_primer_list(assembly_dict['primers'], primes_path)
    assembly_dict['candidate'].write(expect_construct)

    assert os.path.isfile(primes_path)
    assert os.path.isfile(expect_construct)

    return expect_construct, primes_path
