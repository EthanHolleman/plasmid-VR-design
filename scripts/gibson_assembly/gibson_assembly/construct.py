import os
from pandas.core.reshape.tile import cut
import yaml

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

        constructs = []
        for construct, attributes in yaml_dict.items():
            constructs.append(cls(
                name=construct,
                backbone=attributes['backbone'],
                downstream_of=attributes['insert_downstream_of'],
                inserts=attributes['contents']
                )
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
            insert_instances.append(
                Insert(each_insert)
            )
        self._inserts = insert_instances
    
    def make_assembly(self, variable_region):
        pass


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
        unique_cutters_dist_func = lambda loc: loc[0] - int(seqFeature.location.start)
        distances = {enzyme: unique_cutters_dist_func(loc) for enzyme, loc in unique_cutters.items()}
        distances = {enzyme: dist for enzyme, dist in distances.items() if dist * seqFeature.strand > 0}
    
        # get enzyme with RS min distance from seqFeature
        best_enzyme = min(distances, key=lambda e: abs(distances.get(e)))
        cut_details = {
                    'distance': distances[best_enzyme], 
                    'cut_site': unique_cutters[best_enzyme].pop()
                    }
        print(cut_details)
        return best_enzyme, cut_details, seqFeature
        
    
    def insert_fragments(self, inserts, insert_downstream_of):
        cutter = self.get_closest_downstream_unique_RS(insert_downstream_of)
        # pull cutter out of dict
        linear = self.genbank.linearize(cutter)
        amplicons = self._inserts_to_amplicons(inserts)
        frag_list = assembly_fragments(
            [linear] + amplicons + [linear]
        )
        assembly_final = Assembly(frag_list[:-1])

        return assembly_final

    def _inserts_to_amplicons(inserts):
        return [primer_design(each_insert) for each_insert in inserts]
    

class Insert():

    def __init__(self, *args, **kwargs):
        self.__dict__.update(kwargs)


