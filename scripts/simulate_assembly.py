from pydna.readers import readsembly
from pydna.design import primer_design
from pydna.design import assembly_fragments
from pydna.dseqrecord import Dseqrecord
from Bio.Restriction import *
from pydna.amplicon import Amplicon
from Bio.Seq import Seq
from Bio import SeqIO

class Series():

    def __init__(self, name, backbone, insert, linearizers):
        self.name = name
        self.backbone = backbone
        self.insert = insert
        self.linearizers
    
    @property
    def linear_backbone(self):
        self.backbone.cut(self.linearizers)
    
    @property
    def large_fragment(self):
        return max(self.linear_backbone, key=lambda s: len(s))
    
    @property
    def processed_insert(self):
        return self.insert

    def assemble_circular(self):
        return Assembly([self.large_fragment, self.processed_insert])

    

class T7InitiationSeries(Series):

    def __init__(self, name, backbone, insert, linearizers):
        super().__init__(name, backbone, insert, linearizers)
    
    
    def assemble(self):
        # Step 1: digest and linearize backbone, purify large
        # fragment
        lf = self.large_fragment

        # Step 2: add inserts via Gibson assembly


    
    
    
    