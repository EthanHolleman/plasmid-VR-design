from Bio import SeqIO
from Bio.Restriction import *
from Bio.Seq import Seq
from pydna.amplicon import Amplicon
from pydna.assembly import Assembly
from pydna.design import assembly_fragments, primer_design
from pydna.dseqrecord import Dseqrecord
from pydna.readers import read
from pydna.common_sub_strings import terminal_overlap
from pydna.genbankrecord import GenbankRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
import datetime
from pydna.amplify import pcr
from pathlib import Path

MOD_DATE = datetime.date.strftime(datetime.datetime.now(), "%d-%b-%Y").upper()

def read_genbank_record(filepath):
    return GenbankRecord(
        read(filepath)
    )


class Series():

    linearizers = [KpnI, EcoRI]

    @staticmethod
    def select_largest_fragment(frag_list):
        # returns largest fragment by length from a collection
        return max(
            frag_list, key=lambda f: len(f)
        )

    def __init__(self, name, backbone, insert):
        self.name = name
        self.backbone = backbone
        self.insert = insert
    
    @property
    def backbone(self):
        return self._backbone
    
    @backbone.setter
    def backbone(self, new_bone):
        if isinstance(new_bone, str):
            # assume it is a filepath
            self._backbone = read_genbank_record(new_bone)
        elif isinstance(new_bone, Dseqrecord):
            self._backbone = new_bone
        else:
            raise TypeError('Not a filepath or a Dseqrecord!')

    @property
    def linear_backbone(self):
        class_ = type(self)
        cut = self.backbone.cut(class_.linearizers)
        return cut

    @property
    def large_fragment(self):
        return Series.select_largest_fragment(self.linear_backbone)

    @property
    def processed_insert(self):
        return self.insert

    def assemble_circular(self):
        a = Assembly(
            [self.large_fragment, self.processed_insert], limit=15
        ).assemble_circular()
        if len(a) < 1:
            raise Exception('No complete constructs created!')
        else:
            return a[0]

    def write_assembly(self, output_path, name, version, author='your mom'):
        assembly = self.assemble_circular()
        record = GenbankRecord(assembly)

        record.locus = name
        record.stamp()
        record.id = record.seq.seguid()

        record.annotations['data_file_division'] = 'SYN'
        record.annotations['date'] = MOD_DATE
        record.annotations['sequence_version'] = version
        record.annotations['author'] = author

        record.write(str(output_path))



class T7InitiationSeries(Series):

    @classmethod
    def write_assemblies_from_insert_list(cls, insert_paths, backbone_path, 
                                         output_dir):
        output_files = []
        Path(output_dir).mkdir(exist_ok=True, parents=True)
        for each_insert in insert_paths:
            insert_record = read_genbank_record(each_insert)
            insert_name = Path(each_insert).stem

            backbone_record = read_genbank_record(backbone_path)
            backbone_name = Path(backbone_path).stem

            series_name = f'T7_init_{insert_name}'

            series = cls(series_name, backbone_record, insert_record)

            output_path = Path(output_dir).joinpath(
                Path(series_name).with_suffix('.gb')
                )
            output_files.append(output_path)
            series.write_assembly(output_path, series_name, version='1.0')
        
        return output_files

    def __init__(self, name, backbone, insert):
        super().__init__(name, backbone, insert)


class T7TerminationSeries(Series):

    # linearizers for strong initiator insert
    si_linearizers = [HindIII, EcoRI]

    @classmethod
    def write_assemblies_from_insert_list(cls, insert_paths, backbone_path, 
                                        initiator, output_dir):
        output_files = []
        Path(output_dir).mkdir(exist_ok=True, parents=True)
        for each_insert in insert_paths:
            insert_record = read_genbank_record(each_insert)
            insert_name = Path(each_insert).stem

            backbone_record = read_genbank_record(backbone_path)
            backbone_name = Path(backbone_path).stem

            series_name = f'T7_term_{insert_name}'

            series = cls(series_name, backbone_record, insert_record, initiator)

            output_path = Path(output_dir).joinpath(
                Path(series_name).with_suffix('.gb')
                )

            series.write_assembly(output_path, series_name, version='1.0')
            output_files.append(output_path)
        
        return output_files

    def __init__(self, name, backbone, insert, initiator):
        super().__init__(name, backbone, insert)
        self.initiator = initiator

    @property
    def digested_initiator(self):

        di = Series.select_largest_fragment(
            self.initiator.cut(T7TerminationSeries.si_linearizers)
        )
        di.write('digested_record.gb')
        return di
    
    @property
    def initiator(self):
        return self._initiator
    
    @initiator.setter
    def initiator(self, new_init):
        if isinstance(new_init, Dseqrecord):
            self._initiator = new_init
        elif isinstance(new_init, str):
            self._initiator = read_genbank_record(new_init)
        else:
            raise TypeError('Not a filepath or a Dseqrecord!', type(new_init))


    @property
    def backbone_initiator_digest_lf(self):

        # linearize backbone with HindIII and EcoRI to insert strong init
        lf = Series.select_largest_fragment(
            self.backbone.cut(T7TerminationSeries.si_linearizers)
        )
        return lf

    @property
    def backbone_with_initiator(self):
        # see notes/Debug_simulate_assembly.ipynb
        a = Assembly(
            [self.backbone_initiator_digest_lf, self.digested_initiator],
            limit=4, algorithm=terminal_overlap
        ).assemble_circular()
        if len(a) < 1:
            raise Exception('No complete constructs created!')
        else:
            return a[1]

    @property
    def large_fragment(self):
        # want the large fragment produced after insertion of the strong
        # initiator followed by digest

        return Series.select_largest_fragment(
            self.backbone_with_initiator.cut(Series.linearizers)
        )

    @property
    def processed_insert(self):
        return max(self.insert.cut(Series.linearizers), key=lambda f: len(f))

    def assemble_circular(self):
        self.large_fragment.write('t7_large_fragment.gb')
        self.processed_insert.write('t7_processed_insert.gb')
        a = Assembly(
            [self.large_fragment, self.processed_insert], limit=4,
            algorithm=terminal_overlap
        ).assemble_circular()
        if len(a) < 1:
            raise Exception('No complete constructs created!')
        else:
            return a[1]


class TacSeries(Series):

    linearizers = [HindIII, KpnI]
    series_id='Tac'
          
    @classmethod
    def write_assemblies_from_parent_constructs(cls, parent_construct_paths,
                                        backbone_path, primers, output_dir):
        output_files = []
        Path(output_dir).mkdir(exist_ok=True, parents=True)
        backbone_record = read_genbank_record(backbone_path)
        backbone_name = Path(backbone_path).stem

        for each_parent_path in parent_construct_paths:
            
            parent_construct = read_genbank_record(each_parent_path)
            parent_name = Path(each_parent_path).stem

            series_name = f'{cls.series_id}_{parent_name}'

            output_path = Path(output_dir).joinpath(
                Path(series_name).with_suffix('.gb')
                )
            series = cls(series_name, parent_construct, backbone_record, primers)
            # cls.all_primers.append(series.primers)

            series.write_assembly(output_path, series_name, version='1.0')
            output_files.append(output_path)
        
        return output_files


    def __init__(self, name, parent_construct, backbone, primers):
        '''All Tac series constructs (pFC53tacT1T2) are created from existing
        T7 initiation and termination constructs. TacSeries class is the base
        class that specific Tac series (initiation and termination) should
        inherit from.

        Args:
            name (str): Name of series
            parent_construct (str): Path to genbank record of the parent construct
            from which insert is extracted from.
            backbone (str): Path to genbank record describing the backbone into
            which the amplified insert will be inserted into. This really should
            be pFC53T1T2.
        '''
        self.parent_construct = parent_construct
        self.backbone = backbone
        self.primers = primers
        super().__init__(name, backbone, self.pcr_fragment)

    @property
    def parent_construct(self):
        return self._parent_construct
    
    @parent_construct.setter
    def parent_construct(self, new_parent):
        if isinstance(new_parent, GenbankRecord):
            self._parent_construct = new_parent
        elif Path(new_parent).is_file():
            self._parent_construct = GenbankRecord(read(new_parent))
        else:
            raise TypeError
    
    @property
    def pcr_fragment(self):
        '''Defines the region that will be amplified from the parent_construct.
        '''
        pcr_frag = pcr(*self.primers, self.parent_construct)
        print(pcr_frag)
        pcr_frag.write('test_frag.gb')
        return pcr_frag
    
    @property
    def fragment_list(self):
        self.large_fragment.write('test_large_frag.gb')
        return (
            self.large_fragment, self.pcr_fragment
        )
    
    @property
    def primers(self):
        return self._primers
    
    @primers.setter
    def primers(self, new_primers):
        
        if isinstance(new_primers, str) or isinstance(new_primers, Path):
            # is string or Path intrepret as filepath to fasta records
            # of primers
            new_primers = str(new_primers)
            self._primers = list(SeqIO.parse(new_primers, format='fasta'))
            assert len(self._primers) == 2
        else:
            raise TypeError('Primers must be specified as a path to a fasta file')
        
    def assemble_circular(self):
        a = Assembly(self.fragment_list, limit=15).assemble_circular()
        if len(a) < 1:
            raise Exception('No complete constructs created!')
        else:
            return a[0]
    

class TacInitiationSeries(TacSeries):

    series_id='tac_init'

    def __init__(self, name, parent_construct, backbone, primers):
        super().__init__(name, parent_construct, backbone, primers)
    
   
class TacTerminationSeries(TacSeries):

    series_id='tac_term'

    def __init__(self, name, parent_construct, backbone, primers):
        super().__init__(name, parent_construct, backbone, primers)
    
    
def main():

    insert_sequences = snakemake.input['inserts']
    output_dir = str(snakemake.output)
    t7_init_backbone = str(snakemake.input['t7_init_backbone'])
    t7_term_backbone = str(snakemake.input['t7_term_backbone'])
    tac_backbone = str(snakemake.input['tac_backbone'])
    initiator = str(snakemake.input['initiator'])

    initiation_dir = Path(output_dir).joinpath('T7_initiation_series')
    t7_init_assemblies = T7InitiationSeries.write_assemblies_from_insert_list(
        insert_sequences, t7_init_backbone, initiation_dir
    )

    termination_dir = Path(output_dir).joinpath('T7_termination_series')
    t7_term_assemblies = T7TerminationSeries.write_assemblies_from_insert_list(
        insert_sequences, t7_term_backbone, initiator,
        termination_dir
    )

    tac_initiation_dir = Path(output_dir).joinpath('Tac_initiation_series')
    tac_init_primers = str(snakemake.input['tac_init_primers'])
    TacInitiationSeries.write_assemblies_from_parent_constructs(
        t7_init_assemblies, tac_backbone, tac_init_primers, tac_initiation_dir
    )

    tac_termination_dir = Path(output_dir).joinpath('Tac_termination_series')
    tac_term_primers = str(snakemake.input['tac_term_primers'])
    TacTerminationSeries.write_assemblies_from_parent_constructs(
        t7_term_assemblies, tac_backbone, tac_term_primers, tac_termination_dir, 
    )

if __name__ == '__main__':
    main()




