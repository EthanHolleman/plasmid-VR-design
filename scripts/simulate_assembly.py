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
import datetime
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
        cut = self.backbone.cut(Series.linearizers)
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
            return a[1]

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

            series_name = f'T7_term_{backbone_name}_{insert_name}'

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
        return di
    
    @property
    def initiator(self):
        return self._initiator
    
    @initiator.setter
    def initiator(self, new_init):
        if isinstance(new_init, str):
            # assume it is a filepath
            self._initiator = read_genbank_record(new_init)
        elif isinstance(new_init, Dseqrecord):
            self._initiator = new_init
        else:
            raise TypeError('Not a filepath or a Dseqrecord!')


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
        a = Assembly(
            [self.large_fragment, self.processed_insert], limit=4,
            algorithm=terminal_overlap
        ).assemble_circular()
        if len(a) < 1:
            raise Exception('No complete constructs created!')
        else:
            return a[1]
    

def main():

    insert_sequences = snakemake.input['inserts']
    output_dir = str(snakemake.output)
    t7_init_backbone = str(snakemake.input['t7_init_backbone'])
    t7_term_backbone = str(snakemake.input['t7_term_backbone'])
    initiator = str(snakemake.input['initiator'])

    initiation_dir = Path(output_dir).joinpath('T7_initiation_series')
    T7InitiationSeries.write_assemblies_from_insert_list(
        insert_sequences, t7_init_backbone, initiation_dir
    )

    termination_dir = Path(output_dir).joinpath('T7_termination_series')
    T7TerminationSeries.write_assemblies_from_insert_list(
        insert_sequences, t7_term_backbone, initiator,
        termination_dir
    )

if __name__ == '__main__':
    main()




