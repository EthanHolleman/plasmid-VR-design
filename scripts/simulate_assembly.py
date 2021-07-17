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
    all_primers = []

  
    @classmethod
    def write_primers(cls, output_path):
        # get primers with unique sequences
        primer_ids = {}
        for each_primer_pair in cls.all_primers:
            for each_primer in each_primer_pair:
                primer_id = each_primer.seguid()
                if primer_id in primer_ids:
                    primer_ids[primer_id].append(each_primer)
                else:
                    primer_ids[primer_id] = [each_primer]
        
       
        with open(str(output_path), 'w') as handle:
            for each_primer_id, each_primer_list in primer_ids.items():
                final_primer = each_primer_list[0]
                final_primer.description = f"{final_primer.description.split(' ')[0]} {cls.__name__}"
                final_primer.stamp()
                handle.write(final_primer.format('fasta'))
        
        return output_path

        


    @classmethod
    def write_assemblies_from_parent_constructs(cls, parent_construct_paths,
                                        backbone_path, output_dir):
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
            series = cls(series_name, parent_construct, backbone_record)
            cls.all_primers.append(series.primers)

            series.write_assembly(output_path, series_name, version='1.0')
            output_files.append(output_path)
        
        return output_files


    def __init__(self, name, parent_construct, backbone):
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
        super().__init__(name, backbone, self.pcr_target)

    
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
    def pcr_target(self):
        '''Defines the region that will be amplified from the parent_construct.
        Should be overwritten to provide specificity to the particular series
        the child class represents.

        Returns:
            Dseqrecord: Dseqrecord of the region to be amplified from the parent
            construct.
        '''
        raise NotImplementedError

    @property
    def insert_amplicon(self):
        # Amplicon instance version of target_insert. In theory should not
        # be overwritten
        # limit should be length of anchor region
        return primer_design(self.pcr_target, limit=15)
    
    @property
    def fragment_list(self):
       
        return assembly_fragments((
            self.large_fragment, self.insert_amplicon, self.large_fragment
        ))
    
    @property
    def primers(self):
        amplicons = [x for x in self.fragment_list if isinstance(x, Amplicon)]
        primers = [(y.forward_primer, y.reverse_primer) for y in amplicons]
        assert len(primers) == 1 # should only be one set of primers
        primers = primers[0]
        primers[0].description = self.parent_construct.name + '1'
        primers[1].description = self.parent_construct.name + '2'
        return primers
    
    def assemble_circular(self):
        a = Assembly(self.fragment_list, limit=15).assemble_circular()
        if len(a) < 1:
            raise Exception('No complete constructs created!')
        else:
            return a[0]
    


class TacInitiationSeries(TacSeries):

    series_id='tac_init'
    extension_length = 1000  # number of nucleotides downstream of the end of
    # the actual insert sequence to bring along for the ride into the tac
    # backbone
    all_primers = []

    def __init__(self, name, parent_construct, backbone):
        super().__init__(name, parent_construct, backbone)
    
    def _extension_region(self):
        # get 1000 bp before the 3 prime homology arm for extension region
        three_prime_arm = self.parent_construct.features[5]
        start = three_prime_arm.location.start + len(three_prime_arm)
        end = start + TacInitiationSeries.extension_length
        position = FeatureLocation(
            int(start),
            int(end)
        )
        ext_region = SeqFeature(
            position
        ).extract(self.parent_construct)
        assert len(ext_region) == TacInitiationSeries.extension_length
        ext_region.label = 'Extension region'
        return ext_region

    @property
    def pcr_target(self):
        # extract all features from the anchor region to the 3' homology
        # arm
        # add extenstion region as a feature to the parent insert if not already
        # present
        if len(self.parent_construct.features) == 6:
            three_prime_arm = self.parent_construct.features[5]
            start = three_prime_arm.location.start
            end = start + TacInitiationSeries.extension_length
            feature_names = [self.parent_construct.extract_feature(i).name for i in range(len(self.parent_construct.features))]
            if 'Extension region' not in feature_names:
                self.parent_construct.add_feature(
                start, end,
                    name='Extension region', label='Extension region',
                    strand=self.parent_construct.features[5].strand
                )
        # self.parent_construct.write('temp_parent.gb')
        # self.parent_construct = GenbankRecord(read('temp_parent.gb'))
        # pcr_target = self.parent_construct.extract_feature(4)  # anchor
        # feature
        # for i in [2, 3, 4]:
        #     pcr_target += self.parent_construct.extract_feature(i)
        
        # pcr_target = self._orient_pcr_target(pcr_target)
        # return pcr_target
        feature_target_indicies = [3, 4, 5, 6]
        pcr_target  = self.parent_construct.extract_feature(
            feature_target_indicies.pop(0)
        )
        
        while feature_target_indicies:
            pcr_target += self.parent_construct.extract_feature(
                feature_target_indicies.pop(0)
            )
        #pcr_target += self._extension_region()
        #pcr_target.write('target_before_complement.gb')
        #pcr_target = pcr_target.reverse_complement()
        #pcr_target = self._orient_pcr_target(pcr_target)
        return pcr_target
    
    @property
    def insert_amplicon(self):
        # Amplicon instance version of target_insert. In theory should not
        # be overwritten
        # limit should be length of anchor region
        return primer_design(self.pcr_target, limit=15).reverse_complement()

    
    def _orient_pcr_target(self, pcr_target):
        # in order to make sure primers are designed so that insert
        # has correct orientation relative to the tac promotor get strand of
        # the tac promotor and if it is not the same as the variable region
        # take the reverse complement of the variable region.
        variable_region_strand = pcr_target.features[1].strand
        tac_promotor_strand = self.backbone.features[6].strand  # tac promoter index
        assert self.backbone.features[6].type == 'promoter'
        if variable_region_strand != tac_promotor_strand:
            pcr_target = pcr_target.reverse_complement()
        return pcr_target
    

class TacTerminationSeries(TacSeries):

    series_id='tac_term'
    all_primers = []

    def __init__(self, name, parent_construct, backbone):
        super().__init__(name, parent_construct, backbone)
    
    @property
    def pcr_target(self):
        feature_target_indicies = [3, 4, 5]
        # initiator region, variable region and anchor region
        pcr_target = self.parent_construct.extract_feature(
            feature_target_indicies.pop(0)
        )
        while feature_target_indicies:
            pcr_target += self.parent_construct.extract_feature(
                feature_target_indicies.pop(0)
            )
        #pcr_target = self._orient_pcr_target(pcr_target)
        pcr_target = pcr_target.reverse_complement()
        return pcr_target
    
    def _orient_pcr_target(self, pcr_target):
        variable_region_strand = pcr_target.extract_feature(1)
        tac_promoter_strand = self.backbone.features[6].strand
        assert self.backbone.features[6].type == 'promoter'
        if variable_region_strand == tac_promoter_strand:
            pcr_target = pcr_target.reverse_complement()
        return pcr_target
    
    @property
    def fragment_list(self):
       
        return assembly_fragments((
            self.large_fragment, self.insert_amplicon, self.large_fragment
        ))
    
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

    # tac initiation and termination assemblies are having issues that
    # seem to be due to details in pydna that I am working on figuring
    # out. For now creating primers for these assemblies manually and adding
    # script to confirm that primers will produce assemblable amplicons.

    # tac_initiation_dir = Path(output_dir).joinpath('Tac_initiation_series')
    # tac_init_primers = tac_initiation_dir.joinpath('Tac_initiation_primers.fa')
    # TacInitiationSeries.write_assemblies_from_parent_constructs(
    #     t7_init_assemblies, tac_backbone, tac_initiation_dir
    # )
    # TacInitiationSeries.write_primers(tac_init_primers)

    # tac_termination_dir = Path(output_dir).joinpath('Tac_termination_series')
    # tac_term_primers = tac_termination_dir.joinpath('Tac_terminination_primers.fa')
    # TacTerminationSeries.write_assemblies_from_parent_constructs(
    #     t7_term_assemblies, tac_backbone, tac_termination_dir
    # )
    # TacTerminationSeries.write_primers(tac_term_primers)

if __name__ == '__main__':
    main()




