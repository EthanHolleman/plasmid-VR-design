from scripts.gibson_assembly.gibson_assembly.construct import Construct
from scripts.gibson_assembly.gibson_assembly.fa_to_gb import convert_variable_region_fasta_to_genbank
from gibson_assembly.fa_to_gb import *

def main():
    
    vr_region_fasta = str(snakemake.input['vr_fasta'])
    vr_row_def = snakemake.params['vr_row_def']
    construct_yaml= snakemake.params['constructs']

    vr_region_genbank = str(snakemake.output['vr_genbank'])
    assembly_dir = str(snakemake.output['assembly_dir'])

    vr_genbank_file = convert_variable_region_fasta_to_genbank(
        vr_region_fasta
    )
    constructs = Construct.init_from_yaml(construct_yaml)

    construct = constructs[vr_row_def['construct']]  # get the construct this VR uses
    
    # insert the variable region creating new construct
    vr_construct = construct.specify_variable_region(vr_genbank_file)
    vr_construct.write_assembly()





if __name__ == '__main__':
    main()