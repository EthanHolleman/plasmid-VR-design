from scripts.gibson_assembly.gibson_assembly.construct import Construct
from scripts.gibson_assembly.gibson_assembly.fa_to_gb import convert_variable_region_fasta_to_genbank
from gibson_assembly.fa_to_gb import *

def main():
    
    vr_region_fasta = str(snakemake.input['vr_fasta'])
    vr_row_def = snakemake.params['vr_row_def']
    construct_yaml= snakemake.params['constructs']

    vr_region_genbank = str(snakemake.output['vr_genbank'])

    vr_genbank_file = convert_variable_region_fasta_to_genbank(
        vr_region_fasta
    )
    constructs = Construct.init_from_yaml(construct_yaml)
    




if __name__ == '__main__':
    main()