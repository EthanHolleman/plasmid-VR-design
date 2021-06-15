

from varmids.variable_region_maker import *
from varmids.utils import *



def main():
    # input_var_regions = 'variable_defs/termination_plasmids.tsv'
    # output_fasta = 'test.fasta'
    # output_tsv = 'test.tsv'
    input_var_regions = str(snakemake.params['input_file'])
    output_fasta = str(snakemake.output['fasta'])
    output_tsv = str(snakemake.output['tsv'])


    # Create variable regions instances
    variable_regions_dict = read_variable_region_config_file(input_var_regions)
    var_regions = [
        VariableRegion.init_from_csv_row(row) for row in variable_regions_dict
    ] 

    # Generate sequences from variable regions
    seqs = []
    for vr in var_regions:
        vr.calculate_nucleotide_counts()
        seqs.append(vr.generate_sequence())
    
    # Write sequences to output files
    write_sequence_list_to_output_files(seqs, output_fasta, output_tsv)


if __name__ == '__main__':
    main()


