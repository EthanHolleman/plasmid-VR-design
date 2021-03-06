

from varmids.variable_region_maker import *
from varmids.utils import *


def main():
    
    p_record = read_variable_region_dataframe(
        snakemake.params['p_record']).pop()

    vr = VariableRegion.init_from_csv_row(p_record)
    excluded_restriction_sites = snakemake.params['excluded_restriction_sites']
    print(excluded_restriction_sites)
    
    for i in range(0, int(snakemake.params.num_cases)):
        output_fasta = str(snakemake.output['fasta'][i])
        output_tsv = str(snakemake.output['tsv'][i])
        seq = RestrictionSiteChecker(
            vr, excluded_restriction_sites
            ).generate_RE_free_sequence(seq_id_num=i+1)  # make base 1
        assert seq.id_num == (i+1)
        write_sequence_list_to_output_files([seq], output_fasta, output_tsv)


if __name__ == '__main__':
    main()


