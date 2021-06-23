# combine seq metric files for each uniquely defined plasmid
import pandas as pd

def read_file(filepath):
    fp = str(filepath)
    table = pd.read_table(fp, sep='\t')

    return table 


def concatenate_files(all_filepaths):
    tables = [read_file(filepath) for filepath in all_filepaths]
    concat_table = pd.concat(tables)

    return concat_table


def main():
    input_files = snakemake.input
    output_file = snakemake.output

    concat_table = concatenate_files(input_files)
    concat_table.to_csv(str(output_file), sep='\t', index=False)


if __name__ == '__main__':
    main()