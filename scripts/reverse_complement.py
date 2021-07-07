from Bio import SeqIO
import pandas as pd
import os


def get_RC_description(record):
    # take a Seq record and return its description with RC added to name
    # to distinguish reverse complemented seq
    split = record.description.split('_')
    split[1] = f'RC-{split[1]}'
    description = '_'.join(split)
    return description


def RC_record(record):
    # create the revese complement version of a record
    rc_record = record.reverse_complement()
    rc_record.description = get_RC_description(record)
    rc_record.id=''
    return rc_record


def update_tsv_record(table, RC_record):
    # Update the original record's table with values for its reverse complement
    # should only be one record
    table.at[0, 'description'] = RC_record.description
    table.at[0, 'GC_skew'] = table.at[0, 'GC_skew'] * -1
    table.at[0, 'AT_skew'] = table.at[0, 'AT_skew'] * -1
    table.at[0, 'Sequence'] = str(RC_record.seq)
    table.at[0, 'name'] = f'RC-{table.at[0, "name"]}'
    return table


def main():
    input_fasta = str(snakemake.input['fasta'])
    input_tsv = str(snakemake.input['tsv'])
    output_fasta = str(snakemake.output['fasta'])
    output_tsv = str(snakemake.output['tsv'])

    table = pd.read_table(input_tsv)
    record = SeqIO.read(input_fasta, 'fasta')

    rc_record = RC_record(record)
    rc_table = update_tsv_record(table, rc_record)

    SeqIO.write([rc_record], output_fasta, 'fasta')
    assert os.path.isfile(output_fasta)
    rc_table.to_csv(output_tsv, sep='\t', index=False)
    assert os.path.isfile(output_tsv)


if __name__ == '__main__':
    main()




