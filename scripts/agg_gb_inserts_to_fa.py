from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# script for converting a collection of genbank formated records to a single
# large fasta file


def read_genbank_file(filepath):
    with open(filepath) as handle:
        return list(SeqIO.parse(handle, 'genbank'))


def write_fasta(records, filepath):
    with open(filepath, 'w') as handle:
        SeqIO.write(records, handle, 'fasta')


def convert_genbank_records_to_fa(genbank_filepaths, output_fasta):
    records = []
    for each_gf in genbank_filepaths:
        records += read_genbank_file(each_gf)
    write_fasta(records, output_fasta)
    return output_fasta


def main():
    genbank_filepaths = snakemake.input['inserts']
    output_fasta = str(snakemake.output)
    convert_genbank_records_to_fa(genbank_filepaths, output_fasta)


if __name__ == '__main__':
    main()
