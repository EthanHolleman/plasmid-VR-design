from Bio import SeqIO
import os

from pathlib import Path


template_seq_id = "{} range=FIXED:1-{} 5'pad=0 3'pad=0 strand=+ repeatMasking=none"


def main():
    input_fasta = str(snakemake.input['fasta'])
    output_fasta = str(snakemake.output)

    record_name = Path(output_fasta).stem
    # due to snakemake weirdness and me not being able to
    # find a way around not being able to pass a lambda
    # function to output we do this one fasta record at a time
    # so the snakemake output can be a single file instead
    # of a directory
    records = SeqIO.parse(input_fasta, 'fasta')
    for r in records:
        name = r.id.split('_')[1]  # name of seq specified in config tsv
        if name == record_name:
            seq_len = len(r)
            r.description = format_seq_id_for_rlooper(r.id, seq_len)
            with open(output_fasta, 'w') as handle:
                SeqIO.write([r], handle, 'fasta')


def format_seq_id_for_rlooper(seq_id, seq_len):
    return template_seq_id.format(
        seq_id, seq_len
    )


if __name__ == '__main__':
    main()
