from Bio import SeqIO
import pandas as pd

def main():

    fasta = str(snakemake.input)
    seqs = SeqIO.parse(fasta, 'fasta')
    records = []
    for i, each_seq in enumerate(seqs):
        records.append(
            {
            'Number': i+1,
            'Name': each_seq.description,
            'Sequence': str(each_seq.seq)
            }
        )
    table = pd.DataFrame(records)
    table.to_csv(str(snakemake.output), index=False, sep='\t')


if __name__ == '__main__':
    main()
    
