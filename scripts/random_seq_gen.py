import numpy as np

nucs = ['A', 'T', 'G', 'C']

HEADER_TEMPLATE = "RANDOM:{} range=FIXED:1-{} 5'pad=0 3'pad=0 strand=+ repeatMasking=none"

def generate_random_fasta_record(record_name, length):
    seq = ''.join(list(np.random.choice(nucs, length, replace=True)))
    header = HEADER_TEMPLATE.format(record_name, length)
    return f'>{header}\n{seq}'


def main():
    output = str(snakemake.output)
    record_name = str(snakemake.params['record_name'])
    length = int(snakemake.params['length'])
    with open(output, 'w') as handle:
        rand_record = generate_random_fasta_record(record_name, length)
        handle.write(rand_record)


if __name__ == '__main__':
    main()

