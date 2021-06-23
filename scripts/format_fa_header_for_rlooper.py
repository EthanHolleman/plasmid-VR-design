
HEADER_TEMPLATE = ">{} range=FIXED:1-{} 5'pad=0 3'pad=0 strand=+ repeatMasking=none\n"


def read_single_record_fasta(filepath):
    with open(filepath) as handle:
        return handle.readlines()


def format_header(fasta_lines, record_name):
    seq = (''.join([line.strip() for line in fasta_lines[1:]]))
    length = len(seq)
    formated_header = HEADER_TEMPLATE.format(record_name, length)
    return [formated_header, seq]

    return fasta_lines


def write_formated_fasta(formated_fasta_lines, output_path):
    with open(output_path, 'w') as handle:
        for line in formated_fasta_lines:
            handle.write(line)

    return output_path


def main():
    input_fasta = str(snakemake.input)
    output_fasta = str(snakemake.output)
    record_name = str(snakemake.params['record_name'])

    lines = read_single_record_fasta(input_fasta)
    formated_lines = format_header(lines, record_name)
    write_formated_fasta(formated_lines, output_fasta)


if __name__ == '__main__':
    main()