import os
from Bio import SeqIO

TEMPLATE_PATH = os.path.abspath('resources/files/primer3_settings_template.txt')

assert os.path.isfile(TEMPLATE_PATH), 'template file does not exist!'


def read_fasta(fasta_path):
    records = list(SeqIO.parse(fasta_path, "fasta"))
    assert len(fasta_path) == 1, 'Only one record should be passed in fasta!'
    return records.pop()


def format_template(fasta_record):

    param_dict = {
        'sequence': str(fasta_record.seq),
        'seq_id': str(fasta_record.id)
    }

    with open(TEMPLATE_PATH) as handle:
        formated_tempalate = handle.read().format(**param_dict)
        return formated_tempalate


def write_formated_template(formated_template, output_path):
    with open(output_path, 'w') as handle:
        output_path.write(formated_template)
    return output_path


def main():
    fasta = str(snakemake.input['fasta'])
    output_path = str(snakemake.output)
    record = read_fasta(fasta)
    template_str = format_template(record)
    write_formated_template(template_str)


if __name__ == '__main__':
    main()
    


