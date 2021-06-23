import pandas as pd

def format_template(template_path, **kwargs):
    with open(template_path) as handle:
        text = handle.read()
        return text.format(**kwargs)


def write_formated_text(formated_text, output_path):
    with open(output_path, 'w') as handle:
        handle.write(formated_text)
    return output_path 


def extract_format_fields_from_tsv(tsv_path):
    # should be a single row table (1 sequence)
    row = pd.read_table(tsv_path, sep='\t').to_dict(orient='records').pop()
    return {
        'seq_id': f"{row['name']}.{row['id_num']}",
        'sequence': row['Sequence']
    }


def main():
    input_tsv = str(snakemake.input['tsv'])
    template = str(snakemake.input['template'])
    output_path = str(snakemake.output)

    format_dict = extract_format_fields_from_tsv(input_tsv)
    formated_text = format_template(template, **format_dict)
    write_formated_text(formated_text, output_path)


if __name__ == '__main__':
    main()




    