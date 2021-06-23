# One of the options in the variable region definition tsv file is to
# calculate sequences for the reverse complement. While previous seq
# generation was aware of this parameter, currently it is ignored as snakemake
# expects a 1 input -> 1 output relationship. 

# To get around this, this script just adds additional rows representing the
# reverse complement sequences to a temporary version of the user input 

from pathlib import Path
import pandas as pd
import argparse

def get_args():

    parser = argparse.ArgumentParser()
    parser.add_argument('input', help='Path to user specified variable region definition tsv file.')
    parser.add_argument('output', help='Path to output file.')

    return parser.parse_args()


def read_var_def(filepath):
    return pd.read_table(filepath, sep='\t')


def swap_sign(value):
    # swaps the sign of a value (- or +) of a value if it is a valid value
    try:
        value = float(value)
        return value * -1
    except Exception:
        return value


def reverse_complement(row):
    # content metrics are the same but skew takes opposite sign
    temp_row = row
    temp_row['gc_skew'] = swap_sign(row['gc_skew'])
    temp_row['at_skew'] = swap_sign(row['at_skew'])
    temp_row['name'] = f"{row['name']}-RC"

    return temp_row


def apply_reverse_complement(table):
    rows = []
    for index, row in table.iterrows():
        if row['reverse_complement'] == 1:
            rows.append(reverse_complement(row))
    for row in rows:
        table = table.append(row)
    return table


def add_reverse_complement_definitions(input_file):
    table = read_var_def(input_file)
    table_RC = apply_reverse_complement(table)

    output_file = str(Path(input_file).with_suffix('.RC.tsv'))

    table_RC.to_csv(output_file, sep='\t', index=False)

    return output_file


def main():
    args = get_args()
    table = read_var_def(args.input)
    table_RC = apply_reverse_complement(table)
    table_RC.to_csv(args.output, sep='\t', index=False)


if __name__ == '__main__':
    main()




