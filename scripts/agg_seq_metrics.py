# Script for aggregating the output of programs that produce metrics
# used to access the quality of a given variable region into a tsv file which
# is then used to rank all sequences.
import statistics
import csv
import pandas as pd
import os

            
def add_inputs_to_out_dict(out_dict, snakemake):
    # helper function to just transfer over some information from the snakemake
    # input to the final output file
    out_dict['fasta'] = snakemake.input['fasta']
    out_dict['tsv'] = snakemake.input['tsv']
    out_dict['length'] = snakemake.params['length']
    out_dict['p_name'] = snakemake.params['p_name']
    out_dict['id_num'] = snakemake.params['id_num']

    return out_dict


def parse_rlooper_wig(filepath):
    # read numeric values from rlooper produced wigfile 
    values = []
    with open(filepath) as handle:
        lines = handle.readlines()[5:]  # skip first 4 lines
    return [float(v.strip()) for v in lines]


def parse_rlooper_files(bp_prob, lae):
    return {
        'local_average_energy': statistics.mean(parse_rlooper_wig(lae)),
        'bp_prob': statistics.mean(parse_rlooper_wig(bp_prob))
    }


def parse_rna_file(filepath):
    # parse RNA secondary structure file
    # columns that each attribute is located at
    NAMES_DICT = {'prop_hairpin': 4, 
                  'prop_unpaired': 1
                }
    d = {}
    with open(filepath) as handle:
        reader = csv.reader(handle, delimiter='\t')
        record = next(reader)
        for attribute, index in NAMES_DICT.items():
            d[attribute] = record[index]
    return d


def write_out_dict_as_tsv(out_dict, output_path):
    with open(output_path, 'w') as handle:
        fieldnames = list(out_dict.keys())
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        writer.writerow(out_dict)
    
    return output_path

    
def main():
    #read input files 
    input_bp_prob = str(snakemake.input['bp_prob'])
    input_lae = str(snakemake.input['lae'])
    input_rna = str(snakemake.input['RNAss'])

    length = int(snakemake.params['length'])  # length of sequence

    # read metric files
    rlooper_metrics = parse_rlooper_files(input_bp_prob, input_lae)
    rna_metrics = parse_rna_file(input_rna)

    # add some identifiers from the input
    out_dict = {}
    out_dict = add_inputs_to_out_dict(out_dict, snakemake)

    # update with metrics
    out_dict.update(rlooper_metrics)
    out_dict.update(rna_metrics)

    # write as a tsv file, will contain a single record.
    output_path = str(snakemake.output)
    write_out_dict_as_tsv(out_dict, output_path)


if __name__ == '__main__':
    main()





        

    


    
    