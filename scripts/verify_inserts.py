# Final check for all complete insert sequences
from pathlib import Path
import pandas as pd
from Bio.Restriction import *
from pydna.genbankrecord import GenbankRecord
from pydna.readers import read
from Bio import SeqIO

# helper functions


def read_genbank_file(filepath):
    return GenbankRecord(read(str(filepath)))


def read_fasta_records(filepath):
    return SeqIO.parse(filepath)


def get_feature_by_name(record, name):
    # pull a feature out of a record using its name
    # if not present return -1. Send all names to
    # lower case first so check is not case sensitive
    features = [record.extract_feature(i)
                for i in range(len(record.features))
                ]
    name_dict = {f.name.lower(): f for f in features}
    if name.lower() in name_dict:
        return name_dict[name]
    else:
        return -1


def no_cutters_as_string_set(record):
    return set([str(s) for s in set(record.no_cutters())])


def one_cutters_as_string_set(record):
    return set([str(s) for s in set(record.once_cutters())])

def calculate_content(record, nuc_a, nuc_b):
    a_count = str(record.seq).count(nuc_a)
    b_count = str(record.seq).count(nuc_b)
    return (a_count + b_count) / len(record)

def calculate_skew(record, nuc_a, nuc_b):
    a_count = str(record.seq).lower().count(nuc_a.lower())
    b_count = str(record.seq).lower().count(nuc_b.lower())
    
    return (a_count - b_count) / (a_count + b_count)

# checks


def check_insert_for_sequence(insert_record, seq_record):
    assert insert_record.seq.find(str(seq_record.seq)) != -1


def check_homology_arms(insert_record, homology_arm, source_backbone):
    assert insert_record.seq.find(str(homology_arm.seq)) != -1
    assert source_backbone.seq.find(str(homology_arm.seq)) != -1


def check_homology_arms_for_required_cutter(homology_arm, cutter):
    assert cutter in one_cutters_as_string_set(homology_arm)

def check_variable_region_gc_skew(vr_record, attr_dict):
    gc_skew = calculate_skew(vr_record, 'G', 'C')
    gc_skew_expected = float(attr_dict['gc_skew'])
    assert abs(gc_skew - gc_skew_expected) < 0.03

def check_variable_region_gc_content(vr_record, attr_dict):
    gc_content = calculate_content(vr_record, 'G', 'C')
    gc_content_expected = float(attr_dict['gc_content'])
    assert abs(gc_content - gc_content_expected) < 0.03

def write_passed_stamp(insert_record, output_path):
    with open(str(output_path), 'w') as handle:
        handle.write(insert_record.seguid())


def main():

    # read all snakemake input stuff
    anchor_record = read_genbank_file(snakemake.input['anchor_seq'])
    insert_record = read_genbank_file(snakemake.input['insert_record'])
    five_prime_arm = read_genbank_file(snakemake.input['five_prime_arm'])
    three_prime_arm = read_genbank_file(snakemake.input['three_prime_arm'])
    homology_target_backbone = read_genbank_file(
        snakemake.input['homology_target'])
    attr_dict = snakemake.params['attributes'].to_dict(orient='records')
    assert len(attr_dict) == 1
    attr_dict = attr_dict[0]
    print(attr_dict)
    # check insert to make sure contains expected sequences
    check_insert_for_sequence(insert_record, anchor_record)
    check_insert_for_sequence(insert_record, three_prime_arm)
    check_insert_for_sequence(insert_record, five_prime_arm)

    # check the homology arms to make sure they are actually homologous
    check_homology_arms(insert_record, five_prime_arm,
                        homology_target_backbone)
    check_homology_arms(insert_record, five_prime_arm,
                        homology_target_backbone)
    
    check_homology_arms_for_required_cutter(five_prime_arm, 'KpnI')
    check_homology_arms_for_required_cutter(three_prime_arm, 'EcoRI')

    vr_record = get_feature_by_name(insert_record, 'variable_region')
    check_variable_region_gc_content(vr_record, attr_dict)
    check_variable_region_gc_skew(vr_record, attr_dict)

    write_passed_stamp(insert_record, snakemake.output)

if __name__ == '__main__':
    main()


