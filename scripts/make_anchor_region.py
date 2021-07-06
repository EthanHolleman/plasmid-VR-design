from pydna.dseqrecord import Dseqrecord
from pydna.common_sub_strings import common_sub_strings
from pydna.readers import read
from Bio import SeqIO
from Bio.Restriction import *
from Bio.SeqUtils import MeltingTemp as mt
import numpy as np

# read in all sequence files and test randoms
np.random.seed(12311997)

def check_enzymes(user_enzymes):
    # check in make sure all user provided enzymes are in enzymes
    enzymes = set([str(e) for e in list(AllEnzymes)])
    for each_e in user_enzymes:
        if each_e not in enzymes:
            raise TypeError(f'{each_e} is not in list of all restriction enzymes!')


def read_genbank_records(record_list):
    for record in record_list:
        yield read(record)


def random_sequence(length):
    return Dseqrecord(''.join(list(np.random.choice(['A', 'T', 'G', 'C'], length))))


def check_for_common_substrings(seq_x, seq_y, limit=None):
    '''Determines if two records share common substrings longer than
    the limit argument. If limit is None then by default will be set to 3/4
    the length of seq_y. 
    '''
    if not limit:
        limit = int(len(seq_y) / 0.75) 

    seq_x_rc = seq_x.reverse_complement()
    seq_y_rc = seq_y.reverse_complement()

    pairs = [
        (seq_x, seq_y),
        (seq_x, seq_y_rc),
        (seq_x_rc, seq_y),
        (seq_x_rc, seq_y_rc)
    ]
    for seq_x, seq_y in pairs:
        c = common_sub_strings(str(seq_x.seq), str(seq_y.seq), limit=limit)
        if c:
            return True
    return False

    
def check_for_prohibited_restriction_sites(candidate_anchor, prohibited_cutters):
    '''Check to ensure that no prohibited restriction enzyme recognition sites
    are present in the anchor sequence. If a prohibited site is present return
    False, otherwise returns True.

    Args:
        candidate_anchor (Dseqrecord): Candidate anchor sequence to test for 
        recognition sites.
        prohibited_cutters (list): List of names of prohibited restriction enzymes.

    Returns:
        bool: True if no prohibited sites are found, False if any are.
    '''
    no_cutters = set([str(n_cutter) for n_cutter in list(candidate_anchor.no_cutters())])
    # get all non-cutting enzymes as a set of strings
    for each_prohibited_site in prohibited_cutters:
        if not each_prohibited_site in no_cutters:
            return False
        # make sure prohibited enzymes do not cut
    return True


def check_melting_temp(candidate_anchor, min_temp=48):
    '''Check the melting temp using nearest neighbor thermodynamics. If temp is
    below the min_temp return False, otherwise return True. Helps to ensure that
    the returned sequence is not just a bunch of As and Cs.

    Args:
        candidate_anchor (Dseqrecord): Candidate anchor sequence.
        min_temp (int, optional): Min melting temperature. Defaults to 48.

    Returns:
        bool: True if candidate anchor melting temp is above min_temp False otherwise.
    '''
    if mt.Tm_NN(candidate_anchor.seq) < min_temp:
        return False
    else:
        return True
    

def make_anchor_seq(records, anchor_length, prohibited_cutters, 
                    max_attempts=10000, min_melting_temp=48):
    '''Make a suitable anchor sequence that does not share common sub strings
    with any record in records arg, has length of anchor_length, and does
    not have recognition sites for enzymes in prohibited_cutters.

    Args:
        records (list): List of records to check for common sub-strings against.
        anchor_length (int): Length of anchor sequence in.
        prohibited_cutters (list): List of names of prohibited restriction enzymes
        as strings.
        max_attempts (int): Max number of attempts to generate a suitable anchor
        sequence. Throws exception if exceeded in order to avoid endless looping.

    Returns:
        [Dseqrecord]: Dseqrecord of suitable anchor sequence. 
    '''
    assert len(records) > 0
    assert anchor_length > 0
    attempts = 0
    while True:
        print('Number of attempts:', attempts)
        if attempts == max_attempts:
            raise Exception('Max attempts has exceeded allowed value!')
        candidate_anchor = random_sequence(anchor_length)
        sub_str_safe = True
        attempts += 1
        for each_record in records:
            sub_strings = check_for_common_substrings(each_record, candidate_anchor)
            if sub_strings:
                sub_str_safe = False
                break

        if sub_str_safe:  # no common substrings now check for cutters
            additional_checks = [
                check_melting_temp(candidate_anchor, min_melting_temp),
                check_for_prohibited_restriction_sites(candidate_anchor, prohibited_cutters)
            ]
            if False in additional_checks:  # one or more checks failed
                continue
            else:
                return candidate_anchor

def main():
    print('='*10)
    print('Reading snakemake input')

    backbones = snakemake.params['backbones']
    variable_regions = str(snakemake.input['variable_regions'])

    prohibited_cutters = [str(c) for c in snakemake.params['prohibited_cutters']]
    check_enzymes(prohibited_cutters)

    max_attempts = int(snakemake.params['max_attempts'])
    anchor_length = int(snakemake.params['anchor_length'])
    min_melting_temp= int(snakemake.params['min_melting_temp'])

    output_path = str(snakemake.output)

    print('='*10)
    print('Reading records')
    backbone_records = list(read_genbank_records(backbones))
    vr_records = list(SeqIO.parse(variable_regions, 'fasta'))

    all_records = backbone_records + vr_records
    print(all_records)

    print('='*10)
    print('Dropping anchor')
    anchor_seq = make_anchor_seq(
        all_records, anchor_length=15, 
        prohibited_cutters=prohibited_cutters, min_melting_temp=min_melting_temp,
        max_attempts=max_attempts
    
    )
    anchor_seq.id = 'Anchor_sequence'
    anchor_seq.description = ''
    SeqIO.write([anchor_seq], output_path, 'fasta')
    print('Done')


          
if __name__ == '__main__':
    main()