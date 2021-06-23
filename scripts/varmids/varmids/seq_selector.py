from varmids.utils import *

import pandas as pd

class AttributeExpectation():

    def __init__(self, expectation_path):
        self.expectation_path = expectation_path
        self.expect_table = self._read_expect_tsv()
    
    def _read_expect_tsv(self):
        return pd.read_table(self.expectation_path)
    

    def _access_rlooper(self):
        pass


class SeqSelector():


    def __init__(self, variable_region, global_skew_thres, global_content_thres,
                bprpob_max_dist=None, lae_max_dist=None, prop_unpaired_max_dist=None,
                prop_hairpin_max_dist=None, 
                local_skew_thres, local_content_thres, max_attempts=1000):
        '''When sequence objects are generated there is a fair bit of
    randomness that comes into play in terms of where each nucleotide
    is places. In some cases this leads to extreme local differences
    from specified parameters. SeqSelector instances can be used to
    generate sequences until one is produced that passes the
    thresholds of the SeqSelector instance.

        Args:
            variable_region (VariableRegion): VariableRegion instance from which
            Sequences are generated and tested.
            global_skew_thres (float): Amount of allowed difference between
            global skew (AT and GC) parameters and the actual values calculated
            from the sequence.
            global_content_thres (float): Amount of allowed difference between
            global content (AT and GC) parameters and the actual values
            calculated from the sequence.
            local_skew_thres (float): Amount fo allowed difference between skew
            captured by sliding windows across the sequence and specified
            value.
            local_content_thres (float): Amount fo allowed difference between content
            captured by sliding windows across the sequence and specified
            value.
            max_attempts (int, optional): Max number of sequences to generate.
            This prevents endless loops if impossible thresholds or parameters
            are passed. Defaults to 1000.
        '''
        self.variable_region = variable_region
        self.global_skew_thres = skew_thres
        self.global_content_thres = content_thres
        self.local_skew_thres = local_skew_thres
        self.local_content_thres = local_content_thres
        self.max_attempts = max_attempts
    
    def make_sequence(self):
        pass
    

    def _test_sequence(self, sequence):
        nuc_counts = count_nucleotides_in_seq(sequence.nuc_seq)
        seq_len = len(sequence.nuc_seq)

        

    
    def _test_skew(self, x_count, y_count, parameter_val, thres):
        diff = abs(skew(x_count, y_count) - parameter_val)
        if diff < thres:
            return True
        else:
            return False


    def _test_content(self, x_count, y_count, seq_len, parameter_val, thres):
        diff = abs(content(x_count, y_count, seq_len) - parameter_val)
        if diff < thres:
            return True
        else:
            return False
        



    
    
    