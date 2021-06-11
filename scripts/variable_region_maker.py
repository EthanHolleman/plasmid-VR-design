
BRUTE_FORCE_CACHE = {}

def brute_force_gc_count_calculator(gc_skew, gc_content, seq_len, closest=True):
    '''Calculate the number of G and C nucleotides to include in a DNA
    sequence given a sequence length, GC skew and content. 

    Args:
        gc_skew (float): Level of GC skew, between 0 and 1.
        gc_content (float): Level of GC content between 0 and 1.
        seq_len (int): Length of sequence in nucleotides.
        closest (bool, optional): If no exact result returns the closest. Defaults to True.
    '''
    pass



class VairableRegion():

    name_prefix = 'VR'  # variable region

    def __init__(self, length, gc_content, gc_skew, at_skew=None, 
                at_content=None, role=None):
        self.gc_skew = gc_skew
        self.gc_content = gc_content
        self.at_skew = at_skew
        self.at_content = at_content
        self.role = role
        self.gc_count = None
        self._generated_sequences = 1  # base 1 for plebs
    
    
    def calculate_g_c_count(self):
        '''Calculate the number of G and C nucleotides that should be included
        in the sequence given the gc skew and content.
        '''
        self.gc_count = gc_count

        return gc_count
    

    def generate_sequence(self):
        if not self.gc_count:
            self.calculate_g_c_count
        
        self._generated_sequences += 1
        
        return name, sequence

    def _sequence_name(self):
        return f'{VairableRegion.name_prefix}_{self.role}_{self._generated_sequences}_{self.gc_skew}_{self.gc_content}_unclustered'

    def __len__(self):
        return self.length
    

    



class ClusteredVariableRegion():
        
        def __init__(self, cluster_length, cluster_nucleotide, length, 
                    gc_content, gc_skew, at_skew=None, at_content=None):
            super().__init__(length, gc_content, gc_skew, at_skew, at_content)
            self.cluster_length = cluster_length
            self.cluster_nucleotide = cluster_nucleotide
        