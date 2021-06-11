import numpy as np

from utils import find_available_random_range 

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


class Sequence():
    
    def __init__(self, count_dict, cluster_length=None, cluster_nuc=None,
                cluster_dist_func=None):
        self.count_dict = count_dict
        self.cluster_length = cluster_length
        self.cluster_nuc = cluster_nuc
        self.cluster_dist_func = cluster_dist_func
        self.nuc_seq = self._make_sequence()

        if self.cluster_length:
            assert cluster_nuc, 'Must specify cluster nuc if passing cluster len!'
            assert cluster_nuc in count_dict
    

    @cluster_dict_func.setter
    def cluster_dict_func(self, new_func):
        if callable(new_func):
            self._cluster_dict_func = new_func
        else:
            return find_available_random_range


    @property
    def cluster_dict_func(self):
        return self._cluster_dict_func


    @property
    def __len__(self):
        return sum(self.count_dict.values())
    

    def _make_sequence(self):
        '''Private method that actually generates the nucleotide string for
        the Sequence instance. Should be called once when the instance is
        created as everything needed to make it is known during initiation.

        Returns:
            str: Nucleotide sequence as a string.
        '''
        bins = [0 for _ in range(len(self))]
        if cluster_length:
            bins = self._distribute_nucleotides_with_clustering(bins)
            pass
        else:
            bins = self._randomly_distribute_all_nucleotides(bins)
        
        return ''.join(bins)
    

    def _distribute_nucleotides_with_clustering(self, bins):
        '''Fill nucleotide bins, but cluster together the nucleotide
        specified by cluster_nuc in groups of size of cluster_length. 

        Args:
            bins (list): Empty list that represent the slots which are filled
            with nucleotide symbols.

        Returns:
            list: Bins populated with clusters of nucleotide specified by
            cluster_nuc.
        '''
        number_clusters = int(self.count_dict[self.cluster_nuc] / self.cluster_length)
        for i in range(number_clusters):
            rand_range = self.cluster_dist_func(bins, self.cluster_length)
            if rand_range:
                start, end = rand_range
                bins[start, end] = self.cluster_nuc
        
        return bins


    def _randomly_distribute_all_nucleotides(self, bins, ignore_nucs=[]):
        '''Most basic form of sequence creation, uses the nucleotide symbols
        and respective counts in count_dict to randomly fill in "bins"
        (empty slots where nucleotides go).

        Args:
            bins (list): Empty list that represent the slots which are filled
            with nucleotide symbols.
            ignore_nucs (list): Do not add any nucleotides specified by this
            list to bins. This is useful is a particular nucleotide as already
            been added and you want to prevent adding twice. 

        Returns:
            list: List filled with randomly distributed nucleotides.
        '''
        available_bins = np.where(np.any(bins==0, axis=1))

        random_bins = list(np.random.choice(
                available_bins, len(available_bins), replace=False)
                )
        
        temp_count_dict = self.count_dict
        for nuc in ignore_nucs:
            temp_count_dict.pop(nuc)

        for nucleotide, count in self.count_dict.values():
            for _ in range(count):
                bins[random_bins.pop()] = nucleotide
        return bins


            






class VairableRegion():
    # defines the content of a Sequence

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
    

    



