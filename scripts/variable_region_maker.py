import numpy as np

from utils import *

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
    
    def __init__(self, name, count_dict, variable_region, cluster_length=None, cluster_nuc=None,
                cluster_dist_func=None):
        self.name = name
        self.count_dict = count_dict
        self.variable_region = variable_region
        self.cluster_length = cluster_length
        self.cluster_nuc = cluster_nuc
        self.cluster_dist_func = cluster_dist_func
        self.nuc_seq = self._make_sequence()

        if self.cluster_length:
            assert cluster_nuc, 'Must specify cluster nuc if passing cluster len!'
            assert cluster_nuc in count_dict
    
    @property
    def cluster_dict_func(self):
        return self._cluster_dict_func

    @cluster_dict_func.setter
    def cluster_dict_func(self, new_func):
        if callable(new_func):
            self._cluster_dict_func = new_func
        else:
            return find_available_random_range

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

    @classmethod
    def init_from_csv_row(cls, row_dict):
        return cls(**row_dict)


    def __init__(self, name, length, gc_content=None, gc_skew=None, at_skew=None, 
                at_content=None, cluster_length=None, cluster_nuc=None,
                cluster_dist_func=None, role=None, **kwargs):
        
        self.gc_count = (0, 0)
        self.at_count = (0, 0)

        self.name = name
        self.length = length
        self.gc_skew = gc_skew
        self.gc_content = gc_content
        self.at_skew = at_skew
        self.at_content = at_content
        self.cluster_length = cluster_length
        self.cluster_nuc = cluster_nuc
        self.cluster_dist_func = cluster_dist_func
        self.role = role

        self._generated_sequences = 0
    
    @property
    def cluster_nuc(self):
        return self._cluster_nuc
    
    @cluster_nuc.setter
    def cluster_nuc(self, new_nuc):

        if isinstance(new_nuc, float) or new_nuc in self.nuc_dict:
            self._cluster_nuc = new_nuc
        else:
            raise TypeError(
                f'Must pass a cannonical nucleotide or None! Not this shit {new_nuc}')


    @property
    def nuc_dict(self):
        return {
            'A': self.at_count[0],
            'T': self.at_count[1],
            'G': self.gc_count[0],
            'C': self.gc_count[1]
        }


    def calculate_nucleotide_counts(self):
        # gc arbitrary takes priority in calculation because
        # it is my favorite
        self._calculate_gc_count()
        self._calculate_at_count()


    def _calculate_gc_count(self):
        '''Calculate the number of G and C nucleotides that should be included
        in the sequence given the gc skew and content.
        '''
        if self.at_count:  # already calculated AT count
            effective_length = self.length - sum(self.gc_count)
        else:
            effective_length = self.length

        if self.gc_skew and self.gc_content:
            print(self.gc_skew, self.gc_content, 'gc')
            self.gc_count = nuc_count_calculator(
                skew=self.gc_skew,
                content=self.gc_content,
                seq_len=effective_length
            )
            print(self.gc_count, 'gc_count')
        else:
            self.gc_count = get_int_half_length(effective_length)
    

    def _calculate_at_count(self):
        if self.gc_count:  # already calculated GC count
            effective_length = self.length - sum(self.gc_count)
        else:
            effective_length = self.length
        
        if self.at_skew > 0 and self.at_content > 0:
            self.at_count = nuc_count_calculator(
                skew=self.at_skew,
                content=self.at_content,
                seq_len=effective_length
            )
        else:
            self.at_count = get_int_half_length(effective_length)
        
        
    def generate_sequence(self):
        if not self.gc_count:
            self.calculate_nucleotide_counts()
        
        self._generated_sequences += 1
        return Sequence(
            self._sequence_name(),
            self._nuc_count,
            self,
            self.cluster_length,
            self.cluster_nuc,
            self.cluster_dist_func
            )
        
    def _sequence_name(self):
        return f'{VairableRegion.name_prefix}_{self.role}_{self._generated_sequences}_{self.gc_skew}_{self.gc_content}_unclustered'

    def __len__(self):
        return self.length


p = '/home/ethan/Documents/GitHub/plasmid-design/variable_defs/initiation_plamids.csv'
vrs = read_variable_region_config_file(p)

for row_dict in vrs:
    v = VairableRegion.init_from_csv_row(row_dict)
    print(v.nuc_dict)


#v = VairableRegion.init_from_csv('




