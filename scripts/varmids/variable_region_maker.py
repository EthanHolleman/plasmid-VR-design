import numpy as np

from utils import *


class Sequence():
    
    def __init__(self, name, count_dict, variable_region, cluster_length=None, 
                cluster_nuc=None, cluster_dist_func=None):
        self.name = name
        self.count_dict = count_dict
        self.variable_region = variable_region
        self.cluster_length = cluster_length
        self.cluster_nuc = cluster_nuc
        self.cluster_dist_func = cluster_dist_func
        self.nuc_seq = self._make_sequence()

        assert len(self.nuc_seq) == self.variable_region.length, f'Nuc seq: {len(self.nuc_seq)} VR: {self.variable_region.length}'

        if self.cluster_length:
            assert cluster_nuc, 'Must specify cluster nuc if passing cluster len!'
    
    @property
    def cluster_dist_func(self):
        return self._cluster_dist_func

    @cluster_dist_func.setter
    def cluster_dist_func(self, new_func):
        if callable(new_func):
            self._cluster_dist_func = new_func
        else:
            self._cluster_dist_func = find_available_random_range


    def __len__(self):
        return sum(self.count_dict.values())
    

    def _make_sequence(self):
        '''Private method that actually generates the nucleotide string for
        the Sequence instance. Should be called once when the instance is
        created as everything needed to make it is known during initiation.

        Returns:
            str: Nucleotide sequence as a string.
        '''
        bins = np.zeros(len(self))
        if self.cluster_length:
            bins = self._distribute_nucleotides_with_clustering(bins)
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
        if number_clusters == 0:
            # if cluster length is less than nucleotide count reset
            # cluster length to number of available nucleotides
            number_clusters = 1
            self.cluster_length = int(self.count_dict[self.cluster_nuc])

        seq_bins = [None for _ in range(len(bins))]
        for i in range(number_clusters):
            rand_range = self.cluster_dist_func(bins, self.cluster_length)

            if rand_range:
                start, end = rand_range
                for i in range(start, end):
                    bins[i] = 1
            else:
                break  # exit since no more clusters can be added
        
        # set up seq_list as random nucleotide sequence. This is to avoid having
        # empty bins in cases where rounding in order to get the closest possible
        # content / skew counts causes sum of nucleotide count to differ
        # from specified length by one or two nucleotides
        seq_list = [str(np.random.choice(['A', 'T', 'G', 'C'], 1)[0]) for _ in range(len(bins))]
        for i in np.where(bins == 1)[0]:
            seq_list[i] = self.cluster_nuc
        
        seq_list = self._randomly_distribute_all_nucleotides(bins, 
                ignore_nucs=[self.cluster_nuc], seq_list=seq_list)
        return seq_list



    def _randomly_distribute_all_nucleotides(self, bins, ignore_nucs=[], seq_list=[]):
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
        available_bins = np.where(bins==0)[0]

        random_bins = list(np.random.choice(
                available_bins, len(available_bins), replace=False)
                )
        
        if not seq_list:
            seq_list = [None for _ in range(len(bins))]
        for nucleotide, count in self.count_dict.items():
            if nucleotide not in ignore_nucs:
                for _ in range(count):
                    seq_list[random_bins.pop()] = nucleotide
        return seq_list
    

    def as_fasta_entry(self):
        '''Return sequence as fasta formated string.
        '''
        return f'>{self.name}\n{fasta_seq_printer(self.nuc_seq)}\n'
    

    def to_dict(self):
        return {
            'name': self.name,
            'GC_content': self.variable_region.gc_content,
            'GC_skew': self.variable_region.gc_skew,
            'AT_content': self.variable_region.at_content,
            'AT_skew': self.variable_region.at_skew,
            'Cluster_length': self.cluster_length,
            'Clustered_nucleotide': self.cluster_nuc,
            'Clustering method': self.cluster_dist_func.__name__,
            'Sequence': self.nuc_seq,
        }
        



        
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

        self._infer_contents()
        self._check_content_skew()
    
    def _check_content_skew(self):
        '''Private method to check if all content and skew type attributes are
        valid.

        Raises:
            AssertionError: Raise if content or skew is not valid value (> 1).

        Returns:
            
            None: Return None if all checks passed.
        '''
        if self.at_content and self.gc_content:
            assert self.at_content + self.gc_content == 1, f'AT + GC content must equal 1: {self.at_content}, {self.gc_content}'
        for param in (self.gc_skew, self.gc_content, self.at_content, self.at_skew):
            if param:
                assert param <= 1 and param >= -1, f'{param} not valid value'
    

    def _infer_contents(self):
        '''Setting GC content implies an AT content and vice versa. This
        private function detects when one has been set and not the other
        and infers the unset value. GC + AT content must equal 1. 
        '''
        if self.gc_content and not self.at_content:
            self.at_content = 1 - self.gc_content
        elif self.at_content and not self.gc_content:
            self.gc_content = 1 - self.at_content
    
    @property
    def gc_skew(self):
        return self._gc_skew
    
    @gc_skew.setter
    def gc_skew(self, new_skew):
        if type(new_skew) == None:
            new_skew = 0
        elif not isinstance(new_skew, float):
            raise TypeError
        else:
            self._gc_skew = new_skew

    @property
    def at_skew(self):
        return self._gc_skew
    
    @at_skew.setter
    def at_skew(self, new_skew):
        if new_skew == None:
            new_skew = 0
        elif not isinstance(new_skew, float):
            raise TypeError(f'{type(new_skew)}')
        else:
            self._at_skew = new_skew
    
        
    @property
    def cluster_nuc(self):
        return self._cluster_nuc
    
    @cluster_nuc.setter
    def cluster_nuc(self, new_nuc):

        if new_nuc == None or new_nuc in self.nuc_dict:
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
        # The order which counts are calculated matters when only at or gc
        if self.at_content or self.at_skew:
            self._calculate_at_count()
            self._calculate_gc_count()
        else:
            self._calculate_gc_count()
            self._calculate_at_count()



    def _calculate_gc_count(self):
        '''Calculate the number of G and C nucleotides that should be included
        in the sequence given the gc skew and content.
        '''
        if self.gc_skew or self.gc_content:
            self.gc_count = nuc_count_calculator(
                skew=self.gc_skew,
                content=self.gc_content,
                seq_len=self.length
            )
        else:
            self.gc_count = get_int_half_length(self.length - sum(self.at_count))
            

    def _calculate_at_count(self):
        
        if self.at_skew or self.at_content:
            self.at_count = nuc_count_calculator(
                skew=self.at_skew,
                content=self.at_content,
                seq_len=self.length
            )
        else:
            self.at_count = get_int_half_length(self.length - sum(self.gc_count))
        
        
    def generate_sequence(self):
        if sum(self.nuc_dict.values()) == 0:
            self.calculate_nucleotide_counts()
        
        self._generated_sequences += 1
        return Sequence(
            self._sequence_name(),
            self.nuc_dict,
            self,
            self.cluster_length,
            self.cluster_nuc,
            self.cluster_dist_func
            )
    
        
    def _sequence_name(self):
        return f'{VairableRegion.name_prefix}_{self.name}_{self.role}_{self._generated_sequences}_{self.gc_skew}_{self.gc_content}'




