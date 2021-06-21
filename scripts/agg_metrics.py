# Script for aggregating the output of programs that produce metrics
# used to access the quality of a given variable region into a tsv file which
# is then used to rank all sequences.

class Metric():

    def __init__(self, name, filepaths):
        self.name = name
        self.filepaths = filepaths
    
    def convert_to_columns(self):
        return None


class rlooperMetric(Metric):

    def __init__(self, filepath):
        super().__init__(self, filepath)
    
    def parse_rlooper_wig(self):
        values = []
        with open(self.filepath) as handle:
            lines = handle.readlines()[5:]  # skip first 4 lines
        return [float(v.strip()) for v in lines]


class bpProb(rlooperMetric):
    # basepair probability from rlooper
    NAME = 'bp_prob'

    def __init__(self, filepaths):
        super().__init__(self, name, filepath)
    
    def convert_to_columns(self):
        

    


    
    