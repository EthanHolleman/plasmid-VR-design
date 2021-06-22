# Script for aggregating the output of programs that produce metrics
# used to access the quality of a given variable region into a tsv file which
# is then used to rank all sequences.
import statistics
import csv
import pandas as pd
import os

class Expectation():

    expectations = {}


    @classmethod
    def extend_expect_dict_from_table(cls, filepath, length):
        # create list of expectation instances from expectation tsv
        # with colnames length, mean, sd, attribute
        expect_dict = {}
        table = pd.read_table(filepath, sep='\t')
        len_rows = table.loc[table['length']==length]
        for index, row in len_rows.iterrows():
            attribute, mean, sd = str(row['attribute']), float(row['mean']), float(row['sd'])
            cls.expectations[attribute] = cls(attribute, mean, sd)


    def __init__(self, name, expect_mean, expect_sd):
        self.name = str(name)
        self.expect_mean = float(expect_mean)
        self.expect_sd = float(expect_sd)

    def distance(self, metric):
        # distance in standard deviations from the mean
        print(metric, self.expect_mean, 'EXPECTED MEAN')
        z = (metric.value - self.expect_mean) / self.expect_sd
        print(metric.value, self.expect_mean, z, 'VAL, MEAN, Z')
        assert isinstance(z, float)
        
        return f'distance_{self.name}',  z
        
    def __repr__(self):
        return ' '.join([f'{str(key)}:{str(val)}' for key, val in self.__dict__.items()])
    

class Metric():

    NAME = 'metric'

    def __init__(self, filepath, expectation, direction, divergence):
        self.filepath = filepath
        self.expectation= expectation
        self.direction = direction
        self.divergence = divergence

        assert os.path.isfile(filepath)
        assert isinstance(expectation, Expectation)

        print(self, expectation)

    @property
    def distance(self):
        dist_name, val = self.expectation.distance(self)
        val = self._apply_scaling(val)
        return {dist_name: val}
    
    @property
    def value(self):
        return None
    
    @property
    def direction(self):
        return self._direction
    
    @direction.setter
    def direction(self, new_dir):
        if new_dir == 1 or new_dir == -1:
            self._direction = new_dir
        else:
            raise TypeError('Direction must be -1 or 1!')
    
    @property
    def divergence(self):
        return self._divergence
    
    @direction.setter
    def divergence(self, new_dir):
        if new_dir == 1 or new_dir == -1:
            self._divergence = new_dir
        else:
            raise TypeError('Divergence must be -1 or 1!')
    
    def parse(self):
        d = {
            type(self).NAME: self.value,  # original value
        }
        d.update(self.distance)
        return d
    
    def _apply_scaling(self, z):
        return (z * self.direction) ** self.divergence


class rlooperMetric(Metric):

    def __init__(self, filepath, expectation, direction, divergence):
        super().__init__(filepath, expectation, direction, divergence)
    
    def _parse_rlooper_wig(self):
        values = []
        with open(self.filepath) as handle:
            lines = handle.readlines()[5:]  # skip first 4 lines
        return [float(v.strip()) for v in lines]


class bpProb(rlooperMetric):
    # basepair probability from rlooper
    NAME = 'bp_prob'

    def __init__(self, filepath, expectation, direction, divergence):
        super().__init__(filepath, expectation, direction, divergence)

    @property
    def value(self):
        return statistics.mean(self._parse_rlooper_wig())


class localAverageEnergy(rlooperMetric):

    NAME = 'local_average_energy'

    def __init__(self, filepath, expectation, direction, divergence):
        super().__init__(filepath, expectation, direction, divergence)

    @property
    def value(self):
        return statistics.mean(self._parse_rlooper_wig())


class rnaMetric(Metric):

    NAMES_DICT = {'prop_hairpin': 1, 
                  'prop_unpaired': 4
                }

    def __init__(self, filepath, expectation, direction, divergence):
        super().__init__(filepath, expectation, direction, divergence)
    

    def _parse_rna_file(self):
        d = {}
        with open(self.filepath) as handle:
            reader = csv.reader(handle, delimiter='\t')
            record = next(reader)
            for attribute, index in rnaMetric.NAMES_DICT.items():
                d[attribute] = record[index]
        return d


class propHairpin(rnaMetric):

    NAME = 'prop_hairpin'

    def __init__(self, filepath, expectation, direction, divergence):
        super().__init__(filepath, expectation, direction, divergence)

    @property
    def value(self):
        return float(self._parse_rna_file()[propHairpin.NAME])
    

class propUnpaired(rnaMetric):

    NAME = 'prop_unpaired'


    def __init__(self, filepath, expectation, direction, divergence):
        super().__init__(filepath, expectation, direction, divergence)
    
    @property
    def value(self):
        return float(self._parse_rna_file()[propUnpaired.NAME])

            
def add_inputs_to_out_dict(out_dict, snakemake):
    out_dict['fasta'] = snakemake.input['fasta']
    out_dict['tsv'] = snakemake.input['tsv']
    out_dict['length'] = snakemake.params['length']
    out_dict['p_name'] = snakemake.params['p_name']
    out_dict['id_num'] = snakemake.params['id_num']

    # out_dict['fasta'] = 'output/initiation_regions/files/init-1/10/init-1.10.fasta'
    # out_dict['tsv'] = 'output/initiation_regions/files/init-1/10/init-1.10.tsv'
    # out_dict['length'] = 200

    return out_dict


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

    # read metric expectation files 
    rlooper_expect_path = str(snakemake.input['rlooper_expect'])
    rna_expect_path = str(snakemake.input['RNAss_expect'])

    output_path = str(snakemake.output)

    expect_def_dict = snakemake.params['config_dict']['EXPECTATION_DEFS']

    # input_bp_prob = 'output/initiation_regions/files/init-1/10/init-1.1_bpprob.wig'
    # input_lae = 'output/initiation_regions/files/init-1/1/init-1.1_avgG.wig'
    # input_rna = 'output/initiation_regions/files/init-1/1/parsedRNA/init-1.tsv'

    # length = 200  # length of sequence

    # # read metric expectation files 
    # rlooper_expect_path = 'output/expectations/rlooper/rlooper_expect.200.tsv'
    # rna_expect_path = 'output/expectations/SPOTRNA/spotRNA_expect.200.tsv'


    # initialize expectation instances from expect files
    Expectation.extend_expect_dict_from_table(rlooper_expect_path, length)
    Expectation.extend_expect_dict_from_table(rna_expect_path, length)

    
    metrics = [
        bpProb(
            input_bp_prob, Expectation.expectations[bpProb.NAME],
            **expect_def_dict[bpProb.NAME]),
        localAverageEnergy(
            input_lae, Expectation.expectations[localAverageEnergy.NAME],
            **expect_def_dict[localAverageEnergy.NAME]),
        propHairpin(
            input_rna, Expectation.expectations[propHairpin.NAME],
            **expect_def_dict[propHairpin.NAME]),
        propUnpaired(
            input_rna, Expectation.expectations[propUnpaired.NAME],
            **expect_def_dict[propUnpaired.NAME])
    ]


    out_dict = {}

    out_dict = add_inputs_to_out_dict(out_dict, snakemake)

    [out_dict.update(metric.parse()) for metric in metrics]  # add all parsed metrics

    write_out_dict_as_tsv(out_dict, output_path)


if __name__ == '__main__':
    main()





        

    


    
    