import pandas as pd
import csv
from shutil import copyfile


def read_metrics_file(filepath):
    return pd.read_table(str(filepath), sep='\t').to_dict('records')


def copy_top_rank_tsv_fasta(top_ranker, fasta_output, tsv_output):
    print(top_ranker)
    assert top_ranker.rank == 0
    copyfile(top_ranker.seq_data['fasta'], fasta_output)
    copyfile(top_ranker.seq_data['tsv'], tsv_output)


class Sequence():


    def __init__(self, *args, **kwargs):
        self.seq_data = kwargs
        self.rank_dict = {}
        self.rank = None
    
    @property
    def zscore_dict(self):
        # get dictionary of all metrics and zscores for this sequence
        z_dict = {}
        for attribute, val in self.seq_data.items():
            print(attribute)
            if 'zscore' in attribute:
                print('ZSCORE IN ATTRIBUTE')
                attribute_name = attribute.replace('zscore_', '')
                z_dict[attribute_name] = float(val)

        assert z_dict
        return z_dict
    
    @property
    def overall_distance(self):
        # sum of all ranks
        return sum(self.rank_dict.values())
    
    def to_dict(self):
        labeled_rank_dict = {}
        labeled_rank_dict = {f'rank_{attribute}': rank for attribute, rank in self.rank_dict.items()}
        labeled_rank_dict.update(self.seq_data)
        labeled_rank_dict.update({'overall_rank': self.rank})
        labeled_rank_dict.update({'overall_distance': self.overall_distance})
        return labeled_rank_dict
    
    def __lt__(self, other_seq):
        self.overall_distance < other_seq.overall_distance
        

class Ranker():

    def __init__(self, sequences):
        self.sequences = sequences
     
    def rank_sequences_by_metric(self, attribute, direction, divergence):
        assert direction == 1 or direction == -1
        assert divergence == 1 or divergence == -1

        if divergence == 1:  # minimize distance to mean don't care about direction
            sorted_seqs = sorted(
                self.sequences, 
                key=lambda s: abs(s.zscore_dict[attribute])
                )
        else:  # maximize distance to mean
            if direction == 1:  # desired value right of mean
                sorted_seqs = sorted(
                    self.sequences, 
                    key=lambda s: s.zscore_dict[attribute]
                    )
            else:  # left of mean
                sorted_seqs = sorted(
                self.sequences,
                key=lambda s: s.zscore_dict[attribute], revese=True
                )
            
        # assign rank to attribute
        for i, each_seq in enumerate(sorted_seqs):
            sorted_seqs[i].rank_dict[attribute] = i
        
        self.sequences = sorted_seqs


def get_metrics(dict_list):
    metrics = []
    keys = list(dict_list[0].keys())
    for attribute in keys:
        if 'zscore_' in attribute:
            attribute_name = attribute.replace('zscore_', '')
            metrics.append(attribute_name)
    
    return metrics


def main():

    input_path = snakemake.input
    output_path_ranked = snakemake.output['ranked']
    output_path_fasta = str(snakemake.output['fasta'])
    output_path_tsv = str(snakemake.output['tsv'])

    expect_def_dict = snakemake.params['config_dict']['EXPECTATION_DEFS']

    dict_list = read_metrics_file(input_path)

    seqs = [Sequence(**row_dict) for row_dict in dict_list]
    ranker = Ranker(seqs)

    metrics = get_metrics(dict_list)

    for metric in metrics:
        print(metric, 'thIS IS IS METRIC')
        ranker.rank_sequences_by_metric(metric, **expect_def_dict[metric])
    
    ranked_seqs = sorted(ranker.sequences, key=lambda x: x.overall_distance)
    for i, seq in enumerate(ranked_seqs):
        ranked_seqs[i].rank = i
    
    ranked_seq_dicts = [seq.to_dict() for seq in ranked_seqs]

    pd.DataFrame(ranked_seq_dicts).to_csv(output_path_ranked, sep='\t', index=False)

    copy_top_rank_tsv_fasta(ranked_seqs[0], output_path_fasta, output_path_tsv)


if __name__ == '__main__':
    main()