from Bio import SeqIO
import csv
from shutil import copyfile
import pandas as pd


class Sequence():
    # stores variable region data for ranking

    def __init__(self, *args, **kwargs):
        self.seq_data = kwargs
        self.rank_dict = {}

    @property
    def overall_distance(self):
        # sum of all ranks
        return sum(self.rank_dict.values())
    
    def to_dict(self):
        '''Convert instance for dictionary. Used to easily write to tsv output.
        '''
        labeled_rank_dict = {}
        labeled_rank_dict = {f'rank_{attribute}': rank for attribute, rank in self.rank_dict.items()}
        labeled_rank_dict.update(self.seq_data)
        return labeled_rank_dict
            

class Ranker():
    # ranks variable regions (as Sequence instances) based on the metrics
    # specified by metrics_dict which is in turn read from the configfile

    def __init__(self, sequences, metrics_dict):
        self.sequences = sequences
        self.metrics_dict = metrics_dict
     
    def rank_sequences(self):

        for each_metric, op in self.metrics_dict.items():
            print(each_metric, op)
            # sort sequences according to metric definition
            # op is how metric should be interpreted (min or max)
            if op.lower() == 'min':  # lowest values get best rank
                reverse=False
            else:
                reverse=True  # highest vals get best rank

            metric_ranks = sorted(self.sequences, key=lambda x: float(x.seq_data[each_metric]))
             # store the ranks in the sequence instances
            for i in range(len(metric_ranks)):
                metric_ranks[i].rank_dict[each_metric] = i
            self.sequences = metric_ranks
    
    @property
    def overall_ranks(self):
        # Order sequences by overall rank. Best sequence will have lowest
        # average rank. Metrics are therefore all equally weighted currently.
        return sorted(self.sequences, key=lambda s: sum(s.rank_dict.values()))
    
    @property
    def top_seq(self):
        return self.overall_ranks[0]
    
    def _copy_top_ranking_fasta(self, output_path):
        top_record = SeqIO.read(self.top_seq.seq_data['fasta'], 'fasta')
        # read with SeqIO so output gets formated nicely
        SeqIO.write(top_record, output_path, 'fasta')
    
    def _copy_top_ranking_tsv(self, output_path):
        top_tsv = self.top_seq.seq_data['tsv']
        copyfile(top_tsv, output_path)
    
    def _write_all_rankings(self, output_path):
        all_rankings = []
        overall_ranking = self.overall_ranks
        for rank, sequence in enumerate(overall_ranking):
            seq_dict = sequence.to_dict()
            seq_dict['overall_rank'] = rank
            all_rankings.append(seq_dict)
        rankings_df = pd.DataFrame(all_rankings)
        rankings_df.to_csv(output_path, sep='\t', index=False)

    def write_rankings(self, all_ranked_path, top_fasta_path, top_tsv_path):
        self._copy_top_ranking_fasta(top_fasta_path)
        self._copy_top_ranking_tsv(top_tsv_path)
        self._write_all_rankings(all_ranked_path)
    


def read_metrics_file(filepath):
    # reads the file produced by agg_seq_metrics.py. This should be a single
    # record tsv file that contains all metrics that will be used to compare
    # sequences
    return pd.read_table(str(filepath), sep='\t').to_dict('records')
           

def main():

    # read snakemake input
    input_path = snakemake.input  # aggregated seq metrics 

    output_path_ranked = snakemake.output['ranked']
    output_path_fasta = str(snakemake.output['fasta'])
    output_path_tsv = str(snakemake.output['tsv'])

    # defines how metrics should be interpreted
    expect_def_dict = snakemake.params['config_dict']['EXPECTATION_DEFS']
    
    dict_list = read_metrics_file(input_path)
    seqs = [Sequence(**row_dict) for row_dict in dict_list]

    ranker = Ranker(seqs, expect_def_dict)
    ranker.rank_sequences()
    ranker.write_rankings(output_path_ranked, output_path_fasta, output_path_tsv)


if __name__ == '__main__':
    main()