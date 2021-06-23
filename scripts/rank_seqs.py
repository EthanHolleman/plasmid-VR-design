import pandas as pd
import csv
from shutil import copyfile


def read_metrics_file(filepath):
    return pd.read_table(str(filepath), sep='\t').to_dict('records')

def write_ranker_list_to_file(output_path, sorted_ranker_list):

    with open(str(output_path), 'w') as handle:
        headers = list(sorted_ranker_list[0].__dict__.keys())
        writer = csv.DictWriter(handle, delimiter='\t', fieldnames=headers)
        writer.writeheader()
        for rank in sorted_ranker_list:
            rank_dict = rank.__dict__
            print(rank_dict)
            writer.writerow(rank_dict)
        return output_path


def copy_top_rank_tsv_fasta(top_ranker, fasta_output, tsv_output):
    assert top_ranker.rank == 0
    copyfile(top_ranker.fasta, fasta_output)
    copyfile(top_ranker.tsv, tsv_output)



class Ranker():

    def __init__(self, *args, **kwargs):
        self.__dict__.update(kwargs)
        self.rank = None
    
    @property
    def overall_fit(self):
        distance = 0
        for key, val in self.__dict__.items():
            if 'distance' in key:
                distance += float(val)
        return distance
    
    def __lt__(self, other_ranker):
        return self.overall_fit < other_ranker.overall_fit


def main():

    input_path = snakemake.input
    output_path_ranked = snakemake.output['ranked']
    output_path_fasta = str(snakemake.output['fasta'])
    output_path_tsv = str(snakemake.output['tsv'])

    dict_list = read_metrics_file(input_path)
    ranker_list = [Ranker(**row_dict) for row_dict in dict_list]
    ranker_list = sorted(ranker_list, reverse=True)  # highest value first = best seq first
    for i, ranker in enumerate(ranker_list):
        ranker.rank = i
        ranker_list[i] = ranker
    write_ranker_list_to_file(output_path_ranked, ranker_list)
    copy_top_rank_tsv_fasta(ranker_list[0], output_path_fasta, output_path_tsv)


if __name__ == '__main__':
    main()