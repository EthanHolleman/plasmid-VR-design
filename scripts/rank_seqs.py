import pandas as pd


def read_metrics_file(filepath):
    return pd.read_table(str(filepath), sep='\t').to_dict('records')

def write_ranker_list_to_file(output_path, sorted_ranker_list):

    with open(str(output_path)) as handle:
        headers = list(ranker_list[0].__dict__.keys())
        writer = csv.DictWriter(handle, delimiter='\t', fieldnames=headers)
        writer.writeheader()
        for rank in sorted_ranker_list:
            writer.writerow(rank)
        return output_path


class Ranker():

    def __init__(**kwargs):
        self.__dict__.update(kwargs)
    
    @property
    def overall_fit(self):
        distance = 0
        for key, val in self.__dict__.items():
            if 'distance' in key:
                distance += float(val)
        return distance
    
    def __lt__(self, other_ranker):
        return self.overall_fil < other_ranker.overall_fit


def main():

    input_path = snakemake.input
    output_path_ranked = snakemake.output['ranked']

    dict_list = read_metrics_file(input_path)
    ranker_list = [Ranker(**row_dict) for row_dict in dict_list]
    ranker_list = sorted(ranker_list)
    write_ranker_list_to_file(output_path_ranked, ranker_list)

    # still need to copy best fasta and tsv file 
