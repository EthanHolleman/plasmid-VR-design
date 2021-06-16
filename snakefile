

variable_regions = {
    'initiation_regions': 'variable_defs/initiation_plamids.tsv',
    'termination_regions': 'variable_defs/termination_plasmids.tsv',
}

plots = expand(
    'output/{var_name}/plots/{var_name}.pdf',
    var_name=variable_regions.keys()
)

RNA_ss = expand(
    'output/rlooper/{var_name}/fasta',
    var_name=variable_regions.keys()
)

import pandas as pd

vr_tables = {vr_name: pd.read_table(variable_regions[vr_name]).set_index('name', drop=False) for vr_name in variable_regions}

RLOOPER_FILE_SUFFI = [
    'avgG.wig',
    'bpprob.wig',
    'bpprob.bed',
    'mfe.bed',
    'mfe.wig'
]

def rlooper_output(*args, **kwargs):
    output_files = []
    for vr in vr_tables:
        table = vr_tables[vr]
        output_files += expand(
            'output/rlooper/{var_name}/completed_runs/{record}/{record}_{rlooper_suffix}',
            var_name=vr, record=table['name'], rlooper_suffix=RLOOPER_FILE_SUFFI
        )

    return output_files


include: 'rules/make_variable_regions.smk'
include: 'rules/plot_variable_regions.smk'
include: 'rules/RNA_sec_struct.smk'
include: 'rules/rlooper.smk'



rule all:
    input:
        RNA_ss=rlooper_output()

