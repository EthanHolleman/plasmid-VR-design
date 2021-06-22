import pandas as pd


# plots = expand(
#     'output/{var_name}/plots/{var_name}.pdf',
#     var_name=variable_regions.keys()
# )



variable_regions = {
    'initiation_regions': 'variable_defs/initiation_plamids_short.tsv',
}

vr_tables = {
    vr_name: pd.read_table(variable_regions[vr_name]).set_index(
        'name', drop=False) for vr_name in variable_regions
    }


RLOOPER_FILE_SUFFI = [
    'avgG.wig',
    'bpprob.wig',
    'bpprob.bed',
    'mfe.bed',
    'mfe.wig'
]


SPOT_RNA_EXTS = [
    'bpseq',
    'ct',
    'prob'
]

wildcard_constraints:
   var_name = '\w+'

include: 'rules/make_variable_regions.smk'
include: 'rules/calculate_expectations.smk'
include: 'rules/RNA_sec_struct.smk'
include: 'rules/rlooper.smk'



rule all:
    input:
        expand(
            'output/initiation_regions/files/init-1/{id_num}/aggregatedMetrics/init-1.tsv',
            id_num=CASE_RANGE
        )
