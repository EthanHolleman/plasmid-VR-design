import pandas as pd


variable_regions = config['variable_region_definitions']

vr_tables = {
    vr_name: pd.read_table(variable_regions[vr_name]).set_index(
        'name', drop=False) for vr_name in variable_regions
    }

RLOOPER_FILE_SUFFI = config['RLOOPER_FILE_SUFFI']
SPOT_RNA_EXTS = config['SPOT_RNA_EXTS']

NUM_CASES = config['EXPECTATION_SAMPLES']
CASE_RANGE = range(1, NUM_CASES+1)

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
