import pandas as pd

from scripts.expand_RC_seqs import add_reverse_complement_definitions

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
include: 'rules/plot_variable_regions.smk'


rule all:
    input:
        expand(
            'output/{var_regions}/sequences/plasmid_sequences.fasta',
            var_regions=variable_regions.keys()
        ),
        expand(
            'output/{var_regions}/plots/{var_regions}.pdf',
            var_regions=variable_regions.keys()
        )
