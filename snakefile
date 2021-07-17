import pandas as pd

from scripts.expand_RC_seqs import add_reverse_complement_definitions

variable_regions = config['variable_region_definitions']

vr_tables = {
    vr_name: pd.read_table(variable_regions[vr_name]).set_index(
        'name', drop=False) for vr_name in variable_regions
    }

RLOOPER_FILE_SUFFI = config['RLOOPER_FILE_SUFFI']
SPOT_RNA_EXTS = config['SPOT_RNA_EXTS']

NUM_CASES = config['NUMBER_REPLICATES']
CASE_RANGE = range(1, NUM_CASES+1)

NUMBER_REPLICATES = config['NUMBER_REPLICATES']
EXPECT_SAMPLES = config['EXPECT_SAMPLES']
RAND_SEQ_NAMES = [f'RAND-SEQ:{i}' for i in range(0, EXPECT_SAMPLES)]

wildcard_constraints:
   var_name = '\w+'

include: 'rules/make_variable_regions.smk'
include: 'rules/calculate_expectations.smk'
include: 'rules/RNA_sec_struct.smk'
include: 'rules/rlooper.smk'
include: 'rules/plot_variable_regions.smk'
include: 'rules/insert_assembly.smk'
include: 'rules/simulate_constructs.smk'
include: 'rules/primer3.smk'


rule all:
    input:
        expand(  # make all complete insert sequences
            'output/{var_regions}/inserts/complete_inserts.fa',
            var_regions=variable_regions.keys()
        ),
        expand(  # make complete insert checksums 
            'output/{var_name}/inserts/complete_inserts.md5sum',
            var_name=variable_regions.keys()

        ),
        expand(  # verify all complete inserts
            'output/{var_name}/inserts/.all_checks.passed',
            var_name=variable_regions.keys()
        ), # rlooper results all sequences
        expand(
            'output/expectations/{var_name}/rlooper/rlooper_expect.png',
            var_name=variable_regions.keys()
        ),  # same thing but for spot-rna predictions
        expand(
            'output/expectations/{var_name}/SPOT-RNA/rna_secondary_structure.png',
            var_name=variable_regions.keys()
        ),
        expand(  # Simulate construct assembly
            'output/{var_name}/constructs',
            var_name=variable_regions.keys()
        ),
        expand(  # tac series primers
            ['output/{var_name}/sequences/tac_initiation_series_primers.fa',
            'output/{var_name}/sequences/tac_termination_series_primers.fa'],
            var_name=variable_regions.keys()
        )
        
        
