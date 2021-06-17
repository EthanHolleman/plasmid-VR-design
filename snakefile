import pandas as pd

variable_regions = {
    'initiation_regions': 'variable_defs/initiation_plamids.tsv',
    'termination_regions': 'variable_defs/termination_plasmids.tsv',
}

plots = expand(
    'output/{var_name}/plots/{var_name}.pdf',
    var_name=variable_regions.keys()
)

RNA_ss_rnaFold = expand(
    'output/RNA_sec_struct/viennaRNA/{var_name}.out',
    var_name=variable_regions.keys()
)


RNA_ss_SPOTRNA = expand(
    'output/RNA_sec_struct/SPOT-RNA-motifs/{var_name}',
    var_name=variable_regions.keys()
)


vr_tables = {
    vr_name: pd.read_table(variable_regions[vr_name]).set_index(
        'name', drop=False) for vr_name in variable_regions
    }

brRNA = []
for name, table in vr_tables.items():
    for index, row in table.iterrows():
        brRNA.append(
            f'output/RNA_sec_struct/{name}/{row["name"]}.sc'
        )






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
include: 'rules/plot_variable_regions.smk'
include: 'rules/RNA_sec_struct.smk'
include: 'rules/rlooper.smk'



rule all:
    input:
        brRNA=brRNA

