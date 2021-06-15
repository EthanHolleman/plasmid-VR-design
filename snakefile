

variable_regions = {
    'initiation_regions': 'variable_defs/initiation_plamids.tsv',
    'termination_regions': 'variable_defs/termination_plasmids.tsv',
}

plots = expand(
    'output/{var_name}/plots/{var_name}.pdf',
    var_name=variable_regions.keys()
)

RNA_ss = expand(
    'output/RNA_sec_struct/{var_name}',
    var_name=variable_regions.keys()
)

include: 'rules/make_variable_regions.smk'
include: 'rules/plot_variable_regions.smk'
include: 'rules/RNA_sec_struct.smk'


rule all:
    input:
        plots=plots,
        RNA_ss=RNA_ss

