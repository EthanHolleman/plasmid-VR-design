

variable_regions = {
    'initiation_regions': 'variable_defs/initiation_plamids.tsv',
    'termination_regions': 'variable_defs/termination_plasmids.tsv',
}

plots = expand(
    'output/{var_name}/plots/{var_name}.pdf',
    var_name=variable_regions.keys()
)

include: 'rules/make_variable_regions.smk'
include: 'rules/plot_variable_regions.smk'


rule all:
    input:
        plots

