

rule plot_variable_regions:
    conda:
        '../envs/R.yml'
    input:
        variable_regions='output/{var_name}/files/{var_name}.tsv',
        rlooper_calculations=lambda wildcards: expand(
            'output/rlooper/{var_name}/completed_runs/{record}/{record}_{rlooper_suffix}',
            var_name=wildcards.var_name, record=vr_tables[wildcards.var_name]['name'],
            rlooper_suffix='avgG.wig'
        )
    output:
        'output/{var_name}/plots/{var_name}.pdf'
    script:'../scripts/plot.R'
    