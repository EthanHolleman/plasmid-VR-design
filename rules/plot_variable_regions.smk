

rule plot_variable_regions:
    conda:'../envs/R.yml'
    input:
        'output/{var_name}/files/{var_name}.tsv'
    output:
        'output/{var_name}/plots/{var_name}.pdf'
    shell:'../scripts/plot.R'
    