

rule generate_variable_regions:
    output:
        'output/{var_name}/files/{var_name}.fasta',
        'output/{var_name}/files/{var_name}.tsv'
    params:
        input = lambda wildcards: variable_regions[wildcards.var_name]
    script:'../scripts/varmid/varmid.py'
        