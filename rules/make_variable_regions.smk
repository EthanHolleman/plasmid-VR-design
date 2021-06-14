

rule generate_variable_regions:
    conda:
        '../envs/python.yml'
    output:
        fasta='output/{var_name}/files/{var_name}.fasta',
        tsv='output/{var_name}/files/{var_name}.tsv'
    params:
        input_file = lambda wildcards: variable_regions[wildcards.var_name]
    script:'../scripts/varmids/varmids.py'
        