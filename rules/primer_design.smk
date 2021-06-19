# design PCR primers for all generated plasmid variable regions

rule make_primer3_boulder_io:
    conda:
        '../envs/python.yml'
    input:
        ''
    output:
        'output/{var_name}/primers/{record_name}.primer3'
    script:'../scripts/make_primer'


rule design_primer_from_boulder_io:
    conda:
        '../envs/primer3.yml'
    input:
        'output/make_boulder.io'
    output:
        ''
    shell:'''
    primer3 < {input} > {output}
    '''
    