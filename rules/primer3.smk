

rule format_primer3_template_for_variable_region:
    conda:
        '../envs/python.yml'
    input:
        tsv='output/{var_name}/files/{p_name}/rankedSeqs/{p_name}.top_seq.tsv',
        template=config['PRIMER3_TEMPLATE']
    output:
        'output/{var_name}/files/{p_name}/primers/{p_name}.primer3.template.txt'
    script:'../scripts/format_primer3_template.py'


rule make_pcr_primers_for_variable_region:
    conda:
        '../envs/primer3.yml'
    input:
        'output/{var_name}/files/{p_name}/primers/{p_name}.primer3.template.txt'
    output:
        'output/{var_name}/files/{p_name}/primers/{p_name}.primers.primer3'
    shell:'''
    primer3_core {input} > {output}
    '''

rule make_all_variable_region_primers:
    input:
        lambda wildcards: expand(
            'output/{var_name}/files/{p_name}/primers/{p_name}.primers.primer3',
            allow_missing=True, p_name=get_all_p_names(wildcards)
        )
    output:
        'output/{var_name}/.variable_region_primers.done'
    shell:'''
    touch {output}
    '''