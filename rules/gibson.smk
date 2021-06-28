rule make_assembly:
    conda:
        '../envs/pyGibson.yml'
    input:
        vr_fasta='output/{var_name}/files/{p_name}/rankedSeqs/{p_name}.top_seq.fasta',
    output:
        vr_genbank='output/{var_name}/files/{p_name}/constructs/{p_name}.gb',
        assembly_dir=directory('output/{var_name}/files/{p_name}/constructs')
    params:
        vr_row_def = lambda wildcards: get_p_record(wildcards),
        constructs = config['construct_definitions']
    script:'../scripts/gibson_assembly/gibson_assembly.py'


rule assemble_all_plasmid_constructs:
    input:
        lambda wildcards: expand(
            'output/{var_name}/files/{p_name}/constructs',
            var_name=variable_regions.keys(),
            p_name=get_all_p_names(wildcards)
        )
    output:
        'output/{var_name}/files/.assemble_all_plasmid_constructs.done'
    shell:'''
    touch {output}
    '''