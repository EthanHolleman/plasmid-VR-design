
rule make_anchor_sequence:
    conda:
        '../envs/pyGibson.yml'
    input:
        variable_regions='output/{var_name}/sequences/variable_regions.fasta'
    params:
        backbones=lambda wildcards: list(config['backbones'].values()),
        prohibited_cutters=config['prohibited_restriction_enzyme_recognition_sites'],
        max_attempts=config['anchor_sequence']['max_attempts'],
        anchor_length=config['anchor_sequence']['anchor_length'],
        min_melting_temp=config['anchor_sequence']['min_melting_temp']
    output:
        'output/{var_name}/sequences/anchor.gb'
    script:'../scripts/make_anchor_region.py'


rule md5sum_anchor:
    input:
        'output/{var_name}/sequences/anchor.gb'
    output:
        'output/{var_name}/sequences/anchor.md5sum'
    shell:'''
    md5sum {input} > {output}
    '''


rule assemble_inserts:
    conda:
        '../envs/pyGibson.yml'
    input:
        variable_region='output/{var_name}/files/{p_name}/rankedSeqs/{p_name}.top_seq.tsv',
        upstream_sequences=[  # in the order they should appear in final insert
            config['insert_design']['upstream_regions'],
            'output/{var_name}/sequences/anchor.gb'
        ],
        downstream_sequences=config['insert_design']['downstream_regions'],
        anchor_checksum='output/{var_name}/sequences/anchor.md5sum',

    params:
        design_version=config['insert_design']['design_version'],
        feature_annotations=config['insert_design']['feature_annotations']
    output:
        'output/{var_name}/inserts/genbank_files/{p_name}.insert.gb'
    script:'../scripts/assemble_inserts.py'


rule insert_checksum:
    input:
        lambda wildcards: expand(
            'output/{var_name}/inserts/genbank_files/{p_name}.insert.gb',
            p_name=get_all_p_names(wildcards), allow_missing=True
        )
    output:
        'output/{var_name}/inserts/genbank_files/inserts.md5sum'
    shell:'''
    md5sum {input} > {output}
    '''


rule aggregate_inserts_into_fasta:
    conda:
        '../envs/pyGibson.yml'
    input:
        inserts=lambda wildcards: expand(
            'output/{var_name}/inserts/genbank_files/{p_name}.insert.gb',
            p_name=get_all_p_names(wildcards), allow_missing=True
        ),
        checksums='output/{var_name}/inserts/genbank_files/inserts.md5sum'
    output:
        'output/{var_name}/inserts/complete_inserts.fa'
    script:'../scripts/agg_gb_inserts_to_fa.py'


rule fasta_checksum:
    input:
        'output/{var_name}/inserts/complete_inserts.fa'
    output:
        'output/{var_name}/inserts/complete_inserts.md5sum'
    shell:'''
    md5sum {input} > {output}
    '''











