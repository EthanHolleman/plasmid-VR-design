
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


rule final_insert_check:
    conda:
        '../envs/pyGibson.yml'
    input:
        anchor_seq='output/{var_name}/sequences/anchor.gb',
        insert_record='output/{var_name}/inserts/genbank_files/{p_name}.insert.gb',
        homology_target=config['insert_design']['homology_target_backbone'],
        five_prime_arm = config['insert_design']['upstream_regions'],
        three_prime_arm = config['insert_design']['downstream_regions']
    params:
        attributes=lambda wildcards: get_p_record(wildcards),
        prohibited_cutters=config['prohibited_restriction_enzyme_recognition_sites'],
    output:
        'output/{var_name}/inserts/genbank_files/.{p_name}.insert.passed'
    script:'../scripts/verify_inserts.py'


rule check_all_inserts:
    input:
        lambda wildcards: expand(
            'output/{var_name}/inserts/genbank_files/.{p_name}.insert.passed',
             p_name=get_all_p_names(wildcards), allow_missing=True
        )
    output:
        'output/{var_name}/inserts/.all_checks.passed'
    shell:'''
    touch {output}
    '''

rule aggregate_inserts_into_fasta:
    conda:
        '../envs/pyGibson.yml'
    input:
        inserts=lambda wildcards: expand(
            'output/{var_name}/inserts/genbank_files/{p_name}.insert.gb',
            p_name=get_all_p_names(wildcards), allow_missing=True
        ),
        checksums='output/{var_name}/inserts/genbank_files/inserts.md5sum',
        checks='output/{var_name}/inserts/.all_checks.passed'
    output:
        'output/{var_name}/inserts/complete_inserts.fa'
    script:'../scripts/agg_gb_inserts_to_fa.py'

rule insert_fasta_to_tsv:
    conda:
        '../envs/pyGibson.yml'
    input:
        'output/{var_name}/inserts/complete_inserts.fa'
    output:
        'output/{var_name}/inserts/complete_inserts.tsv'
    script:'../scripts/fasta_to_tsv.py'


rule fasta_checksum:
    input:
        'output/{var_name}/inserts/complete_inserts.fa'
    output:
        'output/{var_name}/inserts/complete_inserts.md5sum'
    shell:'''
    md5sum {input} > {output}
    '''

rule summary_notebook:
    conda:
        '../envs/notebook.yml'
    input:
        insert_genbank=lambda wildcards: expand(
            'output/{var_name}/inserts/genbank_files/{p_name}.insert.gb',
            p_name=get_all_p_names(wildcards), allow_missing=True
        ),
        vr_defs_file=lambda wildcards: config['variable_region_definitions'][wildcards.var_name]
    log:
        notebook='output/{var_name}/inserts/insert_summary.ipynb'
    params:
        vr_defs=lambda wildcards: variable_regions[wildcards.var_name],
        cutters=config['prohibited_restriction_enzyme_recognition_sites']
    notebook:
        "../notebooks/insert_summary.py.ipynb"











