

rule simulate_construct_assembly:
    conda:
        '../envs/pyGibson.yml'
    input:
        inserts=lambda wildcards: expand(
            'output/{var_name}/inserts/genbank_files/{p_name}.insert.gb',
            p_name=get_all_p_names(wildcards), allow_missing=True
        ),
        t7_init_backbone=config['backbones']['t7_initiation'],
        t7_term_backbone=config['backbones']['t7_termination'],
        tac_backbone=config['backbones']['tac'],
        initiator=config['termination_series_initiator']
    output:
        directory('output/{var_name}/constructs')
    script:'../scripts/simulate_assembly.py'


rule make_tac_primers:
    conda:
        '../envs/pyGibson.yml'
    input:
        constructs_dir='output/{var_name}/constructs',
        tac_backbone=config['backbones']['tac']
    output:
        init_primers='output/{var_name}/sequences/tac_initiation_series_primers.fa',
        term_primers='output/{var_name}/sequences/tac_termination_series_primers.fa'
    params:
        # TODO: add these names to config
        t7_init_series='T7_initiation_series',
        t7_term_series='T7_termination_series',
    notebook:'../notebooks/tac_series_primers.ipynb'


