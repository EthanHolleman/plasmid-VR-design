

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
        initiator=config['termination_series_initiator']
    output:
        directory('output/{var_name}/constructs')
    script:'../scripts/simulate_assembly.py'
