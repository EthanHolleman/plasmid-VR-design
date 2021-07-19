

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
        initiator=config['termination_series_initiator'],
        tac_init_primers='output/{var_name}/sequences/tac_initiation_series_primers.fa',
        tac_term_primers='output/{var_name}/sequences/tac_termination_series_primers.fa'
    output:
        directory('output/{var_name}/constructs')
    script:'../scripts/simulate_assembly.py'


rule make_tac_primers:
    conda:
        '../envs/pyGibson.yml'
    input:
        constructs_dir='output/{var_name}/constructs',
        tac_backbone=config['backbones']['tac'],
        protocol='output/{var_name}/protocols/t6_init_inserts.tsv'
    output:
        init_primers='output/{var_name}/sequences/tac_initiation_series_primers.fa',
        term_primers='output/{var_name}/sequences/tac_termination_series_primers.fa'
    params:
        # TODO: add these names to config
        t7_init_series='T7_initiation_series',
        t7_term_series='T7_termination_series',
    notebook:'../notebooks/tac_series_primers.ipynb'


rule write_protocol_notebook:
    conda:
        '../envs/pyGibson.yml'
    input:
        constructs_dir='output/{var_name}/constructs',
        t7_init_backbone = config['backbones']['t7_initiation'],
        inserts=lambda wildcards: expand(
            'output/{var_name}/inserts/genbank_files/{p_name}.insert.gb',
            allow_missing=True, p_name=get_all_p_names(wildcards)
        )
    params:
        insert_concentration=None,  # placeholder value currently unknown ng / ul
        backbone_concentration=None  # placeholder value ng / ul
    output:
        library='output/{var_name}/protocols/t7_init_library.tsv',
        inserts='output/{var_name}/protocols/t7_init_inserts.tsv'
    log:
        notebook='output/{var_name}/protocols/protocol_notebook.ipynb'
    notebook: '../notebooks/neb_gibson_protocol.ipynb'

        





