

rule simulate_construct_assembly_t7:
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
        t7_init=directory('output/{var_name}/constructs/T7_initiation_series'),
        t7_term=directory('output/{var_name}/constructs/T7_termination_series')
    script:'../scripts/simulate_assembly.py'


rule simulate_construct_assembly_tac:
    conda:
        '../envs/pyGibson.yml'
    input:
        inserts=lambda wildcards: expand(
            'output/{var_name}/inserts/genbank_files/{p_name}.insert.gb',
            p_name=get_all_p_names(wildcards), allow_missing=True
        ),
        tac_backbone=config['backbones']['tac'],
        tac_init_primers='output/{var_name}/sequences/tac_initiation_series_primers.fa',
        tac_term_primers='output/{var_name}/sequences/tac_termination_series_primers.fa',
        t7_init='output/{var_name}/constructs/T7_initiation_series',
        t7_term='output/{var_name}/constructs/T7_termination_series'
    output:
        tac_init=directory('output/{var_name}/constructs/Tac_initiation_series'),
        tac_term=directory('output/{var_name}/constructs/Tac_termination_series')
    script:'../scripts/simulate_assembly.py'


rule simulate_all_constructs:
    input:
        tac_init='output/{var_name}/constructs/Tac_initiation_series',
        tac_term='output/{var_name}/constructs/Tac_termination_series',
        t7_init='output/{var_name}/constructs/T7_initiation_series',
        t7_term='output/{var_name}/constructs/T7_termination_series'
    output:
        'output/{var_name}/constructs/.all_costructs.done'
    shell:'''
    touch {output}
    '''

rule make_tac_primers:
    conda:
        '../envs/pyGibson.yml'
    input:
        t7_init='output/{var_name}/constructs/T7_initiation_series',
        t7_term='output/{var_name}/constructs/T7_termination_series',
        tac_backbone=config['backbones']['tac'],
        protocol='output/{var_name}/protocols/t7_init_inserts.tsv'
    output:
        init_primers='output/{var_name}/sequences/tac_initiation_series_primers.fa',
        term_primers='output/{var_name}/sequences/tac_termination_series_primers.fa'
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

        





