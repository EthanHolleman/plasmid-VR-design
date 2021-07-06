rule make_anchor_sequence:
    conda:
        '../envs/pyGibson.yml'
    input:
        variable_regions='output/{var_name}/sequences/plasmid_sequences.fasta',
    params:
        backbones=config['backbones'],
        prohibited_cutters=config['prohibited_restriction_enzyme_recognition_sites'],
        max_attempts=config['anchor_sequence']['max_attempts'],
        anchor_length=config['anchor_sequence']['anchor_length'],
        min_melting_temp=config['anchor_sequence']['min_melting_temp']
    output:
        'output/{var_name}/sequences/anchor.fa'
    script:'../scripts/make_anchor_region.py'

