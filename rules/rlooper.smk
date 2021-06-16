
# files produced by one Rlooper run

rule download_rlooper_and_compile_rlooper:
    output:
        'submodules/rlooper/bin/rlooper'
    shell:'''
    mkdir -p submodules
    cd submodules
    rm -rf rlooper
    git clone https://github.com/chedinlab/rlooper.git
    cd rlooper
    make all
    '''

rule seperate_VR_into_individual_fastas_format_header_for_rlooper:
    conda:
        '../envs/python.yml'
    input:
        fasta='output/{var_name}/files/{var_name}.fasta',
        rlooper='submodules/rlooper/bin/rlooper'
    output:
        directory('output/rlooper/{var_name}/fasta')
    script: '../scripts/split_fasta.py'


rule rlooper_sequence:
    input:
        fasta='output/rlooper/{var_name}/fasta/{record}.fa',
        rlooper='submodules/rlooper/bin/rlooper'
    output:
        'output/rlooper/{var_name}/completed_runs/{record}/{record}_{rlooper_suffix}'
    params:
        superhelicity='-0.07',
        domain_size='auto',
        out_dir=lambda wildcards: format(
            'output/rlooper/{}/completed_runs/{}/{}'.format(
                wildcards.var_name, wildcards.record, 
                wildcards.record
                )
        )

    shell:'''
    mkdir -p {output}
    chmod 777 {input.rlooper}
    ./{input.rlooper} {input.fasta} {params.out_dir} --N {params.domain_size} \
    --sigma {params.superhelicity}
    '''



