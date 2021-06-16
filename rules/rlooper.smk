
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
        directory('output/rlooper/{var_name}')
    script: '../scripts/split_fasta.py'

