import os

rule download_SPOT_RNA:
    output:
        directory('submodules/SPOT-RNA/SPOT-RNA-models')
    params:
        url='https://github.com/jaswindersingh2/SPOT-RNA'
    shell:'''
    cd submodules
    git submodule add {params.url}
    cd SPOT-RNA
    wget 'https://www.dropbox.com/s/dsrcf460nbjqpxa/SPOT-RNA-models.tar.gz' || wget -O SPOT-RNA-models.tar.gz 'https://app.nihaocloud.com/f/fbf3315a91d542c0bdc2/?dl=1'
    tar -xvzf SPOT-RNA-models.tar.gz && rm SPOT-RNA-models.tar.gz
    '''


rule transcribe_dna:
    conda:
        '../envs/biopython.yml'
    input:
        'output/{var_name}/files/{var_name}.fasta'
    output:
        'output/{var_name}/files/{var_name}.RNA.fasta'
    shell:'''
    '''


rule run_SPOT_RNA_prediction:
    conda:
        '../envs/SPOT-RNA.yml'
    input:
        spot='submodules/SPOT-RNA/SPOT-RNA-models',
        fasta='output/{var_name}/files/{var_name}.fasta'
    output:
        directory('output/RNA_sec_struct/SPOT-RNA/{var_name}')
    params:
        spot_exe=lambda wildcards: os.path.join('submodules/SPOT-RNA', 'SPOT-RNA.py')
    shell:'''
    mkdir -p {output}
    python3 {params.spot_exe} --inputs {input.fasta} --outputs '{output}'
    '''


rule run_vienna_RNA_prediction:
    conda:
        '../envs/viennaRNA.yml'
    input:
        'output/{var_name}/files/{var_name}.fasta'
    output:
        'output/RNA_sec_struct/viennaRNA/{var_name}.out'
    shell:'''
    RNAfold -p -d2 --noLP < {input} > {output}
    '''


rule plot_vienna_RNA_prediction:
    conda:
        '../envs/viennaRNA.yml'
    input:
        'output/RNA_sec_struct/viennaRNA/{var_name}.out'
    output:
        'output/RNA_sec_struct/viennaRNA/{var_name}.plot'
    shell:'''
    RNAplot -o svg < {input}
    touch {output}
    '''




    