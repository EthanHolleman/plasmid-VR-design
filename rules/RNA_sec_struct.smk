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


rule download_bpRNA:
    output:
        directory('submodules/bpRNA')
    params:
        url='https://github.com/EthanHolleman/bpRNA.git'
    shell:'''
    cd submodules
    git clone {params.url}
    '''


rule download_rna_tools:
    output:
        directory('submodules/rna-tools')
    params:
        url='https://github.com/mmagnus/rna-tools'
    shell:'''
    cd submodules
    git clone {params.url}
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


rule run_vienna_RNA_prediction:
    conda:
        '../envs/viennaRNA.yml'
    input:
        'output/{var_name}/files/{var_name}.fasta'
    output:
        'output/RNA_sec_struct/viennaRNA/{var_name}.out'
    shell:'''
    RNAfold -p -d2 --noLP < {input} > {output}
    rm *.ps  # RNAfold creates a bunch of junk files with .ps extension 
    '''


rule run_SPOT_RNA_prediction:
    # SPOT RNA insists on adding the fasta header to output files
    # so had to write quick python scruipt to rename back to what
    # snakemake is expecting after run is complete
    conda:
        '../envs/SPOT-RNA.yml'
    input:
        spot='submodules/SPOT-RNA/SPOT-RNA-models',
        fasta='output/rlooper/{var_name}/fasta/{record}.fa'
    output:
        expand(
            'output/RNA_sec_struct/SPOT-RNA/{var_name}/{record}/{record}.{SPOT_EXT}',
            SPOT_EXT=SPOT_RNA_EXTS, allow_missing=True)
    params:
        spot_exe=lambda wildcards: os.path.join('submodules/SPOT-RNA', 'SPOT-RNA.py'),
        output_dir=lambda wildcards: f'output/RNA_sec_struct/SPOT-RNA/{wildcards.var_name}/{wildcards.record}',
        rename_script='scripts/truncate_rename.py'
    shell:'''
    mkdir -p {params.output_dir}
    python3 {params.spot_exe} --inputs {input.fasta} --outputs {params.output_dir}
    python {params.rename_script} {params.output_dir} --index 1
    '''


rule run_bpRNA:
    input:
        bpseq='output/RNA_sec_struct/SPOT-RNA/{var_name}/{record}/{record}.bpseq',
        bpRNA='submodules/bpRNA'
    output:
        'output/RNA_sec_struct/bpRNA/{var_name}/st/{record}.st'
    params:
        bp_script='submodules/bpRNA/bpRNA.pl',
        script_output=lambda wildcards: f'{wildcards.record}.st'  # script just spews into CWD 
    shell:'''
    perl {params.bp_script} {input.bpseq}
    mv {params.script_output} {output} && [[ -s {output} ]]
    '''

rule parse_bpRNA_st_files:
    input:
        'output/RNA_sec_struct/bpRNA/{var_name}/{record}.st'
    output:
        'output/RNA_sec_struct/bpRNA/{var_name}/tsv/{record}.tsv'
    script:'../scripts/bpRNA_parser.py'


# Run SPOT-RNA program on random nucleotide sequences of various lengths 
# in order to establish expectations about what "good" sequences might
# look like. This is also done for rlooper, same idea here.
rule run_random_SPOT_RNA_prediction:
    conda:
        '../envs/SPOT-RNA.yml'
    input:
        fasta='testing/rlooper_benchmarking/random_fasta/{rand_fasta}.{length}.fa',
        spot='submodules/SPOT-RNA/SPOT-RNA-models'
    output:
        expand(
            'testing/RNA_sec_struct/SPOT-RNA/{rand_fasta}.{length}/{rand_fasta}.{length}.{SPOT_EXT}',
            SPOT_EXT=SPOT_RNA_EXTS, allow_missing=True
        )
    params:
        spot_exe=lambda wildcards: os.path.join('submodules/SPOT-RNA', 'SPOT-RNA.py'),
        output_dir=lambda wildcards: f'testing/RNA_sec_struct/SPOT-RNA/{wildcards.rand_fasta}.{wildcards.length}',
        rename_script='scripts/truncate_rename.py',
        stem=lambda wildcards: f'{wildcards.rand_fasta}.{wildcards.length}'
    shell:'''
    mkdir -p {params.output_dir}
    python3 {params.spot_exe} --inputs {input.fasta} --outputs {params.output_dir}
    python {params.rename_script} {params.output_dir} --constant_stem {params.stem}
    '''

use rule run_bpRNA as run_random_bpRNA with:
    input:
        bpseq='testing/RNA_sec_struct/SPOT-RNA/{rand_fasta}.{length}/{rand_fasta}.{length}.bpseq',
        bpRNA='submodules/bpRNA'
    output:
        'testing/RNA_sec_struct/bpRNA/st/{rand_fasta}.{length}.st'
    params:
        bp_script='submodules/bpRNA/bpRNA.pl',
        script_output=lambda wildcards: f'{wildcards.rand_fasta}.{wildcards.length}.st'  # script just spews into CWD 


use rule parse_bpRNA_st_files as parse_random_bpRNA with:
    input:
        'testing/RNA_sec_struct/bpRNA/st/{rand_fasta}.{length}.st'
    output:
        'testing/RNA_sec_struct/bpRNA/tsv/{rand_fasta}.{length}.tsv'



    


    




    