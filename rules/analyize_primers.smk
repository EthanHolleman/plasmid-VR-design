from pathlib import Path
MFE_PATH = 'software/mfeprimer-3.2.2-linux-amd64'

rule download_MFEprimer_3_0:
    output:
        MFE_PATH
    params:
        url='https://github.com/quwubin/MFEprimer-3.0/releases/download/v3.2.2/mfeprimer-3.2.2-linux-amd64.gz',
        gz_name='mfeprimer-3.2.2-linux-amd64.gz'
    shell:'''
    mkdir -p software
    wget {params.url}
    gzip -d {params.gz_name}
    chmod 777 {params.gz_name}
    '''

rule analyize_dimers:
    input:
        fasta='input/primers/fasta'
        MFE=MFE_PATH
    output:
        'output/{wildcards}/dimer_analysis.out'
    shell:'''
    ./{input.MFE} dimer -i {input.fasta} -o {output}
    '''

rule analyize_hairpins:
    input:
        fasta='input/primers/fasta',
        MFE=MFE_PATH
    output:
        'output/{wildcards}/dimer_analysis.out'
    shell:'''
    ./{input.MFE} hairpin -i {input.fasta} -o {output}
    '''

    






    