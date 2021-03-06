
rule random_fasta_record_file:
    conda:
        '../envs/python.yml'
    output:
        'testing/rlooper_benchmarking/random_fasta/{rand_fasta}.{length}.fa'
    params:
        record_name=lambda wildcards: wildcards.rand_fasta,
        length=lambda wildcards: wildcards.length
    script:'../scripts/random_seq_gen.py'


rule plot_rlooper_rand_seq_distrabution:
    conda:
        '../envs/R.yml'
    input:
        ale=expand(
            'testing/rlooper_benchmarking/completed_runs/{rand_fasta}.{length}/{rand_fasta}.{length}_avgG.wig',
            rand_fasta=RAND_SEQ_NAMES, length=RAND_SEQ_LENS
        ),
        bpprob=expand(
            'testing/rlooper_benchmarking/completed_runs/{rand_fasta}.{length}/{rand_fasta}.{length}_bpprob.wig',
            rand_fasta=RAND_SEQ_NAMES, length=RAND_SEQ_LENS
        )
    output:
        plot='testing/rlooper_benchmarking/plots/rand_seq_LAE_dist.png',
        expect='testing/rlooper_benchmarking/plots/expectations.tsv'
    script:'../scripts/plot_rlooper_expect.R'


rule plot_spot_rna_rand_seq:
    conda:
        '../envs/R.yml'
    input:
        expand('testing/RNA_sec_struct/bpRNA/tsv/{rand_fasta}.{length}.tsv',
        rand_fasta=RAND_SEQ_NAMES, length=RAND_SEQ_LENS
        )
    output:
        plot='testing/RNA_sec_struct/plots/plot.png',
        expect='testing/RNA_sec_struct/plots/expectations.tsv'
    script:'../scripts/plot_rna_struct_expect.R'

