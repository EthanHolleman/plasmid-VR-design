
rule random_fasta_record_file:
    conda:
        '../envs/python.yml'
    output:
        'testing/rlooper_benchmarking/random_fasta/{rand_fasta}.fa'
    params:
        record_name=lambda wildcards: wildcards.rand_fasta,
        length=200
    script:'../scripts/random_seq_gen.py'


rule plot_rlooper_rand_seq_distrabution:
    conda:
        '../envs/R.yml'
    input:
        expand(
            'testing/rlooper_benchmarking/completed_runs/{rand_fasta}/{rand_fasta}_avgG.wig',
            rand_fasta=RAND_SEQ_NAMES
        )
    output:
        'testing/rlooper_benchmarking/plots/rand_seq_LAE_dist.png'
    script:'../scripts/plot_rlooper_energy_dist.R'

