
def get_all_sequence_lengths(wildcards):
    table = vr_tables[wildcards.var_name]
    return list(set(table['length']))


def all_expectations_by_length(wildcards, expectation_path):
    lengths = get_all_sequence_lengths(wildcards)
    return expand(expectation_path, length=lengths)


def all_expectations(wildcards):
    rna = all_expectations_by_length(
        wildcards, 'output/expectations/SPOTRNA/spotRNA_expect.{length}.tsv'
        )
    rlooper = all_expectations_by_length(
        wildcards, 'output/expectations/rlooper/rlooper_expect.{length}.tsv'
    )
    return rna + rlooper


rule random_fasta_record_file:
    conda:
        '../envs/python.yml'
    output:
        'testing/rlooper_benchmarking/random_fasta/{rand_fasta}.{length}.fa'
    params:
        record_name=lambda wildcards: wildcards.rand_fasta,
        length=lambda wildcards: wildcards.length
    script:'../scripts/random_seq_gen.py'


rule plot_rlooper_rand_seq_distribution_parameter_informed:
    conda:
        '../envs/R.yml'
    input:
        ale = lambda wildcards: expand(
            'output/{var_name}/files/{p_name}/candidate_seqs/{id_num}/{p_name}.{id_num}_avgG.wig',
            p_name=get_all_p_names(wildcards),  id_num=CASE_RANGE, allow_missing=True
        ),
        bpprob=lambda wildcards: expand(
            'output/{var_name}/files/{p_name}/candidate_seqs/{id_num}/{p_name}.{id_num}_bpprob.wig',
            p_name=get_all_p_names(wildcards),  id_num=CASE_RANGE, allow_missing=True
        ),
        tsv_files=lambda wildcards: expand(
            'output/{var_name}/files/{p_name}/candidate_seqs/{id_num}/{p_name}.{id_num}.tsv',
            p_name=get_all_p_names(wildcards), id_num=CASE_RANGE, allow_missing=True
        ),
        random_ale=lambda wildcards: expand(
            expand('testing/rlooper_benchmarking/completed_runs/{rand_fasta}.{length}/{rand_fasta}.{length}_avgG.wig', allow_missing=True, length=get_all_sequence_lengths(wildcards)),
            zip, rand_fasta=RAND_SEQ_NAMES
        ),
        random_bpprob=lambda wildcards: expand(
            expand('testing/rlooper_benchmarking/completed_runs/{rand_fasta}.{length}/{rand_fasta}.{length}_bpprob.wig', allow_missing=True, length=get_all_sequence_lengths(wildcards)),
            zip, rand_fasta=RAND_SEQ_NAMES
        )
    output:
        plot='output/expectations/{var_name}/rlooper/rlooper_expect.png',
    params:
        p_names=lambda wildcards: get_all_p_names(wildcards),
        id_nums=lambda wildcards: CASE_RANGE
    script:'../scripts/plot_rlooper_expect_params.R'


rule plot_rlooper_rand_seq_distribution:
    conda:
        '../envs/R.yml'
    input:
        ale=expand(
            'testing/rlooper_benchmarking/completed_runs/{rand_fasta}.{length}/{rand_fasta}.{length}_avgG.wig',
            rand_fasta=RAND_SEQ_NAMES, allow_missing=True
        ),
        bpprob=expand(
            'testing/rlooper_benchmarking/completed_runs/{rand_fasta}.{length}/{rand_fasta}.{length}_bpprob.wig',
            rand_fasta=RAND_SEQ_NAMES, allow_missing=True
        )
    output:
        plot='output/expectations/rlooper/rlooper_expect.{length}.png',
        expect='output/expectations/rlooper/rlooper_expect.{length}.tsv'
    script:'../scripts/plot_rlooper_expect.R'


rule plot_spot_rna_rand_seq:
    conda:
        '../envs/R.yml'
    input:
        expand('testing/RNA_sec_struct/bpRNA/tsv/{rand_fasta}.{length}.tsv',
        rand_fasta=RAND_SEQ_NAMES, allow_missing=True
        )
    output:
        plot='output/expectations/SPOTRNA/spotRNA_expect.{length}.png',
        expect='output/expectations/SPOTRNA/spotRNA_expect.{length}.tsv'
    script:'../scripts/plot_rna_struct_expect.R'

