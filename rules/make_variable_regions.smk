

def get_p_record(wildcards):
    return vr_tables[wildcards.var_name].loc[
        vr_tables[wildcards.var_name]['name'] == wildcards.p_name]


def get_all_p_names(wildcards):
    return list(vr_tables[wildcards.var_name]['name'])

def get_all_p_names_RC(wildcards):
    # get names of variable regions with reverse complement arg true
    table = vr_tables[wildcards.var_name]
    l = list(table.loc[table['reverse_complement'] == 1]['name'])
    return l


# this is where snakemake needs to become aware of revsere somplement files 
# currently I do not think they are being produced. 
# one option would be to just expand the user input to include reverse
# complemented seqs. Maybe better because slighly more visible in this way
rule generate_random_seq:
    conda:
        '../envs/python.yml'
    output:
        fasta=expand(
            'output/{var_name}/files/{p_name}/candidate_seqs/{id_num}/{p_name}.{id_num}.fasta',
            id_num=CASE_RANGE, allow_missing=True
        ),
        tsv=expand(
            'output/{var_name}/files/{p_name}/candidate_seqs/{id_num}/{p_name}.{id_num}.tsv',
            id_num=CASE_RANGE, allow_missing=True
        )
    params:
        p_record=lambda wildcards: get_p_record(wildcards),
        num_cases=NUM_CASES  # change to config param
    script:
        '../scripts/varmids/varmids.py'


rule SPOT_RNA_prediction_on_variable_regions:
    # SPOT RNA insists on adding the fasta header to output files
    # so had to write quick python scruipt to rename back to what
    # snakemake is expecting after run is complete
    conda:
        '../envs/SPOT-RNA.yml'
    input:
        spot='submodules/SPOT-RNA/SPOT-RNA-models',
        fasta='output/{var_name}/files/{p_name}/candidate_seqs/{id_num}/{p_name}.{id_num}.fasta'
    output:
        expand(
            'output/{var_name}/files/{p_name}/candidate_seqs/{id_num}/SPOTRNA/{p_name}-{id_num}.{SPOT_EXT}',
            SPOT_EXT=SPOT_RNA_EXTS, allow_missing=True)
    params:
        spot_exe=lambda wildcards: os.path.join('submodules/SPOT-RNA', 'SPOT-RNA.py'),
        output_dir=lambda wildcards: f'output/{wildcards.var_name}/files/{wildcards.p_name}/candidate_seqs/{wildcards.id_num}/SPOTRNA',
        rename_script='scripts/truncate_rename.py'
    shell:'''
    mkdir -p {params.output_dir}
    python3 {params.spot_exe} --inputs {input.fasta} --outputs {params.output_dir}
    python3 scripts/truncate_rename.py {params.output_dir} --index 1
    '''


rule annotate_RNAss_predictions:
    input:
        bpseq='output/{var_name}/files/{p_name}/candidate_seqs/{id_num}/SPOTRNA/{p_name}-{id_num}.bpseq',
        bpRNA='submodules/bpRNA'
    output:
        'output/{var_name}/files/{p_name}/candidate_seqs/{id_num}/bpRNA/{p_name}-{id_num}.st'
    params:
        output_dir=lambda wildcards: f'output/{wildcards.var_name}/files/{wildcards.p_name}/{wildcards.id_num}/bpRNA',
        bp_script='submodules/bpRNA/bpRNA.pl',
        script_output=lambda wildcards: f'{wildcards.p_name}-{wildcards.id_num}.st'  # script just spews into CWD 
    shell:'''
    perl {params.bp_script} {input.bpseq}
    mkdir -p {params.output_dir}
    mv {params.script_output} {output} && [[ -s {output} ]]
    '''


rule parse_bpRNA_annotations:
    input:
        'output/{var_name}/files/{p_name}/candidate_seqs/{id_num}/bpRNA/{p_name}-{id_num}.st'
    output:
        tsv_summary='output/{var_name}/files/{p_name}/candidate_seqs/{id_num}/parsedRNA/{p_name}.tsv',
    script:'../scripts/bpRNA_parser.py'


rule format_fasta_header_for_rlooper:
    input:
        'output/{var_name}/files/{p_name}/candidate_seqs/{id_num}/{p_name}.{id_num}.fasta'
    output:
        'output/{var_name}/files/{p_name}/candidate_seqs/{id_num}/{p_name}.{id_num}.rlooper_ready.fasta'
    params:
        record_name = lambda wildcards: f'{wildcards.p_name}.{wildcards.id_num}'
    script:'../scripts/format_fa_header_for_rlooper.py'


rule rlooper_variable_region:
    input:
        fasta='output/{var_name}/files/{p_name}/candidate_seqs/{id_num}/{p_name}.{id_num}.rlooper_ready.fasta',
        rlooper='submodules/rlooper/bin/rlooper'
    output:
        expand(
            'output/{var_name}/files/{p_name}/candidate_seqs/{id_num}/{p_name}.{id_num}_{rlooper_suffix}',
            rlooper_suffix=RLOOPER_FILE_SUFFI, allow_missing=True
            )
    params:
        superhelicity='-0.07',
        domain_size='auto',
        minlength='30',  # value used in R-looper paper
        out_path=lambda wildcards: f'output/{wildcards.var_name}/files/{wildcards.p_name}/candidate_seqs/{wildcards.id_num}/{wildcards.p_name}.{wildcards.id_num}',
        out_dir=lambda wildcards: f'output/{wildcards.var_name}/files/{wildcards.p_name}/candidate_seqs/{wildcards.id_num}'
    shell:'''
    mkdir -p {params.out_dir}
    chmod 777 {input.rlooper}
    ./{input.rlooper} {input.fasta} {params.out_path} --N {params.domain_size} \
    --sigma {params.superhelicity} --localaverageenergy --minlength {params.minlength}
    '''


def sequence_length(wildcards):
    # calculate variable region sequence length from wilcard input
    table = vr_tables[wildcards.var_name]
    length = int(table.loc[table['name'] == wildcards.p_name]['length'])
    return length


rule aggregate_sequence_metrics:
    conda:
        '../envs/python.yml'
    input:
        bp_prob='output/{var_name}/files/{p_name}/candidate_seqs/{id_num}/{p_name}.{id_num}_bpprob.wig',
        lae='output/{var_name}/files/{p_name}/candidate_seqs/{id_num}/{p_name}.{id_num}_avgG.wig',
        RNAss='output/{var_name}/files/{p_name}/candidate_seqs/{id_num}/parsedRNA/{p_name}.tsv',
        fasta='output/{var_name}/files/{p_name}/candidate_seqs/{id_num}/{p_name}.{id_num}.fasta',
        tsv='output/{var_name}/files/{p_name}/candidate_seqs/{id_num}/{p_name}.{id_num}.tsv',
        rlooper_expect = lambda wildcards: f'output/expectations/rlooper/rlooper_expect.{sequence_length(wildcards)}.tsv',
        RNAss_expect =lambda wildcards: f'output/expectations/SPOTRNA/spotRNA_expect.{sequence_length(wildcards)}.tsv'
    output:
        'output/{var_name}/files/{p_name}/candidate_seqs/{id_num}/aggregatedMetrics/{p_name}.tsv'
    params:
        length = lambda wildcards: sequence_length(wildcards),
        p_name = lambda wildcards: wildcards.p_name,
        id_num = lambda wildcards: wildcards.id_num,
    script:'../scripts/agg_seq_metrics.py'


rule concatenate_seq_metrics:
    conda:
        '../envs/python.yml'
    input:
        expand(
            'output/{var_name}/files/{p_name}/candidate_seqs/{id_num}/aggregatedMetrics/{p_name}.tsv',
            id_num=CASE_RANGE, allow_missing=True
        )
    output:
        'output/{var_name}/files/{p_name}/concatAggMetrics/{p_name}.tsv'
    script:'../scripts/concate_tsvs.py'


rule rank_and_select_sampled_sequences:
    conda:
        '../envs/python.yml'
    input:
        'output/{var_name}/files/{p_name}/concatAggMetrics/{p_name}.tsv'
    output:
        fasta='output/{var_name}/files/{p_name}/rankedSeqs/{p_name}.top_seq.fasta',
        tsv='output/{var_name}/files/{p_name}/rankedSeqs/{p_name}.top_seq.tsv',
        ranked='output/{var_name}/files/{p_name}/rankedSeqs/{p_name}.all_ranked.tsv'
    params:
        config_dict = config
    script:'../scripts/rank_seqs.py'


rule reverse_complement_top_sequences:
    conda:
        '../envs/python.yml'
    input:
        fasta='output/{var_name}/files/{p_name}/rankedSeqs/{p_name}.top_seq.fasta',
        tsv='output/{var_name}/files/{p_name}/rankedSeqs/{p_name}.top_seq.tsv'
    output:
        fasta='output/{var_name}/files/{p_name}/rankedSeqs/{p_name}-RC.top_seq.fasta',
        tsv='output/{var_name}/files/{p_name}/rankedSeqs/{p_name}-RC.top_seq.tsv'
    script:'../scripts/reverse_complement.py'


rule concatenate_top_ranked_sequences_fasta:
    input:
        forward_seqs=lambda wildcards: expand(
            'output/{var_name}/files/{p_name}/rankedSeqs/{p_name}.top_seq.fasta',
            p_name=get_all_p_names(wildcards), allow_missing=True
        ),
        reverse_complement=lambda wildcards: expand(
        'output/{var_name}/files/{p_name}/rankedSeqs/{p_name}-RC.top_seq.fasta',
        p_name=get_all_p_names_RC(wildcards), allow_missing=True
        )    
    output:
        'output/{var_name}/sequences/plasmid_sequences.fasta'
    shell:'''
    cat {input} > {output}
    '''


rule concatenate_to_ranked_sequences_tsv:
    conda:
        '../envs/python.yml'
    input:
        forward_seqs=lambda wildcards: expand(
            'output/{var_name}/files/{p_name}/rankedSeqs/{p_name}.top_seq.tsv',
            p_name=get_all_p_names(wildcards), allow_missing=True
        ),
        reverse_complement=lambda wildcards: expand(
        'output/{var_name}/files/{p_name}/rankedSeqs/{p_name}-RC.top_seq.fasta',
        p_name=get_all_p_names_RC(wildcards), allow_missing=True
        )
    output:
        'output/{var_name}/sequences/plasmid_sequences.tsv'
    script:'../scripts/concate_tsvs.py'

    




