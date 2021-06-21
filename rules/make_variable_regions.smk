

# rule generate_variable_regions:
#     conda:
#         '../envs/python.yml'
#     output:
#         fasta='output/{var_name}/files/{var_name}.fasta',
#         tsv='output/{var_name}/files/{var_name}.tsv'
#     params:
#         input_file = lambda wildcards: variable_regions[wildcards.var_name]
#     script:'../scripts/varmids/varmids.py'
NUM_CASES = 10
CASE_RANGE = range(1, NUM_CASES+1)


def get_p_record(wildcards):
    return vr_tables[wildcards.var_name].loc[
        vr_tables[wildcards.var_name]['name'] == wildcards.p_name]


rule generate_random_seq:
    conda:
        '../envs/python.yml'
    output:
        fasta=expand(
            'output/{var_name}/files/{p_name}/{id_num}/{p_name}.{id_num}.fasta',
            id_num=CASE_RANGE, allow_missing=True
        ),
        tsv=expand(
            'output/{var_name}/files/{p_name}/{id_num}/{p_name}.{id_num}.tsv',
            id_num=CASE_RANGE, allow_missing=True
        )
    params:
        p_record=lambda wildcards: get_p_record(wildcards),
        num_cases=NUM_CASES  # change to config param
    script:
        '../scripts/varmids/varmids/variable_region_maker.py'



rule run_SPOT_RNA_prediction_on_variable_regions:
    # SPOT RNA insists on adding the fasta header to output files
    # so had to write quick python scruipt to rename back to what
    # snakemake is expecting after run is complete
    conda:
        '../envs/SPOT-RNA.yml'
    input:
        spot='submodules/SPOT-RNA/SPOT-RNA-models',
        fasta='output/{var_name}/files/{p_name}/{id_num}/{p_name}.{id_num}.fasta'
    output:
        expand(
            'output/{var_name}/files/{p_name}/{id_num}/SPOTRNA/{p_name}-{id_num}.{SPOT_EXT}',
            SPOT_EXT=SPOT_RNA_EXTS, allow_missing=True)
    params:
        spot_exe=lambda wildcards: os.path.join('submodules/SPOT-RNA', 'SPOT-RNA.py'),
        output_dir=lambda wildcards: f'output/{wildcards.var_name}/files/{wildcards.p_name}/{wildcards.id_num}/SPOTRNA',
        rename_script='scripts/truncate_rename.py'
    shell:'''
    mkdir -p {params.output_dir}
    python3 {params.spot_exe} --inputs {input.fasta} --outputs {params.output_dir}
    python3 scripts/truncate_rename.py {params.output_dir} --index 1
    '''


rule annotate_RNAss_predictions:
    input:
        bpseq='output/{var_name}/files/{p_name}/{id_num}/SPOTRNA/{p_name}-{id_num}.bpseq',
        bpRNA='submodules/bpRNA'
    output:
        'output/{var_name}/files/{p_name}/{id_num}/bpRNA/{p_name}-{id_num}.st'
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
        'output/{var_name}/files/{p_name}/{id_num}/bpRNA/{p_name}-{id_num}.st'
    output:
        tsv_summary='output/{var_name}/files/{p_name}/{id_num}/parsedRNA/{p_name}.tsv',
    script:'../scripts/bpRNA_parser.py'


rule rlooper_variable_region:
    input:
        fasta='output/{var_name}/files/{p_name}/{id_num}/{p_name}.{id_num}.fasta',
        rlooper='submodules/rlooper/bin/rlooper'
    output:
        expand(
            'output/{var_name}/files/{p_name}/{id_num}/{p_name}.{id_num}.{rlooper_suffix}',
            rlooper_suffix=RLOOPER_FILE_SUFFI, allow_missing=True
            )
    params:
        superhelicity='-0.07',
        domain_size='auto',
        minlength='30',  # value used in R-looper paper
        out_path=lambda wildcards: 
            f'output/{wildcards.var_name}/files/{wildcards.p_name}/{wildcards.id_num}/{wildcards.p_name}.{wildcards.id_num}'
        ),
        out_dir=lambda wildcards:
            f'output/{wildcards.var_name}/files/{wildcards.p_name}/{wildcards.id_num}'
        )
    shell:'''
    mkdir -p {params.out_dir}
    chmod 777 {input.rlooper}
    ./{input.rlooper} {input.fasta} {params.out_path} --N {params.domain_size} \
    --sigma {params.superhelicity} --localaverageenergy --minlength {params.minlength}
    '''


rule predict_Rlooper_statistics:
    input:
        expand(
            'output/{var_name}/files/{p_name}/{id_num}/{p_name}.{id_num}.{rlooper_suffix}',
            id_num=CASE_RANGE, allow_missing=True, rlooper_suffix=RLOOPER_FILE_SUFFI
            )
    output:
        'output/{var_name}/{p_name}.predict_rlooper.done'
    shell:'''
    touch {output}
    '''


rule predict_RNA_secondary_structure:
    input:
        expand(
            'output/{var_name}/files/{p_name}/{id_num}/parsedRNA/{p_name}.tsv',
            id_num=CASE_RANGE,
            allow_missing=True
            )
    output:
        'output/{var_name}/{p_name}.predict_RNAss.done'
    shell:'''
    touch {output}
    '''


# then need a script to rank all tested sequences and get the best
# ones back out

def sequece_length(wildcards):
    # calculate variable region sequence length from wilcard input
    table = vr_tables[wildcards.var_name]
    length = table.loc[table['name'] == wildcards.p_name]['length']
    return length


rule aggregate_sequence_metrics:
    input:
        rlooper_bpprop='output/{var_name}/files/{p_name}/{id_num}/{p_name}.{id_num}.avgG.wig',
        rlooper_lae='output/{var_name}/files/{p_name}/{id_num}/{p_name}.{id_num}.avgG.wig',
        RNAss='output/{var_name}/files/{p_name}/{id_num}/parsedRNA/{p_name}.tsv'
    output:
        'output/{var_name}/files/{p_name}/{id_num}/aggregatedMetrics/{p_name}.tsv'
    script:''


rule combine_aggregated_sequence_metrics:
    input:
        expand(
            'output/{var_name}/files/{p_name}/{id_num}/aggregatedMetrics/{p_name}.tsv',
             id_num=CASE_RANGE, allow_missing=True
        ),
    output:
        'output/{var_name}/files/{p_name}/aggregations/{p_name}.all_aggregated_metrics.tsv'
    script:''


rule rank_sequences:
    input:
        'output/{var_name}/files/{p_name}/aggregations/{p_name}.all_aggregated_metrics.tsv'
    params:
        # length from pandas dataframe
        rlooper= lambda wildcards: f'output/expectations/rlooper/rlooper_expect.{sequence_length(wildcards)}.tsv',
        SPOTRNA=lambda wildcards: f'output/expectations/SPOTRNA/spotRNA_expect.{sequece_length(wildcards)}.tsv'
    script:'tello.py'









