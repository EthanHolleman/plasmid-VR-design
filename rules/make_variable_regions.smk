

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


rule agg_seqs:
    input:
        expand(
            'output/{var_name}/files/{p_name}/{id_num}/parsedRNA/{p_name}.tsv',
            id_num=CASE_RANGE,
            allow_missing=True
            )
    output:
        'output/test/{var_name}.{p_name}.done'
    shell:'''
    touch {output}
    '''
# then need a script to rank all tested sequences and get the best
# ones back out

def sequece_length(wildcards):
    table = vr_tables[wildcards.var_name]
    length = table.loc[table['name'] == wildcards.p_name]['length']
    return length


rule rank_sequences:
    input:
        
    output:
        ranked='output/{var_name}/files/{p_name}/{id_num}/ranked/{p_name}.fasta',
        #best='output/{var_name}/files/{p_name}/{id_num}/ranked/{p_name}.fasta'
    params:
        # length from pandas dataframe
        rlooper= lambda wildcards: f'output/expectations/rlooper/rlooper_expect.{sequence_length(wildcards)}.tsv',
        SPOTRNA=lambda wildcards: f'output/expectations/SPOTRNA/spotRNA_expect.{sequece_length(wildcards)}.tsv'
    script:'tello.py'









