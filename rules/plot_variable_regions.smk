


rule plot_variable_regions:
    conda:
        '../envs/R.yml'
    input:
        variable_regions='output/{var_name}/sequences/plasmid_sequences.tsv',
        rlooper_lae=lambda wildcards: expand(  # collect local average energy files
            'output/{var_name}/files/{p_name}/{id_num}/{p_name}.{id_num}_{rlooper_suffix}',
            var_name=wildcards.var_name, 
            p_name=get_all_p_names(wildcards),
            rlooper_suffix='avgG.wig',
            id_num=CASE_RANGE
        ),
        rlooper_bprob=lambda wildcards: expand(  # collect probability files
            'output/{var_name}/files/{p_name}/{id_num}/{p_name}.{id_num}_{rlooper_suffix}',
            var_name=wildcards.var_name, 
            p_name=get_all_p_names(wildcards),
            rlooper_suffix='bpprob.wig',
            id_num=CASE_RANGE
        ),

        parsed_RNA=lambda wildcards: expand(
            'output/{var_name}/files/{p_name}/{id_num}/parsedRNA/{p_name}.tsv',
            var_name=wildcards.var_name, p_name=get_all_p_names(wildcards),
            id_num=CASE_RANGE
        ),
        expectation_files = lambda wildcards: all_expectations(wildcards) 
    output:
        'output/{var_name}/plots/{var_name}.pdf'
    script:'../scripts/plot.R'
    