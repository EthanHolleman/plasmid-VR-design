
# if reverse complement is specified snakemake does not actually know that
# revese complement calcs will be included in the output. So for specified
# variable regions where RC is calculated need to explicitly tell snakemake
# to expect reverse complement input from R-looper so all plots can be
# # created
# def rlooper_calculations(wildcards, *args, **kwargs):
#     table = vr_tables[wildcards.var_name]
#     records = []
#     for index, row in table.iterrows():
#         records.append(row['name'])
#         if int(row['reverse_complement']) == 1:
#             records.append(f'{row["name"]}-RC')
#     return records

rule plot_variable_regions:
    conda:
        '../envs/R.yml'
    input:
        variable_regions='output/{var_name}/sequences/plasmid_sequences.tsv',
        rlooper_calculations=lambda wildcards: expand(
            'output/{var_name}/files/{p_name}/{id_num}/{p_name}.{id_num}_{rlooper_suffix}',
            var_name=wildcards.var_name, 
            p_name=get_all_p_names(wildcards),
            rlooper_suffix='avgG.wig',
            id_num=CASE_RANGE
        )
    output:
        'output/{var_name}/plots/{var_name}.pdf'
    script:'../scripts/plot.R'
    