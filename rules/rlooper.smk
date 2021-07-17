
# files produced by one Rlooper run

rule download_rlooper_and_compile_rlooper:
    output:
        'submodules/rlooper/bin/rlooper'
    shell:'''
    mkdir -p submodules
    cd submodules
    rm -rf rlooper
    git clone https://github.com/chedinlab/rlooper.git
    cd rlooper
    make all
    '''

rule seperate_VR_into_individual_fastas_format_header_for_rlooper:
    conda:
        '../envs/python.yml'
    input:
        fasta='output/{var_name}/files/{var_name}.fasta',
        rlooper='submodules/rlooper/bin/rlooper'
    output:
        'output/rlooper/{var_name}/fasta/{record}.fa'
    script: '../scripts/split_fasta.py'


rule rlooper_sequence:
    input:
        fasta='output/rlooper/{var_name}/fasta/{record}.fa',
        rlooper='submodules/rlooper/bin/rlooper'
    output:
        expand(
            'output/rlooper/{var_name}/completed_runs/{record}/{record}_{rlooper_suffix}',
            rlooper_suffix=RLOOPER_FILE_SUFFI, allow_missing=True
            )
        # snakemake thinks this output is a directory for some reason
        # throws error if not marked as such
    params:
        superhelicity='-0.07',
        domain_size='auto',
        minlength='30',  # value used in R-looper paper
        out_path=lambda wildcards: format(
            'output/rlooper/{}/completed_runs/{}/{}'.format(
                wildcards.var_name, wildcards.record, 
                wildcards.record
                )
        ),
        out_dir=lambda wildcards: format(
            'output/rlooper/{}/completed_runs/{}'.format(
                wildcards.var_name, wildcards.record, 
                )
        )

    shell:'''
    mkdir -p {params.out_dir}
    chmod 777 {input.rlooper}
    ./{input.rlooper} {input.fasta} {params.out_path} --N {params.domain_size} \
    --sigma {params.superhelicity} --localaverageenergy --minlength {params.minlength}
    '''


# use rule rlooper_sequence as random_rlooper with:
#     input:
#         fasta='testing/rlooper_benchmarking/random_fasta/{rand_fasta}.{length}.fa',
#         rlooper='submodules/rlooper/bin/rlooper'
#     output:
#         expand(
#             'testing/rlooper_benchmarking/completed_runs/{rand_fasta}.{length}/{rand_fasta}.{length}_{rlooper_suffix}',
#             rlooper_suffix=RLOOPER_FILE_SUFFI, allow_missing=True
#             )
#     params:
#         superhelicity='-0.07',
#         domain_size='auto',
#         minlength='30',  # value used in R-looper paper
#         out_path=lambda wildcards: format(
#             'testing/rlooper_benchmarking/completed_runs/{}.{}/{}.{}'.format(
#                 wildcards.rand_fasta, wildcards.length, wildcards.rand_fasta,
#                 wildcards.length
#             )
#         ),
#         out_dir=lambda wildcards: format(
#             'testing/rlooper_benchmarking/completed_runs/{}.{}'.format(
#                 wildcards.rand_fasta, wildcards.length
#                 )
#         )




