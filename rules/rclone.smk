from pathlib import Path

rule backup_results:
    input:
        constructs='output/{var_name}/constructs'
    output:
        done_marker='output/{var_name}/.synced.done'
    params:
        rclone_exe=config['rclone']['exe'],
        remote_target=lambda wildcards: str(
            Path(
                config['rclone']['remote']
                ).joinpath(wildcards.var_name)
            ),
        local_target=lambda wildcards: str(Path('output').joinpath(wildcards.var_name))
    shell:'''
    {params.rclone_exe} sync {params.local_target} {params.remote_target} 
    touch {output.done_marker}
    '''

rule dont_clone:
    output:
        'output/{var_name}/.no_sync'
    shell:'''
    touch {output}
    '''
    