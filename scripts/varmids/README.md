# Varmids

Varmids scripts are used to generate variable region sequences from user
parameters. These parameters should be defined in a `tsv` and include
all the fields documented in the main `README.md`.

Varmids scripts, like almost all other scripts in the repo, are intended
to be executed from a snakemake rule. `varmids.py` should be called from
a snakemake rule. See rule [generate_variable_region_from_user_params](rules/make_variable_regions.smk) for how this script is executed from snakemake.

## Output

Running `varmids.py` with correct inputs will produce two files. First is a `fasta` formated file that contains the generated sequence
and the second is a `tsv` that also contains additional information
about the generated sequence including GC skew, content, clustering method, etc which is used for plotting as it is easier to parse.


