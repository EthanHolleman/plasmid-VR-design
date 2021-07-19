# Plasmid design for investigating sequence effect on R-loop initiation and termination

![example workflow](https://github.com/ethanholleman/plasmid-design/actions/workflows/tests.yml/badge.svg)


## Dependencies

The workflow is intended to be run on a unix system, I ran on a Ubuntu server.
The vast majority of software dependencies for the workflow are handled by
snakemake and conda. That being said you will need both snakemake and conda
installed.

However there are some programs that will need to be configured / installed
before running. 

### Perl and required libraries

The workflow currently assumes you have the perl language installed and that
the `Graph.pm` module is available for local import using `local::lib`. You can
install [Graph.pm](https://metacpan.org/dist/Graph/view/lib/Graph.pod)
using the command `cpanm Graph`. 


## Running the workflow

The workflow is executed via snakemake, the most basic command to do so it below.

`snakemake -j 1 --use-conda`

Edit the `run.sh` and `cluster.yml` files for your computing environment.
Then make `run.sh` executable and run the workflow by calling `./run.sh`.

### Specifying variable region parameters

Create a csv file in the `variable_defs` directory with the fields listed below.
| Field name                |       Description                                                                                                                                       |
| ------------------------- | -----------------------------------------------------------------------------------------------------------------------------------------------------   |
| name                      | Name of the variable region, no spaces and not underscores!                                                                                                                 |
| length                    | Length ofregion in nucleotides.  |
| gc_content                | GC content, float between 0 and 1. |
| gc_skew                   | GC skew, float between -1 and 1. |
| at_skew                   | AT skew, float between -1 and 1. |
| at_content                | AT content, float between 0 and 1. |
| cluster_length            | If the variable region is to have clusters of nucleotides set to the length of each cluster in nucleotides. Also then requires specifying the cluster_nuc field. |
| cluster_nuc               | Nucleotide that will compose clusters. |
| role                      | Short description of role of variable region. This will appear in fasta headers generated from this region. No spaces!                                  |

If a field is not specified is should be set to `NA`.
