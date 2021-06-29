# Scripts

The scripts dir contains all scripts that are used within the workflow and
that are not submodules.

## Descriptions

- `agg_seq_metrics.py`: General helper to aggregate the output of programs that produce metrics used to access the quality of a given variable region into a tsv file.
- `bpRNA_parser`: Parse the output produced by the `bpRNA` program in order to extract easy to access metrics regarding the predicted RNA secondary structure of a sequence. Produces a tsv file which includes the proportion of the sequence predicted to be in hairpin and proportion of sequence predicted to be unpaired.
- `concat_tsvs.py`: General helper script for concatenating `tsv` files with the same headers.
- `expand_RC_seqs.py`: TODO
- `format_header_for_rlooper.py`: TODO
- `plot_rlooper_expect.R`: TODO
- `plot_rna_struct_expect.R`: TODO
- `primer3_config.py`: TODO
- `random_seq_gen.py`: Generate a random nucleotide sequence fasta file with R-looper ready header.
- `rank_seqs.py`: TODO
- `split_fasta.py`: Splits a multi-record fasta file into individual fasta files, each with one record. This was part of getting around the fact that some programs are ok with multi-record input and others are not and because of this it is easier to just do everything as single records.
- `truncate_rename.py`: TODO.


