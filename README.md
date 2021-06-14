# Plasmid design for investigating sequence effect on R-loop initiation and termination

![example workflow](https://github.com/ethanholleman/plasmid-design/actions/workflows/tests.yml/badge.svg)

## Grant language

Language taken directly from grant describing the work.

### Analyzing initiation sequences

```
To identify characteristics of R-loop initiation regions, we will syn-thesize 20 versions of a 200 bp initiation fragment cloned
immediately downstream of a promoter (T7 or the E. coli tac promoter). These 20 sequences will represent all combinations resulting from varying GC content (40, 50, 60 and 70%) and GC skew (0, 0.1, 0.2, 0.4, and 0.6). Each sequence will be followed by a constant 300 bp R-loop extension region (50% GC and 0.2 GC skew)
```

and


```
We will next design substrates to systematically test the notion, reported in the literature, that G clustering is essential for initiation [27]. For this, we will use favorable initiation regions (60% GC content with GC skew of 0.4) and arrange the ﬁxed number of guanines on the displaced strand in G clusters of sizes 1, 2, 3, or 4. Given reports that GA-rich sequences favor R-loop formation [22], we will test the eﬀect of AT skews for the ﬁrst time by introducing AT skews of 0, 0.2, and 0.4 in addition to G clustering.
```

### Analyzing termination sequences

```
We will add a constant extension region of 100 bp, and synthesize and clone a series of 200 bp sequences immediately thereafter. These potential termination regions will possess decreasing GC content (50, 40, 30%) and decreasing GC skew (0, -0.2, -0.4), thus producing increasingly unfavorable contexts for RNA:DNA base-pairing. We will test the eﬀect of C clustering and negative AT skews, as described above for initiation
```


## Specifing variable region parameters

Create a csv file in the `variable_defs` directory with the fields listed below.

- name: Name of the variable region, no spaces!
- length: Length ofregion in nucleotides.
- gc_content: GC content, float between 0 and 1.
- gc_skew: GC skew, float between 0 and 1.
- at_skew: AT skew, float between 0 and 1.
- at_content: AT content, float between 0 and 1.
- cluster_length: If the variable region is to have clusters of nucleotides set to the length of each cluster in nucleotides. Also then requires specifying the cluster_nuc field.
- cluster_nuc: Nucleotide that will compose clusters.
- clusering_mode: Way in which clusters should be placed in the variable regions. More details coming soon here.
- Role: Short description of role of variable region. This will appear in fasta headers generated from this region.

If a field is not specified is should be set to `NA`.