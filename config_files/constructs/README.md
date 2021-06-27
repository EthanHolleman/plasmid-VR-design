# Construct definitions

Variable regions that are designed via the workflow need to be incorporated
into a plasmid construct in order to be useful. In order to define how and
where variable regions should be inserted into existing plasmid backbones
a specific `yaml` syntax is used and is described below.

## Yaml construct definition

Each construct is described as a `yaml` type dictionary. A construct refers
to a plasmid backbone with all inserted regions. 

An example definition for a construct called `construct_a` is shown below.

```
construct_a:
    backbone: "backbone.gb" 
    insert_downstream_of: "promotor_a"
    contents:
        - "VAR_REGION"  
        - "extension_region.gb"  
```

### Fields

Fields are case sensitive and should be entered exactly as they appear in the
example above.

#### `backbone`

Path to a genbank formatted file that described the plasmid backbone inserts will be inserted into via Gibson assembly.

#### `insert_downstream_of`

Label of a feature within the `backbone` genbank record that inserts should
be inserted downstream of (with respect to the sense of the feature) via
Gibson assembly. This was written with the intent of this feature being a
promotor. The workflow will identify this feature from the `backbone` record
and select the closest downstream restriction site to linearize the backbone
and insert the fragments. 

#### `contents`

This is an ordered list of the fragments to be included in the insert. The
order should be the intended order of the fragments in the 5' -> 3' direction.
Items within the `contents` list can either be keywords or genbank records
which are further described below.

**Note**: The sequence that is actually inserted will be the reverse complement
of what is passed to the workflow if the `insert_downstream_of` feature is in
the antisense orientation. This is because at the time of writing it was
assumed that the inserted region is to be transcribed and in order to produce
an mRNA strand equivalent to the sequence passed to the workflow (assumed to
be in the 5' -> 3' orientation) the reverse complement of that sequence must
be inserted if it is to be transcribed from the opposite direction. 

#### Keywords

- `VAR_REGION`: Specifying this keyword indicates to the workflow that a
variable region should be inserted at the specified position in the `contents`
list. 
- Path to genbank record: Fragment to be inserted at this position is constant
and is specified within the record located at this filepath.






