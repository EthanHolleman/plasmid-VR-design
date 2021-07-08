import datetime
from pydna.readers import read
from pandas import pd


def read_variable_region_tsv(filepath):
    return pd.read_table(filepath, sep='\t', header=True).to_dict(orient='record')


def variable_region_to_labeled_record(vr_dict, id, **kwargs):
    mod_date = datetime.date.strftime(datetime.datetime.now(), "%m/%d/%Y")
    record = Genbankrecord(
        vr_dict['Sequence'].strip()
    )
    record.locus = vr_dict['name'].strip()  # check this keyword
    record.description = vr_dict['role']
    record.stamp()

    record.add_feature(
        0, len(record), name='variable_region',
        gc_skew=vr_dict['GC_skew'], gc_content=vr_dict['gc_content'],
        at_skew=vr_dict['AT_skew'], at_content=vr_dict['AT_content'],
        clustered_nucleotide=vr_dict['Clustered_nucleotide'],
        date=mod_date,
        type='CDS',
        **kwargs
    )
    return record


def assemble_insert(*args):
    insert = Genbankrecord(args[0])
    for each_record in args[1:]:
        insert += each_record
    return insert


def read_records(records):
    return [read(r) for r in records]


def main():

    vr_tsv = str(snakemake.input['variable_region'])
    upstream_sequences = snakemake.input['upstream_sequences']
    downstream_sequences = snakemake.input['downstream_sequences']

    upstream_records = read_records(upstream_sequences)
    downstream_records = read_records(downstream_records)

    output_file = str(snakemake.output)

    vr_dict = read_variable_region_tsv(vr_tsv)
    vr_record = variable_region_to_labeled_record(vr_dict)
    insert = assemble_insert(*upstream_records, vr_record, *downstream_records)
    insert.locus = f'insert_{vr_dict["name"]}'
    insert.stamp()
    insert.write(output_file)


if __name__ == '__main__':
    main()
