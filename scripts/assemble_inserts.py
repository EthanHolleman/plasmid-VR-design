import datetime
from pydna.readers import read
import pandas as pd
from pydna.genbankrecord import GenbankRecord

MOD_DATE = datetime.date.strftime(datetime.datetime.now(), "%d-%b-%Y").upper()



def read_variable_region_tsv(filepath):
    records = pd.read_table(filepath, sep='\t').to_dict(orient='records')
    assert len(records) == 1
    return records[0]


def variable_region_to_labeled_record(vr_dict, id, **kwargs):
    record = GenbankRecord(
        vr_dict['Sequence'].strip()
    )
    record.locus = vr_dict['name'].strip()  # check this keyword
    record.stamp()

    record.add_feature(
        0, len(record), name='variable_region',
        gc_skew=vr_dict['GC_skew'], gc_content=vr_dict['GC_content'],
        at_skew=vr_dict['AT_skew'], at_content=vr_dict['AT_content'],
        clustered_nucleotide=vr_dict['Clustered_nucleotide'],
        label='Variable region',
        date=MOD_DATE,
        type='CDS',
        **kwargs
    )
    return record


def assemble_insert(*args):
    insert = GenbankRecord(args[0])
    for each_record in args[1:]:
        insert += each_record
    return insert


def read_records(records):
    records = [GenbankRecord(read(r)) for r in records]
    return records


def set_insert_attributes(insert, vr_dict, version):
    # cleanup insert genbank attributes
    name = f'insert_{vr_dict["name"]}'
    insert.locus = name
    insert.id = insert.seq.seguid()

    insert.description = ''  # clear junk
    insert.definition = ''

    insert.definition = name
    insert.annotations['data_file_division'] = 'SYN'
    insert.annotations['date'] = MOD_DATE
    insert.annotations['sequence_version'] = version

    return insert


def main():

    vr_tsv = str(snakemake.input['variable_region'])
    upstream_sequences = snakemake.input['upstream_sequences']
    downstream_sequences = snakemake.input['downstream_sequences']

    upstream_records = read_records(upstream_sequences)
    downstream_records = read_records(downstream_sequences)

    version = str(snakemake.params['design_version'])

    output_file = str(snakemake.output)

    vr_dict = read_variable_region_tsv(vr_tsv)
    vr_record = variable_region_to_labeled_record(
        vr_dict, version,
        **snakemake.params['feature_annotations']
    )
    insert = assemble_insert(*upstream_records, vr_record, *downstream_records)
    insert = set_insert_attributes(insert, vr_dict, version)
    insert.write(output_file)


if __name__ == '__main__':
    main()
