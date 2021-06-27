from pathlib import Path
import datetime
from pydna.genbankrecord import GenbankRecord
from pydna.readers import read
from Bio import SeqIO

def safe_dict_access(d, key):
    if key in d:
        return d[key]
    else:
        return None


class fastaToGenbank():
    '''Class to be used to convert fasta records to Genbank files and do basic
    operations like adding features and changing label names along the way
    if need be.
    '''
    
    def __init__(self, path, data={}):
        # record_kwargs are added to GenBankRecord object
        self.path = path
        self.data = data
        self.record = GenbankRecord(SeqIO.read(path, 'fasta'))
        self._update_record_with_parsed_data()
    
    @property
    def label(self):  # overwrite
        label = safe_dict_access(self.data, 'label')
        if label:
            return label
        else:
            return self.record.id
    
    @property
    def locus(self):  # overwrite
        locus = safe_dict_access(self.data, 'locus')
        if locus:
            return locus
        else:
            return self.record.id
    
    @property
    def definition(self):  # overwrite
        definition = safe_dict_access(self.data, 'definition')
        if definition:
            return definition
        else:
            return self.record.description
    
    def _update_record_with_parsed_data(self):
        self.record.label = self.label
        self.record.locus = self.locus
        self.record.definition = self.definition

    
    def write_record(self, output_path=None):
        if not output_path:
            output_path = str(Path(output_path).with_suffix('.gb'))
        self.record.write(output_path)
    
    
    def add_label_feature(self, label=None, **kwargs):
        mod_date = datetime.date.strftime(datetime.datetime.now(), "%m/%d/%Y")
        kwargs_dict = {'modification date': mod_date}
        kwargs_dict.update(kwargs)

        if not label:
            label = self.label
        self.record.add_feature(
            x=0,
            y=len(self.record), label=label,
            type='CDS',
            **kwargs_dict
            )
