from pathlib import Path
import datetime

def safe_dict_accesss(d, key):
    if key in d:
        return d[key]
    else:
        return None


class fastaToGenbank():
    
    def __init__(self, path, data={}):
        # record_kwargs are added to GenBankRecord object
        self.path = path
        self.data = data
        self.record = GenbankRecord(SeqIO.read(path))
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
    def defintion(self):  # overwrite
        definition = safe_dict_access(self.definition, 'definition')
        if definition:
            return definition
        else:
            return self.record.description
    
    def _update_record_with_parsed_data(self):
        self.record.__dict__.update(
            {
            'label': self.label,
            'locus': self.locus,
            'definition': self.definition
            }
        )

    
    def write_record(self, output_path=None):
        if not output_path:
            output_path = str(Path(output_path).with_suffix('.gb'))
        self.record.write(output_path)
    
    
    def add_label_feature(self, label=None):
        mod_date = datetime.date.strftime(datetime.datetime.now(), "%m/%d/%Y")
        if not label:
            label = self.label
        self.record.add_feature(
            x=0, 
            y=len(self.record), label=label)