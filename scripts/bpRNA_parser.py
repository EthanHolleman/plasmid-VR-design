# script for parsing output of the bpRNA program into easily digestable
# metrics for sequence comparison


class StructDescrip():

    structure_codes = {
        'E': '',
        'S': 'stem',
        'M': 'multiloop',
        'I': 'internal loop',
        'B': 'bulge',
        'H': 'hairpin loop',
    }

    @classmethod
    def init_from_line(cls, line):
        split_line = line.split(' ')  # space delim
        code = split_line[0][0]
        if code in StructDescrip.structure_codes:  # ignore if not in code dict
            number = split_line[0][1:]  # everything from beyond the code
            numbers = number.split('.')  # seperate primary and seconday id
            assert len(numbers) <= 2
            id = int(numbers[0])
            if len(numbers) == 2:
                sec_id = int(numbers[1])
            else:
                sec_id = None
            return cls(code, id, sec_id)
        else:
            return None

    def __init__(self, code, prim_id, sec_id=None):
        self.code = code
        self.prim_id = prim_id
        self.sec_id = sec_id
    

    def __repr__(self):
        return ' '.join([f'{key}: {val}' for key, val in self.__dict__.items()])

class StructureFile():


    def __init__(self, filepath):
        self.filepath = filepath
        self.dot_bracket = None
        self.struct_str = None
        self.seq = None
        self.header = None
        self.struct_descrips = []

        self._read_file()
        self._run_checks()
    
    def _struct_descrip_parser(self, line):
        items = line.split(' ')


    
    def _read_file(self):
        headers = []
        body_lines = 0
        with open(self.filepath) as handle:
            while handle:
                cur_line = handle.readline()
                if cur_line[0] == '#':
                    headers.append(cur_line[1:].strip())
                else:
                    if body_lines == 0:
                        # read the nucleotide sequence
                        self.seq = cur_line.strip()
                    elif body_lines == 1:
                        self.dot_bracket = cur_line.strip()
                    elif body_lines == 2:
                        self.struct_str = cur_line.strip()
                    else:
                        self.struct_descrips.append(cur_line.strip())
                    body_lines += 1
        
        self.headers = headers
    
    def _run_checks(self):
        assert len(self.seq) == len(self.dot_bracket) == len(self.struct_str)
    
    def __len__(self):
        return len(self.seq)

    @property
    def prop_paired(self):
        unpaired = self.dot_bracket.count('.') / len(self)
    
    @property
    def number_of_hairpins(self):
        count = 0
        for each_sd in self.struct_descrips:
            if each_sd.code == 'H':
                count += 1
            
        return count

