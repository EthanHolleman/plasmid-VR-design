# script for parsing output of the bpRNA program into easily digestable
# metrics for sequence comparison
import os

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
            
            # try to get the range of structure (length)
            struct_range = split_line[1].split('..')
            # if not in this format then ignore
            if len(struct_range) == 2:
                struct_range = int(struct_range[0]), int(struct_range[1])
            else:
                struct_range = None

            return cls(code, id, sec_id, struct_range)
        else:
            return None

    def __init__(self, code, prim_id, sec_id=None, struct_range=None):
        self.code = code
        self.prim_id = prim_id
        self.sec_id = sec_id
        self.struct_range = struct_range
    

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

    @property
    def prop_unpaired(self):
        return self.dot_bracket.count('.') / len(self)
    
    @property
    def num_hairpins(self):
        return len(list(self._get_sd_of_type('H')))
            
        return count
    
    @property
    def prop_hairpin(self):
        total_hairpin_span = 0
        for each_sd in self._get_sd_of_type('H'):
            if each_sd.struct_range:
                s, e = each_sd.struct_range
                total_hairpin_span += (e - s)
        return total_hairpin_span / len(self)


    def _read_file(self):
        headers = []
        body_lines = 0
        with open(self.filepath) as handle:
            lines = handle.readlines()
            for cur_line in lines:
                if cur_line:
                    if cur_line[0] == '#':  # header line
                        headers.append(cur_line[1:].strip())
                    else:  # body
                        if body_lines == 0:
                            self.header = headers
                            self.seq = cur_line.strip()
                        elif body_lines == 1:
                            self.dot_bracket = cur_line.strip()
                        elif body_lines == 2:
                            self.struct_str = cur_line.strip()
                        elif body_lines == 3:
                            pass  # some weird string maybe related to knots?
                        else:
                            descrip =  StructDescrip.init_from_line(cur_line.strip())
                            if descrip:
                                self.struct_descrips.append(descrip)

                        body_lines += 1

        
    def _run_checks(self):
        assert len(self.seq) == len(self.dot_bracket) == len(self.struct_str)
    
    def __len__(self):
        return len(self.seq)
    
    def _get_sd_of_type(self, code):
        # return all StructureDescrip instances with the code specified
        # by code argument
        for each_sd in self.struct_descrips:
            if each_sd.code == code:
                yield each_sd
    
    def to_tsv(self, output_path):
        # headers would be
        # filepath, prop_unpaired, length, num_hairpins, prop_hairpin
        with open(output_path, 'w') as handle:
            handle.write(
                '\t'.join([str(item) for item in 
                [self.filepath, self.prop_unpaired, len(self),
                self.num_hairpins, self.prop_hairpin]])
            )
        return output_path


def main():
    input_st = str(snakemake.input)
    output_tsv = str(snakemake.output)
    sf = StructureFile(input_st)
    sf.to_tsv(output_tsv)


if __name__ == '__main__':
    main()
        
