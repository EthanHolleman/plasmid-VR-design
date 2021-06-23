from dot_bracket import Parser


def read_rnaFold_output(filepath):
    structures = {}
    with open(filepath) as handle:
        while handle:
            current_line = handle.readline()
            if current_line[0] == '>': # start of record
                current_header = current_line
                current_seq = handle.readline()
                current_struct = handle.readline()
                structures[current_header.strip()] = Parser(current_struct.strip())

    return structures


def write_structures_as_tsv(structures):
    rows = []
    for header, struct in structures.items():
        row_dict = {}
        row_dict['header'] = header
        struct.parse
        for struct_type in struct.


