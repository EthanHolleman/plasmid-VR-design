import os
import argparse
from pathlib import Path


# rename all files in a given directory by splitting filenames by a delimiter
# and then selecting one of the remaining regions by index (base 0)
# keeps the same file extension. 

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('DIR', help='Path to directory with files to truncate.')
    parser.add_argument('--delim', default='_', help='Delimiter to split files by.')
    parser.add_argument('--index', default=0, type=int, help='Index to keep after split')
    parser.add_argument('--constant_stem', default='', help='Set all file stems to this name but keep file extensions the same.')

    return parser.parse_args()



def main():
    args = get_args()
    for each_file in Path(args.DIR).iterdir():
        ext = each_file.suffix

        if args.constant_stem:
            new_stem = args.constant_stem
        else:
            new_stem = each_file.stem.split(args.delim)[args.index]
        new_name = new_stem + ext
        each_file.rename(each_file.with_name(new_name))


if __name__ == '__main__':
    main()


