"""
sorting tool for BED's and VCF's that contain information on regions of the human genome
sorts by chromosome (1-22, X, Y, mt) then by start coord
assumes all files use utf-8 encoding
"""

import getopt
import sys
import os.path


DEFAULT_FILE_IN = './in.bed'
DEFAULT_FILE_OUT = './out.bed'


class SortableFile:
    def __init__(self, name):
        self.filename = name
        self.validate_file()
        self.type = os.path.splitext(self.filename)[1]

        self.header = []
        self.entries = []

        self.process()

    def __str__(self):
        return f"Input file was of type {self.type.upper()} and contained {len(self.entries)} entries."

    def process(self):
        with open(self.filename, 'r') as fin:
            lines = fin.readlines()

            # specific processing for bed files
            if self.type == ".bed":
                for line in lines:
                    if not line.startswith('chr'):
                        self.header.append(line)
                    elif line != '':
                        self.entries.append(GenericEntry(line))

            # specific processing for vcf files
            elif self.type == ".vcf":
                for line in lines:
                    if line.startswith('#'):
                        self.header.append(line)
                    elif line != '':
                        self.entries.append(GenericEntry(line))

            fin.close()

    def validate_file(self):
        # confirm that the input file exists
        if not os.path.isfile(self.filename):
            invalid_arg_file = self.filename
            self.filename = DEFAULT_FILE_IN

            # if the fallback file doesn't exist, exit
            if not os.path.isfile(self.filename):
                print(f'Neither arg file nor default file found. Try again.\nArg: {invalid_arg_file}\n'
                      f'Def: {self.filename}')
                exit(1)

    def sort(self):
        print(f'Attempting to sort input file: {self.filename}')
        self.entries.sort(key=lambda entry: (int(entry.chr), entry.start))

    def print_to_file(self, name):
        with open(name, 'w') as fout:
            fout.writelines(self.header)

            for entry in self.entries:
                fout.writelines(entry.original)
                if entry.original != self.entries[-1].original:
                    fout.writelines('\n')

            fout.close()

            print(f'Sorted file saved as: {name}')

    def check(self):
        print("Attempting to check if input file is sorted")
        comparison_list = sorted(self.entries, key=lambda entry: (int(entry.chr), entry.start))
        for i in range(len(self.entries)):
            if comparison_list[i].original != self.entries[i].original:
                return False
        return True


class GenericEntry:
    def __init__(self, line):
        self.original = line.strip()
        self.chr, self.start, self.id_or_end = line.split('\t', 2)

        if self.chr.startswith('chr'):
            self.chr = self.chr[3:]

        if self.chr.upper() == 'X':
            self.chr = '23'
        elif self.chr.upper() == 'Y':
            self.chr = '24'
        elif self.chr.upper() == 'MT' or self.chr.upper() == 'M':
            self.chr = '25'

    def __str__(self):
        return f'{self.chr}\t{self.start}\t{self.id_or_end}'


if __name__ == '__main__':
    file_in = DEFAULT_FILE_IN
    file_out = DEFAULT_FILE_OUT
    b_sort = False
    b_check = False

    # opts
    try:
        opts, args = getopt.getopt(sys.argv[1:], "f:o:sc", ["file=", "out=", "sort", "check"])

        for o, a in opts:
            if o in ("-f", "--file"):
                file_in = a

            elif o in ("-o", "--out"):
                file_out = a

            elif o in ("-s", "--sort"):
                b_sort = True

            elif o in ("-c", "--check"):
                b_check = True

            else:
                print(f"Unrecognized opt ({o}, {a})")

    except getopt.GetoptError as err:
        print(err)

    # try to process file
    f = SortableFile(file_in)

    if b_sort:
        f.sort()
        f.print_to_file(file_out)

    if b_check:
        if f.check():
            print(f'File is sorted')
        else:
            print (f'File is not sorted')
