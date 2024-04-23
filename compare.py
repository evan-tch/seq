"""
vcf/vcf file comparison tool
will check if variants in a first vcf are present in a second vcf
"""

import getopt
import sys
import os.path
from pathlib import Path


DEFAULT_FILE_IN_VCF = './in.vcf'
DEFAULT_FILE_IN_VCF_REF = './ref.vcf'
DEFAULT_FILE_OUT = './out.vcf'


class VariantFile:
    def __init__(self, file_path, default):
        self.path = file_path
        self.header = []
        self.variants = []
        self.var_keys = []
        self.concordant_variants = []
        self.unique_to_input_variants = []
        self.unique_to_other_variants = []

        self.process(default)

    def process(self, default):
        self.path = validate_file(self.path, default)

        with open(self.path, 'r') as fin:
            lines = fin.readlines()

            for line in lines:
                if line.startswith('#'):
                    self.header.append(line)
                elif line != '':
                    self.variants.append(Variant(line))
                    self.var_keys.append(str(self.variants[-1]))

    def write_to_file(self, other):
        concordant_out = f'{Path(self.path).parent}/{Path(self.path).stem}_concordant{Path(self.path).suffix}'
        unique_to_input_out = f'{Path(self.path).parent}/{Path(self.path).stem}_unique{Path(self.path).suffix}'
        unique_to_other_out = f'{Path(self.path).parent}/{Path(other.path).stem}_unique{Path(self.path).suffix}'

        print(f"Concordant: {concordant_out}\n"
              f"Unique to input VCF: {unique_to_input_out}\n"
              f"Unique to reference VCF: {unique_to_other_out}")

        with open(concordant_out, 'w') as fout:
            fout.writelines(self.header)

            for variant in self.concordant_variants:
                fout.writelines(variant.original)
                if variant.original != self.concordant_variants[-1].original:
                    fout.writelines('\n')

            fout.close()

            print(f'Concordant results file saved as: {concordant_out}')

        with open(unique_to_input_out, 'w') as fout:
            fout.writelines(self.header)

            for variant in self.unique_to_input_variants:
                fout.writelines(variant.original)
                if variant.original != self.unique_to_input_variants[-1].original:
                    fout.writelines('\n')

            fout.close()

            print(f'Unique to input VCF file saved as: {unique_to_input_out}')

        with open(unique_to_other_out, 'w') as fout:
            fout.writelines(self.header)

            for variant in self.unique_to_other_variants:
                fout.writelines(variant.original)
                if variant.original != self.unique_to_other_variants[-1].original:
                    fout.writelines('\n')

            fout.close()

            print(f'Unique to referemce VCF file saved as: {unique_to_other_out}')

    def compare_with(self, other):
        # simple comparison
        for var_key in self.var_keys:
            if var_key in other.var_keys:
                self.concordant_variants.append(self.variants[self.var_keys.index(var_key)])
            else:
                self.unique_to_input_variants.append(self.variants[self.var_keys.index(var_key)])

        for var_key in other.var_keys:
            if var_key not in self.var_keys:
                self.unique_to_other_variants.append(other.variants[other.var_keys.index(var_key)])


class Variant:
    def __init__(self, line):
        self.original = line.strip()
        split_line = self.original.split('\t')
        self.chr = split_line[0]
        self.pos = split_line[1]
        self.ref = split_line[3]
        self.alt = split_line[4]
        split_line.clear()

        if self.chr.upper() == 'X':
            self.chr = '23'
        elif self.chr.upper() == 'Y':
            self.chr = '24'
        elif self.chr.upper() == 'MT' or self.chr.upper() == 'M':
            self.chr = '25'

    def __str__(self):
        return f'{self.chr}:{self.pos}{self.ref}_{self.alt}'


# return a default file if the input file doesn't exist; exits if the default doesn't exist
def validate_file(input_path, default_path):
    path = input_path

    # confirm that the file exists
    if not os.path.isfile(path):
        path = default_path

        # if the fallback file doesn't exist, exit
        if not os.path.isfile(path):
            print(f'Neither arg file nor default file found. Try again.\nArg: {input_path}\n'
                  f'Def: {path}')
            exit(1)

    return path


if __name__ == '__main__':
    file_in_vcf = DEFAULT_FILE_IN_VCF
    file_in_ref = DEFAULT_FILE_IN_VCF_REF
    file_out = DEFAULT_FILE_OUT

    # opts
    try:
        opts, args = getopt.getopt(sys.argv[1:], "v:r:s", ["variant_file=", "reference_file=", "sorted"])

        for o, a in opts:
            if o in ("-v", "--variant_file"):
                file_in_vcf = a

            elif o in ("-r", "--reference_file"):
                file_in_ref = a

            elif o in ("-s", "--sorted"):
                print('unimplemented')

            else:
                print(f"Unrecognized opt ({o}, {a})")

    except getopt.GetoptError as err:
        print(err)

    vcf = VariantFile(file_in_vcf, DEFAULT_FILE_IN_VCF)
    ref = VariantFile(file_in_ref, DEFAULT_FILE_IN_VCF_REF)

    vcf.compare_with(ref)
    vcf.write_to_file(ref)
