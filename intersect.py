"""
vcf/bed file intersection tool
will attempt to filter a file containing variants (VCF) by coverage regions defined in a second file (BED)
offers two intersection methods; if both files are sorted (run through sort.py), the more efficient method can be used
by providing the --sorted flag
assumes that regions in BED don't overlap (run through merge.py)
the region file must contain at least 3 tab-separated columns: chrom, chromStart, chromEnd
the variant file must contain at least 2 tab-separated columns: chrom, pos
assumes all files use utf-8 encoding
"""

import getopt
import sys
import os.path


DEFAULT_FILE_IN_VCF = './in.vcf'
DEFAULT_FILE_IN_BED = './in.bed'
DEFAULT_FILE_OUT = './out.vcf'


class VariantFile:
    def __init__(self, file_path):
        self.path = file_path
        self.header = []
        self.variants = []

        self.process()

    def process(self):
        self.path = validate_file(self.path, DEFAULT_FILE_IN_VCF)

        with open(self.path, 'r') as fin:
            lines = fin.readlines()

            for line in lines:
                if line.startswith('#'):
                    self.header.append(line)
                elif line != '':
                    self.variants.append(Variant(line))

    def intersect_with(self, region_file, intersect_type='simple'):
        if intersect_type == 'simple':
            # simple, time consuming intersection
            for variant in reversed(self.variants):
                b_intersect = False
                for region in region_file.regions:
                    if int(variant.chr) == int(region.chrom) and \
                       int(region.chromStart) <= int(variant.pos) <= int(region.chromEnd):
                        b_intersect = True
                        break

                # remove the variant if it doesn't fall in any region defined in the bed file
                if not b_intersect:
                    self.variants.remove(variant)
        elif intersect_type == 'complex':
            # more complex, faster intersection
            # assumes both files are sorted
            iarr = []
            current_region = 0
            for variant in self.variants:

                for region in region_file.regions[current_region:]:
                    if int(variant.chr) == int(region.chrom):
                        if int(region.chromStart) <= int(variant.pos) <= int(region.chromEnd):
                            break
                        elif int(variant.pos) < int(region.chromStart):
                            iarr.append(self.variants.index(variant))
                            break
                        else:
                            current_region += 1
                    else:
                        current_region += 1

            # remove by index array
            for i in reversed(iarr):
                self.variants.remove(self.variants[i])

    def write_to_file(self, output_path):
        with open(output_path, 'w') as fout:
            fout.writelines(self.header)

            for variant in self.variants:
                fout.writelines(variant.original)
                if variant.original != self.variants[-1].original:
                    fout.writelines('\n')

            fout.close()

            print(f'Intersect file saved as: {output_path}')


class RegionFile:
    def __init__(self, file_path):
        self.path = file_path
        self.header = []
        self.regions = []

        self.process()

    def process(self):
        self.path = validate_file(self.path, DEFAULT_FILE_IN_BED)

        with open(self.path, 'r') as fin:
            lines = fin.readlines()

            for line in lines:
                if not line.startswith('chr'):
                    self.header.append(line)
                elif line != '':
                    self.regions.append(Region(line))

            fin.close()


class Region:
    def __init__(self, line):
        split_line = (line.strip()).split('\t')
        self.chrom = split_line[0]
        self.chromStart = split_line[1]
        self.chromEnd = split_line[2]
        split_line.clear()

        self.chrom = self.chrom[3:]

        if self.chrom.upper() == 'X':
            self.chrom = '23'
        elif self.chrom.upper() == 'Y':
            self.chrom = '24'
        elif self.chrom.upper() == 'MT' or self.chrom.upper() == 'M':
            self.chrom = '25'

    def __str__(self):
        return f'{self.chrom}\t{self.chromStart}\t{self.chromEnd}'


class Variant:
    def __init__(self, line):
        self.original = line.strip()
        split_line = self.original.split('\t')
        self.chr = split_line[0]
        self.pos = split_line[1]
        split_line.clear()

        if self.chr.upper() == 'X':
            self.chr = '23'
        elif self.chr.upper() == 'Y':
            self.chr = '24'
        elif self.chr.upper() == 'MT' or self.chr.upper() == 'M':
            self.chr = '25'

    def __str__(self):
        return f'{self.chr}\t{self.pos}'


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
    file_in_vcf = ''
    file_in_bed = ''
    file_out = DEFAULT_FILE_OUT
    intersect_type = 'simple'

    # opts
    try:
        opts, args = getopt.getopt(sys.argv[1:], "v:r:o:s", ["variant_file=", "region_file=", "out=", "sorted"])

        for o, a in opts:
            if o in ("-v", "--variant_file"):
                file_in_vcf = a

            elif o in ("-r", "--region_file"):
                file_in_bed = a

            elif o in ("-o", "--out"):
                file_out = a

            elif o in ("-s", "--sorted"):
                intersect_type = 'complex'

            else:
                print(f"Unrecognized opt ({o}, {a})")

    except getopt.GetoptError as err:
        print(err)

    # try to process bed file
    bed = RegionFile(file_in_bed)

    # try to process vcf file
    vcf = VariantFile(file_in_vcf)
    vcf.intersect_with(bed, intersect_type)
    vcf.write_to_file(file_out)