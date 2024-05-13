import re


EXPORT_PATH = 'C:\\Users\\WDAGUtilityAccount\\Desktop\\custom.txt'
OUT_PATH = 'C:\\Users\\WDAGUtilityAccount\\Desktop\\out.txt'
OUT_HEADER = 'HGNC ID\tApproved symbol\tPredicted parent\tApproved name'

INDEX_HGNC_ID = None
INDEX_APPROVED_SYMBOL = None
INDEX_APPROVED_NAME = None
INDEX_STATUS = None
INDEX_LOCUS_TYPE = None

class HGNCGene:
    def __init__(self, raw):
        self.hgnc_id = raw[INDEX_HGNC_ID]
        self.approved_symbol = raw[INDEX_APPROVED_SYMBOL]
        self.approved_name = raw[INDEX_APPROVED_NAME]
        self.status = raw[INDEX_STATUS]
        self.locus_type = raw[INDEX_LOCUS_TYPE]
        self.predicted_parent_symbol = None
        self.in_panel = None

        self.process()

    def process(self):
        self.predict_parent()
        self.check_in_panel()

    def predict_parent(self):
        # symbol-based parent gene prediction
        symbol_match = re.search('(.+?)(P\\d*)$', self.approved_symbol)
        if symbol_match is not None:
            self.predicted_parent_symbol = symbol_match.group(1)

        # name-based parent gene prediction
        name_match = re.search('^(.+?),{0,1}\\s(.*)(pseudogene.*)$', self.approved_name)
        if name_match is not None:
            if name_match.group(1).upper() == self.predicted_parent_symbol:
                pass

            # exception for mitochondrial gene naming
            elif self.approved_symbol.startswith('MT') and self.predicted_parent_symbol is not None:
                if name_match.group(1).upper() == (self.predicted_parent_symbol[0:1] + '-' + self.predicted_parent_symbol[2:]):
                    self.predicted_parent_symbol = name_match.group(1).upper()

            # catchall for other cases
            # these are pseudogenes without predicted parents, might look at in the future
            else:
                if self.predicted_parent_symbol is None:
                    pass
                    # print(f'LOOKHERE\t{name_match.group(1)}\t{name_match.group(3)}\t{self.approved_symbol}\t{self.approved_name}')

    def check_in_panel(self):
        pass

    def __str__(self):
        if self.predicted_parent_symbol is not None:
            return str(f'{self.hgnc_id}\t{self.approved_symbol}\t{self.predicted_parent_symbol}\t{self.approved_name}')
        else:
            return str(f'{self.hgnc_id}\t{self.approved_symbol}\tUnable to predict parent\t{self.approved_name}')

    def clear(self):
        self.hgnc_id = None
        self.approved_symbol = None
        self.status = None
        self.locus_type = None
        self.predicted_parent_symbol = None
        self.in_panel = None


if __name__ == '__main__':
    file_in = open(EXPORT_PATH, 'r')
    lines = file_in.readlines()

    # identify indices of fields of interest
    fields = lines.pop(0).split('\t')
    INDEX_HGNC_ID = fields.index('HGNC ID')
    INDEX_APPROVED_SYMBOL = fields.index('Approved symbol')
    INDEX_APPROVED_NAME = fields.index('Approved name')
    INDEX_STATUS = fields.index('Status')
    INDEX_LOCUS_TYPE = fields.index('Locus type')

    pseudogenes = []

    for line in lines:
        gene = HGNCGene(line.split('\t'))
        if gene.locus_type == 'pseudogene':
            pseudogenes.append(gene)

    file_in.close()
    file_out = open(OUT_PATH, 'w')

    file_out.write(OUT_HEADER + '\n')
    for pseudogene in pseudogenes:
        print(pseudogene)
        file_out.write(str(pseudogene) + '\n')

    file_out.close()
