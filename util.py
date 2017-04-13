"""
Utility functions useful in the Arabidopsis gene, phenotype, and relationship
extractor.
"""

import csv

def read_file_lines(filename):
    with open(filename) as file_handle:
        return [line for line in file_handle]

def read_tsv_flat(filename, delimiter="\t"):
    return_list = []
    with open(filename) as file_handle:
        for line in csv.reader(file_handle, delimiter=delimiter):
            for value in line:
                return_list.append(value)
    return return_list
