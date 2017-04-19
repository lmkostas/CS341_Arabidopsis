import csv
import os
import util
from obo_parser import parseGOOBO

BLACKLIST = "dicts/blacklist_words.txt"

GENE_LIST = "dicts/list_genes.txt"
ALLELE_PATTERNS = [
    "_ {gene} _ - {allele}",
    "_ {gene} _ -{allele}",
    "_ {gene} _ - _ {allele} _",
    "{gene} - {allele}"
]

PHENO_LIST = "dicts/list_phenotypes_arabidopsis_filtered.txt"
PHENO_EQ_LIST = "dicts/phenotypes_all_eq_dict.txt"
PHENO_MANUAL = "dicts/phenotypes_manual.txt"
ONTOLOGIES = ["dicts/po.obo", "dicts/chebi.obo", "dicts/go-basic.obo"]
PATO_ONTOLOGY = "dicts/pato.obo"

# Dictionary of linkwords to be added to the 'quality' dictionnary for the phenotype extraction.
LINKWORDS = ['of', 'over', 'in', 'the']

'''
def main():
    """
    Load full gene and phenotype list.
    """
    from snorkel import SnorkelSession
    session = SnorkelSession()
    from snorkel.matchers import DictionaryMatch
    genes = load_gene_list()
    GM = DictionaryMatch(d=genes)

    phenos = load_pheno_list()
    PM = DictionaryMatch(d=phenos, attrib='lemmas')
'''

def read_blacklist():
    """
    Read the blacklist words from disk.
    """
    with open(BLACKLIST) as blacklist:
        # Blacklist words are in the second column, starting on the second row.
        return [line.rstrip().split('\t')[1].lower() for line in blacklist][1:]


def load_gene_list():
    """
    Loads a list of genes from the genes dictionary, sans those present in the
    gene blacklist.
    """
    blacklist = read_blacklist()
    gene_blacklist = [item.lower() for item in blacklist]

    genes = util.read_tsv_flat(GENE_LIST)

    for g in genes:
        if g in gene_blacklist: print 'b', g
        if g.lower()=='also': print 'ALSO'

    genes_filtered = [gene.lower() for gene in genes if gene.lower() not in gene_blacklist]

    genes_filtered = [gene for gene in genes_filtered if len(gene) >= 1]

    genes_filtered.extend(enumerate_allele_extensions(genes_filtered))

    return genes_filtered


def enumerate_allele_extensions(genes):
    """
    Each gene could appear in a particular allele form. These are usually
    designated with a postfixed number. e.g. AGP1233 - 8
    """
    result = []
    allele_range = range(10)
    for allele_pattern in ALLELE_PATTERNS:
        genes_with_alleles = [allele_pattern.format(gene=gene, allele=allele)
                                 for gene in genes for allele in allele_range]
        result.extend(genes_with_alleles)
    return result


def load_pheno_list():
    """
    Loads a list of phenotypes from multiple files.
    """
    result = []
    result.extend(util.read_file_lines(PHENO_LIST))
    result.extend(util.read_tsv_flat(PHENO_EQ_LIST, delimiter=";"))
    result.extend(util.read_file_lines(PHENO_MANUAL))

    result.extend(load_pheno_ontology())

    # Filter by blacklist
    blacklist = read_blacklist()
    return [pheno.lower() for pheno in result if pheno.lower() not in blacklist]

def load_pheno_ontology():
    """
    Load chebi, pato, and go ontologies.
    """
    ontology_terms = []
    for ontology_file in ONTOLOGIES:
        ontology_terms.extend(parse_ontology(ontology_file))

    return ontology_terms


def parse_ontology(ontology_file):
    terms = []
    for elt in parseGOOBO(ontology_file):
        terms.append(elt["name"])
        if 'synonym' in elt:
            if isinstance(elt['synonym'], list):
                for syn in elt['synonym']:
                    try:
                        terms.append(syn.split('"')[1])
                    except:
                        print 'error parsing ontology synonym'
            else:
                try:
                    terms.append(elt['synonym'].split('"')[1])
                except:
                    print 'error parsing ontology synonym non list'
    return terms

def parse_pato(ontology_file):
    terms = []
    for elt in parseGOOBO(ontology_file):
        terms.append(elt["name"])
        if 'synonym' in elt:
            if isinstance(elt['synonym'], list):
                for syn in elt['synonym']:
                    try:
                        terms.append(syn.split('"')[1])
                    except:
                        print 'error parsing ontology synonym'
            else:
                try:
                    terms.append(elt['synonym'].split('"')[1])
                except:
                    print 'error parsing ontology synonym non list'
    return terms

#if __name__ == "__main__":
#    main()

from snorkel import SnorkelSession
session = SnorkelSession()
from snorkel.matchers import DictionaryMatch, Concat, RegexMatchSpan, SlotFillMatch
dict_linkwords = ['of', 'over', 'in', 'the', 'with']
patos = parse_pato(PATO_ONTOLOGY)+ dict_linkwords
'''
PATO_MATCH = Concat(DictionaryMatch(d=patos, attrib='lemmas', longest_match_only=True), DictionaryMatch(d=patos, attrib='lemmas', longest_match_only=True), longest_match_only=True, left_required=True, right_required=False)
END_LIST_MATCH = RegexMatchSpan(rgx=r'^[,;)/\-]')
START_LIST_MATCH = RegexMatchSpan(rgx=r'^[(]')
genes = load_gene_list()
#GENE_MATCH = Concat(DictionaryMatch(d=genes, longest_match_only=True), DictionaryMatch(d=genes, longest_match_only=True))
FULL_GENE = Concat(Concat(PATO_MATCH, DictionaryMatch(d=genes, longest_match_only=True), longest_match_only=True, left_required=False, right_required=True), PATO_MATCH, longest_match_only=True, left_required=True, right_required=False)
GENE_LIST = Concat(Concat(START_LIST_MATCH, FULL_GENE, longest_match_only=True, left_required=False), END_LIST_MATCH, longest_match_only=True, right_required=False)
GM = Concat(GENE_LIST, GENE_LIST, longest_match_only=True, right_required=False)
phenos = load_pheno_list()
FULL_PHENO = Concat(Concat(PATO_MATCH,DictionaryMatch(d=phenos, attrib='lemmas', longest_match_only=True), longest_match_only=True, left_required=False, right_required=True, permutations=False), PATO_MATCH, longest_match_only=True, left_required=True, right_required=False)
LIST_PHENO = Concat(Concat(START_LIST_MATCH, FULL_PHENO, longest_match_only=True, left_required=False), END_LIST_MATCH, longest_match_only=True, right_required=False)
PM = Concat(LIST_PHENO, LIST_PHENO, longest_match_only=True, right_required=False)
'''
load_gene_list()
PATO_MATCH = Concat(DictionaryMatch(d=patos, attrib='lemmas', longest_match_only=True), DictionaryMatch(d=patos, attrib='lemmas', longest_match_only=True), longest_match_only=True, left_required=False, right_required=False)
END_LIST_MATCH = RegexMatchSpan(rgx=r'^[,;/\-]')
#START_LIST_MATCH = RegexMatchSpan(rgx=r'^[(]')
genes = load_gene_list()
GM = DictionaryMatch(d=genes, longest_match_only=True)#Concat(DictionaryMatch(d=genes, longest_match_only=True), DictionaryMatch(d=genes, longest_match_only=True))
SlotFillMatch(GM, pattern = '{0} and {0}')
SlotFillMatch(GM, pattern = '({0})')
RegexMatchSpan(rgx=r'^[/\-][^\s]+')

#FULL_GENE_BEFORE = Concat(PATO_MATCH, DictionaryMatch(d=genes, longest_match_only=True), longest_match_only=True, left_required=False, right_required=True)
#FULL_GENE_AFTER = Concat(DictionaryMatch(d=genes, longest_match_only=True), PATO_MATCH, longest_match_only=True, left_required=True, right_required=False)
#FULL_GENE = Concat(FULL_GENE_BEFORE, FULL_GENE_AFTER, longest_match_only=True, left_required=False, right_required=False, permutations=True)
#GENE_LIST = Concat(Concat(START_LIST_MATCH, FULL_GENE, longest_match_only=True, left_required=False), END_LIST_MATCH, longest_match_only=True, right_required=False)
#GENE_LIST = Concat(FULL_GENE, END_LIST_MATCH, longest_match_only=True, right_required=False)
#FULL_GENE_LIST = Concat(GENE_LIST, GENE_LIST, longest_match_only=True, left_required=False, right_required=False)
#GM = Concat(FULL_GENE_LIST, SlotFillMatch(FULL_GENE_LIST, pattern = '({0})'), permutations=True, left_required=False, right_required=False)

phenos = load_pheno_list()
FULL_PHENO_BEFORE = Concat(PATO_MATCH, DictionaryMatch(d=phenos, longest_match_only=True), longest_match_only=True, left_required=False, right_required=True)
FULL_PHENO_AFTER = Concat(DictionaryMatch(d=phenos, longest_match_only=True), PATO_MATCH, longest_match_only=True, left_required=True, right_required=False)
FULL_PHENO = Concat(FULL_PHENO_BEFORE, FULL_PHENO_AFTER, longest_match_only=True, left_required=False, right_required=False, permutations=True)
LIST_PHENO = Concat(FULL_PHENO, END_LIST_MATCH, longest_match_only=True, right_required=False)
FULL_PHENO_LIST = Concat(LIST_PHENO, LIST_PHENO, longest_match_only=True, left_required=False, right_required=False)
PM = Concat(FULL_PHENO_LIST, SlotFillMatch(FULL_PHENO_LIST, pattern = '({0})'), permutations=True, left_required=False, right_required=False)
