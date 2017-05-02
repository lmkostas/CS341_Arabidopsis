import csv
import os
import util
from obo_parser import parseGOOBO
import re

BLACKLIST = "dicts/blacklist_words.txt"

GENE_LIST = "dicts/list_genes.txt"
ALLELE_PATTERNS = [
    "_ {gene} _ - {allele}",
    "_ {gene} _ -{allele}",
    "_ {gene} _ - _ {allele} _",
    "{gene} - {allele}",
    "{gene}-{allele}"
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

    genes = util.read_tsv_flat(GENE_LIST) + ['st5.1']

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

    full_result = []
    for p in result:
        p = re.sub(r'\([^)]*\)',' ',p)
        full_result.extend(p.split())

    result.extend(full_result)
    result = list(set(result))

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

    full_result = []
    for p in terms:
        p = re.sub(r'\([^)]*\)',' ',p)
        full_result.extend(p.split())

    terms.extend(full_result)
    terms = list(set(terms))
    blacklist = read_blacklist()
    return [pheno.lower() for pheno in terms if pheno.lower() not in blacklist]

#if __name__ == "__main__":
#    main()

from snorkel import SnorkelSession
session = SnorkelSession()
from snorkel.matchers import Sequence, DictionaryMatch, Concat, RegexMatchEach, RegexMatchSpan, SlotFillMatch, Union
genes = load_gene_list()
#dict_linkwords = ['of', 'over', 'in', 'the', 'with', 'to', 'a']
adjs = ['advanced', 'reduced', 'greater', 'less', 'small', 'large', 'short', 'tall', 'increased', 'decreased']
patos = parse_pato(PATO_ONTOLOGY)+adjs
phenos = load_pheno_list()
print phenos
GM = Union(Sequence(DictionaryMatch(d=genes, longest_match_only=True)), DictionaryMatch(d=genes, longest_match_only=True), RegexMatchEach(rgx=r'([A-Za-z]{1,4}\d+(.\d+)?(-\d+)?){2,}'))
PM = Concat(DictionaryMatch(d=phenos, stemmer='porter', longest_match_only=True), DictionaryMatch(d=patos, stemmer='porter', longest_match_only=True), permutations=True, longest_match_only=True)
#PM = Sequence(DictionaryMatch(d=phenos, stemmer='porter', longest_match_only=True), DictionaryMatch(d=patos, stemmer='porter', longest_match_only=True), DictionaryMatch(d=['a', 'the'], longest_match_only=True), longest_match_only=True, required=[1, 0, 0])
'''
PM=Concat(DictionaryMatch(d=phenos, stemmer='porter', longest_match_only=True), DictionaryMatch(d=['a', 'the']), longest_match_only=True, permutations=True)
PM_NN = Concat(RegexMatchEach(rgx=r'(\w+ion(, (and)?)?)+ (of|in)', longest_match_only=True), PM, longest_match_only=True)
PM_PART = Concat(RegexMatchEach(rgx=r'VBN+( CC VBN)?.*', longest_match_only=True, attrib='pos_tags'), PM, longest_match_only=True)
PM_ADJ = Concat(RegexMatchEach(rgx=r'(JJ|JJR)+( CC (JJ|JJR))?.*', longest_match_only=True, attrib='pos_tags'), PM, longest_match_only=True)
PM_PART_POST = SlotFillMatch(PM, DictionaryMatch(d=['was', 'is', 'were', 'are']), RegexMatchEach(rgx=r'(JJ|JJR|VBN)+( CC (JJ|JJR|VBN))?.*', longest_match_only=True, attrib='pos_tags'), pattern='{0} {1} {2}', longest_match_only=True)
PM_PATO = Concat(PM, DictionaryMatch(d=patos, stemmer='porter', longest_match_only=True), permutations=True, longest_match_only=True)

PM = Union(PM_PATO, PM_ADJ, PM_PART, PM_PART_POST, longest_match_only=True)
#PM = Sequence(PM, DictionaryMatch(d=patos, stemmer='porter', longest_match_only=True), RegexMatchEach(rgx=r'IN', longest_match_only=True, attrib='pos_tags'), required = [1, 1, 0])
'''