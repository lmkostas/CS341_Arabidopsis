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
    "{gene} - {allele}",
    "_{gene}_-{allele}",
    "_{gene}_-{allele}",
    "_{gene}_-_{allele}_",
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

    genes = util.read_tsv_flat(GENE_LIST)

    genes_filtered = [gene.lower() for gene in genes if gene.lower() not in gene_blacklist]

    genes_filtered = [gene for gene in genes_filtered if len(gene) >= 1]

    genes_filtered.extend(enumerate_allele_extensions(genes_filtered))

    genes_filtered = [gene.lower() for gene in genes if gene.lower() not in gene_blacklist]

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

    #result.extend(load_pheno_ontology())

    # Filter by blacklist
    blacklist = read_blacklist()
    return [pheno.lower() for pheno in result if pheno.lower() not in blacklist and len(pheno)>1]

def load_pheno_ontology():
    """
    Load chebi, pato, and go ontologies.
    """
    ontology_terms = []
    for ontology_file in ONTOLOGIES:
        ontology_terms.extend(parse_ontology(ontology_file))

    blacklist = read_blacklist()
    return [pheno.lower() for pheno in ontology_terms if pheno.lower() not in blacklist and len(pheno)>1]

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
from snorkel.matchers import Contains, Sequence, DictionaryMatch, Concat, RegexMatchSpan, SlotFillMatch, Union, RegexMatchEach

genes = load_gene_list()
GENE = DictionaryMatch(d=genes, longest_match_only=True)
#GM = Sequence(GENE, longest_match_only=True)
#GENE = Concat(GENE, SlotFillMatch(GENE, pattern=r'{0}-\d+'), longest_match_only=True, right_required=False)
GENE = Union(GENE, SlotFillMatch(GENE, pattern=r'{0}-\d+'), RegexMatchEach(rgx=r'([A-Za-z]{1,4}\d+(\.\d+)?(-\d+)?){2,}'))
GM = Sequence(GENE, longest_match_only=True)

dict_linkwords = ['of', 'over', 'in', 'the', 'with', 'to', 'a']
#adjs = ['advanced', 'reduced', 'greater', 'less', 'small', 'large', 'short', 'tall', 'increased', 'decreased']
patos = parse_pato(PATO_ONTOLOGY) + ['increase', 'decrease', 'level', 'enhance', 'reduce', 'sensitive', 'resistant']
phenos = load_pheno_list() #+ ['root']
obos = load_pheno_ontology() + ['root', 'glucose'] + phenos

QUALIFIERS = RegexMatchEach(rgx=r'JJ.*|JJS.*|JJR.*', attrib='pos_tags')
#QUALIFIERS = RegexMatchSpan(rgx=r'JJ.*|JJS.*|JJR.*|VBN.*|RB.*|RBS.*|RBR.*|(JJ |JJS |JJR |VBN )?(NN|NNS|NNP|NNPS).*', attrib='pos_tags')
LINKS = DictionaryMatch(d=dict_linkwords, longest_match_only=True)
PATOS = DictionaryMatch(d=patos, attrib='lemmas', stemmer='porter', longest_match_only=True)
OBOS = DictionaryMatch(d=obos, attrib='lemmas', stemmer='porter', longest_match_only=True)
PHENOS = DictionaryMatch(d=phenos, attrib='lemmas', stemmer='porter', longest_match_only=True)

QUANT = RegexMatchSpan(rgx=r'(\d+(.\d+)?-)?\d+(.\d+)?%', longest_match_only=True)
ADJS = Sequence(RegexMatchEach(rgx=r'^(JJ|JJS|JJR|VBN)$', attrib='pos_tags', longest_match_only=True), longest_match_only=True)
PREPS = RegexMatchEach(rgx=r'(IN|TO).*', attrib='pos_tags')#DictionaryMatch(d=['of', 'in', 'to', 'over'])
DETS = DictionaryMatch(d=['a', 'the', 'its', 'their'], longest_match_only=True) 
NOUNS = RegexMatchEach(rgx=r'(NN|NNS)+ (IN|TO).*', attrib='pos_tags')
#ADJ_PHRASE = RegexMatchEach(rgx=r'(JJ|JJS|JJR|VBN)(( IN|TO)?( DT)?( NN| NNS)+( IN|TO)*)?.*', attrib='pos_tags')
#PM = Union (OBOS, PM, SlotFillMatch(ADJ_PHRASE, OBOS, pattern = '{0} {1}'))
VB = Union(RegexMatchEach(rgx=r'VBD|VBN.*', attrib='pos_tags', longest_match_only=True), DictionaryMatch(d=['is', 'are', 'was', 'were'], longest_match_only=True))
NN = Sequence(RegexMatchEach(rgx=r'(NN|NNS).*', attrib='pos_tags', longest_match_only=True), longest_match_only=True)
NNS = Sequence(RegexMatchSpan(rgx=r'[A-Za-z]+(ion|ity|ant|ent)', longest_match_only=True), longest_match_only=True)
ADJ_NN = Concat(ADJS, NN, longest_match_only=True)
VB_PHRASE = Concat(DictionaryMatch(d=['is', 'are', 'was', 'were'], longest_match_only=True), RegexMatchEach(rgx=r'VBD|VBN|JJ|JJR|JJS.*', attrib='pos_tags', longest_match_only=True), longest_match_only=True)
#PM = Sequence(OBOS, SlotFillMatch(OBOS, VB, ADJS, pattern='{0} {1} {2}'), SlotFillMatch(ADJS, OBOS, pattern='{0} {1}'), PATOS, PREPS, DETS, longest_match_only=True, required=[1, 0, 0, 0, 0, 0])

#ADJ_PHRASE = RegexMatchEach(rgx = r'(JJ|JJR|JJS) (NN|NNS)( IN)?.*', attrib='pos_tags')
PM = Sequence(OBOS, PATOS, VB, PREPS, DETS, longest_match_only=True, required = [1, 0, 0, 0, 0])
PM = Union(Contains(PM, rgxs=[r'.*(JJ|JJS|JJR).*',], attrib=['pos_tags']), Concat(ADJS, PM, longest_match_only=True, permutations=True), Concat(NNS, PM, longest_match_only=True, permutations=True), Concat(ADJ_NN, PM, longest_match_only=True), Concat(PM, VB_PHRASE,longest_match_only=True), longest_match_only=True)
'''
#PM = SlotFillMatch(PHENOS, DictionaryMatch(d=['is', 'are', 'was', 'were']), ADJ_PHRASE, pattern='{0} {1} {2}')
#OBO_QUAL = Sequence(QUALIFIERS, PATOS, OBOS, LINKS, longest_match_only=True, required=[1, 0, 1, 0])
#OBO_PATO = Sequence(QUALIFIERS, PATOS, OBOS, LINKS, longest_match_only=True, required = [1, 1, 1, 0])
#PM = Union(OBO_QUAL, PHENOS, longest_match_only=True)
#PM = RegexMatchEach(rgx=r'VBN.*', attrib='pos_tags')
#PM = RegexMatchSpan(rgx=r'JJ.*|JJS.*|JJR.*', attrib='pos_tags')
#PM = Sequence(DictionaryMatch(d=phenos, longest_match_only=True, attrib='lemmas'), DictionaryMatch(d=patos, attrib='lemmas', longest_match_only=True), LINKS, QUALIFIERS, longest_match_only=True, required=[1, 0, 0, 0])

PRE_ADJ = RegexMatchEach(rgx = r'(JJ|JJS|JJR|VBN) (NN|NNS) (IN|TO).*', attrib='pos_tags')
ADJ_SEQ = RegexMatchEach(rgx = r'((JJ|JJS|JJR|VBN), |(JJ|JJS|JJR|VBN))+and (JJ|JJS|JJR|VBN).*', attrib='pos_tags')
PRE_VBN = RegexMatchEach(rgx = r'VBN.*', attrib='pos_tags')
POST_ADJ = RegexMatchEach(rgx = r'(IN|TO) (JJ|JJS|JJR|VBN) (NN|NNS).*', attrib='pos_tags')
PRE_VBD = RegexMatchEach(rgx = r'(VBD|VBN) (TO|IN).*', attrib='pos_tags')
RegexMatchEach(rgx = r'(VBD|VBN) IN.*', attrib='pos_tags')
#PM = Concat(Sequence(RegexMatchEach(rgx = r'NN.*', attrib='pos_tags')), PM, longest_match_only=True, left_required=False)
'''
#Sequence(RegexMatchEach(rgx=r'(JJ|JJS|JJR).*', attrib='pos_tags', longest_match_only=True), longest_match_only=True)
#SlotFillMatch(RegexMatchEach(rgx=r'(JJ|JJS|JJR).*', attrib='pos_tags', longest_match_only=True), RegexMatchEach(rgx=r'(NN|NNS).*', attrib='pos_tags', longest_match_only=True), pattern='{0} {1}', longest_match_only=True)
#PM = Concat(ADJS, NN, longest_match_only=True)