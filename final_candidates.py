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

    result.extend(load_pheno_ontology())

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

    #blacklist = read_blacklist()
    #return [pheno.lower() for pheno in ontology_terms if pheno.lower() not in blacklist and len(pheno)>1]

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
GENE_REGEX = RegexMatchEach(rgx=r'([A-Za-z0-9\/\.:-]*[A-Za-z]{2,4}\d+(\.\d+)?(-\d+(\.\d+)?)?[A-Za-z0-9\/:\.-]* ?)+', longest_match_only=True)
GENE_REGEX = RegexMatchEach(rgx=r'([A-Za-z0-9\/\.:-]*[A-Za-z]{2,4}\d+[A-Za-z0-9\/:\.-]*)+', longest_match_only=True)
GENE_SLOTFILL = SlotFillMatch(GENE, pattern=r'(([A-Za-z\.-]+(::|-|\/))*{0}((::|-|\/)[A-Za-z\.-]+)+)|(([A-Za-z\.-]+(::|-|\/))+{0}((::|-|\/)[A-Za-z\.-]+)*)')
PRE_GENE = RegexMatchEach(rgx=r'[A-Za-z0-9\.\/:-]+', longest_match_only=True)
GENE_SLOTFILL= SlotFillMatch(PRE_GENE, GENE, pattern=r'({0}(-|::|\/))?{1}((-|::|\/){0})?')
GM = Sequence(Union(GENE_REGEX, GENE_SLOTFILL, GENE, longest_match_only=True), longest_match_only=True)


phenos = load_pheno_list()
patos = parse_pato(PATO_ONTOLOGY)
NN_ADJ = RegexMatchSpan(rgx=r'([A-Za-z-]+(ion|ment|ance|ence|ity|ive) )+(of|to|in|on|over|against|ic)')
RegexMatchSpan(rgx=r'(IN|TO)', longest_match_only=True)
#PM = RegexMatchSpan(rgx=r'VBN?( NN| NNS| CC)+ IN( NN| NNS| CC)+', attrib='pos_tags', longest_match_only=True)
#PM = Concat(DictionaryMatch(d=['increased', 'decreased', 'reduced'], longest_match_only=True), Union(DictionaryMatch(d=patos, longest_match_only=True), DictionaryMatch(d=phenos, longest_match_only=True)), longest_match_only=True)
#PM = RegexMatchSpan(rgx=r'([A-Za-z-]+(ion|ment|ance|ence|ity|ive)( and)?)+ (of|to|in|on)', longest_match_only=True)
#RegexMatchSpan(rgx=r'VBN')
#PM = 

HELPER_VBS = DictionaryMatch(d=['was', 'is', 'are', 'were', 'become', 'becomes', 'has', 'had', 'have'], longest_match_only=True)
ADJS = RegexMatchSpan(rgx=r'((JJR|JJ|CC) ?)', attrib='pos_tags', longest_match_only=True)
NN = RegexMatchSpan(rgx=r'NN|NNS|NNP|NNPS', attrib='pos_tags', longest_match_only=True)
#PM = SlotFillMatch(RegexMatchSpan(rgx=r'((NN|NNS|CC) ?)+', attrib='pos_tags', longest_match_only=True), DictionaryMatch(d=['was', 'is', 'are', 'were', 'became', 'become'], longest_match_only=True), RegexMatchSpan(rgx=r'((JJ|JJR|CC) ?)+', attrib='pos_tags', longest_match_only=True), pattern='{0} {1} {2}', longest_match_only=True)
PM = SlotFillMatch(HELPER_VBS, ADJS, pattern='{0} {1}', longest_match_only=True)

#PM = SlotFillMatch(VB, DictionaryMatch(d=phenos, longest_match_only=True), pattern = '{0} {1}')
NN_PHRASE = RegexMatchSpan(rgx=r'(VB |VBD |VBZ |VBP |VBG |VBN )?(DT )?(NN ?| NNS ?| NNP ?| NNPS ?|JJ ?)+', attrib='pos_tags', longest_match_only=True)
PHENOS = DictionaryMatch(d=['phenotype', 'phenotypes'], longest_match_only=True)
PHENO_PHRASE = SlotFillMatch(NN_PHRASE, PHENOS, pattern='{0} {1}')
PM_ADJS = RegexMatchSpan(rgx=r'(RBR JJ|JJR).*(NN|NNS|NNP|NNPS)', attrib='pos_tags', longest_match_only=True)
PM_VBS = RegexMatchSpan(rgx=r'(VBD|VBN)( DT)?( NN|NNS|NNP|NNPS)+ (IN|TO)( JJ| JJR| NN| NNS| NNP| NNPS)*', attrib='pos_tags', longest_match_only=True)
PM = Union(PM_ADJS, PM_VBS, PHENO_PHRASE, longest_match_only=True)