import csv
import os
import util
import re
from obo_parser import parseGOOBO
if 'CI' not in os.environ:
    try:
        from nltk.stem.porter import PorterStemmer
    except ImportError:
        warnings.warn("nltk not installed- some default functionality may be absent.")

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
    split = []
    for r in result:
        split.extend(r.split(r'[\.,;]'))
    result.extend(split)
    result.extend(util.read_tsv_flat(PHENO_EQ_LIST, delimiter=";"))
    result.extend(util.read_file_lines(PHENO_MANUAL))

    #result.extend(load_pheno_ontology())

    # Filter by blacklist
    blacklist = read_blacklist()
    #blacklist = [stem_word(b) for b in blacklist]

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
    blacklist = read_blacklist()
    #blacklist = [stem_word(b) for b in blacklist]

    return [pheno.lower() for pheno in terms if pheno.lower() not in blacklist and len(pheno)>1]

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

    blacklist = read_blacklist()
    blacklist.extend([stem(b) for b in blacklist])
    #blacklist = [stem_word(b) for b in blacklist]
    return [pheno.lower() for pheno in terms if pheno.lower() not in blacklist and len(pheno)>1]

    return terms

def stem(w):
        """Apply stemmer, handling encoding errors"""
        stemmer = PorterStemmer()
        try:
            return stemmer.stem(w)
        except UnicodeDecodeError:
            return w

#if __name__ == "__main__":
#    main()
from snorkel import SnorkelSession
session = SnorkelSession()
from snorkel.matchers import Seq, Contains, Sequence, DictionaryMatch, Concat, RegexMatchSpan, SlotFillMatch, Union, RegexMatchEach

genes = load_gene_list()
GENE = DictionaryMatch(d=genes, longest_match_only=True)
GENE_REGEX = RegexMatchSpan(rgx=r'((\S+)?[A-Za-z]{2,4}\d(\S+)?)', longest_match_only=True)
PRE_GENE = RegexMatchEach(rgx=r'[A-Za-z0-9\.\/:-]+', longest_match_only=True)
GENE_SLOTFILL= SlotFillMatch(PRE_GENE, GENE, pattern=r'({0}(-|::|\/))?{1}((-|::|\/){0})?')
GM = Sequence(Union(GENE_REGEX, GENE_SLOTFILL, GENE, longest_match_only=True), longest_match_only=True)

blacklist = read_blacklist()
blacklist.extend([stem(b) for b in blacklist])
phenos = load_pheno_list()
obos = load_pheno_ontology() + ['stem', 'leaves', 'phenotype', 'carpel']
full_phenos = load_pheno_list()
patos = parse_pato(PATO_ONTOLOGY) + ['alter']
#obos.extend(patos)

PATO = DictionaryMatch(d=patos, longest_match_only=False, stemmer='porter', blacklist=blacklist)
OBO = DictionaryMatch(d=obos, longest_match_only=False, stemmer='porter', blacklist=blacklist)
PHENOS = DictionaryMatch(d=phenos, longest_match_only=True)
NN_ADJ = RegexMatchSpan(rgx=r'^[A-Za-z-]+(ion|ment|ance|ence|ity|ive)$', longest_match_only=True)
ADJS = DictionaryMatch(d=['JJ', 'JJR', 'VBN', 'RB', 'RBR'], longest_match_only=True, attrib='pos_tags')
LINK_VBS = DictionaryMatch(d=['is', 'are', 'was', 'were', 'has', 'had', 'became', 'become', 'no', 'not'], longest_match_only=True, stemmer='porter')
LINKS = RegexMatchEach(rgx=r'(IN|TO|DT|CC)', longest_match_only=True, attrib='pos_tags')
NUMS = SlotFillMatch(RegexMatchEach(rgx=r'CD', longest_match_only=True, attrib='pos_tags'), pattern= r'{0}(%| percent|-fold)?', longest_match_only=True)
PM = Sequence(PATO, OBO, PHENOS, NN_ADJ, ADJS, GENE, GENE_REGEX, GENE_SLOTFILL, links = [LINKS, LINK_VBS, NUMS], required = [1, 1, 0, 0, 0, 0, 0, 0], longest_match_only=True)
'''
LINKS = RegexMatchEach(rgx=r'(IN|TO|DT|CC)', longest_match_only=True, attrib='pos_tags')
NUMS = SlotFillMatch(RegexMatchEach(rgx=r'CD', longest_match_only=True, attrib='pos_tags'), pattern= r'{0}(%| percent|-fold)?', longest_match_only=True)
ADJS = DictionaryMatch(d=['JJ', 'JJR', 'VBN'], longest_match_only=True, attrib='pos_tags')
VBS = DictionaryMatch(d=['VB', 'VBZ', 'VBD', 'RB', 'RBR'], longest_match_only=True, attrib='pos_tags')
VBS = RegexMatchEach(rgx=r'(VB|VBZ|VBD|RB|RBR)', longest_match_only=True, attrib='pos_tags')
NNS = DictionaryMatch(d=['NN', 'NNS', 'NNP', 'NNPS'], longest_match_only=True, attrib='pos_tags')
LINK_VBS = ['is', 'are', 'was', 'were', 'has', 'had', 'became', 'become']

#adj + phrase
PRE_NN = RegexMatchSpan(rgx=r'^(JJ |JJR |VBN )?(NN|NNS|NNP|NNPS)(( IN| TO)( IN| TO| DT)+)?$', longest_match_only=True, attrib='pos_tags')
POST_NN = RegexMatchSpan(rgx=r'^((IN |TO )(IN |TO |DT )+ )?(JJ |JJR |VBN )*(NN|NNS|NNP|NNPS)$', longest_match_only=True, attrib='pos_tags')


#PM = Union(PATO, OBO, PHENOS, ADJS)
#PM = Seq(d=patos + obos + phenos, longest_match_only=True)
PM1 = Union(SlotFillMatch(PATO, LINKS, pattern = '{0} {1} {0}', longest_match_only=True), SlotFillMatch(PATO, LINKS, pattern = '{1} {0} {1}', longest_match_only=True), SlotFillMatch(PATO, LINKS, pattern = '{1} {0}', longest_match_only=True), SlotFillMatch(PATO, LINKS, pattern = '{0} {1}', longest_match_only=True), PATO, longest_match_only=True)
PM2 = Concat(PM1, PM1, longest_match_only=True, right_required=False)
PM3 = Union(SlotFillMatch(OBO, pattern = r'([a-zA-Z]+-)?{0}(-[a-zA-Z]+)?', longest_match_only=True), SlotFillMatch(OBO, pattern = r'{0},? (and|or) {0}', longest_match_only=True), OBO, longest_match_only=True) 
PM4 = Union(SlotFillMatch(PM3, LINKS, pattern = '{0} {1} {0}', longest_match_only=True), SlotFillMatch(PM3, pattern = '{0} {0}', longest_match_only=True), PM3, longest_match_only=True)#Concat(PM3, PM3, longest_match_only=True, right_required=False)
PM = Union(PM2, PM4, longest_match_only=True)
#PM = Concat(PM2, PM4, longest_match_only=True, permutations=True)
#PM = Union(SlotFillMatch(PM2, PM4, pattern = '{1} {0} {1}', longest_match_only=True), SlotFillMatch(PM2, PM4, pattern = '{0} {1} {0}', longest_match_only=True), SlotFillMatch(PM2, PM4, pattern = '{1} {0}', longest_match_only=True), SlotFillMatch(PM2, PM4, pattern = '{0} {1}', longest_match_only=True), longest_match_only=True) 
#PM = Union(Concat(ADJS, PM, longest_match_only=True), PM, longest_match_only=True)
#PM = Union(SlotFillMatch(PRE_NN, POST_NN, PM, pattern = '{0} {2} {1}', longest_match_only=True), SlotFillMatch(PRE_NN, PM, pattern = '{0} {1}', longest_match_only=True), SlotFillMatch(POST_NN, PM, pattern = '{1} {0}', longest_match_only=True), PM, longest_match_only=True)
#PM = Union(PM1, OBO, longest_match_only=True) 
'''
'''
PM1 = Union(SlotFillMatch(PATO, LINKS, pattern = '{0} {1} {0}', longest_match_only=True), SlotFillMatch(PATO, pattern = '{0} {0}', longest_match_only=True), PATO, longest_match_only=True)
PM2 = Union(SlotFillMatch(OBO, pattern = r'([a-zA-Z]+-)?{0}(-[a-zA-Z]+)?', longest_match_only=True), SlotFillMatch(OBO, pattern = r'{0},? (and|or) {0}', longest_match_only=True), OBO, longest_match_only=True) 
PM2 = Union(SlotFillMatch(PM2, LINKS, pattern = '{0} {1} {0}', longest_match_only=True), SlotFillMatch(PM2, pattern = '{0} {0}', longest_match_only=True), PM2, longest_match_only=True)
#PM = Concat(PM2, PM4, longest_match_only=True, permutations=True)
PM = Union(SlotFillMatch(PM1, PM2, LINKS, pattern = '{1} {2} {0} {2} {1}', longest_match_only=True), SlotFillMatch(PM1, PM2, LINKS, pattern = '{1} {2} {0} {1}', longest_match_only=True), SlotFillMatch(PM1, PM2, LINKS, pattern = '{1} {0} {2} {1}', longest_match_only=True), SlotFillMatch(PM1, PM2, pattern = '{1} {0} {1}', longest_match_only=True), SlotFillMatch(PM1, PM2, LINKS, pattern = '{0} {2} {1} {2} {0}', longest_match_only=True), SlotFillMatch(PM1, PM2, LINKS, pattern = '{0} {2} {1} {0}', longest_match_only=True), SlotFillMatch(PM1, PM2, LINKS, pattern = '{0} {1} {2} {0}', longest_match_only=True), SlotFillMatch(PM1, PM2, pattern = '{0} {1} {0}', longest_match_only=True), SlotFillMatch(PM1, PM2, LINKS, pattern = '{1} {2} {0}', longest_match_only=True), SlotFillMatch(PM1, PM2, pattern = '{1} {0}', longest_match_only=True), SlotFillMatch(PM1, PM2, LINKS, pattern = '{0} {2} {1}', longest_match_only=True), SlotFillMatch(PM1, PM2, pattern = '{0} {1}', longest_match_only=True), longest_match_only=True) 
'''