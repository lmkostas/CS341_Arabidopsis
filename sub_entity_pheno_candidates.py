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

PHENO_LIST = "dicts/list_phenotypes_arabidopsis_filtered.txt"  #Leave out for now 9/28/17
PHENO_EQ_LIST = "dicts/phenotypes_all_eq_dict.txt"
PHENO_MANUAL = "dicts/phenotypes_manual.txt"
ONTOLOGIES = ["dicts/po.obo", "dicts/chebi.obo", "dicts/go-basic.obo", "dicts/pato.obo"]
PATO_ONTOLOGY = "dicts/pato.obo"


def read_blacklist():
    """
    Read the blacklist words from disk.
    """
    with open(BLACKLIST) as blacklist:
        # Blacklist words are in the second column, starting on the second row.
        return [line.rstrip().split('\t')[1].lower() for line in blacklist][1:]

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

    # Filter by blacklist
    blacklist = read_blacklist()
    #blacklist = [stem_word(b) for b in blacklist]

    return [pheno.lower() for pheno in result if pheno.lower() not in blacklist and len(pheno)>1]

def load_pheno_ontology(O):
    """
    Load chebi, pato, and go ontologies.
    """
    ontology_terms = []
    for ontology_file in O:
        if ontology_file == "dicts/po.obo" or ontology_file == "dicts/pato.obo":
            ontology_terms.extend(parse_ontology(ontology_file, split=True))
        else:
            ontology_terms.extend(parse_ontology(ontology_file))

    #blacklist = read_blacklist()
    #return [pheno.lower() for pheno in ontology_terms if pheno.lower() not in blacklist and len(pheno)>1]

    return ontology_terms


def parse_ontology(ontology_file, split=False):
    terms = []
    for elt in parseGOOBO(ontology_file):
        terms.append(elt["name"])
        if 'synonym' in elt:
            if isinstance(elt['synonym'], list):
                for syn in elt['synonym']:
                    try:
                        term = syn.split('"')[1]
                        terms.append(term)
                        if ' ' in term: terms.extend(term.split())
                    except:
                        print 'error parsing ontology synonym'
            else:
                try:
                    term = elt['synonym'].split('"')[1]
                    terms.append(term)
                    if ' ' in term: terms.extend(ter.split())
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

from snorkel import SnorkelSession
session = SnorkelSession()
from extended_matchers import DictionaryMatch, Concat, RegexMatchSpan, SlotFillMatch, Union, RegexMatchEach, Intersection, SetDiff

blacklist = read_blacklist()
blacklist.extend([stem(b) for b in blacklist])
#phenos = load_pheno_list()
entity = load_pheno_ontology(ONTOLOGIES) + ['stem', 'leaves', 'phenotype', 'carpel', 'tip', 'type' ] #formerly named obos
descriptor = parse_pato(PATO_ONTOLOGY) + ['alter', 'growth', 'develop', 'affect', 'display', 'twice', 'inhibit', 'type'] #formerly named patos


DESC = DictionaryMatch(d=patos, longest_match_only=True, blacklist=blacklist)
ALL_ADJ_LIKE = DictionaryMatch(d=['JJ','JJR','JJS','RB','RBR','RBS'], longest_match_only=True, attrib='pos_tags')

DESC_ENT = SetDiff(DESC, ALL_ADJ_LIKE)
ENT = DictionaryMatch(d=obos, longest_match_only=True, blacklist=blacklist)
NN_ADJ = RegexMatchSpan(rgx=r'^[A-Za-z-]+(ion|ment|ance|ence|ity)s?$', longest_match_only=True)
ALL_ENTS = Union(Union(ENT, DESC_ENT, longest_match_only=True), NN_ADJ, longest_match_only=True)

DESC = Intersection(DESC, ALL_ADJ_LIKE)
ALL_DESCS = Union(DESC,DictionaryMatch(d=['JJ','JJR','RB','RBR','VB', 'VBD'], longest_match_only=True, attrib='pos_tags'), longest_match_only=True)

PERC = RegexMatchSpan(rgx=r'\d+.?\d*%', longest_match_only=True)

PREP = DictionaryMatch(d=['IN', 'TO'], longest_match_only=True, attrib='pos_tags')

ALL_DESCS = Union(Concat(ALL_DESCS, PERC, right_required=False, permutations=True), SlotFillMatch(ALL_DESCS, PREP, PERC, '{0} {1} {2}', longest_match_only=True))

# COMP_ADJS = DictionaryMatch(d=['JJR'], longest_match_only=True, attrib='pos_tags')
# ADJS = DictionaryMatch(d=['JJ'], longest_match_only=True, attrib='pos_tags')
# PAST_PART = DictionaryMatch(d=['VBN'], longest_match_only=True, attrib='pos_tags')

# ADVBS = DictionaryMatch(d=['RB'], longest_match_only=True, attrib='pos_tags')
# COMP_ADVBS = DictionaryMatch(d=['RBR'], longest_match_only=True, attrib='pos_tags')

# VRB = DictionaryMatch(d=['VB', 'VBD'], longest_match_only=True, attrib='pos_tags')
#TODO Concat(VB, Union(ADVBS,COMP_ADVBS))