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

PHENO_LIST = "dicts/list_phenotypes_arabidopsis_filtered.txt"
PHENO_EQ_LIST = "dicts/phenotypes_all_eq_dict.txt"
PHENO_MANUAL = "dicts/phenotypes_manual.txt"
ONTOLOGIES = ["dicts/po.obo", "dicts/chebi.obo", "dicts/go-basic.obo"]
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

from snorkel import SnorkelSession
session = SnorkelSession()
from extended_matchers import DictionaryMatch, Concat, RegexMatchSpan, SlotFillMatch, Union, RegexMatchEach, Sequence

blacklist = read_blacklist()
blacklist.extend([stem(b) for b in blacklist])
#phenos = load_pheno_list()
obos = load_pheno_ontology() + ['stem', 'leaves', 'phenotype', 'carpel', 'tip', 'type' ]
patos = parse_pato(PATO_ONTOLOGY) + ['alter', 'growth', 'develop', 'affect', 'display', 'twice', 'inhibit', 'type']


PATO = DictionaryMatch(d=patos, longest_match_only=True, stemmer='porter', blacklist=blacklist)
OBO = DictionaryMatch(d=obos, longest_match_only=True, stemmer='porter', blacklist=blacklist)

G = RegexMatchSpan(rgx=r'^[A-Z]{1,4}\d+$', ignore_case=False, longest_match_only=True)

NN_ADJ = RegexMatchSpan(rgx=r'^[A-Za-z-]+(ion|ment|ance|ence|ity)s?$', longest_match_only=True)
ADJS = DictionaryMatch(d=['JJR', 'RBR'], longest_match_only=True, attrib='pos_tags')
NUMS = RegexMatchSpan(rgx=r'^CD$', longest_match_only=True, attrib='pos_tags')
NNS = RegexMatchSpan(rgx=r'^(NN|NNS|NNP|NNPS)$', attrib='pos_tags', longest_match_only=False)
LINK_VBS = DictionaryMatch(d=['is', 'are', 'was', 'were', 'has', 'had', 'became', 'become', 'no', 'not'], longest_match_only=True, stemmer='porter')
LINKS = DictionaryMatch(d=['IN', 'DT', 'TO', 'CC', 'WDT', 'PRP', 'WRB'], longest_match_only=True, attrib='pos_tags')
ADVBS = DictionaryMatch(d=['RB', 'RBR', 'RBS', 'VBD', 'VBG', 'VBN', 'JJ', 'JJR', 'JJS'], longest_match_only=True, attrib='pos_tags')
NN = DictionaryMatch(d=['NN', 'NNS', 'NNP', 'NNPS'], longest_match_only=True, attrib='pos_tags')
COMPS = DictionaryMatch(d=['compared', 'contrast', 'relative', 'than', 'similar', 'different', 'same'])
NOS = DictionaryMatch(d=['no', 'not'], longest_match_only=True)
PREPS = DictionaryMatch(d=['IN', 'DT', 'TO', 'PRP', 'WRB'], longest_match_only=True, attrib='pos_tags')
PREP_PHRASE = RegexMatchEach(rgx=r'(IN|TO|DT|WDT|PRP)', attrib='pos_tags', longest_match_only=True)
DETS = DictionaryMatch(d='DT|WDT|PRP', longest_match_only=True, attrib='pos_tags')
MODS = DictionaryMatch(d=['RB', 'RBR', 'RBS', 'VBN', 'JJ', 'JJR', 'JJS'], longest_match_only=True, attrib='pos_tags')
NN_PHRASE = Sequence(NN, MODS, required=[1,0], longest_match_only=True)
OBO = Sequence(Union(OBO,G, longest_match_only=True), NNS, COMPS, NOS, NUMS, MODS, required=[1,0,0,0,0,0], punct='punct', links = [PREPS], longest_match_only=True)

PATO1 = Sequence(Union(PATO, ADJS, NN_ADJ, longest_match_only = True), NUMS, ADVBS, NOS, COMPS, required=[1,0,0,0,0], punct='punct', links = [LINKS, LINK_VBS], longest_match_only=True)
PATO2 = Sequence(DictionaryMatch(d=['JJ', 'VBN'], longest_match_only=True, attrib='pos_tags'), DictionaryMatch(d=['RB'], longest_match_only=True, attrib='pos_tags'), LINK_VBS, PATO, NN_ADJ, NUMS, ADVBS, NOS, COMPS, required=[1,1,1,0,0,0,0,0,0], punct='punct', links = [LINKS], longest_match_only=True)

PATO = Union(PATO1, PATO2)
OBO = Concat(OBO, PATO, right_required=False, longest_match_only=True, permutations=True)