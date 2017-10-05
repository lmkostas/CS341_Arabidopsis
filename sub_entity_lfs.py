import re
import os
from snorkel.lf_helpers import (
    get_left_tokens,
    get_between_tokens,
    get_right_tokens,
    contains_token,
    get_text_between,
    get_text_splits,
    get_tagged_text,
    is_inverted,
    get_tagged_text,
    rule_regex_search_tagged_text,
    rule_regex_search_btw_AB,
    rule_regex_search_btw_BA,
    rule_regex_search_before_A,
    rule_regex_search_before_B,
    
)

#HELPERS:
def inverted(c):
    return 1 if is_inverted(c) else 0

def distance_btwn(c):
    span0 = c[0]
    span1 = c[1]
    indices0 = set(np.arange(span0.get_word_start(), span0.get_word_end() + 1))
    indices1 = set(np.arange(span1.get_word_start(), span1.get_word_end() + 1))
    if len(indices0.intersection(indices1)) > 0: return 0
    if span0.get_word_start() < span1.get_word_start():
        return span1.get_word_start() - span0.get_word_end() - 1
    else:
        left_span = span1
        return span0.get_word_start() - span1.get_word_end() - 1
    
def overlap(c):
    span0 = c[0]
    span1 = c[1]
    indices0 = set(np.arange(span0.get_word_start(), span0.get_word_end() + 1))
    indices1 = set(np.arange(span1.get_word_start(), span1.get_word_end() + 1))
    if len(indices0.intersection(indices1)) > 0: return 1
    return 0

def ends_in(ci, val, attrib):
    return val == ci.get_attrib_tokens(attrib)[-1]
    
def starts_with(ci, val, attrib):
    return val == ci.get_attrib_tokens(attrib)[0]

# WORD GROUPS:

action_link_words = set(['affect', 'lead', 'led', 'show', 'display', 'exhibit', 'cause', 'result in'])
mutant_words = set(['mutant', 'mutation', 'plant', 'line', 'phenotype', 'seedlings', 'variant'])
helper_vbs = set(['is', 'was', 'are', 'were', 'become', 'became'])
tester_words = set(['sequence', 'published', 'diagram', 'hypothesis', 'hypothesize', 'aim', 'goal', 'understand', 'examine', 'we', 'our', 'experiment', 'test', 'study', 'design', 'analyze', 'analysis', 'results', 'research'])
neg_words = set(['strategy', 'public', 'examine', 'measure', 'subject', 'statistic', 'instance'])


#LABELING_FUNCTIONS:

#"DESC ENT" OR #"ENT DESC"
def DIST_BTWN_0(c):
    return 1 if distance_btwn(c) == 0 else 0

def DESC_ENT(c):
    idx0 = c[0].get_word_end() #desc
    idx1 = c[1].get_word_start() #ent
    return 1 if idx1 == (idx0 + 1) else 0

def ENT_DESC(c):
    idx0 = c[1].get_word_end() #ent
    idx1 = c[0].get_word_start() #desc
    return 1 if idx1 == (idx0 + 1) else 0

#"ENT DESC(if verb past tense)"
def DESC_VERB_PAST(c):
    return 1 if contains_token(c[0], 'VBD', attrib='pos_tags') else 0

#"ENT DESC(if verb present tense)"
def DESC_VERB_PRESENT(c):
    return 1 if contains_token(c[0], 'VB', attrib='pos_tags') else 0

#"ENT DESC(if past participle)"
def DESC_PAST_PART(c):
    return 1 if contains_token(c[0], 'VBN', attrib='pos_tags') else 0    

#"ENT <helper_vb> DESC" or "DESC <helper_vb> ENT" #NOT DONE YET
def HELPER_VERB_BTWN(c):
    return 1 if len(helper_vbs.intersection(set(get_between_tokens(c, attrib='lemmas', n_max=3)))) > 0 else 0

#"ENT <prep phrase> DESC"  #NOT DONE YET
def PREP_PHRASE_BTWN(c):
    return 1 if len(helper_vbs.intersection(set(get_between_tokens(c, attrib='lemmas', n_max=3)))) > 0 else 0

