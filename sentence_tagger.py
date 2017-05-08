import multiprocessing
import os
import sys
import cPickle
import multiprocessing

sys.path.insert(1, '../snorkel')

from snorkel import SnorkelSession
from snorkel.matchers import DictionaryMatch
from final_candidates import GM, PM
from snorkel.candidates import Ngrams, CandidateSpace, CandidateExtractor
from snorkel.models import Document, Sentence, candidate_subclass
from snorkel.viewer import SentenceNgramViewer

SPLIT_ON_DOCS = False
ALL_DOCS = True # if true, create train dev and test. if false, push everything to dev cands.

session = SnorkelSession()

GenePhenoPair = candidate_subclass('GenePhenoPair',['gene', 'pheno'])

gene_ngrams = Ngrams(n_max=5)
pheno_ngrams = Ngrams(n_max=10)
cand_extractor = CandidateExtractor(GenePhenoPair, 
                                    [gene_ngrams, pheno_ngrams], [GM, PM],
                                    symmetric_relations=True)

print "Splitting Docs..."
pathname = 'small_data/' if os.environ['AGP_DATA_SIZE'] == 'small-data' else 'data/'
with open(pathname+'pmcids_400.pkl', 'rb') as f:
    sent_dicts = cPickle.load(f)
train_ids, dev_ids, test_ids = set(sent_dicts['train']), set(sent_dicts['dev']), set(sent_dicts['test'])
all_ids = train_ids.union(dev_ids).union(test_ids)
# 40, 10, 10
train_sents, dev_sents, test_sents, all_sents = set(), set(), set(), set()
train_docs, dev_docs, test_docs = set(), set(), set()
docs = session.query(Document).order_by(Document.name).all()
doc_sents = dict()
for doc_num, doc in enumerate(docs):
    if len(train_docs) >= 40 and len(dev_docs) >= 10 and len(test_docs) >= 10:break
    doc_sents[doc_num] = set()
    for s in doc.sentences:
	all_sents.add(s)
	doc_sents[doc_num].add(s)
	name = doc.name.split('-')[0]
	if name in train_ids:
            train_docs.add(name)
	    train_sents.add(s)
	elif name in dev_ids:
            dev_docs.add(name)
	    dev_sents.add(s)
	elif name in test_ids:
            test_docs.add(name)
	    test_sents.add(s)
	else:
	    raise Exception('ID <{0}> not found in any id set'.format(doc.name))

print "Docs Split"
print "Extracting Candidates..."


if SPLIT_ON_DOCS:
    for split, sents in doc_sents.iteritems():
        cand_extractor.apply(sents, split=split, parallelism=multiprocessing.cpu_count())
    all_cands = session.query(GenePhenoPair).filter(GenePhenoPair.split < len(doc_sents)).all()
    print "Number of candidates:", len(all_cands)
else:
    if ALL_DOCS:
        cand_extractor.apply(train_sents, split=0, parallelism=multiprocessing.cpu_count())
        train_cands = session.query(GenePhenoPair).filter(GenePhenoPair.split==0).all()
        cand_extractor.apply(dev_sents, split=1, parallelism=multiprocessing.cpu_count())
        dev_cands = session.query(GenePhenoPair).filter(GenePhenoPair.split==1).all()
        cand_extractor.apply(test_sents, split=2, parallelism=multiprocessing.cpu_count())
        test_cands = session.query(GenePhenoPair).filter(GenePhenoPair.split==2).all()
        print "Number of train candidates:", len(train_cands)
        print "Number of dev candidates:", len(dev_cands)
        print "Number of test candidates:", len(test_cands)
    else:
        cand_extractor.apply(all_sents, split=1, parallelism=multiprocessing.cpu_count())
        dev_cands = session.query(GenePhenoPair).filter(GenePhenoPair.split==1).all()
        print "Number of dev candidates:", len(dev_cands)
        
# NOTE: This if-then statement is only to avoid opening the viewer during automated testing of this notebook
# You should ignore this!
#sv = SentenceNgramViewer(dev_cands, session, annotator_name = 'gold')
