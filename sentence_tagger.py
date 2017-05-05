import multiprocessing
import os
import sys
import cPickle
import multiprocessing

sys.path.insert(1, '../snorkel')

from snorkel import SnorkelSession
from snorkel.matchers import DictionaryMatch
from revised_extraction import GM, PM
from snorkel.candidates import Ngrams, CandidateSpace, CandidateExtractor
from snorkel.models import Document, Sentence, candidate_subclass
from snorkel.viewer import SentenceNgramViewer

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

train_sents, dev_sents, test_sents, all_sents = set(), set(), set(), set()
docs = session.query(Document).order_by(Document.name).all()
doc_sents = dict()
for split, doc in enumerate(docs):
    if len(doc_sents) >= 50:break
    doc_sents[split] = set()
    for s in doc.sentences:
    	all_sents.add(s)
        doc_sents[split].add(s)
    	name = doc.name.split('-')[0]
        if name in train_ids:
            train_sents.add(s)
        elif name in dev_ids:
            dev_sents.add(s)
        elif name in test_ids:
            test_sents.add(s)
        else:
            raise Exception('ID <{0}> not found in any id set'.format(doc.name))

print "Docs Split"
print "Extracting Candidates..."

#cand_extractor.apply(train_sents, split=0)#, parallelism=multiprocessing.cpu_count())
#train_cands = session.query(GenePhenoPair).filter(GenePhenoPair.split==0).all()
#cand_extractor.apply(dev_sents, split=1)#, parallelism=8)
#dev_cands = session.query(GenePhenoPair).filter(GenePhenoPair.split==1).all()
#cand_extractor.apply(test_sents, split=2)#, parallelism=8)
#test_cands = session.query(GenePhenoPair).filter(GenePhenoPair.split==2).all()

#print "Number of train candidates:", len(train_cands)
#print "Number of dev candidates:", len(dev_cands)
#print "Number of test candidates:", len(test_cands)

for split, sents in doc_sents.iteritems():
    cand_extractor.apply(sents, split=split, parallelism=multiprocessing.cpu_count())
all_cands = session.query(GenePhenoPair).filter(GenePhenoPair.split < len(doc_sents)).all()
print "Number of candidates:", len(all_cands)

# NOTE: This if-then statement is only to avoid opening the viewer during automated testing of this notebook
# You should ignore this!
#sv = SentenceNgramViewer(dev_cands, session, annotator_name = 'gold')
