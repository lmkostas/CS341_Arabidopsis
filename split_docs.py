'''
This file takes a tsv file containing a series of pmc papers in where
each paper is represented on a single line in the format of [pmc_id <tab> document text]
and stores each doc-id, senetence tuple in each document in a postgress db, snorkel-small-data
or snorkel-full-data, depending on the dataset being used

The script takes one argument:
1)'small' or 'full' to flag whether a subset or the full set of articles, respectively, is being processed.
e.g. run split_docs.py <small/full>
''' 

import os
import sys
#top-level snorkel folder must be in the same directory as the project folder for this to run
sys.path.insert(1, '../snorkel')

def docs_to_sentences(data_size):
	db_name = 'full-data'
	if data_size == 'small':
		db_name = 'small-data'
	#need this line for parallel processing, cannot use sqlite db for this
	os.environ['SNORKELDB'] = 'postgres:///snorkel-'+db_name

	from snorkel import SnorkelSession

	session = SnorkelSession()

	from snorkel.parser import TSVDocPreprocessor
	from snorkel.parser import CorpusParser
	from snorkel.models import Document, Sentence

	pathname = 'data/full_pp.tsv'
	if data_size == 'small':
		pathname = 'small_data_pp/small_pp.tsv'
	doc_preprocessor = TSVDocPreprocessor(pathname)

	corpus_parser = CorpusParser()
	corpus_parser.apply(doc_preprocessor, parallelism=50)

	print "Documents:", session.query(Document).count()
	print "Sentences:", session.query(Sentence).count()

if __name__ == '__main__':
	docs_to_sentences(sys.argv[1])
