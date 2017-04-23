"""
This file takes a tsv file containing a series of pmc papers in where
each paper is represented on a single line in the format of [pmc_id <tab> document text]
and stores each doc-id, senetence tuple in each document in a postgress db, snorkel-small-data
or snorkel-full-data, depending on the dataset being used

The script takes one argument:
1)'small' or 'full' to flag whether a subset or the full set of articles, respectively, is being processed.
e.g. run split_docs.py <small/full>
""" 

import os
import sys

def docs_to_sentences():
	# Must set SNORKELDB before importing SnorkelSession
	from snorkel import SnorkelSession
	from snorkel.parser import TSVDocPreprocessor
	from snorkel.parser import CorpusParser
	from snorkel.models import Document, Sentence
	session = SnorkelSession()

	pathname = 'small_data_pp/small_pp.tsv' if os.environ['AGP_DATA_SIZE'] == 'small-data' else 'data/full_pp.tsv'
	doc_preprocessor = TSVDocPreprocessor(pathname)

	corpus_parser = CorpusParser()
	corpus_parser.apply(doc_preprocessor, parallelism=8)

	print "Documents:", session.query(Document).count()
	print "Sentences:", session.query(Sentence).count()

if __name__ == '__main__':
	docs_to_sentences()
