import os
import sys
sys.path.insert(1, '../snorkel')

def docs_to_sentences(data_size):
	# TO USE A DATABASE OTHER THAN SQLITE, USE THIS LINE
	# Note that this is necessary for parallel execution amongst other things...
	db_name = 'full-data'
	if data_size == 'small':
		db_name = 'small-data'
	#os.environ['SNORKELDB'] = 'postgres:///snorkel.db'#+db_name
	#print os.environ['SNORKELDB']

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
	corpus_parser.apply(doc_preprocessor)#, parallelism=20)

	print "Documents:", session.query(Document).count()
	print "Sentences:", session.query(Sentence).count()

if __name__ == '__main__':
	docs_to_sentences(sys.argv[1])
