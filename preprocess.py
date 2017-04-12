'''
This file takes an xml file containing a series of pmc papers in html format and
saves them as a tsv file where each paper is represented on a single line in the format of 
[pmc_id <tab> document text].

The script takes three arguments:
1) raw xml data file name,
2) name of tsv file to which processed data will be saved, and 
3)'small' or 'full' to flag whether a subset or the full set of articles, respectively, is being processed.
e.g. run preprocess.py <small_data.xml/full_xml_data.xml> <small_pp.tsv/full_data_pp.tsv> <small/full>
''' 


from bs4 import BeautifulSoup
import sys
import cPickle
import random
import math

TRAIN_PCT = 0.8
DEV_PCT = 0.1

def parseAggrHTML(raw_xml_file, tsv_file, id_file):
	full_body_parsed = 0
	abstract_parsed = 0
	no_data_parsed = 0
	total_files = 0
	pmc_ids = []
	with open(tsv_file, 'w') as wf:

		with open(raw_xml_file, 'rb') as f:
			print "Parsing XML file..."
			html = BeautifulSoup(f, 'xml')
			for i, doc in enumerate(html.find_all('article')):
				#if i%1000 ==0: print i
				print i
				total_files += 1
				pmc_id = (doc.find('article-id', {'pub-id-type':'pmc'}).text).encode('utf-8')
				text = ''
				if doc.find('body'):
					full_body_parsed += 1
					#remove all tabs and newlines from article text and convert to unicode
					#conversion to unicode removes characters that cause later errors in snorkel preprocessing
					text = unicode(str(''.join(doc.find('body').find_all(text=True)).encode('utf-8').strip())+'.', errors='ignore')
					#replace newliens with period for later processing by sentence
					text = text.replace("\n", ". ")
					text = text.replace("\t", " ")

					split_count = 1
					if len(text) > 100000:
						sent_start_idx = 0
						while sent_start_idx < len(text):
							end_idx = min(sent_start_idx+100000, len(text))
							last_sent = text[sent_start_idx: end_idx].rindex('.')+1+sent_start_idx
							wf.write("%s\t%s\n" % (pmc_id+'-'+str(split_count), text[sent_start_idx: last_sent]))
							sent_start_idx = last_sent
							split_count += 1
					else:
						wf.write("%s\t%s\n" % (pmc_id, text))
					pmc_ids.append(pmc_id)
                		#if the article has no body, check for abstract
                		elif doc.find('abstract'):
					abstract_parsed += 1
					print pmc_id, ' - full article text not found, parsing abstract instead'
					#remove all tabs and newlines from article text and convert to unicode
					#conversion to unicode removes characters that cause later errors in snorkel preprocessing
					text = unicode(str(''.join(doc.find('abstract').find_all(text=True)).encode('utf-8').strip())+'.', errors='ignore')
					#replace newlines with period for later processing by sentence
					text = text.replace("\n", ". ")
					text = text.replace("\t", " ")
					
					#split docs that contain more than 100K characters into multiple lines in the tsv (ie sub-docs)
					#snorkel doc to sentence preprocessing can only handle at most 100k docs
					#due to db unique key constraints, doc ids for split docs are stored as: <pmc_id>-<split_count>
					split_count = 1
					if len(text) > 100000:
						start_idx = 0
						while start_idx < len(text):
							end_idx = min(start_idx+100000, len(text))
							last_sent = text[start_idx: end_idx].rindex('.')+1+sent_start_idx
							wf.write("%s\t%s\n" % (pmc_id+'-'+str(split_count), text[start_idx: last_sent]))
							if last_sent == len(text): break
							start_idx = last_sent
					else:
						wf.write("%s\t%s\n" % (pmc_id, text))
					pmc_ids.append(pmc_id)
				#if the article has no body or abstract, there is no text to parse so we skip this article
				else:
					no_data_parsed += 1
					print pmc_id, ' - Error: no text data' 
		f.close()
	wf.close()
	print "File Parsed!"
	print 'full data parsed: %d / %d files' % (full_body_parsed, total_files)
	print 'abstract parsed: %d / %d files' % (abstract_parsed, total_files)
	print 'no data parsed: %d / %d files' % (no_data_parsed, total_files)

	#get pkl file of document ids for train, dev, and test set
	print "Segmenting and Writing PMC_ids..."
	random.shuffle(pmc_ids)

	parsed_files = full_body_parsed + abstract_parsed
	train_end_idx = int(math.ceil(TRAIN_PCT*parsed_files))
	dev_end_idx = train_end_idx + int(math.ceil(DEV_PCT*parsed_files))

	with open(id_file, 'wb') as id_f:
		cPickle.dump({'train': pmc_ids[0:train_end_idx], 'dev': pmc_ids[train_end_idx:dev_end_idx], 'test': pmc_ids[dev_end_idx:]}, id_f)
	
	print "PMC_ids written!"


if __name__ == '__main__':
	raw_xml_file = sys.argv[1]
	tsv_file = sys.argv[2]
	data_type = sys.argv[3]
	path = './data/'
	if data_type == 'small':
		path = './small_data_pp/'
	parseAggrHTML(path+raw_xml_file, path+tsv_file, path+'pmc_ids.pkl')
