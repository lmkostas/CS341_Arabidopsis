'''
This file takes an xml file containing a series of pmc papers in html format and saves them as a tsv file where each paper is represented on a single line in the format of [pmc_id <tab> document text].

The script takes three arguments: 1) raw xml data file name, 2) name of tsv file to which processed data will be saved, 3)'small' or 'full' to flag whether a subset or the full set of articles, respectively, is being processed.
'''


from bs4 import BeautifulSoup
import sys

def parseAggrHTML(raw_xml_file, tsv_file):
	full_body_parsed = 0
	abstract_parsed = 0
	no_data_parsed = 0
	total_files = 0
	with open(tsv_file, 'w') as wf:
		with open(raw_xml_file, 'rb') as f:
			print "Parsing XML file..."
			html = BeautifulSoup(f, 'xml')
			for i, doc in enumerate(html.find_all('article')):
				if i%1000 ==0: print i
				total_files += 1
				pmc_id = (doc.find('article-id', {'pub-id-type':'pmc'}).text).encode('utf-8')
				text = ''
				if doc.find('body'):
					full_body_parsed += 1
					#remove all tabs and newlines from article text
					text = ''.join(doc.find('body').find_all(text=True)).encode('utf-8').strip()
					#replace newliens with period for later processing by sentence
					text = text.replace("\n", ". ")
                                	text = text.replace("\t", " ")
                                	wf.write(pmc_id + '\t' + text + '\n')
				#if the article has no body, check for abstract
				elif doc.find('abstract'):
					abstract_parsed += 1
					print 'no body', pmc_id
					#remove all tabs and newlines from article text
					text = ''.join(doc.find('abstract').find_all(text=True)).encode('utf-8').strip()
					#replace newlines with period for later processing by sentence
					text = text.replace("\n", ". ")
                                	text = text.replace("\t", " ")
                                	wf.write(pmc_id + '\t' + text + '\n')
				#if the article has no body or abstract, there is no text to parse so we skip this article
				else:
					no_data_parsed += 1
					print 'no data', pmc_id
		f.close()
	wf.close()
	print "File Parsed!"
	print 'full data parsed: %d / %d files' % (full_body_parsed, total_files)
	print 'abstract parsed: %d / %d files' % (abstract_parsed, total_files)
	print 'no data parsed: %d / %d files' % (no_data_parsed, total_files)

if __name__ == '__main__':
	raw_xml_file = sys.argv[1]
	tsv_file = sys.argv[2]
	data_type = sys.argv[3]
	path = './data/'
	if data_type == 'small':
		path = './small_data_pp/'
	print path+raw_xml_file, path+tsv_file
	parseAggrHTML(path+raw_xml_file, path+tsv_file)
