"""
This file takes an xml file containing a series of pmc papers in html format and
saves them as a tsv file where each paper is represented on a single line in the format of 
[pmc_id <tab> document text].

*Note - this script returns the entire contents of each article
*Note - all new lines are treated as new sentences

The script takes three arguments:
1) raw xml data file name,
2) name of tsv file to which processed data will be saved, and 
3) name of pkl file to which pmc_ids will be saved
e.g. run preprocess.py <small_data.xml/full_xml_data.xml> <small_pp.tsv/full_data_pp.tsv> <pmc_ids.pkl>
"""

import argparse
from bs4 import BeautifulSoup
import sys
import cPickle
import random
import math

TRAIN_PCT = 0.8
DEV_PCT = 0.1
ESCAPE_CHARS = {
    u"\u2212": "-",
    u"\u2010": "-",
    u"\u2011": "-",
    u"\u2012": "-",
    u"\u2013": "-",
    u"\u2014": "-",
    u"\u0391": "alpha",
    u"\u03B1": "alpha",
    u"\u0392": "beta",
    u"\u03B2": "beta",
    u"\u0393": "gamma",
    u"\u03B3": "gamma",
    u"\u0394": "delta",
    u"\u03B4": "delta",
    u"\u0395": "epsilon",
    u"\u03B5": "epsilon",
    u"\u0396": "zeta",
    u"\u03B6": "zeta",
    u"\u0397": "eta",
    u"\u03B7": "eta",
    u"\u0398": "theta",
    u"\u03B8": "theta",
    u"\u0399": "iota",
    u"\u03B9": "iota",
    u"\u039A": "kappa",
    u"\u03BA": "kappa",
    u"\u039B": "lambda",
    u"\u03BB": "lambda",
    u"\u039C": "mu",
    u"\u03BC": "mu",
    u"\u039D": "nu",
    u"\u03BD": "nu",
    u"\u039E": "xi",
    u"\u03BE": "xi",
    u"\u039F": "omicron",
    u"\u03BF": "omicron",
    u"\u03A0": "pi",
    u"\u03C0": "pi",
    u"\u03A1": "rho",
    u"\u03C1": "rho",
    #u"\u03A2": "sigma",
    u"\u03C2": "sigma",
    u"\u03C3": "sigma",
    u"\u03A4": "tau",
    u"\u03C4": "tau",
    u"\u03A5": "upsilon",
    u"\u03C5": "upsilon",
    u"\u03A6": "phi",
    u"\u03C6": "phi",
    u"\u03A7": "chi",
    u"\u03C7": "chi",
    u"\u03A8": "psi",
    u"\u03C8": "psi",
    u"\u03A9": "omega",
    u"\u03C9": "omega"
}

"""
Function to replace escape characters their unicode representable equivalent and
replace greek letter symbols with their english spelling equivalents
"""
def replace_all_escaped(text):
    for char in ESCAPE_CHARS:
        text = text.replace(char, ESCAPE_CHARS[char])
    return text

def parseAggrHTML(raw_xml_file, tsv_file, id_file, data_size_flag, required_pmcids):
    full_body_parsed = 0
    abstract_parsed = 0
    no_data_parsed = 0
    total_files = 0
    pmc_ids = []
    
    small_file_count = 0
    nonreq_count = 0
    req_file_count = 0
    random.seed(5)

    if (data_size_flag == 'small'):
        print "running small version...."

    with open(tsv_file, 'w') as wf:

        with open(raw_xml_file, 'rb') as f:
            print "Parsing XML file..."
            html = BeautifulSoup(f, 'xml')
            for i, doc in enumerate(html.find_all('article')):
                if i%100 ==0: print i
                # print i
                total_files += 1
                pmc_id = 'PMC'+(doc.find('article-id', {'pub-id-type':'pmc'}).text).encode('utf-8')
                idISSelected = (pmc_id in required_pmcids)
                if (data_size_flag == 'small') and (not idISSelected):
                    if nonreq_count >= 174 or random.random() > .5:
                    	continue

                print "parsing following id: ", pmc_id
                    
                text = ''
                if doc.find('body'):
                    small_file_count += 1
                    if idISSelected:
                        req_file_count += 1
                    else:
                    	nonreq_count += 1

                    full_body_parsed += 1
                    #remove all tabs and newlines from article text and convert to unicode
                    #conversion to unicode removes characters that cause later errors in snorkel preprocessing
                    text = replace_all_escaped(''.join(doc.find('body').find_all(text=True)))
                    text = unicode(str(text.encode('utf-8').strip())+'.', errors='ignore')
                    # text = unicode(str(''.join(doc.find('body').find_all(text=True)).encode('utf-8').strip())+'.', errors='ignore')
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
                elif doc.find('abstract'):
                    #if the article has no body, check for abstract

                    abstract_parsed += 1
                    print pmc_id, ' - full article text not found, parsing abstract instead'
                    #remove all tabs and newlines from article text and convert to unicode
                    #conversion to unicode removes characters that cause later errors in snorkel preprocessing
                    text = replace_all_escaped(''.join(doc.find('body').find_all(text=True)))
                    text = unicode(str(text.encode('utf-8').strip())+'.', errors='ignore')                    
                    # text = unicode(str(''.join(doc.find('abstract').find_all(text=True)).encode('utf-8').strip())+'.', errors='ignore')
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
    print "num total files used: ", small_file_count
    print "not req files used: ", nonreq_count
    print "req files used: ", req_file_count

def buildSelectedPMCIDSet(file):
    selected_set = set(line.strip() for line in open(file))
    return selected_set

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("raw_xml_file",
                        help="Name of raw xml file from PubMed which contains all articles. e.g. small_data_pp/small_data.xml")
    parser.add_argument("tsv_file",
                        help="Name of tsv output file to which processed data will be saved.")
    parser.add_argument("pmc_ids_file",
                        help="Name of pmc_ids pkl file to which ids will be saved")
    parser.add_argument('-d', '--data_size', choices=["small","all"], default="small", help="what size dataset (small or all)")
    args = parser.parse_args()
    # print args.data_size
    required_pmcids = set()
    if (args.data_size == 'small'):
        required_pmcids = buildSelectedPMCIDSet("small_data/selected_pmcid.txt")
        # print required_pmcids, len(required_pmcids)
    parseAggrHTML(args.raw_xml_file, args.tsv_file, args.pmc_ids_file, args.data_size, required_pmcids)
