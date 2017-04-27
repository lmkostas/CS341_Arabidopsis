#Usage: python download_data.py pmcid_all.txt
import sys

list_PMC = sys.argv[1]

file_hash={}

with open(list_PMC, 'r') as file_PMC_read:
	with open('PMC_url_list.txt', 'w') as write_file:
		for line in file_PMC_read:
			write_file.write('ftp://ftp.ncbi.nlm.nih.gov/pub/pmc/' + str(line) + '\n')
			# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3563272/
