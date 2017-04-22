import cPickle


def main():
	with open('pmcids_1.pkl', 'rb') as f_1:
	    data_1 = cPickle.load(f_1)

	with open('pmcids_2.pkl', 'rb') as f_2:
	    data_2 = cPickle.load(f_2)

	full_data = {}

	full_data['test'] = data_1['test'] + data_2['test']
	full_data['train'] = data_1['train'] + data_2['train']
	full_data['dev'] = data_1['dev'] + data_2['dev']

	with open('pmcids_400.pkl', 'wb') as id_f:
		cPickle.dump(full_data, id_f)

main()
