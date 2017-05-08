with open ("data_400.tsv") as data:
    for line in data:
        pmc_id,text = line.split("\t")
        with open("400_files/{}.txt".format(pmc_id), 'w') as doc_file:
            doc_file.write(text.replace(". ", ". \n"))
