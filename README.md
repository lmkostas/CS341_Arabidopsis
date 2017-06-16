# Arabidopsis
cs341 Spring 2017 project. Learning relationships between genes and phenotypes using snorkel.

## Environment
Clone the snorkel repository. Checkout the brat branch.
Create a virtualenv.
Install snorkel python requirements.
Run the following to install arabidopsis requirements:
```
pip install -r requirements.txt
```
Set up snorkel env:
```
source ../snorkel/set_env.sh
source agp_env.sh small-data
```
Modify snorkel files:
```
cp brat/brat_tools.py ../snorkel/snorkel/contrib/brat/tools.py
cp extended_matchers.py ../snorkel/snorkel/matchers.py
```
Run Jupyter Notebook
```
jupyter notebook
```

## Files
```
Follow the instructions in the section below to download the necessary data.
Then run the following python scripts to process the data:
1. preprocess_results.py (to obtain only th results sections from the selected papers - use preprocess.py to get full text)
2. split_docs.py (transform the processed data into a form that is compatible with snorkel)

Finally, proceed through the following notebooks in the specified order:
1. small_matcher_test-complexPheno.ipynb - Candidate extraction
2. BRAT Import-complexPheno.ipynb - import labels from brat annotations
3. LabelingFunctions-complexPheno.ipynb - Run through generative and discriminative models
```

## Data
Follow the steps below to retrive the data. This will take a while.
```
go to webpage https://www.ncbi.nlm.nih.gov/pmc/?term=Arabidopsis[Text%20Word]
click on the 'send to' dropdown in the top right corner
select file
select XML
click 'create file'
```

A pmc_result.xml file will be downloaded. 
Add the file in the CS341_Arabidopsis github under the following name and path

```
CS341_Arabidopsis/data/full_xml_data.xml
```

The large data file has been added to the .gitignore already.
