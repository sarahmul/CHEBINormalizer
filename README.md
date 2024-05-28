# CHEBINormalizer
This is a repository corresponding to the paper 'Chemical Entity Normalization for Successful Translational Development of Alzheimer’s Disease and Dementia Therapeutics' by Mullin et al.

## Dependencies:
Python 3 (running PubMED+Chebi after mapping file is created):Pytorch, sentence transformers, scikit-learns, scipy, pandas, nltk, pybrat, abbreviations, Levenshtein, fuzzywuzzy

R (for creation of mapping file and data visualization): tidyverse, dplyr, ggplot2, treemapify, caret
## Running PubMED+Chebi
The BERT model with ChEBI incorporated is available on [hugging face](https://huggingface.co/sarahmul/PubmedBERTChEBI).

## Gold Standard
The gold standard file for alzheimer and dementia is available here: '/Data_files/Gold_standard_alzdem.csv.'

## Additional Bioportal mappings
The Bioportal mappings are available here: '/Data_files/MeSHChebi_mapping.csv'



