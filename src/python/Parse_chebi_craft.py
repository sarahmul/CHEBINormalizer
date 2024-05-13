#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script cleans CRAFT 3.0
@author: sarahmullin
"""
import json
from pymongo import MongoClient
import pandas as pd
import nltk
import pickle as pkl
import os
import numpy as np
import dataclasses
from pybrat.parser import BratParser, Entity, Event, Example, Relation
from nltk.stem import WordNetLemmatizer
from nltk.corpus import wordnet
from abbreviations import schwartz_hearst


# convert greek alphabet symbols into word form 

greek_alphabet = {
    u'\u0391': 'Alpha',
    u'\u0392': 'Beta',
    u'\u0393': 'Gamma',
    u'\u0394': 'Delta',
    u'\u0395': 'Epsilon',
    u'\u0396': 'Zeta',
    u'\u0397': 'Eta',
    u'\u0398': 'Theta',
    u'\u0399': 'Iota',
    u'\u039A': 'Kappa',
    u'\u039B': 'Lamda',
    u'\u039C': 'Mu',
    u'\u039D': 'Nu',
    u'\u039E': 'Xi',
    u'\u039F': 'Omicron',
    u'\u03A0': 'Pi',
    u'\u03A1': 'Rho',
    u'\u03A3': 'Sigma',
    u'\u03A4': 'Tau',
    u'\u03A5': 'Upsilon',
    u'\u03A6': 'Phi',
    u'\u03A7': 'Chi',
    u'\u03A8': 'Psi',
    u'\u03A9': 'Omega',
    u'\u03B1': 'alpha',
    u'\u03B2': 'beta',
    u'\u03B3': 'gamma',
    u'\u03B4': 'delta',
    u'\u03B5': 'epsilon',
    u'\u03B6': 'zeta',
    u'\u03B7': 'eta',
    u'\u03B8': 'theta',
    u'\u03B9': 'iota',
    u'\u03BA': 'kappa',
    u'\u03BB': 'lamda',
    u'\u03BC': 'mu',
    u'\u03BD': 'nu',
    u'\u03BE': 'xi',
    u'\u03BF': 'omicron',
    u'\u03C0': 'pi',
    u'\u03C1': 'rho',
    u'\u03C3': 'sigma',
    u'\u03C4': 'tau',
    u'\u03C5': 'upsilon',
    u'\u03C6': 'phi',
    u'\u03C7': 'chi',
    u'\u03C8': 'psi',
    u'\u03C9': 'omega',
    u'\u0391': 'Alpha',
    u'\u0392': 'Beta',
    u'\u0393': 'Gamma',
    u'\u0394': 'Delta',
    u'\u0395': 'Epsilon',
    u'\u0396': 'Zeta',
    u'\u0397': 'Eta',
    u'\u0398': 'Theta',
    u'\u0399': 'Iota',
    u'\u039A': 'Kappa',
    u'\u039B': 'Lamda',
    u'\u039C': 'Mu',
    u'\u039D': 'Nu',
    u'\u039E': 'Xi',
    u'\u039F': 'Omicron',
    u'\u03A0': 'Pi',
    u'\u03A1': 'Rho',
    u'\u03A3': 'Sigma',
    u'\u03A4': 'Tau',
    u'\u03A5': 'Upsilon',
    u'\u03A6': 'Phi',
    u'\u03A7': 'Chi',
    u'\u03A8': 'Psi',
    u'\u03A9': 'Omega',
    u'\u03B1': 'alpha',
    u'\u03B2': 'beta',
    u'\u03B3': 'gamma',
    u'\u03B4': 'delta',
    u'\u03B5': 'epsilon',
    u'\u03B6': 'zeta',
    u'\u03B7': 'eta',
    u'\u03B8': 'theta',
    u'\u03B9': 'iota',
    u'\u03BA': 'kappa',
    u'\u03BB': 'lamda',
    u'\u03BC': 'mu',
    u'\u03BD': 'nu',
    u'\u03BE': 'xi',
    u'\u03BF': 'omicron',
    u'\u03C0': 'pi',
    u'\u03C1': 'rho',
    u'\u03C3': 'sigma',
    u'\u03C4': 'tau',
    u'\u03C5': 'upsilon',
    u'\u03C6': 'phi',
    u'\u03C7': 'chi',
    u'\u03C8': 'psi',
    u'\u03C9': 'omega',
    u'\u03B1':  'degree',
    u'\u03B1'  : 'alpha',
u'\u03B2':'beta',
u'\u03B3':'gamma',
u'\u03B4':'delta',
u'\u03B5':'epsilon',
u'\u03B6':'zeta',
u'\u03B7':'eta',
u'\u03B8':'theta',
u'\u03B9':'iota',
u'\u03BA':'kappa',
u'\u03BB':'lamda',
u'\u03BC':'mu',
u'\u03BD':'nu',
u'\u03BE':'xi',
u'\u03BF':'omicron',
u'\u03C0':'pi',
u'\u03C1':'rho',
u'\u03C2':'final sigma',
u'\u03C3':'sigma',
u'\u03C4':'tau',
u'\u03C5':'upsilon',
u'\u03C6':'phi',
u'\u03C7':'chi',
u'\u03C8':'psi',
u'\u03C9':'omega',
'α' : 'alpha',
'β' : 'beta',
'γ' : 'gamma',
'δ' : 'delta',
'ε' : 'epsilon',
'ζ' : 'zeta',
'η' : 'eta',
'θ' : 'theta',
'ι' : 'iota',
'κ' : 'kappa',
'λ' : 'lamda',
'μ' : 'mu',
'ν' : 'nu',
'ξ' : 'xi',
'ο' : 'omicron',
'π' : 'pi',
'ρ' : 'rho',
'ς' : 'final sigma',
'σ' : 'sigma',
'τ' : 'tau',
'υ' : 'upsilon',
'φ' : 'phi',
'χ' : 'chi',
'ψ' : 'psi',
'ω' : 'omega'
}

os.chdir('../../Data_files/craft-3.0/ontology-concepts/CHEBI/CHEBI/')

brat = BratParser(error="ignore")
examples = brat.parse("brat/")

examples2 = [*map(dataclasses.asdict, examples)]

list_data=[]    
for i in range(len(examples2)):
    id_=examples2[i]['id']
    text_=examples2[i]['text']
    pairs = schwartz_hearst.extract_abbreviation_definition_pairs(doc_text=text_)
    for j in range(len(examples2[i]['entities'])):
        mention2=examples2[i]['entities'][j]['mention'].replace('µ', 'mu').replace('°',' degrees ').replace('-','-') #.replace('', 'alpha').replace('', 'beta').replace('', 'gamma').replace('', 'delta').replace('','epsilon').replace('','tau')
        abbrev=''
        for key,value in pairs.items():
            if key==examples2[i]['entities'][j]['mention']:
                abbrev=value
        for emote, replacement in greek_alphabet.items():
            mention2 = mention2.replace(emote, replacement)
        chebi_id=examples2[i]['entities'][j]['type']
        start=examples2[i]['entities'][j]['spans'][0]['start']
        end=examples2[i]['entities'][j]['spans'][0]['end']
        id_2=str(id_)+mention2+str(start)+str(end)
        text2=text_[max(start-300,0):min(end+300, len(text_))].replace('\n', ' ').replace('\r', '').replace('Œ|¬|∞|~|°','').replace('µ', 'mu').replace('°',' degrees ')
        for emote, replacement in greek_alphabet.items():
            text2 = text2.replace(emote, replacement)
        text3=text_[max(start-200,0):min(end+200, len(text_))].replace('\n', ' ').replace('\r', '').replace('Œ|¬|∞|º|°','').replace('µ', 'mu').replace('°',' degrees ')
        for emote, replacement in greek_alphabet.items():
            text3 = text3.replace(emote, replacement)
        text22=nltk.sent_tokenize(text2)
        text23='\n'.join(s for s in text22 if mention2 in s)
        print(abbrev)
        list_data.append([id_,id_2,mention2, abbrev,chebi_id, start, end, text23, text3])

        
df_chebi=pd.DataFrame(list_data, columns=["id_","id_2","mention", "abbreviation","chebi_id", "start", "end", "sentence", "surroung"])        
df_chebi.to_csv('../../Data_files/CRAFT_chebi.csv', index=False)        


##### chebi extension
os.chdir('../../Data_files/craft-3.0/ontology-concepts/CHEBI/CHEBI+extensions')
examples_ext = brat.parse("brat/")
examples2_ext = [*map(dataclasses.asdict, examples_ext)]


list_data_ext=[]    
for i in range(len(examples2_ext)):
    id_=examples2_ext[i]['id']
    text_=examples2_ext[i]['text']
    #text_sent=nltk.sent_tokenize(text_)
    #text_sent
    for j in range(len(examples2_ext[i]['entities'])):
        mention2=examples2_ext[i]['entities'][j]['mention']
        chebi_id=examples2_ext[i]['entities'][j]['type']
        start=examples2_ext[i]['entities'][j]['spans'][0]['start']
        end=examples2_ext[i]['entities'][j]['spans'][0]['end']
        
        text2=text_[max(start-300,0):min(end+300, len(text_))].replace('\n', ' ').replace('\r', '').replace('Œ|¬|∞|~|°','').replace('µ', 'mu').replace('°',' degrees ')
        text3=text_[max(start-200,0):min(end+200, len(text_))].replace('\n', ' ').replace('\r', '').replace('Œ|¬|∞|º|°','').replace('µ', 'mu').replace('°',' degrees ')
        text22=nltk.sent_tokenize(text2)
        text23='\n'.join(s for s in text22 if mention2 in s)
        list_data_ext.append([id_,mention2, chebi_id, start, end, text23, text3])
 
        
df_chebi=pd.DataFrame(list_data_ext, columns=["id_","mention", "chebi_id", "start", "end", "sentence", "surroung"])        
df_chebi.to_csv('../../Data_files/CRAFT_chebiext.csv', index=False)        


###cleaning, lemmatization, normalization

import pandas as pd
title_df = pd.DataFrame(columns = ['pmid', 'title_proc', 'title_lem'])

for paper in progress_bar(alz.pubmed.find()):
    pmid=paper['pmid']
    title=paper['derived']['processed_title']
    title_lem=paper['derived']['lemmatized_title']
    title_df = title_df.append({'pmid' : pmid, 'title_proc' : title, 'title_lem' :title_lem}, 
                ignore_index = True)
    
title_df.to_csv('/Users/sarahmullin/Desktop/Menagerie/Chemical_normalization/chemical_titles_12621.csv', index=False)
    
#####make dataset for our annotation

def normalize_term(term):
    return " ".join(nltk.word_tokenize(term))

# function to convert nltk tag to wordnet tag
def nltk_tag_to_wordnet_tag(nltk_tag):
    if nltk_tag.startswith("J"):
        return wordnet.ADJ
    elif nltk_tag.startswith("V"):
        return wordnet.VERB
    elif nltk_tag.startswith("N"):
        return wordnet.NOUN
    elif nltk_tag.startswith("R"):
        return wordnet.ADV
    else:
        return None


def lemmatize_multiword(term):
    lemmatizer = WordNetLemmatizer()
    # tokenize the sentence and find the POS tag for each token
    nltk_tagged = nltk.pos_tag(nltk.word_tokenize(term))
    # tuple of (token, wordnet_tag)
    wordnet_tagged = map(lambda x: (x[0], nltk_tag_to_wordnet_tag(x[1])), nltk_tagged)
    lemmatized_term = []
    for word, tag in wordnet_tagged:
        if tag is None:
            # if there is no available tag, append the token as is
            lemmatized_term.append(word)
        else:
            # else use the tag to lemmatize the token
            lemmatized_term.append(lemmatizer.lemmatize(word, tag))
    return " ".join(lemmatized_term)

chem_df = pd.DataFrame(columns = ['pmid', 'id_2',  'term','abbrev', 'term_norm','term_lem','term_lem_normal','chebi_id', 'title_proc', 'title_lem', 'title_proc_sur', 'title_lem_sur' ])

for index, row in df_chebi.iterrows():
    term=row['mention']
    chem_df=chem_df.append({'pmid': row['id_'], 
                            'id_2':row['id_2'], 
                            'term':term,
                            'abbrev':row['abbreviation'],
                            'term_norm':normalize_term(term),
                            'term_lem':lemmatize_multiword(term),
                            'term_lem_normal':lemmatize_multiword(normalize_term(term)),
                            'chebi_id':row['chebi_id'], 
                            'title_proc':row['sentence'], 
                            'title_lem':lemmatize_multiword(row['sentence']), 
                            'title_proc_sur':row['surroung'], 
                            'title_lem_sur':lemmatize_multiword(row['surroung'])}, 
                           ignore_index=True)

### remove Bad Chebi ids (aka. keep only 2/3 star chebi)


Chebi_compounds2 = pd.read_csv('../ChEBI/compounds.tsv', sep='\t', nrows=183324)
Chebi_compounds_filt2 = Chebi_compounds2.loc[(Chebi_compounds2['STATUS'].isin(['C', 'E'])) & (Chebi_compounds2['STAR'].isin([2, 3]))]
Chebi_compounds_filt2['new'] = 1

Chebi_check_across = pd.merge(Chebi_compounds_filt2[['ID', 'new']], Chebi_compounds_filt[['ID', 'old']], how='outer')
Chebi_check_across['KIND'] = (Chebi_check_across['old'] == 1) & (Chebi_check_across['new'] == 1)

Chebi_compounds_filt_badIDs = Chebi_compounds2.loc[(Chebi_compounds2['STAR'].isin([1, 0])) | ((Chebi_compounds2['STAR'].isin([2, 3])) & (Chebi_compounds2['STATUS'].isin(['D', 'F', 'O', 'S'])))]

#### 


Chebi_compounds_filt_badIDs = Chebi_compounds_filt_badIDs.drop_duplicates(subset=['ID'])
Chebi_compounds_filt_badIDs['BadChebi'] = 1
Chebi_compounds_filt_badIDs = Chebi_compounds_filt_badIDs.rename(columns={'ID': 'chebi_id'})

chem_df['chebi_id'] = chem_df['chebi_id'].str.replace('CHEBI_', '').astype(float)

chem_df_bad = pd.merge(chem_df, Chebi_compounds_filt_badIDs, how='left')
chem_df_bad = chem_df_bad[chem_df_bad['BadChebi'] == 1]

chem_df = pd.merge(chem_df, Chebi_compounds_filt_badIDs, how='left')
chem_df = chem_df[chem_df['BadChebi'].isna()]

chem_df['term_lower'] = chem_df['term'].str.lower()
chem_df = chem_df.drop_duplicates(subset=['id_2', 'term'], keep='first')

chem_df['term_norm'] = chem_df['term_norm'].str.strip()
chem_df['term_lem_normal'] = chem_df['term_lem_normal'].str.strip()

chem_df['term_lower_lemma'] = chem_df['term_lem_normal'].str.lower()
chem_df['punc'] = chem_df['term_lower_lemma'].str.replace(r'[^a-zA-Z\s]', ' ')

##filter out the gene and protein related mentions

chem_df = chem_df[~chem_df['term'].isin(['amyloid', 'Amyloid', 'aqueous', 
                                         'Aqueous', 'chow', 'Chow', 'Feed', 
                                         'm.', 'message', 'messages', 'messenger', 
                                         'messengers', 'compound', 'compounds', 
                                         'water', 'Water', 'H20'])]

            
chem_df.to_csv('../../Data_files/craft_anno.csv', index=False) 

    

#### create craft-chebi json for REEl

###
f = lambda x: x.tolist() if len(x) > 1 else x
df_REEL = chem_df.groupby(['pmid'])['term_lower'].agg(f).reset_index().reindex(chem_df.columns, axis=1)

df_REEL_json=df_REEL[['pmid','term_lower']]
df_REEL_json.set_index('pmid', inplace=True)
df_dict=df_REEL_json.to_dict(orient='dict')

json_object = json.dumps(df_dict, indent=0)
 
# Writing to sample.json
with open("../../Data_files/craft_REEL_slim2.json", "w") as outfile:
    outfile.write(json_object)
#df_REEL_json.to_json(orient='values', index=True)



#####
import json
 
# Opening JSON file
with open('/Users/sarahmullin/Desktop/Menagerie/Chemical_normalization/REEL_alz/craft_kb_corpus_linkREEL_craft_results.json') as json_file:
    data_ = json.load(json_file)
data_set=[]
for k,v in data_.items():
    pmid=k
    for m,l in v.items():
        term=m
        ID_level1=l
        data_set.append([pmid, term,ID_level1])
        

df_REEL=pd.DataFrame(data_set, columns=["pmid","term", 'ID_level1'])
df_REEL.to_csv('/Users/sarahmullin/Desktop/Menagerie/Chemical_normalization/CRAFT_REEL_kbcorpus_link.csv', index=False)        


