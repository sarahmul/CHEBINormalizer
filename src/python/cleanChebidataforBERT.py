#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 13 22:06:11 2021

@author: sarahmullin
"""
import pandas as pd

aa=pd.read_csv('/Users/sarahmullin/Desktop/Menagerie/Chemical_normalization/Chebi_relations_BERT.csv')

with open('drugs_onto_input.txt', 'w') as a:
    for index, row in aa.iterrows():
        group=row['grouped2']
        print(group)
        a.write(group+'\n')
        a.write('\n')
