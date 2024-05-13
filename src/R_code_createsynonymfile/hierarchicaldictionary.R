#######################################################################################
#ChEBI normalization using BERT: hierarchical matching and candidate generation
#######################################################################################


#####
library(dplyr)
library(readr)
library(stringr)
library(curl)

### pull in things to match
data_to_match_CHEBI<-read_csv('/Users/sarahmullin/Desktop/Menagerie/Chemical_normalization/chemical_anno_12921.csv') 
data_to_match_CHEBI$term_lower<-tolower(data_to_match_CHEBI$term)
data_to_match_CHEBI<-data_to_match_CHEBI %>% distinct(pmid, term_lower, .keep_all=T) %>% select(-term_lower)
data_to_match_CHEBI$term_norm<-str_trim(data_to_match_CHEBI$term_norm, side='both')
data_to_match_CHEBI$term_lem_normal<-str_trim(data_to_match_CHEBI$term_lem_normal, side='both')
data_to_match_CHEBI$term_lower<-tolower(str_trim(data_to_match_CHEBI$term, side='both'))

data_to_match_CHEBI$term_lower<-tolower(data_to_match_CHEBI$term_norm)
data_to_match_CHEBI$term_lower_lemma<-tolower(data_to_match_CHEBI$term_lem_normal)
data_to_match_CHEBI$punc<-gsub("[^[:alnum:][:space:]]"," ", data_to_match_CHEBI$term_lower_lemma)

## exact match
data_to_match_CHEBI_out<-left_join(data_to_match_CHEBI %>% rename(text_norm_lower=term_lower), Chebi_compounds_filt2 %>% filter(TYPE=='ChebiName') %>% select(ID_level1, text_norm_lower) ) %>% filter(is.na(ID_level1))
data_to_match_CHEBI_exact<-left_join(data_to_match_CHEBI %>% rename(text_norm_lower=term_lower), Chebi_compounds_filt2 %>% filter(TYPE=='ChebiName') %>% select(ID_level1, text_norm_lower) ) %>% filter(!is.na(ID_level1))
data_to_match_CHEBI_exact$kind<-'exact'

###exact relaxed 
data_to_match_CHEBI_exactlem<-left_join(data_to_match_CHEBI_out %>% select(-ID_level1) %>% rename(text_lower_lemma=term_lower_lemma) , Chebi_compounds_filt2 %>% filter(TYPE=='ChebiName') %>% select(ID_level1, text_lower_lemma) ) %>% filter(!is.na(ID_level1))
data_to_match_CHEBI_out<-left_join(data_to_match_CHEBI_out %>% select(-ID_level1)  %>% rename(text_lower_lemma=term_lem) , Chebi_compounds_filt2 %>% filter(TYPE=='ChebiName') %>% select(ID_level1, text_lower_lemma) ) %>% filter(is.na(ID_level1))
data_to_match_CHEBI_exactlem$kind<-'exactlem'
###exact relaxed punc
data_to_match_CHEBI_exactlempunc<-left_join(data_to_match_CHEBI_out %>% select(-ID_level1)  , Chebi_compounds_filt2 %>% filter(TYPE=='ChebiName') %>% select(ID_level1, punc) ) %>% filter(!is.na(ID_level1))
data_to_match_CHEBI_out<-left_join(data_to_match_CHEBI_out %>% select(-ID_level1)  , Chebi_compounds_filt2 %>% filter(TYPE=='ChebiName') %>% select(ID_level1, punc) ) %>% filter(is.na(ID_level1))
data_to_match_CHEBI_exactlempunc$kind<-'exactlempunc'
###exact relaxed punc
data_to_match_CHEBI_exactlempunc2<-left_join(data_to_match_CHEBI_out %>% select(-ID_level1) %>% rename(punc2=punc) , Chebi_compounds_filt2 %>% filter(TYPE=='ChebiName') %>% select(ID_level1, punc2) ) %>% filter(!is.na(ID_level1))
data_to_match_CHEBI_out<-left_join(data_to_match_CHEBI_out %>% select(-ID_level1) %>% rename(punc2=punc) , Chebi_compounds_filt2 %>% filter(TYPE=='ChebiName') %>% select(ID_level1, punc2) ) %>% filter(is.na(ID_level1))
data_to_match_CHEBI_out<-data_to_match_CHEBI_out %>% rename(punc=punc2)
data_to_match_CHEBI_exactlempunc2$kind<-'exactlempunc'
## exact match to syn
data_to_match_CHEBI_syn<-left_join(data_to_match_CHEBI_out %>% select(-ID_level1) , Chebi_compounds_filt2 %>% filter(TYPE!='ChebiName') %>% select(ID_level1, text_norm_lower) ) %>% filter(!is.na(ID_level1))
data_to_match_CHEBI_out<-left_join(data_to_match_CHEBI_out %>% select(-ID_level1)   , Chebi_compounds_filt2 %>% filter(TYPE!='ChebiName') %>% select(ID_level1, text_norm_lower) ) %>% filter(is.na(ID_level1))
data_to_match_CHEBI_syn$kind<-'syn'
###exact relaxed  to syn
data_to_match_CHEBI_synlem<-left_join(data_to_match_CHEBI_out %>% select(-ID_level1)  , Chebi_compounds_filt2 %>% filter(TYPE!='ChebiName' ) %>% select(ID_level1, text_lower_lemma) ) %>% filter(!is.na(ID_level1))
data_to_match_CHEBI_out<-left_join(data_to_match_CHEBI_out %>% select(-ID_level1)  , Chebi_compounds_filt2 %>% filter(TYPE!='ChebiName') %>% select(ID_level1, text_lower_lemma) ) %>% filter(is.na(ID_level1))
data_to_match_CHEBI_synlem$kind<-'synlem'
###exact relaxed punc to syn
data_to_match_CHEBI_synlempunc<-left_join(data_to_match_CHEBI_out %>% select(-ID_level1)  , Chebi_compounds_filt2 %>% filter(TYPE!='ChebiName') %>% select(ID_level1, punc) ) %>% filter(!is.na(ID_level1))
data_to_match_CHEBI_out<-left_join(data_to_match_CHEBI_out %>% select(-ID_level1)  , Chebi_compounds_filt2 %>% filter(TYPE!='ChebiName') %>% select(ID_level1, punc) ) %>% filter(is.na(ID_level1))
data_to_match_CHEBI_synlempunc$kind<-'synlempunc'
###exact relaxed punc to syn
data_to_match_CHEBI_synlempunc2<-left_join(data_to_match_CHEBI_out %>% select(-ID_level1) %>% rename(punc2=punc)  , Chebi_compounds_filt2 %>% filter(TYPE!='ChebiName') %>% select(ID_level1, punc2) ) %>% filter(!is.na(ID_level1))
data_to_match_CHEBI_out<-left_join(data_to_match_CHEBI_out %>% select(-ID_level1) %>% rename(punc2=punc) , Chebi_compounds_filt2 %>% filter(TYPE!='ChebiName') %>% select(ID_level1, punc2) ) %>% filter(is.na(ID_level1))
data_to_match_CHEBI_out<-data_to_match_CHEBI_out %>% rename(punc=punc2)
data_to_match_CHEBI_synlempunc2$kind<-'synlempunc'
###pubchem matching
###pull in pubchem
Pubchem2<-read_csv('Pubchem.csv')
Pubchem2$text_norm<-str_trim(Pubchem2$text_norm, side='both')
Pubchem2$text_lower_lemma<-str_trim(Pubchem2$text_lower_lemma, side='both')
Pubchem2$text_norm_lower<-tolower(Pubchem2$text_norm)
########### 
Pubchem2$punc<-gsub("[^[:alnum:][:space:]]"," ", Pubchem2$text_lower_lemma)
Pubchem2$punc2<-gsub("[^[:alnum:][:space:]]","", Pubchem2$text_lower_lemma)


## exact match to pub
data_to_match_CHEBI_pubchem<-left_join(data_to_match_CHEBI_out %>% select(-ID_level1) , Pubchem2 %>% select(ID_level1, text_norm_lower) ) %>% filter(!is.na(ID_level1))
data_to_match_CHEBI_out<-left_join(data_to_match_CHEBI_out %>% select(-ID_level1)   , Pubchem2 %>% select(ID_level1, text_norm_lower) ) %>% filter(is.na(ID_level1))
data_to_match_CHEBI_pubchem$kind<-'pubchem'
###exact relaxed  to pub
data_to_match_CHEBI_publem<-left_join(data_to_match_CHEBI_out %>% select(-ID_level1)  , Pubchem2 %>% select(ID_level1, text_lower_lemma) ) %>% filter(!is.na(ID_level1))
data_to_match_CHEBI_out<-left_join(data_to_match_CHEBI_out %>% select(-ID_level1)  , Pubchem2 %>% select(ID_level1, text_lower_lemma) ) %>% filter(is.na(ID_level1))
data_to_match_CHEBI_publem$kind<-'pubchemlem'
###exact relaxed punc to pub
data_to_match_CHEBI_publempunc<-left_join(data_to_match_CHEBI_out %>% select(-ID_level1)  , Pubchem2 %>% select(ID_level1, punc) ) %>% filter(!is.na(ID_level1))
data_to_match_CHEBI_out<-left_join(data_to_match_CHEBI_out %>% select(-ID_level1)  , Pubchem2 %>% select(ID_level1, punc) ) %>% filter(is.na(ID_level1))
data_to_match_CHEBI_publempunc$kind<-'pubchempunc'
###exact relaxed punc to pub
data_to_match_CHEBI_publempunc2<-left_join(data_to_match_CHEBI_out %>% select(-ID_level1)  %>% rename(punc2=punc) , Pubchem2 %>% select(ID_level1, punc2) ) %>% filter(!is.na(ID_level1))
data_to_match_CHEBI_out<-left_join(data_to_match_CHEBI_out %>% select(-ID_level1)  %>% rename(punc2=punc) , Pubchem2 %>% select(ID_level1, punc2) ) %>% filter(is.na(ID_level1))
data_to_match_CHEBI_publempunc2$kind<-'pubchempunc'
#fluorodeoxyglucose
#f-fdg==Fluorodeoxyglucose F18
#write.csv(data_to_match_CHEBI_out, 'leftovers_121021.csv', row.names = F)

#leftovers<-read_csv('/Users/sarahmullin/Desktop/Menagerie/Chemical_normalization/leftovers_121021.csv')

###fuzzymatching in file chemical_normalization2.py
#leftovers<-read_csv('/Users/sarahmullin/Desktop/Menagerie/Chemical_normalization/leftovers_matched.csv')
exact_match_dict_fuzzy<-read_csv('/Users/sarahmullin/Desktop/Menagerie/Chemical_normalization/dictionary_match_chebi.csv')
exact_match_dict_fuzzy<-left_join(exact_match_dict_fuzzy %>% rename(text_lower_lemma=possible), Chebi_compounds_filt2)
exact_match_dict_fuzzy<-exact_match_dict_fuzzy %>% arrange(ID_level1, leftover,desc(score)) %>% distinct(ID_level1, leftover,text_lower_lemma, .keep_all=T) %>% rename(possible=text_lower_lemma)

exact_match_dict_fuzzy2= exact_match_dict_fuzzy %>% 
  group_by(possible) %>% 
  mutate(duplicate.flag = n() > 1) 

### fuzzy match Chebi

data_to_match_CHEBI_fuzzycheb<-left_join(data_to_match_CHEBI_out %>% select(-ID_level1) %>% rename(leftover=term_lower_lemma) , exact_match_dict_fuzzy %>% select(ID_level1, leftover, score, possible) ) %>% filter(!is.na(ID_level1))
data_to_match_CHEBI_out<-left_join(data_to_match_CHEBI_out %>% select(-ID_level1) %>% rename(leftover=term_lower_lemma) , exact_match_dict_fuzzy %>% select(ID_level1, leftover) ) %>% filter(is.na(ID_level1))
data_to_match_CHEBI_fuzzycheb$kind<-'fuzzychebi1'
###fuzzy match pubchem
exact_match_dict_fuzzypub<-read_csv('/Users/sarahmullin/Desktop/Menagerie/Chemical_normalization/dictionary_match_pubchem.csv')
exact_match_dict_fuzzypub<-left_join(exact_match_dict_fuzzypub %>% rename(text_lower_lemma=possible), Pubchem2)
exact_match_dict_fuzzypub<-exact_match_dict_fuzzypub %>% arrange(ID_level1, leftover,desc(score)) %>% distinct(ID_level1, leftover,text_lower_lemma,.keep_all=T) %>% rename(possible=text_lower_lemma)

data_to_match_CHEBI_fuzzypub<-left_join(data_to_match_CHEBI_out %>% select(-ID_level1)  , exact_match_dict_fuzzypub %>% select(ID_level1, leftover, score, possible) ) %>% filter(!is.na(ID_level1))
data_to_match_CHEBI_out<-left_join(data_to_match_CHEBI_out %>% select(-ID_level1) , exact_match_dict_fuzzypub %>% select(ID_level1, leftover) ) %>% filter(is.na(ID_level1))
data_to_match_CHEBI_fuzzypub$kind<-'fuzzypub'
#### now try it with something less matchy
exact_match_dict_fuzzychem2<-read_csv('/Users/sarahmullin/Desktop/Menagerie/Chemical_normalization/dictionary_match_chebi2.csv')

exact_match_dict_fuzzychem2<-left_join(exact_match_dict_fuzzychem2 %>% rename(text_lower_lemma=possible), Chebi_compounds_filt2)
exact_match_dict_fuzzychem2<-exact_match_dict_fuzzychem2 %>% arrange(ID_level1, leftover,desc(score)) 
a<-exact_match_dict_fuzzychem2 %>% filter(leftover=='omega-6'| leftover=='omega-3') %>% distinct(leftover, .keep_all=T)
a$ID_level1[a$leftover=='omega-3']<-25681
a$text_lower_lemma[a$leftover=='omega-3']<-'omega-3 fatty acid'
a$ID_level1[a$leftover=='omega-6']<-36009
a$text_lower_lemma[a$leftover=='omega-6']<-'omega-6 fatty acid'
exact_match_dict_fuzzychem2<-full_join(exact_match_dict_fuzzychem2,a)
exact_match_dict_fuzzychem2<-exact_match_dict_fuzzychem2 %>% distinct(leftover, text_lower_lemma, score, .keep_all=T) %>% rename(possible=text_lower_lemma)
exact_match_dict_fuzzychem22 <- exact_match_dict_fuzzychem2 %>%                                      # Top N highest values by group
  arrange(desc(score)) %>% 
  group_by(leftover) %>%
  slice(1:10)

data_to_match_CHEBI_fuzzychem2<-left_join(data_to_match_CHEBI_out %>% select(-ID_level1)  , exact_match_dict_fuzzychem22 %>% select(ID_level1, leftover, score, possible) ) %>% filter(!is.na(ID_level1))
data_to_match_CHEBI_out<-left_join(data_to_match_CHEBI_out %>% select(-ID_level1) , exact_match_dict_fuzzychem22 %>% select(ID_level1, leftover) ) %>% filter(is.na(ID_level1))
data_to_match_CHEBI_fuzzychem2$kind<-'fuzzychebi3'
data_to_match_CHEBI_fuzzychem2<-data_to_match_CHEBI_fuzzychem2 %>% arrange(ID_level1, desc(score))

#### combining all options together
number<-c(27213, 989, 251)
complete<-full_join(data_to_match_CHEBI_exact %>% select(pmid,term, kind, ID_level1, term_norm) %>% distinct(),data_to_match_CHEBI_exactlem %>% select(pmid,term, kind, ID_level1)  %>% distinct())
complete<-full_join(complete, data_to_match_CHEBI_exactlempunc %>% select(pmid,term, kind, ID_level1, term_norm)  %>% distinct() )
complete<-full_join(complete, data_to_match_CHEBI_exactlempunc2 %>% select(pmid,term, kind, ID_level1, term_norm)  %>% distinct() )
complete<-full_join(complete, data_to_match_CHEBI_syn %>% select(pmid,term, kind, ID_level1, term_norm)  %>% distinct() )
complete<-full_join(complete, data_to_match_CHEBI_synlem %>% select(pmid,term, kind, ID_level1)  %>% distinct())
complete<-full_join(complete, data_to_match_CHEBI_synlempunc %>% select(pmid,term, kind, ID_level1, term_norm) %>% distinct() )
complete<-full_join(complete, data_to_match_CHEBI_synlempunc2 %>% select(pmid,term, kind, ID_level1, term_norm) %>% distinct() )
complete<-full_join(complete, data_to_match_CHEBI_pubchem %>% select(pmid,term, kind, ID_level1, term_norm) %>% distinct() )
complete<-full_join(complete, data_to_match_CHEBI_publem %>% select(pmid,term, kind, ID_level1, term_norm) %>% distinct() )
complete<-full_join(complete, data_to_match_CHEBI_publempunc %>% select(pmid,term, kind, ID_level1, term_norm) %>% distinct() )
complete<-full_join(complete, data_to_match_CHEBI_publempunc2 %>% select(pmid,term, kind, ID_level1, term_norm) %>% distinct() )
complete<-full_join(complete, data_to_match_CHEBI_fuzzycheb %>% select(pmid,term, kind, ID_level1, leftover, score, possible, term_norm) %>% distinct() )
complete<-full_join(complete, data_to_match_CHEBI_fuzzypub %>% select(pmid,term, kind, ID_level1, leftover, score, possible, term_norm)  %>% distinct())
complete<-full_join(complete, data_to_match_CHEBI_fuzzychem2 %>% select(pmid,term, kind, ID_level1, leftover, score, possible, term_norm) %>% distinct() )
complete<-complete %>% arrange(pmid, term)
library(campfin)
complete$term_low<-tolower(complete$term)
complete2= complete %>% 
  group_by(pmid, term_low) %>% 
  mutate(duplicate.flag = n() > 1)

complete_only1<-complete2 %>% filter(duplicate.flag==FALSE) %>% distinct(pmid, kind, ID_level1, term_low, .keep_all=T)

cantmatch<-complete_only1 %>% filter(score<83)
cantmatch<-full_join(cantmatch, data_to_match_CHEBI_out %>% select(pmid,term, ID_level1) %>% distinct())
complete_only1<-complete_only1%>% filter(score>=83|  is.na(score))

#### grab some more so that the samples are even
complete_only1_check<-complete_only1 %>% filter(!(kind=='exact'|kind=='synlem'|kind=='syn'|kind=='exactlem'))
complete_only1_check<-left_join(complete_only1_check,Chebi_compounds_filt_ID_name %>% distinct(ID_level1, .keep_all=T))
complete_only1_check$possible<-with(complete_only1_check, ifelse(is.na(possible), NAME, possible))
#complete_only1_check<-complete_only1 %>% filter((kind=='fuzzychebi1'|kind=='fuzzychebi3'|kind=='fuzzypub'|kind=='exactlem'|kind=='exactlempunc'))

titles<-read_csv('/Users/sarahmullin/Desktop/Menagerie/Chemical_normalization/chemical_titles_12621.csv')
complete_only1_check<-left_join(complete_only1_check, titles)
complete_only1_check<-complete_only1_check %>% distinct(pmid, term, ID_level1, leftover, score, possible, title_proc, .keep_all=T)
complete_only1_check$Sentence2<-complete_only1_check$title_proc
#write.csv(complete_only1_check, 'complete_only1_check.csv', row.names = F)

complete_only1_final<-complete_only1 %>% filter((kind=='exact'|kind=='synlem'|kind=='syn'|kind=='exactlem'))
complete_only1_final2<-complete_only1 %>% filter((kind=='exact'|kind=='synlem'|kind=='syn'|kind=='exactlem'|kind=='synlempunc'|kind=='exactlempunc'))
