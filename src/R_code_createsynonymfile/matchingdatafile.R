#######################################################################################
#ChEBI normalization using BERT: making synonym matching datafile from multiple sources
#######################################################################################


#####
library(dplyr)
library(readr)
library(stringr)
library(curl)

setwd('~/CHEBINormalizer/Data_files')

url = "https://ftp.ebi.ac.uk/pub/databases/chebi/Flat_file_tab_delimited/compounds.tsv.gz"
download.file(url, )

Chebi_compounds<-read_delim('CHEBI/compounds.tsv', delim='\t', guess_max=163324)
Names_CHEBI<-read_delim('/Users/sarahmullin/Desktop/Menagerie/Chemical_normalization/names.tsv', delim='\t')
Chebi_compounds_filt<-Chebi_compounds %>% filter(STATUS=='C'|STATUS=='E') %>% filter(STAR==2|STAR==3)
Chebi_compounds_filt$old<-1
###compounds_again
Chebi_compounds2<-read_delim('/Users/sarahmullin/Downloads/compounds.tsv', delim='\t', guess_max=183324 )
Chebi_compounds_filt2<-Chebi_compounds2 %>% filter(STATUS=='C'|STATUS=='E') %>% filter(STAR==2|STAR==3)
Chebi_compounds_filt2$new<-1
Chebi_check_across<-full_join(Chebi_compounds_filt %>% select(ID,old),Chebi_compounds_filt2 %>% select(ID, new))
Chebi_check_across$KIND<-with(Chebi_check_across, ifelse(old==1 & new==1, 1,0))

A<-Chebi_check_across %>% filter(is.na(KIND))
#### parent_id has subsumed other ids-- switch all to parent
Chebi_compounds_filt$ID_level1<-Chebi_compounds_filt$ID
Chebi_compounds_filt$ID_level1<-with(Chebi_compounds_filt, ifelse(PARENT_ID!='null', PARENT_ID, ID_level1))
Chebi_compounds_filt$TYPE<-'ChebiName'

Chebi_compounds_filt_ID<-Chebi_compounds_filt %>% distinct(ID, ID_level1) %>% rename(COMPOUND_ID=ID)
Names_CHEBI<- left_join(Names_CHEBI, Chebi_compounds_filt_ID) %>% filter(!is.na(ID_level1))
Chebi_compounds_filt<-full_join(Chebi_compounds_filt, Names_CHEBI %>% select(COMPOUND_ID,ID_level1, TYPE, SOURCE, NAME)) %>% arrange(ID_level1) 
#%>% distinct(NAME, .keep_all=T)

Chebi_compounds_filt_ID_name<-Chebi_compounds_filt %>% filter(TYPE=='ChebiName') %>% distinct(ID_level1, NAME)
Chebi_compounds_filt_ID_name$ID_level1<-as.integer(Chebi_compounds_filt_ID_name$ID_level1)
##notes: there are null names in this file-- remove
Chebi_compounds_filt_ID_name2<-Chebi_compounds_filt_ID_name %>% filter(NAME!='null')

All_mappings<-readRDS('/Users/sarahmullin/Desktop/Menagerie/All_mappings.rds')
PubChem<-All_mappings %>% filter(TYPE=='Pubchem')
Chebi_compounds_filt<-full_join(Chebi_compounds_filt, Names_CHEBI %>% select(COMPOUND_ID,ID_level1, 'TYPE', 'SOURCE', 'NAME')) %>% arrange(ID_level1) 
PubChem<- left_join(PubChem %>% rename(COMPOUND_ID=ID_level1, NAME=text_lower_lemma) %>% select(COMPOUND_ID, TYPE, Importance, NAME), Chebi_compounds_filt_ID) %>% filter(!is.na(ID_level1))
Chebi_compounds_filt<-full_join(Chebi_compounds_filt, PubChem) %>% arrange(ID_level1) %>% distinct(NAME, .keep_all=T)
#write.csv(PubChem,'PubChem_out2.csv', row.names=F)
#write.csv(Chebi_compounds_filt,'All_mappings3.csv', row.names=F)
#Chebi_compounds_filt333<-read_csv('All_mappings2.csv')
Chebi_compounds_filt2<-read_csv('All_mappings3.csv')
Chebi_compounds_filt2$text_norm<-str_trim(Chebi_compounds_filt2$text_norm, side='both')
Chebi_compounds_filt2$text_lower_lemma<-str_trim(Chebi_compounds_filt2$text_lower_lemma, side='both')
Chebi_compounds_filt2$text_norm_lower<-tolower(Chebi_compounds_filt2$text_norm)
########### 
Chebi_compounds_filt2$punc<-gsub("[^[:alnum:][:space:]]"," ", Chebi_compounds_filt2$text_lower_lemma)
Chebi_compounds_filt2$punc2<-gsub("[^[:alnum:][:space:]]",'', Chebi_compounds_filt2$text_lower_lemma)
Chebi_compounds_filt2$name_lower<-tolower(str_trim(Chebi_compounds_filt2$NAME, side='both'))
