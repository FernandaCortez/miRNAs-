#load packages

require(readr)
require(dbplyr)

setwd("~/Documentos/Proyecto/Nuevos datos/Diccionario de fármacos")

#load data 

drug_atc<- read.table(file = 'drug_atc.tsv', sep = '\t', header = F)
drug_names <- read.table(file = 'drug_names.tsv', sep = '\t', header = F)
     
colnames(drug_atc)[1] <- "CID"  
colnames(drug_atc)[2] <- "ATC"
colnames(drug_names)[1] <- "CID"  
colnames(drug_names)[2] <- "DRUG"

#Generate dictionary with ATC code and drug name

drugs <- merge(x=drug_atc , y=drug_names, by="CID" ) 
atc_code <- as.data.frame( substr(drugs$ATC, start = 1, stop = 3))
drugs_table <-cbind(drugs,atc_code)
drugs_table<- drugs_table[-2]
colnames(drugs_table)[3] <- "ATC"
drugs_table<- unique(drugs_table)

#Add ATC Hierarchy information 

ATC_hierarchy <- read.table(file = 'classification.txt', sep = '\t')
colnames(ATC_hierarchy)[1] <- "ATC"
colnames(ATC_hierarchy)[2] <- "Classification"

drug_dictionary <- merge(x= drugs_table ,  y= ATC_hierarchy , by= "ATC")


#Save table
write.table(drug_dictionary, file = "drug_dictionary(complete).csv", sep = "\t", 
            row.names = F, col.names = TRUE)
            





           
