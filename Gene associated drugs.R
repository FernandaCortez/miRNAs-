##Identify drugs associated with genes in a network inferred by mutual information##

#Load Package
require(readr)

#repeat by subtype
#load mutual information  networks

genes_mi <- read.table(file = 'miRNAs_gene_basal_.txt', sep = '\t', header = T)
genes_mi <- genes_mi[-3]
colnames(genes_mi)[1] <- "Gene"
colnames(genes_mi)[2] <- "miRNA"

#load genes GO
genes_GO <- read.table(file = 'genes_basal.txt', sep = '\t', header = F)
colnames(genes_GO) [1] <- "Gene"

#load pharmacomiR drug list 
new_drugs_list <- read_delim("new_drugs_list.csv", 
                             delim = "\t", escape_double = FALSE, 
                             trim_ws = TRUE)
                             
 new_drugs_list <- data.frame(Drug= new_drugs_list$Drug,
                                 Gene = new_drugs_list$Gene)
                                 
new_drugs_list<- unique(new_drugs_list)

#Filter genes from mutual information with genes from curated lists with gene ontology

new_genes <- merge(x=genes_GO , y= genes_mi , by = "Gene")
new_genes <- unique(new_genes)

#filter drug-associated genes 

basal_set <- merge(x=new_genes , y= new_drugs_list, by="Gene")  
basal_set <- unique(basal_set)

write.table(conjuntos_basal, file = "conjuntos_basal.csv", sep = ",", 
            row.names = F, col.names = TRUE)


####repeat by subtype#####
