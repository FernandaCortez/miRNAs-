###############load package#################
require(readr)
require(dplyr)
require(igraph)
require(tidygraph)
require(ggraph)
require(ggplot2)

setwd("~/Documentos/Proyecto/Selección de genes/Búsqueda_fármacos/Conjuntos")

# Load data
#repeat for each subtype
basal_data <- read_csv("conjuntos_basal.csv")

# basal_data <- conjuntos_basal[, -4]

#Generate bipartite networks
gene_drug <- cbind(basal_data[,1], basal_data[,2])
mir_gene <- cbind(basal_data[,3], basal_data[,1])

# Remove duplicate values
gene_drug <- unique(gene_drug)
mir_gene<-unique(mir_gene)

#Create an igraph object
g1 <- graph_from_data_frame(gene_drug)
g2 <- graph_from_data_frame(mir_gene)

# Generate a tripartite network
g3 = g2 +  g1
as_tbl_graph(g3)

#Create dictionary for the network.  
g3 <- g3 %>% as_tbl_graph %>% mutate(tipo = case_when(name %in% unique(gene_drug$Genes)~ "Genes",
                                                      name %in% unique(gene_drug$Farmaco)~ "Farmacos",
                                                      name %in% unique(mir_gene$miRNAs)~ "miRNAs"))


#Count the number of elements in my network
g3 %>%  as_tibble %>% group_by(tipo) %>% count()



# filter first and second neighbors
g3 <- g3 %>% mutate(first_neighbors = local_size(order = 1, mode = "all", mindist = 1))
g3 <- g3 %>% mutate(second_neighbors = local_size(order = 2, mode = "all", mindist = 2))
g3 %>% filter(tipo== "miRNAs")

# Generate scatter plot
g3 %>% filter(tipo== "miRNAs") %>% 
  as_tibble() %>% 
  ggplot() + geom_point(aes(first_neighbors, second_neighbors)) 

g3 %>% filter(tipo== "miRNAs") %>% 
  as_tibble() %>% 
  ggplot() + geom_point(aes(first_neighbors, second_neighbors))+
  geom_label(aes(vecinos_uno, vecinos_dos, label= name))

#Filter miRNAs with higher number of second neighbors
miRNAS <- g3 %>% filter(tipo== "miRNAs") %>% 
  filter(vecinos_dos >= 40) %>% 
  pull(name)

miRNAS
# [1] "hsa-mir-3942" "hsa-mir-141"  "hsa-mir-99b"  "hsa-mir-6824" "hsa-mir-561"
# "hsa-mir-4715" [7] "hsa-mir-133b"
names(miRNAS) <- miRNAS

# Get number of first and second neighbors
neighborhood(g3, order = 2, nodes = miRNAS, mindist = 2)
lapply(miRNAS, FUN = function(i){neighborhood(g3, order = 2
                                                  , nodes = i, mindist = 2)})

