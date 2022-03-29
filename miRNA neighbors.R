#Generation of networks and search for second neighbors associated with a miRNA.

###############load package#################
require(readr)
require(dplyr)
require(igraph)
require(tidygraph)
require(ggraph)
require(ggplot2)

#setwd("~/Documentos/Proyecto/Nuevos datos/miRNAs Segundos vecinos")

######Repeat by subtype#####

# Load data
basal_data <- read_csv("conjuntos_basal.csv")

#Generate bipartite networks
mir_gene <- cbind(basal_data[,2], basal_data[,1])
gene_drug <- cbind(basal_data[,1], basal_data[,3])

# Remove duplicate values
gene_drug <- unique(gene_drug)
mir_gene<-unique(mir_gene)

# Generate a tripartite network
g3 = g1 + g2
as_tbl_graph(g3)

#Save complex networks

write_graph(g1, "miR_gene_basal.graphml", "graphml")
write_graph(g2, "gene_drug_basal.graphml", "graphml")
write_graph(g3, "tripartite_basal.graphml", "graphml")


#Create an igraph object
g1 <- graph_from_data_frame(d = mir_gene, directed = TRUE, vertices = NULL)
g2 <- graph_from_data_frame(d = gene_drug, directed = TRUE, vertices = NULL)


#Create dictionary for the network.  
g3 <- g3 %>% as_tbl_graph %>% mutate(tipo = case_when(name %in% unique(gene_drug$Gene)~ "Gene",
                                                      name %in% unique(gene_drug$Drug)~ "Drug",
                                                      name %in% unique(mir_gene$miRNA)~ "miRNA"))


#Count the number of elements in my network
g3 %>%  as_tibble %>% group_by(tipo) %>% count()

# filter first and second neighbors for each miRNA in the network
g3 <- g3 %>% mutate(first_neighbors = local_size(order = 1, mode = "out", mindist = 1))
g3 <- g3 %>% mutate(second_neighbors = local_size(order = 2, mode = "out", mindist = 2))
g3 %>% filter(tipo== "miRNA")

#add percentiles 
g3 <- g3 %>% mutate(p1 = ntile(first_neighbors, 100))
g3 <- g3 %>% mutate(p2 = ntile(second_neighbors, 100))

#Filter miRNAs with higher number of second neighbors
miRNAS <- g3 %>% filter(tipo== "miRNA",p1 >=75, p2 >=75) %>% 
  filter(second_neighbors >= 20) %>% 
  pull(name)

names(miRNAS) <- miRNAS

#Identify only drugs associated with filtered miRNAs

neighborhood(g3, order = 2, nodes = miRNAS, mindist = 2)
lapply(miRNAS, FUN = function(i){neighborhood(g3, order = 2
                                              , nodes = i, mindist = 2,mode = "out")})


#Generate scatter plots of the data

g3 %>% filter(tipo== "miRNA", p1 >=75, p2 >=75) %>% 
  as_tibble() %>% 
  ggplot() + geom_point(aes(first_neighbors, second_neighbors, color = second_neighbors)) +
  scale_color_gradient(low = "#FFB6C1", high = "#8B5F65")+
  labs(title = "Basal subtype") +
  coord_cartesian(xlim =c(NA, 16), ylim = c(NA, 120))

#with labels 
g3 %>% filter(tipo== "miRNA", p1 >=75, p2 >=75) %>% 
  as_tibble() %>% 
  ggplot() + geom_point(aes(first_neighbors, second_neighbors))+
  geom_label(aes(first_neighbors, second_neighbors, label= name)) +
  labs(title = "Basal subtype")+
  xlim(NA,16)+
  ylim(NA,120)

####repeat by subtype#####

