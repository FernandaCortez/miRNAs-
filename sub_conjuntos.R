library(TCGAbiolinks)
library(dplyr)
library(readr)
#Cargar Matriz 
RNA <- read.table(file = 'RNAseqnormalized.tsv', sep = '\t', header = TRUE)
#cargar lista con subtipos 
subtipos = read.table(file = 'subtype.tsv', sep = '\t', header = TRUE)
#cargar miRNAs 
mirnas = read.table(file = 'miRNAseqNormi.tsv', sep = '\t', header = TRUE)
#generar un objeto con el nombre de todas las columnas. 
sample_id_lista <- data.frame (subtipos[,1])
matriz_exp_id <-data.frame( colnames(RNA))
#cambiar el nombre 
colnames(matriz_exp_id)[1] <- "sample_id_matriz"
colnames(sample_id_lista)[1] <- "sample_id_lista"
#generar una lista con los subtipos de mi lista
tipo_tumor <- data.frame (subtipos[,4])
colnames(tipo_tumor)[1] <- "tipo_tumor"
unidas <- cbind(matriz_exp_id, sample_id_lista, tipo_tumor)


lumB<- subset(unidas, tipo_tumor== "LumB")
lumA <- subset(unidas,tipo_tumor=="LumA")
basal <- subset(unidas, tipo_tumor == "Basal")
her2 <- subset(unidas, tipo_tumor == "Her2")
normal <- subset(unidas, tipo_tumor=="Normal")
metastatic <-subset(unidas, tipo_tumor == "Metastatic")

length(lumB$tipo_tumor)
#140
length(lumA$tipo_tumor)
# 416
length(her2$tipo_tumor)
#46 
length(basal$tipo_tumor)
#128
length(normal$tipo_tumor)
#75
# length(metastatic$tipo_tumor)
# #0

#140+ 416 + 46 + 128 + 75 = 805 

# her2 <-RNA[26:90]
# 
# Her2 <- names(her2$sample_id_lista)[which(colnames(RNA) == "Her2")]


her2_prueba <- RNA[colnames(RNA) %in% her2$sample_id_matriz ]
#17077    46
basal_prueba <- RNA[colnames(RNA) %in% basal$sample_id_matriz ]
#17077   128
lumA_prueba <- RNA[colnames(RNA) %in% lumA$sample_id_matriz ]
# 17077   416
lumB_prueba <- RNA[colnames(RNA) %in% lumB$sample_id_matriz ]
#17077   140
normal_prueba <- RNA[colnames(RNA) %in% normal$sample_id_matriz ]
# 17077    75
dim(normal_prueba)


#######mirnas 

id_mirnas <- data.frame(colnames(mirnas))
colnames(id_mirnas)[1] <- "sample_id_mirnas"
unidas <- cbind(unidas,id_mirnas)

mirnas
# 604 805
prueba_lumA_mir <- mirnas[colnames(mirnas) %in% lumA$sample_id_matriz ]
# 604 416
prueba_lumB_mir <- mirnas[colnames(mirnas) %in% lumB$sample_id_matriz ]
# 604 140
prueba_basal_mir <- mirnas[colnames(mirnas) %in% basal$sample_id_matriz ]
# 604 128
prueba_her2_mir <- mirnas[colnames(mirnas) %in% her2$sample_id_matriz ]
#604  46
prueba_normal_mir <- mirnas[colnames(mirnas) %in% normal$sample_id_matriz ]
#  604  75
# 416+ 140+ 128 + 46+ 75 =805

#her2
her2.1 <- data.frame (colnames(prueba_her2_mir))
her2.2 <- data.frame (colnames(her2_prueba))
her2_table <- cbind(her2.1, her2.2)
#aparecen en el mismo orden 

identical(colnames(prueba_her2_mir), colnames(her2_prueba))
#TRUE 

exp_table_Her2 <- rbind(her2_prueba, prueba_her2_mir)
write.table(exp_table_Her2, file='exp_her2.tsv', quote=FALSE, sep='\t', col.names = T)

#===========================================================================================#
#BASAL 
basal1 <- data.frame (colnames(prueba_basal_mir))
basal2 <- data.frame (colnames(basal_prueba))

basal_table <- cbind(basal1, basal2)
#aparecen en el mismo orden 

identical(colnames(prueba_basal_mir), colnames(basal_prueba))
#FALSE  


# exp_table_Her2 <- rbind(her2_prueba, prueba_her2_mir)
# write.table(exp_table_Her2, file='exp_her2.tsv', quote=FALSE, sep='\t', col.names = T)



#=======================================================
#LUMA

luma1 <- data.frame (colnames(prueba_lumA_mir))
luma2 <- data.frame (colnames(lumA_prueba))
luma_table <- cbind(luma1, luma2)
#aparecen en el mismo orden 

identical(colnames(luma1), colnames(luma2))
#FALSE 

# exp_table_Her2 <- rbind(her2_prueba, prueba_her2_mir)
# write.table(exp_table_Her2, file='exp_her2.tsv', quote=FALSE, sep='\t', col.names = T)

#=========================================================
#LUMB

lumb1 <- data.frame (colnames(prueba_lumB_mir))
lumb2 <- data.frame (colnames(lumB_prueba))
lumb_table<- cbind(lumb1, lumb2)
#aparecen en el mismo orden 

identical(colnames(lumb1), colnames(lumb2))
#TRUE 

# exp_table_Her2 <- rbind(her2_prueba, prueba_her2_mir)
# write.table(exp_table_Her2, file='exp_her2.tsv', quote=FALSE, sep='\t', col.names = T)



#=========================================================
# normal


normal1 <- data.frame (colnames(prueba_normal_mir))
normal1 <-normal1[, ncol(normal1):1]
normal2 <- data.frame (colnames(normal_prueba))
normal2[, ncol(normal2):1]
normal_table<- cbind(normal1, normal2)
#aparecen en el mismo orden 
identical(colnames(normal1), colnames(normal2))
#TRUE 



# exp_table_Her2 <- rbind(her2_prueba, prueba_her2_mir)
# write.table(exp_table_Her2, file='exp_her2.tsv', quote=FALSE, sep='\t', col.names = T)

#Verificar si el orden de mi lista coincide con mi matriz 
# write.table(unidas, file='columnas_lista.tsv', quote=FALSE, sep='\t', col.names = T)
# unidas = read.table(file = 'columnas_lista.tsv', sep = '\t', header = TRUE)
# identical(rownames(unidas$sample_id_matriz), rownames(unidas$sample_id_lista))
# #TRUE 
# # sample_id_lista <- data.frame (subtipos[,1])
# # matriz_exp_id <-data.frame( colnames(RNA))
# # #cambiar el nombre 
# # colnames(matriz_exp_id)[1] <- "sample_id_matriz"
# # colnames(sample_id_lista)[1] <- "sample_id_lista"
# # #generar una lista con los subtipos de mi lista
# # tipo_tumor <- data.frame (subtipos[,4])
# # colnames(tipo_tumor)[1] <- "tipo_tumor"
# # unidas <- cbind(matriz_exp_id, sample_id_lista, tipo_tumor)
# # id_mirnas <- data.frame(colnames(mirnas))
# # colnames(id_mirnas)[1] <- "sample_id_mirnas"
# # unidas <- cbind(unidas,id_mirnas)
# #¿ son identicos los rownames de mirnas con los de genes?
# identical(rownames(unidas$sample_id_lista), rownames(unidas$sample_id_mirnas))
# #TRUE
# write.table(unidas, file='columnas_lista.tsv', quote=FALSE, sep='\t', col.names = T)


#Verificar si el orden de mi lista coincide con mi matriz 
write.table(unidas, file='columnas_lista.tsv', quote=FALSE, sep='\t', col.names = T)
unidas = read.table(file = 'columnas_lista.tsv', sep = '\t', header = TRUE)
identical(rownames(unidas$sample_id_matriz), rownames(unidas$sample_id_lista))
#TRUE 
# sample_id_lista <- data.frame (subtipos[,1])
# matriz_exp_id <-data.frame( colnames(RNA))
# #cambiar el nombre 
# colnames(matriz_exp_id)[1] <- "sample_id_matriz"
# colnames(sample_id_lista)[1] <- "sample_id_lista"
# #generar una lista con los subtipos de mi lista
# tipo_tumor <- data.frame (subtipos[,4])
# colnames(tipo_tumor)[1] <- "tipo_tumor"
# unidas <- cbind(matriz_exp_id, sample_id_lista, tipo_tumor)
# id_mirnas <- data.frame(colnames(mirnas))
# colnames(id_mirnas)[1] <- "sample_id_mirnas"
# unidas <- cbind(unidas,id_mirnas)
#¿ son identicos los rownames de mirnas con los de genes?
identical(rownames(unidas$sample_id_lista), rownames(unidas$sample_id_mirnas))
#TRUE
write.table(unidas, file='columnas_lista.tsv', quote=FALSE, sep='\t', col.names = T)


