speaqeasy_data <- file.path(tempdir(), "rse_speaqeasy.RData")

download.file("https://github.com/LieberInstitute/SPEAQeasy-example/blob/master/rse_speaqeasy.RData?raw=true", speaqeasy_data, mode = "wb")

library("SummarizedExperiment")

load(speaqeasy_data, verbose = TRUE)


#sort.by en toptable()
#para no usar el match (se puede pero es mas complicado)

##¿Hay diferencias en totalAssignedGene o mitoRate entre los grupos de diagnosis (PrimaryDx)?
with(colData(rse_gene), tapply(mitoRate))

#ejerciciooooooooooooooooooooooooo

table(rse_gene$PrimaryDx)
##eliminar diagnosis "other" porque no tiene info
rse$PrimaryDx <- droplevels(rse_gene$PrimaryDx)
table(rse_gene$PrimaryDx)

#explorar diferencias enter grupos de diagnostico
with(colData(rse_gene), tapply(totalAssignedGene, PrimaryDx,summary))
with(colData(rse_gene), tapply(mitoRate, PrimaryDx, summary))


#aver
with(colData)

##ggplot no sabe usar objetos DataFrame
#convertir a dataframe para que pueda manejarlo
ggplot(
  as.data.frame(colData(rse_gene)),
  aes(y = mitoRate, group = BrainRegion, x = BrainRegion)
) +
  geom_bowplot() +
  theme_bw(base:size = 20)+
  xlab("Brain Region")

#hacer dataframe
df <- data.frame(
  expression = assay(rse_gene)[i,]

)

#scater pa visualizar dato de single cell


#despues de usar droplevels!!!!
#
mod <- with(
    colData(rse_gene,),
    model.matrix(~ PrimaryDx + totalAssignedGene + mitoRate + RNA_rate + BrainRegion + Sex + AgeDeath)
)

#vamo a visualizarlo
library(iSEE)
iSEE::iSEE(rse_gene)
colData(rse_gene)[,c("PrimaryDx","totalAssignedGene","rRNA_rate")]
