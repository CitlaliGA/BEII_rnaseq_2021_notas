## Lets build our first SummarizedExperiment object
library("SummarizedExperiment")
## ?SummarizedExperiment

## De los ejemplos en la ayuda oficial

## Creamos los datos para nuestro objeto de tipo SummarizedExperiment
## para 200 genes a lo largo de 6 muestras
nrows <- 200
ncols <- 6
## Números al azar de cuentas
set.seed(20210223)
counts <- matrix(runif(nrows * ncols, 1, 1e4), nrows)
## Información de nuestros genes
rowRanges <- GRanges(
  rep(c("chr1", "chr2"), c(50, 150)),
  IRanges(floor(runif(200, 1e5, 1e6)), width = 100),
  strand = sample(c("+", "-"), 200, TRUE),
  feature_id = sprintf("ID%03d", 1:200)
)
names(rowRanges) <- paste0("gene_", seq_len(length(rowRanges)))
## Información de nuestras muestras
colData <- DataFrame(
  Treatment = rep(c("ChIP", "Input"), 3),
  row.names = LETTERS[1:6]
)
## Juntamos ahora toda la información en un solo objeto de R
rse <- SummarizedExperiment(
  assays = SimpleList(counts = counts),
  rowRanges = rowRanges,
  colData = colData
)

## Exploremos el objeto resultante
rse


dim(rse)

dimnames(rse)

assayNames(rse)

head(assay(rse))

#objeto G-ranges
rowRanges(rse)

#para ver todos los cromosomas
seqlevels(rse)
#vector comprimido para que sea mas eficiente de memoria

seqnames(rowRanges(rse))
#que tambien seria lo mismo que lo de abajo
#pero mas chido el de arriba
unique(as.vector(seqnames(rowRanges(rse))))

##PA VER LA MEMORIA QUE ABARCA
 pryr::object_size(rse)

# QUE ES LO QUE ESTA PASANDO EN ESTOS DOS COMANDOS
 rse[1:2, ]
 rse[ , c("A","B","C")]

 ## Explora el objeto rse de forma interactiva
 library("iSEE")
 iSEE::iSEE(rse)
###########
 #########
 #
 #
 #
 #
 #
 #
 #
 ## Descarguemos unos datos de spatialLIBD
 sce_layer <- spatialLIBD::fetch_data("sce_layer")
 sce_layer
 library("iSEE")
 iSEE::iSEE(sce_layer)

