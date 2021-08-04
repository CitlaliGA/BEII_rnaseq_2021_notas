```{r download_SRP045638}
library("recount3")
human_projects <- available_projects()
rse_gene_SRP211963 <- create_rse(
  subset(
    human_projects,
    project == "SRP211963" & project_type == "data_sources"
  )
)
assay(rse_gene_SRP211963, "counts") <- compute_read_counts(rse_gene_SRP211963)
```

Una vez descargados y con los números de lecturas podemos usar `expand_sra_attributes()`. Sin embargo, tenemos un problema con estos datos.

```{r describe_issue}
rse_gene_SRP211963$sra.sample_attributes[1:3]
```

Vamos a intentar resolverlo eliminando información que está presente solo en ciertas muestras.

```{r solve_issue}
rse_gene_SRP211963$sra.sample_attributes <- gsub("dev_stage;;Fetal\\|", "", rse_gene_SRP211963$sra.sample_attributes)
rse_gene_SRP211963$sra.sample_attributes[1:3]
```

Ahora si podemos continuar con el mismo código de ayer.

```{r attributes}
rse_gene_SRP211963 <- expand_sra_attributes(rse_gene_SRP211963)
colData(rse_gene_SRP211963)[
  ,
  grepl("^sra_attribute", colnames(colData(rse_gene_SRP211963)))
]
```

Como ahora si vamos a usar esta información para un modelo estadístico, será importante que tengamos en el formato correcto de R a la información que vamos a usar.

```{r re_cast}
## Pasar de character a nuemric o factor
rse_gene_SRP211963$sra_attribute.age <- as.numeric(rse_gene_SRP211963$sra_attribute.age)
rse_gene_SRP211963$sra_attribute.disease <- factor(rse_gene_SRP211963$sra_attribute.disease)
rse_gene_SRP211963$sra_attribute.RIN <- as.numeric(rse_gene_SRP211963$sra_attribute.RIN)
rse_gene_SRP211963$sra_attribute.sex <- factor(rse_gene_SRP211963$sra_attribute.sex)
## Resumen de las variables de interés
summary(as.data.frame(colData(rse_gene_SRP211963)[
  ,
  grepl("^sra_attribute.[age|disease|RIN|sex]", colnames(colData(rse_gene_SRP211963)))
]))
```

<a href="https://r4ds.had.co.nz/explore-intro.html"><img src="images/data-science-explore.png" width="800px" /></a>

  Ahora crearemos un par de variables para que las podamos usar en nuestro análisis.

```{r new_variables}
## Encontraremos diferencias entre muestra prenatalas vs postnatales
rse_gene_SRP211963$prenatal <- factor(ifelse(rse_gene_SRP211963$sra_attribute.age < 0, "prenatal", "postnatal"))
table(rse_gene_SRP211963$prenatal)
## http://research.libd.org/recount3-docs/docs/quality-check-fields.html
rse_gene_SRP211963$assigned_gene_prop <- rse_gene_SRP211963$recount_qc.gene_fc_count_all.assigned / rse_gene_SRP211963$recount_qc.gene_fc_count_all.total
summary(rse_gene_SRP211963$assigned_gene_prop)
with(colData(rse_gene_SRP211963), plot(assigned_gene_prop, sra_attribute.RIN))
## Hm... veamos si hay una diferencia entre los grupos
with(colData(rse_gene_SRP211963), tapply(assigned_gene_prop, prenatal, summary))
```

A continuación podemos eliminar algunas muestras que consideremos de baja calidad y genes con niveles de expresión muy bajos.


```{r filter_rse}
## Guardemos nuestro objeto entero por si luego cambiamos de opinión
rse_gene_SRP211963_unfiltered <- rse_gene_SRP211963
## Eliminemos a muestras malas
hist(rse_gene_SRP211963$assigned_gene_prop)
table(rse_gene_SRP211963$assigned_gene_prop < 0.77)
rse_gene_SRP211963 <- rse_gene_SRP211963[, rse_gene_SRP211963$assigned_gene_prop > 0.77]
## Calculemos los niveles medios de expresión de los genes en nuestras
## muestras.
## Ojo: en un análisis real probablemente haríamos esto con los RPKMs o CPMs
## en vez de las cuentas.
gene_means <- rowMeans(assay(rse_gene_SRP211963, "counts"))
summary(gene_means)
## Eliminamos genes
rse_gene_SRP211963 <- rse_gene_SRP211963[gene_means > 0.1, ]
## Dimensiones finales
dim(rse_gene_SRP211963)
## Porcentaje de genes que retuvimos
round(nrow(rse_gene_SRP211963) / nrow(rse_gene_SRP211963_unfiltered) * 100, 2)
```
Ahora ya estamos listos para continuar con el análisis de expresión diferencial, bueno, casi.


## Normalización de datos

* Lean _A hypothetical scenario_ en uno de los artículos sobre `edgeR` https://genomebiology.biomedcentral.com/articles/10.1186/gb-2010-11-3-r25#Sec2 para entender un poco sobre el concepto de _composition bias_.
* Sigue siendo relevante con datos de scRNA-seq como pueden ver en http://bioconductor.org/books/release/OSCA/normalization.html#normalization-by-deconvolution. Ahí descubren una serie de pasos para usar métodos desarrollados para bulk RNA-seq y como se pueden usar en scRNA-seq.

###############VAMOS AKIIIIIIIIII

```{r normalize}
library("edgeR") # BiocManager::install("edgeR", update = FALSE)
dge <- DGEList(
  counts = assay(rse_gene_SRP211963, "counts"),
  genes = rowData(rse_gene_SRP211963)
)
dge <- calcNormFactors(dge)
```

## Expresión diferencial

Primero que nada, definamos nuestro modelo estadístico. Típicamente, exploraríamos más los datos para revisar que no haya otros problemas con las muestras y para explorar la relación entre nuestras variables.

```{r explore_gene_prop_by_age}
library("ggplot2")
ggplot(as.data.frame(colData(rse_gene_SRP211963)), aes(y = assigned_gene_prop, x = prenatal)) +
  geom_boxplot() +
  theme_bw(base_size = 20) +
  ylab("Assigned Gene Prop") +
  xlab("Age Group")
```
Por ejemplo, usando el paquete de [`variancePartition`](https://bioconductor.org/packages/variancePartition) y [`scater`](https://bioconductor.org/packages/scater) entre otros tal como exploramos en el siguiente video del club de R de LIBD (_[notes in English](https://docs.google.com/document/d/1hil3zwPN6BW6HlwldLbM1FdlLIBWKFNXRUqJJEK_-eY/edit)_)/

  <iframe width="560" height="315" src="https://www.youtube.com/embed/OdNU5LUOHng" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>


  Por ahora continuaremos con el siguiente modelo estadístico.

```{r statiscal_model}
mod <- model.matrix(~ prenatal + sra_attribute.RIN + sra_attribute.sex + assigned_gene_prop,
                    data = colData(rse_gene_SRP211963)
)
colnames(mod)
```


Ya teniendo el modelo estadístico, podemos usar `limma` para realizar el análisis de expresión diferencial como tal.

```{r run_limma}
library("limma")
vGene <- voom(dge, mod, plot = TRUE)
eb_results <- eBayes(lmFit(vGene))
de_results <- topTable(
  eb_results,
  #tiene que ser el coeficiente que se pasó anteriormetnte( en este caso el que sigue despues de ~ el prenatal,
  #no puede ser nombre por lo tanto estar pendiente!!!!!!!!!!!!!!!!!IMPORTANTE)
  coef = 2,
  number = nrow(rse_gene_SRP211963),
  sort.by = "none"
)
dim(de_results)
head(de_results)

#####
identical(sign(de_results$logFC), sign(de_results$t))

## Genes diferencialmente expresados entre pre y post natal con FDR < 5%
table(de_results$adj.P.Val < 0.05)
## Visualicemos los resultados estadísticos
plotMA(eb_results, coef = 2)
volcanoplot(eb_results, coef = 2, highlight = 3, names = de_results$gene_name)
de_results[de_results$gene_name %in% c("ZSCAN2", "VASH2", "KIAA0922"), ]
```

* https://www.genecards.org/cgi-bin/carddisp.pl?gene=ZSCAN2
* https://www.genecards.org/cgi-bin/carddisp.pl?gene=VASH2
* https://www.genecards.org/cgi-bin/carddisp.pl?gene=KIAA0922

## Visualizando genes DE


De `vGene$E` podemos extraer los datos normalizados por `limma-voom`. Revisemos los top 50 genes diferencialmente expresados.

```{r pheatmap}
## Extraer valores de los genes de interés
##$E los datos ya normalizads
##de los 50 primeros, ver cuales son y sus valores de expresion normalizados!!!!!!!!!!!!!!!!
exprs_heatmap <- vGene$E[rank(de_results$adj.P.Val) <= 50, ]
## Creemos una tabla con información de las muestras
## y con nombres de columnas más amigables
df <- as.data.frame(colData(rse_gene_SRP211963)[, c("prenatal", "sra_attribute.RIN", "sra_attribute.sex")])
colnames(df) <- c("AgeGroup", "RIN", "Sex")
## Hagamos un heatmap
library("pheatmap")
pheatmap(
  exprs_heatmap,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = FALSE,
  show_colnames = FALSE,
  annotation_col = df
)
```
