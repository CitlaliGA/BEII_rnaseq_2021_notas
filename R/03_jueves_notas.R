# Modelos estadísticos

* Revisión de regresión lineal https://lcolladotor.github.io/bioc_team_ds/helping-others.html#linear-regression-example
* Con R, usamos mucho la función `model.matrix()` y la sintáxis de fórmula `Y ~ X1 + X2` tal como en el siguiente ejemplo.

```{r model.matrix}
## ?model.matrix
mat <- with(trees, model.matrix(log(Volume) ~ log(Height) + log(Girth)))
mat
colnames(mat)
```

* ¿Cómo interpretamos los nombres de las columnas de `mat`?

  ```{r lm_example}
summary(lm(log(Volume) ~ log(Height) + log(Girth), data = trees))
```

## ExploreModelMatrix

* Es un paquete de Bioconductor que nos ayuda a entender los modelos estadísticos que estamos usando gracias a visualizaciones http://www.bioconductor.org/packages/ExploreModelMatrix/ que está descrito en el siguiente artículo
* Revisaremos los ejemplos en http://www.bioconductor.org/packages/release/bioc/vignettes/ExploreModelMatrix/inst/doc/ExploreModelMatrix.html


### Ejemplo 1

```{r EMM_example1}
## Datos de ejemplo
(sampleData <- data.frame(
  genotype = rep(c("A", "B"), each = 4),
  treatment = rep(c("ctrl", "trt"), 4)
))
## Creemos las imágenes usando ExploreModelMatrix
vd <- ExploreModelMatrix::VisualizeDesign(
  sampleData = sampleData,
  designFormula = ~ genotype + treatment,
  textSizeFitted = 4
)
## Veamos las imágenes
cowplot::plot_grid(plotlist = vd$plotlist)
```

De forma interactiva podemos correr el siguiente código:

  ```{r EMM_example1_interactive, eval = FALSE}
## Usaremos shiny otra ves
app <- ExploreModelMatrix(
  sampleData = sampleData,
  designFormula = ~ genotype + treatment
)
if (interactive()) shiny::runApp(app)
```

### Ejemplo 2

http://bioconductor.org/packages/release/bioc/vignettes/ExploreModelMatrix/inst/doc/ExploreModelMatrix.html#example-2

### Ejemplo 3

http://bioconductor.org/packages/release/bioc/vignettes/ExploreModelMatrix/inst/doc/ExploreModelMatrix.html#example-3

### Ejercicio

* Interpreta `ResponseResistant.Treatmentpre` del ejercicio 2. Puede ser útil tomar un _screenshot_ (captura de pantalla) y anotarla con líneas de colores. Si haces eso, puedes incluir la imagen en tus notas.
* ¿Por qué es clave el `0` al inicio de la fórmula en el ejercicio 3?

  ### Para aprender más

  _A guide to creating design matrices for gene expression experiments_:

  * http://bioconductor.org/packages/release/workflows/vignettes/RNAseq123/inst/doc/designmatrices.html
* https://f1000research.com/articles/9-1444

_“Model matrix not full rank”_

* http://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#model-matrix-not-full-rank

## Datos de SRP045638

Vamos a usar datos de https://www.ncbi.nlm.nih.gov/sra/?term=SRP045638 procesados con `recount3`. Primero hay que descargar los datos con los comandos que vimos ayer.

```{r download_SRP045638}
library("recount3")
human_projects <- available_projects()
rse_gene_SRP045638 <- create_rse(
  subset(
    human_projects,
    project == "SRP045638" & project_type == "data_sources"
  )
)
assay(rse_gene_SRP045638, "counts") <- compute_read_counts(rse_gene_SRP045638)
```

Una vez descargados y con los números de lecturas podemos usar `expand_sra_attributes()`. Sin embargo, tenemos un problema con estos datos.

```{r describe_issue}
rse_gene_SRP045638$sra.sample_attributes[1:3]
```

Vamos a intentar resolverlo eliminando información que está presente solo en ciertas muestras.

```{r solve_issue}
rse_gene_SRP045638$sra.sample_attributes <- gsub("dev_stage;;Fetal\\|", "", rse_gene_SRP045638$sra.sample_attributes)
rse_gene_SRP045638$sra.sample_attributes[1:3]
```

Ahora si podemos continuar con el mismo código de ayer.

```{r attributes}
rse_gene_SRP045638 <- expand_sra_attributes(rse_gene_SRP045638)
colData(rse_gene_SRP045638)[
  ,
  grepl("^sra_attribute", colnames(colData(rse_gene_SRP045638)))
]
```

Como ahora si vamos a usar esta información para un modelo estadístico, será importante que tengamos en el formato correcto de R a la información que vamos a usar.

```{r re_cast}
## Pasar de character a nuemric o factor
rse_gene_SRP045638$sra_attribute.age <- as.numeric(rse_gene_SRP045638$sra_attribute.age)
rse_gene_SRP045638$sra_attribute.disease <- factor(rse_gene_SRP045638$sra_attribute.disease)
rse_gene_SRP045638$sra_attribute.RIN <- as.numeric(rse_gene_SRP045638$sra_attribute.RIN)
rse_gene_SRP045638$sra_attribute.sex <- factor(rse_gene_SRP045638$sra_attribute.sex)
## Resumen de las variables de interés
summary(as.data.frame(colData(rse_gene_SRP045638)[
  ,
  grepl("^sra_attribute.[age|disease|RIN|sex]", colnames(colData(rse_gene_SRP045638)))
]))
```

<a href="https://r4ds.had.co.nz/explore-intro.html"><img src="images/data-science-explore.png" width="800px" /></a>

  Ahora crearemos un par de variables para que las podamos usar en nuestro análisis.

```{r new_variables}
## Encontraremos diferencias entre muestra prenatalas vs postnatales
rse_gene_SRP045638$prenatal <- factor(ifelse(rse_gene_SRP045638$sra_attribute.age < 0, "prenatal", "postnatal"))
table(rse_gene_SRP045638$prenatal)
## http://research.libd.org/recount3-docs/docs/quality-check-fields.html
rse_gene_SRP045638$assigned_gene_prop <- rse_gene_SRP045638$recount_qc.gene_fc_count_all.assigned / rse_gene_SRP045638$recount_qc.gene_fc_count_all.total
summary(rse_gene_SRP045638$assigned_gene_prop)
with(colData(rse_gene_SRP045638), plot(assigned_gene_prop, sra_attribute.RIN))
## Hm... veamos si hay una diferencia entre los grupos
with(colData(rse_gene_SRP045638), tapply(assigned_gene_prop, prenatal, summary))
```

A continuación podemos eliminar algunas muestras que consideremos de baja calidad y genes con niveles de expresión muy bajos.


```{r filter_rse}
## Guardemos nuestro objeto entero por si luego cambiamos de opinión
rse_gene_SRP045638_unfiltered <- rse_gene_SRP045638
## Eliminemos a muestras malas
hist(rse_gene_SRP045638$assigned_gene_prop)
table(rse_gene_SRP045638$assigned_gene_prop < 0.3)
rse_gene_SRP045638 <- rse_gene_SRP045638[, rse_gene_SRP045638$assigned_gene_prop > 0.3]
## Calculemos los niveles medios de expresión de los genes en nuestras
## muestras.
## Ojo: en un análisis real probablemente haríamos esto con los RPKMs o CPMs
## en vez de las cuentas.
gene_means <- rowMeans(assay(rse_gene_SRP045638, "counts"))
summary(gene_means)
## Eliminamos genes
rse_gene_SRP045638 <- rse_gene_SRP045638[gene_means > 0.1, ]
## Dimensiones finales
dim(rse_gene_SRP045638)
## Porcentaje de genes que retuvimos
round(nrow(rse_gene_SRP045638) / nrow(rse_gene_SRP045638_unfiltered) * 100, 2)
```
Ahora ya estamos listos para continuar con el análisis de expresión diferencial, bueno, casi.


## Normalización de datos

* Lean _A hypothetical scenario_ en uno de los artículos sobre `edgeR` https://genomebiology.biomedcentral.com/articles/10.1186/gb-2010-11-3-r25#Sec2 para entender un poco sobre el concepto de _composition bias_.
* Sigue siendo relevante con datos de scRNA-seq como pueden ver en http://bioconductor.org/books/release/OSCA/normalization.html#normalization-by-deconvolution. Ahí descubren una serie de pasos para usar métodos desarrollados para bulk RNA-seq y como se pueden usar en scRNA-seq.

###############VAMOS AKIIIIIIIIII

```{r normalize}
library("edgeR") # BiocManager::install("edgeR", update = FALSE)
dge <- DGEList(
  counts = assay(rse_gene_SRP045638, "counts"),
  genes = rowData(rse_gene_SRP045638)
)
dge <- calcNormFactors(dge)
```

## Expresión diferencial

Primero que nada, definamos nuestro modelo estadístico. Típicamente, exploraríamos más los datos para revisar que no haya otros problemas con las muestras y para explorar la relación entre nuestras variables.

```{r explore_gene_prop_by_age}
library("ggplot2")
ggplot(as.data.frame(colData(rse_gene_SRP045638)), aes(y = assigned_gene_prop, x = prenatal)) +
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
                    data = colData(rse_gene_SRP045638)
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
  number = nrow(rse_gene_SRP045638),
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
df <- as.data.frame(colData(rse_gene_SRP045638)[, c("prenatal", "sra_attribute.RIN", "sra_attribute.sex")])
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
####en grafica de -log10(pvalue)~ log2 fold ..valores mas extremos de pvalue ::
###
#####genecards



Los resultados que tenemos no son tan sorprendentes porque hay una diferencia enorme en los perfiles de expresión en el DLPFC entre muestra pre y post-natales. Eso lo podemos ver con MDS (multidimensional scaling) tal como describen en [este workflow](http://bioconductor.org/packages/release/workflows/vignettes/RNAseq123/inst/doc/limmaWorkflow.html#unsupervised-clustering-of-samples).

                                                                                                                                                                                                                                                            ```{r plot_mds}
                                                                                                                                                                                                                                                            ## Para colores
                                                                                                                                                                                                                                                            library("RColorBrewer")
                                                                                                                                                                                                                                                            ## Conviertiendo los grupos de edad a colores
                                                                                                                                                                                                                                                            col.group <- df$AgeGroup
                                                                                                                                                                                                                                                            levels(col.group) <- brewer.pal(nlevels(col.group), "Set1")
                                                                                                                                                                                                                                                            col.group <- as.character(col.group)
                                                                                                                                                                                                                                                            ## MDS por grupos de edad
                                                                                                                                                                                                                                                            plotMDS(vGene$E, labels = df$AgeGroup, col = col.group)
                                                                                                                                                                                                                                                            ## Conviertiendo los valores de Sex a colores
                                                                                                                                                                                                                                                            col.sex <- df$Sex
                                                                                                                                                                                                                                                            levels(col.sex) <- brewer.pal(nlevels(col.sex), "Dark2")
                                                                                                                                                                                                                                                            col.sex <- as.character(col.sex)
                                                                                                                                                                                                                                                            ## MDS por sexo
                                                                                                                                                                                                                                                            plotMDS(vGene$E, labels = df$Sex, col = col.sex)
                                                                                                                                                                                                                                                            ```

                                                                                                                                                                                                                                                            ## Ejercicio

                                                                                                                                                                                                                                                            Agreguen los nombres de los genes a nuestro `pheatmap`.

                                                                                                                                                                                                                                                            Pistas:

                                                                                                                                                                                                                                                              * Revisen la información de `rowRanges(rse_gene_SRP045638)`.
                                                                                                                                                                                                                                                            * Exploren que hace la función `match()`.


                                                                                                                                                                                                                                                            ## Comunidad

                                                                                                                                                                                                                                                            Algunxs de lxs autores de `ExploreModelMatrix`:

                                                                                                                                                                                                                                                              * https://twitter.com/CSoneson
                                                                                                                                                                                                                                                            * https://twitter.com/FedeBioinfo
                                                                                                                                                                                                                                                            * https://twitter.com/mikelove

                                                                                                                                                                                                                                                            Algunxs de lxs autores de `edgeR` y `limma`:

                                                                                                                                                                                                                                                              * https://twitter.com/mritchieau
                                                                                                                                                                                                                                                            * https://twitter.com/davisjmcc
                                                                                                                                                                                                                                                            * https://twitter.com/markrobinsonca
                                                                                                                                                                                                                                                            * https://twitter.com/AliciaOshlack

                                                                                                                                                                                                                                                            <blockquote class="twitter-tweet"><p lang="en" dir="ltr">If you&#39;ve ever been dazed by design matrices or confused by contrasts when performing gene expression analysis in limma, the new article by Charity Law is for you <a href="https://t.co/ZSMOA20tdm">https://t.co/ZSMOA20tdm</a> <a href="https://twitter.com/hashtag/bioconductor?src=hash&amp;ref_src=twsrc%5Etfw">#bioconductor</a> <a href="https://twitter.com/hashtag/rstats?src=hash&amp;ref_src=twsrc%5Etfw">#rstats</a> (1/2)</p>&mdash; Matt Ritchie (@mritchieau) <a href="https://twitter.com/mritchieau/status/1338639551128952832?ref_src=twsrc%5Etfw">December 15, 2020</a></blockquote> <script async src="https://platform.twitter.com/widgets.js" charset="utf-8"></script>


                                                                                                                                                                                                                                                              ## Ejercicio: respuesta


                                                                                                                                                                                                                                                              ```{r respuesta, out.height="800px"}
                                                                                                                                                                                                                                                            ## Tenemos que usar gene_id y gene_name
                                                                                                                                                                                                                                                            rowRanges(rse_gene_SRP045638)
                                                                                                                                                                                                                                                            ## Con match() podemos encontrar cual es cual
                                                                                                                                                                                                                                                            rownames(exprs_heatmap) <- rowRanges(rse_gene_SRP045638)$gene_name[match(rownames(exprs_heatmap), rowRanges(rse_gene_SRP045638)$gene_id)]
                                                                                                                                                                                                                                                            ## Y luego podemos cambiar el valor de show_rownames de FALSE a TRUE
                                                                                                                                                                                                                                                            pheatmap(
                                                                                                                                                                                                                                                              exprs_heatmap,
                                                                                                                                                                                                                                                              cluster_rows = TRUE,
                                                                                                                                                                                                                                                              cluster_cols = TRUE,
                                                                                                                                                                                                                                                              show_rownames = TRUE,
                                                                                                                                                                                                                                                              show_colnames = FALSE,
                                                                                                                                                                                                                                                              annotation_col = df
                                                                                                                                                                                                                                                            )
                                                                                                                                                                                                                                                            ```

                                                                                                                                                                                                                                                            ## Específicaciones del proyecto

                                                                                                                                                                                                                                                            * Con datos de algún estudio disponible vía `recount3`, hagan un análisis de expresión diferencial.
                                                                                                                                                                                                                                                            * Incluyan al menos 3 gráficas en su reporte.
                                                                                                                                                                                                                                                            * Su reporte debe ser público y estar listado en el [Google Sheet](https://docs.google.com/spreadsheets/d/1sOBAnPkN_mP_Tq6-a8TyO7T4ii_hRPLFlKXj4qJdfUs/edit?usp=sharing) del curso.

                                                                                                                                                                                                                                                            Suena fácil, pero cada estudio tiene sus complejidades.

                                                                                                                                                                                                                                                            Hay muchos paquetes que no vimos que les pueden llamar la atención, tal como `ideal`. En http://research.libd.org/SPEAQeasy-example/bootcamp_intro pueden encontrar varias gráficas que tal vez les quieran reproducir. En fin, ¡esto solo es el inicio!

                                                                                                                                                                                                                                                              <blockquote class="twitter-tweet"><p lang="en" dir="ltr">🎉🎉🎉Our new MS is finally out! Given the timing, Santa had an early round with us 🎅<br>💡<a href="https://t.co/a0dHFGWN7V">https://t.co/a0dHFGWN7V</a>, &quot;ideal: an R/Bioconductor package for interactive differential expression analysis&quot;.<br><br>I promise a proper <a href="https://twitter.com/hashtag/rstats?src=hash&amp;ref_src=twsrc%5Etfw">#rstats</a> hexsticker will follow, for now enjoy the package 😉</p>&mdash; Federico Marini (@FedeBioinfo) <a href="https://twitter.com/FedeBioinfo/status/1336944561592078336?ref_src=twsrc%5Etfw">December 10, 2020</a></blockquote> <script async src="https://platform.twitter.com/widgets.js" charset="utf-8"></script>
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            © 2021 GitHub, Inc.
