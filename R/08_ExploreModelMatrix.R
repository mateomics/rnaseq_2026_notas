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

BiocManager::install("ExploreModelMatrix")

## Usaremos shiny otra ves
app <- ExploreModelMatrix(
    sampleData = sampleData,
    designFormula = ~ genotype + treatment
)
if (interactive()) shiny::runApp(app)


## Vamos a usar datos de https://www.ncbi.nlm.nih.gov/sra/?term=SRP045638 procesados con recount3.
## Primero hay que descargar los datos con los comandos que vimos ayer.

library("recount3")
human_projects <- available_projects()

rse_gene_SRP045638 <- create_rse(
    subset(
        human_projects,
        project == "SRP045638" & project_type == "data_sources"
    )
)

assay(rse_gene_SRP045638, "counts") <- compute_read_counts(rse_gene_SRP045638)

## Una vez descargados y con los números de lecturas podemos usar expand_sra_attributes(). 
# Sin embargo, tenemos un problema con estos datos.
rse_gene_SRP045638$sra.sample_attributes[1:3]

## Vamos a intentar resolverlo eliminando información que está presente solo en ciertas muestras.
rse_gene_SRP045638$sra.sample_attributes <- gsub("dev_stage;;Fetal\\|", "", rse_gene_SRP045638$sra.sample_attributes)
rse_gene_SRP045638$sra.sample_attributes[1:3]

## Y ahora sí podemos continuar...
rse_gene_SRP045638 <- expand_sra_attributes(rse_gene_SRP045638)

colData(rse_gene_SRP045638)[
    ,
    grepl("^sra_attribute", colnames(colData(rse_gene_SRP045638)))
]

## Como ahora si vamos a usar esta información para un modelo estadístico, 
# será importante que tengamos en el formato correcto de R a la información que vamos a usar.

## Pasar de character a numeric o factor
rse_gene_SRP045638$sra_attribute.age <- as.numeric(rse_gene_SRP045638$sra_attribute.age)
rse_gene_SRP045638$sra_attribute.disease <- factor(tolower(rse_gene_SRP045638$sra_attribute.disease))
rse_gene_SRP045638$sra_attribute.RIN <- as.numeric(rse_gene_SRP045638$sra_attribute.RIN)
rse_gene_SRP045638$sra_attribute.sex <- factor(rse_gene_SRP045638$sra_attribute.sex)

## Resumen de las variables de interés
summary(as.data.frame(colData(rse_gene_SRP045638)[
    ,
    grepl("^sra_attribute.[age|disease|RIN|sex]", colnames(colData(rse_gene_SRP045638)))
]))

## Ahora crearemos un par de variables para que las podamos usar en nuestro análisis.

## Encontraremos diferencias entre muestra prenatalas vs postnatales
rse_gene_SRP045638$prenatal <- factor(ifelse(rse_gene_SRP045638$sra_attribute.age < 0, "prenatal", "postnatal"))
table(rse_gene_SRP045638$prenatal)

## http://rna.recount.bio/docs/quality-check-fields.html
rse_gene_SRP045638$assigned_gene_prop <- rse_gene_SRP045638$recount_qc.gene_fc_count_all.assigned / rse_gene_SRP045638$recount_qc.gene_fc_count_all.total
summary(rse_gene_SRP045638$assigned_gene_prop)

with(colData(rse_gene_SRP045638), plot(assigned_gene_prop, sra_attribute.RIN))

## Hm... veamos si hay una diferencia entre los grupos
with(colData(rse_gene_SRP045638), tapply(assigned_gene_prop, prenatal, summary))

## A continuación podemos eliminar algunas muestras que consideremos de baja calidad y genes con niveles de expresión muy bajos.
## Guardemos nuestro objeto entero de igual forma por cualquier cosa...
rse_gene_SRP045638_unfiltered <- rse_gene_SRP045638

## Eliminemos a muestras malas
hist(rse_gene_SRP045638$assigned_gene_prop)

table(rse_gene_SRP045638$assigned_gene_prop < 0.3)

rse_gene_SRP045638 <- rse_gene_SRP045638[, rse_gene_SRP045638$assigned_gene_prop > 0.3]

## Calculemos los niveles medios de expresión de los genes en nuestras
## muestras.
## Ojo: en un análisis real probablemente haríamos esto con los RPKMs o CPMs
## en vez de las cuentas.
## En realidad usariamos:
# edgeR::filterByExpr() https://bioconductor.org/packages/edgeR/ https://rdrr.io/bioc/edgeR/man/filterByExpr.html
# genefilter::genefilter() https://bioconductor.org/packages/genefilter/ https://rdrr.io/bioc/genefilter/man/genefilter.html
# jaffelab::expression_cutoff() http://research.libd.org/jaffelab/reference/expression_cutoff.html
#
gene_means <- rowMeans(assay(rse_gene_SRP045638, "counts"))

## Eliminamos genes
rse_gene_SRP045638 <- rse_gene_SRP045638[gene_means > 0.1, ]

## Dimensiones finales
dim(rse_gene_SRP045638)

## Porcentaje de genes que retuvimos
round(nrow(rse_gene_SRP045638) / nrow(rse_gene_SRP045638_unfiltered) * 100, 2)


## Este paquete se usa por razones históricas (en lugar de SummarizedExperiment) para el análisis de datos de RNA-seq, 
# pero también tiene funciones útiles para la normalización y filtrado de datos.
library("edgeR") # BiocManager::install("edgeR", update = FALSE)
dge <- DGEList(
    counts = assay(rse_gene_SRP045638, "counts"),
    genes = rowData(rse_gene_SRP045638)
)
dge <- calcNormFactors(dge)

library("ggplot2")
ggplot(as.data.frame(colData(rse_gene_SRP045638)), aes(y = assigned_gene_prop, x = prenatal)) +
    geom_boxplot() +
    theme_bw(base_size = 20) +
    ylab("Assigned Gene Prop") +
    xlab("Age Group")


## Por ahora continuaremos con el siguiente modelo estadístico.
mod <- model.matrix(~ prenatal + sra_attribute.RIN + sra_attribute.sex + assigned_gene_prop,
    data = colData(rse_gene_SRP045638)
)
colnames(mod)

## Ya teniendo el modelo estadístico, podemos usar limma para realizar el análisis de expresión diferencial como tal.
library("limma")
vGene <- voom(dge, mod, plot = TRUE)

eb_results <- eBayes(lmFit(vGene))

de_results <- topTable(
    eb_results,
    coef = 2,
    number = nrow(rse_gene_SRP045638),
    sort.by = "none"
)
dim(de_results)
head(de_results)

## El p-value ajustado es el que se usa para determinar si un gen es diferencialmente expresado o no
## También es importante observar el logFC para ver la magnitud del cambio.

## Genes diferencialmente expresados entre pre y post natal con FDR < 5%
table(de_results$adj.P.Val < 0.05)
## Visualicemos los resultados estadísticos
plotMA(eb_results, coef = 2)

# Y en este caso el eje Y es el logFC y el eje X es la media expresión de los genes.
# Cuando tenemos valores positivos en el eje Y, gen está más expresado en el grupo prenatal y viceversa

# De vGene$E podemos extraer los datos normalizados por limma-voom. Revisemos los top 50 genes diferencialmente expresados.

## Extraer valores de los genes de interés
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
