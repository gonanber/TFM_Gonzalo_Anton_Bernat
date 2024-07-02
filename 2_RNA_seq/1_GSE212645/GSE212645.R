###############################################################################
# Autor: Gonzalo Anton Bernat
# Fecha: 16/11/2023
# Descripcion: 
#   Este programa permite realizar un análisis exploratorio de los
#   datos de RNA seq y el análisis de expresión diferencial.
###############################################################################


# Se cargan las librerías
library(dplyr)
library(ggplot2) 
library(gridExtra)
library(DESeq2)
library(edgeR)
library(pheatmap)
library(reshape2)
library(DT)
library(AnnotationDbi)
library(org.Hs.eg.db)

## Se cargan los datos
load("se_GSE21645.rda")
conteos <- assay(se)
metadatos <- colData(se)
class(se)
dim(se)
colData(se)
rowData(se)


################################################################################
# PREPARACIÓN DE LOS DATOS                                                     #
################################################################################

## Se crea una nueva columna 'sex_diagnosis' 
datos <- colData(se)

datos$sex_diagnosis <- ifelse(
  grepl("F", datos$sex) & grepl("ASD", datos$status), "female_ASD",
  ifelse(
    grepl("F", datos$sex) & grepl("SIB", datos$status), "female_control",
    ifelse(
      grepl("M", datos$sex) & grepl("ASD", datos$status), "male_ASD",
      ifelse(
        grepl("M", datos$sex) & grepl("SIB", datos$status), "male_control",
        NA
      )
    )
  )
)

colData(se)$sex_diagnosis <- datos$sex_diagnosis

# Se cuentan las diferentes categorías de sex_diagnosis
nrow(datos) # Totales
table(datos$sex_diagnosis) # Por categoria

## Se obtiene la matriz de conteos en la que las filas son los genes y las columnas las muestras
conteos <- assay(se)
conteos <- as.data.frame(conteos)
dimensiones <- dim(conteos)
nombres_columnas <- paste("sample", 1:dimensiones[2], sep="_")
colnames(conteos) <- nombres_columnas

## Se obtienen los metadatos en los que cada fila es una muestra
metadata <- colData(se)
rownames(metadata) <- nombres_columnas
metadata$Sample_ID <- nombres_columnas
metadata

## Preparacion de los datos
conteos_tidy <- log(conteos + 1)
conteos_tidy$genes <- rownames(conteos_tidy)
conteos_tidy <- reshape2::melt(conteos_tidy)
colnames(conteos_tidy) <- c("genes", "Sample_ID", "value")
conteos_tidy <- merge(conteos_tidy, metadata, by.x = "Sample_ID", by.y = "Sample_ID")
conteos_tidy <- as.data.frame(conteos_tidy)


################################################################################
# ANÁLISIS EXPLORATORIO DE LOS GENES                                           #
################################################################################

## Summary de cada muestra
summary(conteos)

## Boxplot de conteos por muestra para ver la distribución de los conteos
ggplot(conteos_tidy,
       aes(x=Sample_ID,
           y=value,
           fill=sex_diagnosis)) +
  geom_boxplot() +
  scale_fill_manual(values=c(male_ASD="#DE369D",
                             male_control="#ffd166",
                             female_ASD="#06d6a0",
                             female_control="#6AD8F6")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  ggtitle("Boxplot") + 
  xlab("Samples") + 
  ylab("Log10 Read counts")


####  Alrededor de la mitad de los valores en cada muestra son ceros. Este conjunto 
####  de datos se beneficiaría de un filtrado de recuento bajo.

## Diagrama de barras de genes detectados en cada muestra
{barplot(colSums(conteos > 0),
         ylab="Number of detected genes",
         las=2)
  abline(h=median(colSums(conteos > 0)))}

## Histograma de en cuantas muestras tiene expresión cada gen
n_samples <- rowSums(conteos > 0)
hist(n_samples,
     xlab="Number of samples",
     breaks=12)


################################################################################
# Filtrado de genes con conteos bajos                                          #
################################################################################

# Mimimo 1 cpm en 2 o más muestras
keep_genes <- rowSums(edgeR::cpm(conteos) > 1) >= 2
table(keep_genes)
conteos_filtrados <- conteos[keep_genes,]
dim(conteos_filtrados)

# Se vuelven a preparar los datos
conteos_tidy = log(conteos_filtrados + 1)
conteos_tidy$genes = rownames(conteos_tidy)
conteos_tidy = reshape2::melt(conteos_tidy)
colnames(conteos_tidy) = c("genes", "Sample_ID", "value")
conteos_tidy = merge(conteos_tidy,
                     metadata,
                     by.x = "Sample_ID",
                     by.y = "Sample_ID")
conteos_tidy <- as.data.frame(conteos_tidy)


################################################################################
# ANÁLISIS EXPLORATORIO PRE-NORMALIZACIÓN                                                     #
################################################################################

# Boxplot
ggplot(conteos_tidy,
       aes(x=Sample_ID,
           y=value,
           fill=sex_diagnosis)) +
  geom_boxplot() +
  scale_fill_manual(values=c(male_ASD="#DE369D",
                             male_control="#ffd166",
                             female_ASD="#06d6a0",
                             female_control="#6AD8F6")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  ggtitle("Boxplot") + 
  xlab("Samples") + 
  ylab("Log10 Read counts")


####  Es una buena idea comprobar la correlación entre las muestras. Las muestras 
####  de RNA-Seq suelen tener una correlación muy alta (R^2 > 0,9). Los valores de 
####  R^2 inferiores a 0,8 pueden indicar que se trata de una muestra atípica. 

## Análisis de correlación
cor_matrix <- as.matrix(cor(conteos_filtrados, method="spearman"))
annot_col <- metadata[,"sex_diagnosis", drop=F]
annot_row <- metadata[,"sex_diagnosis", drop=F]
annot_colors <-list(sex_diagnosis=c(male_ASD="#DE369D",
                                    male_control="#ffd166",
                                    female_ASD="#06d6a0",
                                    female_control="#6AD8F6"))
heatmap_color <- colorRampPalette(c("navy", "white", "red"))(50)

pheatmap(cor_matrix,
         annotation_colors = annot_col,
         border_color="white",
         annotation_legend=T,
         color=heatmap_color)


## PCA
pca.res <- PCAtools::pca(mat=conteos_filtrados,
                         metadata=metadata,
                         scale=TRUE)
PCAtools::screeplot(pca.res)
PCAtools::biplot(pca.res,
                 colby='sex_diagnosis',
                 colkey = c(male_ASD="#DE369D",
                            male_control="#ffd166",
                            female_ASD="#06d6a0",
                            female_control="#6AD8F6"),
                 legendPosition = "right")


## Clustering Jerarquico
distancia = dist(1 - cor(conteos_filtrados,
                         method="pearson"))
nnot_col = metadata[,"sex_diagnosis", drop=F]
annot_row = metadata[,"sex_diagnosis", drop=F]
heatmap_color = colorRampPalette(c("navy", "white", "red"))(50)

pheatmap(cor_matrix,
         border_color="white",
         annotation_legend=T, 
         annotation_colors=annot_colors, 
         color=heatmap_color)


################################################################################
# NORMALIZACIÓN DE LOS DATOS                                                   #
################################################################################

#### DESeq2 corrige internamente los recuentos para la profundidad de secuenciación 
#### y el sesgo de composición del ARN utilizando el método de "median ratios". 
#### Para ejecutar este método, se crea un objeto DESeq2

# Creacion del objeto
metadata$sex_diagnosis <- factor(metadata$sex_diagnosis)
datos <- DESeqDataSetFromMatrix(countData=conteos_filtrados,
                                colData=metadata,
                                design=~sex_diagnosis)

# Normalizacion
datos <- DESeq2::estimateSizeFactors(datos, type="ratio")
head(counts(datos, normalized=FALSE))
head(counts(datos, normalized=TRUE))

# Se guardan los conteos en el slot assay
assay(datos, "conteos_norm") <- counts(datos, normalized=TRUE)
assay(datos, "conteos_norm_log") <- log(assay(datos, "conteos_norm") + 1)
conteos_norm_log <- assay(datos, "conteos_norm_log")


################################################################################
# ANÁLISIS EXPLORATORIO POST NORMALIZACIÓN                                     #
################################################################################

# Boxplot
conteos_tidy = as.data.frame(conteos_norm_log)
conteos_tidy$genes = rownames(conteos_tidy)
conteos_tidy = reshape2::melt(conteos_tidy)
colnames(conteos_tidy) = c("genes", "Sample_ID", "value")
conteos_tidy = merge(conteos_tidy,
                     metadata,
                     by.x = "Sample_ID",
                     by.y = "Sample_ID")
conteos_tidy <- as.data.frame(conteos_tidy)
dir.create("images")
jpeg(filename = "images/boxplot.jpg",  width = 800, height = 400, quality = 100)
ggplot(conteos_tidy,
       aes(x=Sample_ID,
           y=value,
           fill=sex_diagnosis)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  ggtitle("Comparación de la expresión genética entre los grupos experimentales") + 
  labs(fill = "Grupos") +
  xlab("Muestras") + 
  ylab("Log10 conteos normalizados")

dev.off()


## Análisis de correlación
cor_matrix <- as.matrix(cor(conteos_norm_log, method="spearman"))
annot_col <- metadata[,"sex_diagnosis", drop=F]
annot_row <- metadata[,"sex_diagnosis", drop=F]
annot_colors <-list(sex_diagnosis=c(male_ASD="#DE369D",
                                    male_control="#ffd166",
                                    female_ASD="#06d6a0",
                                    female_control="#6AD8F6"))
heatmap_color <- colorRampPalette(c("navy", "white", "red"))(50)

jpeg(filename = "images/clustering.jpg", width = 1000, height = 800, quality = 100)
pheatmap(cor_matrix,
         annotation_colors = annot_col,
         border_color="white",
         annotation_legend=T,
         color=heatmap_color,
         main = "Clustering jerárquico y mapa de calor")
dev.off()

# PCA
pca.res <- PCAtools::pca(mat=conteos_norm_log,
                         metadata=metadata,
                         scale=TRUE)
PCAtools::screeplot(pca.res)
jpeg(filename = "images/PCA.jpg", width = 800, height = 400, quality = 100)
biplot <- PCAtools::biplot(pca.res,
                           colby = 'sex_diagnosis',
                           legendPosition = "right",
                           drawConnectors = T,
                           colLegendTitle = "Grupos",
                           borderWidth = 0,
                           max.overlaps = 0 
)

biplot <- biplot +
  theme(panel.grid.major = element_line(linewidth = 0.5),
        panel.grid.minor = element_blank(), 
        plot.title = element_text(size = 12),  
        axis.title.x = element_text(size = 12),  
        axis.title.y = element_text(size = 12),  
        legend.title = element_text(size = 12),  
        legend.text = element_text(size = 12)  
  ) +
  ggtitle("Análisis de componentes principales") 
biplot
dev.off()


################################################################################
# ANÁLISIS DE EXPRESIÓN DIFERENCIAL CON LIMMA-VOOM                             #
################################################################################

# Se construye la matriz del diseño
condicion <- metadata$sex_diagnosis
experiment <- metadata$experiment
design <- model.matrix(~0 + condicion + experiment)
columnas <- c(levels(factor(metadata$sex_diagnosis)), "batch_2")
colnames(design) <- columnas

# Se transforma los datos de RNA-seq para hacer modelos lineales
v <- limma::voom(conteos_filtrados, design)

# Se genera la matriz de contrastes
cont_matrix <- limma::makeContrasts( 
  dif1 = (female_ASD) - (female_control), # Contraste mujeres
  dif2 = (male_ASD) - (male_control),  # Contraste hombres
  dif3 = (male_ASD + female_ASD) - (male_control + female_control), # Enfermedad sin sexo
  dif4 = (female_ASD - female_control) - (male_ASD - male_control), # Diferencias de sexo
  
  levels=design
)

# Se realiza el ajuste a un modelo lineal
fit <- limma::lmFit(v, design)
fit2 <- contrasts.fit(fit, cont_matrix)
fit2 <- eBayes(fit2)
#save(fit2, file = "../../Microarrays/fits/GSE212645.rda")

# Se extraen los resultados. Cambiar coef para cada contraste
resultados <- topTable(fit2,
                       coef=4,
                       adjust="BH",
                       number = nrow(conteos_norm_log))

table(resultados[, "adj.P.Val"] < 0.05)


################################################################################
# ANÁLISIS DE EXPRESIÓN DIFERENCIAL CON EDGE-R                                 #
################################################################################

# Generamos objeto DGE con los datos crudos (estos ya están filtrados de antes,
# pero también podrías construir primero el objeto con el objeto "conteos" y
# luego filtrar con filterByExprs)
dge = DGEList(counts=conteos_filtrados, samples = metadata, group = metadata$sex_diagnosis)

# Normalización TMM
dge <- calcNormFactors(dge, method = "TMM")

# DGE
design = model.matrix(~ 0 + dge$samples$sex_diagnosis + dge$samples$experiment) 

dge = estimateDisp(dge, design)
fit = glmQLFit(dge, design)

colnames(design) # Comprueba el orden de las columnas para establecer los contrastes

dif1 = c(1,-1, 0, 0, 0) # Contraste mujeres
dif2 = c(0, 0, 1, -1, 0) # Contraste hombres
dif3 = c(1, -1, 1, -1, 0) # Enfermedad sin sexo
dif4 = c(1, -1, -1, 1, 0) # Diferencias de sexo

test.dif = glmQLFTest(fit, contrast=dif4) 
results = topTags(test.dif,
                  n = nrow(dge))

results = results$table
table(results$PValue < 0.05)

# Los resultados obtenidos están sin ajustar el p valor, por lo que es 
# necesaro ajustarlos mediante BH
x <- p.adjust(results$PValue, "BH")
results$p.adj <- x
table(results$p.adj < 0.05)


################################################################################
# ANÁLISIS DE EXPRESIÓN DIFERENCIAL CON DESEQ2                                 #
################################################################################

# Se construye el objeto DESeqDataSet incluyendo el diseño

metadata$experiment = factor(metadata$experiment)
metadata$sex_diagnosis = factor(metadata$sex_diagnosis)

datos <- DESeqDataSetFromMatrix(countData=conteos_filtrados,
                                colData=metadata,
                                design=~ 0 + sex_diagnosis + experiment)

# Normalización median ratios
datos <- DESeq2::estimateSizeFactors(datos, type="ratio")

# DGE
datos.deseq = DESeq(datos)
resultsNames(datos.deseq) # Se comprueba el orden de las columnas para establecer los contrastes

dif1 = c(1,-1, 0, 0, 0) # Contraste mujeres
dif1 = list(c("sex_diagnosisfemale_ASD"), c("sex_diagnosisfemale_control"))

dif2 = c(0, 0, 1, -1, 0) # Contraste hombres
dif2 = list(c("sex_diagnosismale_ASD"), c("sex_diagnosismale_control"))

dif3 = c(1, -1, 1, -1, 0) # Enfermedad sin sexo
dif3 = list(c("sex_diagnosisfemale_ASD", "sex_diagnosismale_ASD"), c("sex_diagnosisfemale_control", "sex_diagnosismale_control"))

dif4 = c(1, -1, -1, 1, 0) # Diferencias de sexo
dif4 = list(c("sex_diagnosisfemale_ASD", "sex_diagnosismale_control"), c("sex_diagnosisfemale_control",  "sex_diagnosismale_ASD"))


res1 = as.data.frame(results(datos.deseq, contrast = dif1,                 
                                    alpha=0.05,
                                    pAdjustMethod = "BH"))

res2 = as.data.frame(results(datos.deseq, contrast = dif2,                
                             alpha=0.05,
                             pAdjustMethod = "BH"))

res3 = as.data.frame(results(datos.deseq, contrast = dif3,               
                             alpha=0.05,
                             pAdjustMethod = "BH"))

res4 = as.data.frame(results(datos.deseq, contrast = dif4,              
                             alpha=0.05,
                             pAdjustMethod = "BH"))

table(res1$padj < 0.05)
table(res2$padj < 0.05)
table(res3$padj < 0.05)
table(res4$padj < 0.05)


# Distribución de p-valores
hist(res1$pvalue[res1$baseMean > 1],
     main="res1 Pval distribution",
     xlab="P-values")

hist(res2$pvalue[res2$baseMean > 1],
     main="res2 Pval distribution",
     xlab="P-values")

hist(res3$pvalue[res3$baseMean > 1],
     main="res3 Pval distribution",
     xlab="P-values")

hist(res4$pvalue[res4$baseMean > 1],
     main="res4 Pval distribution",
     xlab="P-values")



#######################################################################################################################################
################################################################################
# ANOTACIÓN                                                                    #
################################################################################

# Se crea una columna con el nombre de los genes
df_res <- as.data.frame(res1)
df_res$external_gene_name <- rownames(df_res)

# Se lee el fichero de anotación de humano
hg <- read.delim("human_genes.txt",
                 header=T,
                 sep="\t",
                 stringsAsFactors=F)
hg <- hg[!duplicated(hg$ensembl_gene_id),]

# Se realiza el merge de los dataframes
df_res1 <- merge(df_res,
                 hg,
                 by="external_gene_name")

# Se muestra en formato de tabla bonita
DT::datatable(df_res1, filter="top", rownames = FALSE) %>%
  DT::formatRound(columns = c(2:7), digits=4) %>%
  DT::formatStyle(3, background = DT::styleInterval(c(0), c("orange", "yellow")))

# Se corrigen los valores de log2 fold change
lres1 <- lfcShrink(dge_results,
                   coef="sex_diagnosis_female_control_vs_female_ASD",
                   res=res1,
                   type="normal")

lres2 <- lfcShrink(dge_results,
                   coef="sex_diagnosis_male_ASD_vs_female_ASD",
                   res=res2,
                   type="normal")
head(lres1)
head(res1)


################################################################################
# VISUALIZACIÓN                                                                #
################################################################################

## MA-plot
DESeq2::plotMA(res1,
               alpha=0.05,
               colNonSig = "gray60",
               colSig = "blue")

## Volcano plot
sig_limit <- 0.05
lfc_limit <- 2
res1[,"Significative"] <- rep("none", length(rownames(res1)))

# Se cambian las etiquetas
res1[which(res1$padj < sig_limit & abs(res1$log2FoldChange) < lfc_limit ),"Significative"] <- "p.adjusted"
res1[which(res1$padj < sig_limit & abs(res1$log2FoldChange) > lfc_limit ), "Significative"] <- "p.adjusted + logFC"

ggplot(data = as.data.frame(res1),
       aes(y = -log(padj),
           x = log2FoldChange)) +
  geom_point(aes(color = Significative), na.rm = TRUE) +
  scale_color_manual(values = c("#A19D9F", "#E6007B", "#009E73")) + #https://stackoverflow.com/q/57153428
  xlab("logFC") + ylab("-log10(p.adjusted)")

## Counts Plot
plotCounts(datos,
           gene=rownames(res1)[1],
           intgroup="sex_diagnosis",
           normalized=T)

sessionInfo()
