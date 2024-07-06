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
library(SummarizedExperiment)

# Se cargan los datos
load("../../00_Datos/GSE140702.rda")
se <- seGSE140702
class(se)
dim(se)
colData(se)
rowData(se)


################################################################################
# PREPARACIÓN DE LOS DATOS                                                     #
################################################################################

datos <- colData(se)

# Filtrado de los metadatos para quedarnos con los de interés
dim(datos)
datos <- datos[datos$diagnosis.ch1 != "Pervasive Developmental Disorder Not Otherwise Specified (PDDNOS) or Asperger", ]
dim(datos)
datos <- datos[!(datos$treatment.ch1 %in% c("LPS", "LTA")), ]
dim(datos)


# Filtrado de los conteos para quedarnos con los de interés
muestras <- rownames(datos)

# Filtrado de la matriz de conteos para quedarnos con las muestras de interés
dim(assay(se))
conteos <- subset(assay(se), select = muestras)
dim(conteos)

# Se crea el nuevo Summarized Experiment con los datos de interés
se <- SummarizedExperiment(assays = list(conteos),
                           colData = datos)


## Se crea una nueva columna 'sex_diagnosis' 
datos$sex_diagnosis <- ifelse(
  grepl("^Female", datos$subject.sex.ch1) & grepl("Full Autistic Disorder", datos$diagnosis.ch1), "female_ASD",
  ifelse(
    grepl("^Female", datos$subject.sex.ch1) & grepl("Typical", datos$diagnosis.ch1), "female_control",
    ifelse(
      grepl("^Male", datos$subject.sex.ch1) & grepl("Full Autistic Disorder", datos$diagnosis.ch1), "male_ASD",
      ifelse(
        grepl("^Male", datos$subject.sex.ch1) & grepl("Typical", datos$diagnosis.ch1), "male_control",
        NA
      )
    )
  )
)

colData(se)$sex_diagnosis <- datos$sex_diagnosis


## Se obtiene la matriz de conteos en la que las filas son los genes y las columnas las muestras
conteos <- assay(se)
conteos <- as.data.frame(conteos)
dimensiones <- dim(conteos)
nombres_muestras <- colnames(se)
nombres_columnas <- paste("sample", 1:dimensiones[2], sep="_")
colnames(conteos) <- nombres_columnas
equivalencias <- data.frame(nombres_columnas, nombres_muestras)

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
jpeg(filename = "images/GSE140702.jpg",  width = 800, height = 400, quality = 100)
ggplot(conteos_tidy,
       aes(x=Sample_ID,
           y=value,
           fill=sex_diagnosis)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  ggtitle("Comparación de la expresión genética entre los  grupos experimentales") +
  labs(fill = "Groups") +
  xlab("Muestras") + 
  ylab("Log10 Read counts")

dev.off()
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


## PCA ##########################################################################################
pca.res <- PCAtools::pca(mat=conteos_filtrados,
                         metadata=metadata,
                         scale=TRUE)
PCAtools::screeplot(pca.res)

jpeg(filename = "images/GSE140702_PCA.jpg", width = 800, height = 400, quality = 100)
biplot <- PCAtools::biplot(pca.res,
                           colby = 'sex_diagnosis',
                           legendPosition = "right",
                           drawConnectors = T,  
                           borderWidth = 0,
                           max.overlaps = 0 
)

biplot <- biplot +
  theme(panel.grid.major = element_line(linewidth = 0.5), # Grosor uniforme para las líneas principales de la cuadrícula
        panel.grid.minor = element_blank() # Eliminar líneas menores de la cuadrícula
  )
biplot
dev.off()





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
  theme(axis.text.x = element_blank()) +
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


# Permite eliminar los outliers. Volver a ejecutar desde los boxplots
result <- result$data
datos_x <- subset(result, x < -50)
datos_y <- subset(result, y < -1000)
outliers_x <- unique(rownames(datos_x))
outliers_y <- unique(rownames(datos_y))
outliers <- union(outliers_x, outliers_y)

condiciones <- colnames(conteos_norm_log)
outliers <- paste0("sample_", outliers)
eliminar <- subset(equivalencias, nombres_columnas %in% outliers)
eliminar <- c(eliminar$nombres_muestras)
eliminar

################################################################################
# Construccion del nuevo se
################################################################################

cont <- assay(se)
guardar <- setdiff(colnames(se), eliminar)
cont <- subset(cont, select = guardar)

se2 <- se
met <- colData(se)
met <- met[!(met$geo_accession %in% eliminar), ]

se <- SummarizedExperiment(assays = list(cont),
                                    colData = met)


save(se, file = "../../00_Datos/seGSE140702_sin_outliers.rda")

################################################################################
# ANOTACION                                                                    #
################################################################################

load("../../00_Datos/seGSE140702_sin_outliers.rda")

# Obtener los Entrez IDs de los 'rownames' del objeto 'se'
ensemble_ids <- rownames(se)

# Anotar los Entrez IDs a símbolos
genedata <- select(org.Hs.eg.db, keys = ensemble_ids, keytype = "ENTREZID", column = "SYMBOL")
genedata$Symbol <- sub(" ///.*","",genedata$SYMBOL)
genedata[genedata == ""] <- "eliminar"

medianReps <- function(matriz){
  ID <- as.character(rownames(matriz))
  ID <- factor(ID, levels = unique(ID))
  df <- by(matriz, ID, function(x) apply(x, 2, stats::median))
  mat <- do.call("rbind", df)
  return(mat)
}

exp2 <- assay(se)
rownames(exp2) <- genedata$Symbol
exp2 <- medianReps(exp2)
exp3 <- exp2[rownames(exp2) != "eliminar",]
genedata <- genedata[match(rownames(exp3),genedata$Symbol),c("Symbol", "ENTREZID" )]
rownames(genedata) <- genedata$Symbol

metadatos <- colData(se)
se2 <- SummarizedExperiment(assays = list(exp3),
                            colData = metadatos,
                            rowData = genedata)

save(se2, file = "../../00_Datos/se2GSE140702.rda")


################################################################################
# ANÁLISIS DE EXPRESIÓN DIFERENCIAL CON LIMMA-VOOM                             #
################################################################################

load("../../00_Datos/se2GSE140702.rda")

metadata <- colData(se2)
# Mimimo 1 cpm en 2 o más muestras
conteos <- assay(se2)
conteos <- as.data.frame(conteos)
keep_genes <- rowSums(edgeR::cpm(conteos) > 1) >= 2
table(keep_genes)
conteos_filtrados <- conteos[keep_genes,]
dim(conteos_filtrados)


# Se construye la matriz del diseño
condicion <- metadata$sex_diagnosis
design <- model.matrix(~0 + condicion)
columnas <- levels(factor(metadata$sex_diagnosis))
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
head(fit2)
save(fit2, file = "../../3_Genera_todos_fit/fit_GSE140702.rda")

# Se extraen los resultados
res1 <- topTable(fit2,
                 coef=1,
                 adjust="BH",
                 number = nrow(conteos_filtrados))

res2 <- topTable(fit2,
                 coef=2,
                 adjust="BH",
                 number = nrow(conteos_filtrados))

res3 <- topTable(fit2,
                 coef=3,
                 adjust="BH",
                 number = nrow(conteos_filtrados))

res4 <- topTable(fit2,
                 coef=4,
                 adjust="BH",
                 number = nrow(conteos_filtrados))

table(res1[, "adj.P.Val"] < 0.05)
table(res2[, "adj.P.Val"] < 0.05)
table(res3[, "adj.P.Val"] < 0.05)
table(res4[, "adj.P.Val"] < 0.05)


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
design = model.matrix(~ 0 + dge$samples$sex_diagnosis) 

dge = estimateDisp(dge, design)
fit = glmQLFit(dge, design)

colnames(design) # Comprueba el orden de las columnas para establecer los contrastes

dif1 = c(1,-1, 0, 0) # Contraste mujeres
dif2 = c(0, 0, 1, -1) # Contraste hombres
dif3 = c(1, -1, 1, -1) # Enfermedad sin sexo
dif4 = c(1, -1, -1, 1) # Diferencias de sexo

res1 = glmQLFTest(fit, contrast=dif1) # Cambiarías el valor de contrast para hacer cada contraste
res1 = topTags(res1,
                  n = nrow(dge))

res2 = glmQLFTest(fit, contrast=dif2) # Cambiarías el valor de contrast para hacer cada contraste
res2 = topTags(res2,
               n = nrow(dge))
res3 = glmQLFTest(fit, contrast=dif3) # Cambiarías el valor de contrast para hacer cada contraste
res3 = topTags(res3,
               n = nrow(dge))
res4 = glmQLFTest(fit, contrast=dif4) # Cambiarías el valor de contrast para hacer cada contraste
res4 = topTags(res4,
               n = nrow(dge))

table(res1$table$PValue < 0.05)
table(res2$table$PValue < 0.05)
table(res3$table$PValue < 0.05)
table(res4$table$PValue < 0.05)

# Los resultados obtenidos están sin ajustar el p valor, por lo que es 
# necesaro ajustarlos mediante BH
x <- p.adjust(res1$table$PValue, "BH")
res1$p.adj <- x
table(res1$p.adj < 0.05)

x <- p.adjust(res2$table$PValue, "BH")
res2$p.adj <- x
table(res2$p.adj < 0.05)

x <- p.adjust(res3$table$PValue, "BH")
res3$p.adj <- x
table(res3$p.adj < 0.05)

x <- p.adjust(res4$table$PValue, "BH")
res4$p.adj <- x
table(res4$p.adj < 0.05)


################################################################################
# ANÁLISIS DE EXPRESIÓN DIFERENCIAL CON DESEQ2                                 #
################################################################################

# Se construye el objeto DESeqDataSet incluyendo el diseño
save.image()
load(".RData")

metadata$sex_diagnosis = factor(metadata$sex_diagnosis)

datos <- DESeqDataSetFromMatrix(countData=conteos_filtrados,
                                colData=metadata,
                                design=~ 0 + sex_diagnosis)

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

