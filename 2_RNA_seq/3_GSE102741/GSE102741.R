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

## Se cargan los datos
load("../../00_Datos/GSE102741.rda")
se <- seGSE102741
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
  grepl("^Female", datos$Sex.ch1) & grepl("Autism", datos$disease.status.ch1), "female_ASD",
  ifelse(
    grepl("^Female", datos$Sex.ch1) & grepl("Healthy", datos$disease.status.ch1), "female_control",
    ifelse(
      grepl("^Male", datos$Sex.ch1) & grepl("Autism", datos$disease.status.ch1), "male_ASD",
      ifelse(
        grepl("^Male", datos$Sex.ch1) & grepl("Healthy", datos$disease.status.ch1), "male_control",
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
# ANÁLISIS EXPLORATORIO POST-NORMALIZACIÓN                                                     #
################################################################################

# Boxplot
dir.create("images")
jpeg(filename = "images/boxplot.jpg",  width = 800, height = 400, quality = 100)
ggplot(conteos_tidy, aes(x = Sample_ID, y = value, fill = sex_diagnosis)) +
  geom_boxplot() +
  theme(axis.text.x = element_blank()) +
  ggtitle("Comparación de la expresión genética entre los grupos experimentales") + 
  xlab("Muestras") + 
  ylab("Log10 conteos normalizados")

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


jpeg(filename = "images/clustering.jpg", width = 1000, height = 800, quality = 100)
pheatmap(cor_matrix,
         annotation_colors = annot_col,
         border_color="white",
         annotation_legend=T,
         color=heatmap_color,
         main = "Clustering jerárquico y mapa de calor")
dev.off()


## PCA ##########################################################################################
pca.res <- PCAtools::pca(mat=conteos_filtrados,
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

result <- biplot +
  theme(panel.grid.major = element_line(linewidth = 0.5),
        panel.grid.minor = element_blank(), 
        plot.title = element_text(size = 12),  
        axis.title.x = element_text(size = 12),  
        axis.title.y = element_text(size = 12),  
        legend.title = element_text(size = 12),  
        legend.text = element_text(size = 12)  
  ) +
  ggtitle("Análisis de componentes principales") 
result
dev.off()


# Permite eliminar los outliers. Volver a ejecutar desde los boxplots
result <- result$data
datos_x <- subset(result, x > 100)
datos_y <- subset(result, y < -200)
outliers_x <- unique(rownames(datos_x))
outliers_y <- unique(rownames(datos_y))
outliers <- union(outliers_x, outliers_y)

## Se obtiene la matriz de conteos en la que las filas son los genes y las columnas las muestras
conteos <- assay(se)
conteos <- as.data.frame(conteos)
dimensiones <- dim(conteos)
nombres_muestras <- colnames(se)
nombres_columnas <- paste("sample", 1:dimensiones[2], sep="_")
colnames(conteos) <- nombres_columnas
equivalencias <- data.frame(nombres_columnas, nombres_muestras)

condiciones <- colnames(conteos_filtrados)
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

save(se, file = "../../00_Datos/seGSE102741.rda")

################################################################################
# ANOTACION                                                                    #
################################################################################

load("../../00_Datos/seGSE102741.rda")
# Obtener los Entrez IDs de los 'rownames' del objeto 'se'
ensemble_ids <- rownames(se)

# Anotar los Entrez IDs a símbolos
genedata <- select(org.Hs.eg.db, keys = ensemble_ids, keytype = "ENSEMBL", column = "SYMBOL")
genedata$Symbol <- sub(" ///.*","",genedata$SYMBOL)
genedata[genedata == ""] <- "eliminar"


# Definir el número deseado de filas
numero_filas_deseadas <- 57659

# Contador de filas eliminadas
filas_eliminadas <- 0

# Iniciar el bucle
while (nrow(genedata) > numero_filas_deseadas) {
  # Buscar la primera fila con NA y eliminarla
  fila_con_na <- which(rowSums(is.na(genedata)) > 0)[1]
  if (!is.na(fila_con_na)) {
    genedata <- genedata[-fila_con_na, ]
    filas_eliminadas <- filas_eliminadas + 1
  } else {
    break  # Salir del bucle si no hay más filas con NA
  }
}


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
genedata <- genedata[match(rownames(exp3),genedata$Symbol),c("Symbol", "ENSEMBL" )]
rownames(genedata) <- genedata$Symbol

metadatos <- colData(se)
se2 <- SummarizedExperiment(assays = list(exp3),
                           colData = metadatos,
                           rowData = genedata)

save(se2, file = "../../00_Datos/se2GSE102741.rda")

################################################################################
# ANÁLISIS DE EXPRESIÓN DIFERENCIAL CON LIMMA-VOOM                             #
################################################################################
load("../../00_Datos/se2GSE102741.rda")
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
save(fit2, file = "../../3_Genera_todos_fit/fit_GSE102741.rda")

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

# Los resultados obtenidos están sin ajustar el p valor, por lo que es 
# necesaro ajustarlos mediante BH
x <- p.adjust(res1$table$PValue, "BH")
res1$p.adj <- x
table(res1$p.adj < 0.05)

x <- p.adjust(res2$table$PValue, "BH")
res2$p.adj <- x
res2 <- as.data.frame((res2))
table(res2$p.adj < 0.05)

x <- p.adjust(res3$table$PValue, "BH")
res3$p.adj <- x
res3 <- as.data.frame((res3))
table(res3$p.adj < 0.05)

x <- p.adjust(res4$table$PValue, "BH")
res4$p.adj <- x
table(res4$p.adj < 0.05)
