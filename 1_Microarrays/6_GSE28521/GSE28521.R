###############################################################################
# Autor: Gonzalo Anton Bernat
# Fecha: 03/09/2023
# Descripcion: 
#   Este programa permite realizar un análisis exploratorio de los
#   datos previamente descargados y el análisis de expresión diferencial.
###############################################################################


## 0. Se cargan las librerías

library(GEOquery)
library(affy)
library(Biobase)
library(geneplotter)
library(ggplot2)
library(reshape2)
library(limma)
library(statmod)
library(ggfortify)
library(factoextra)
library(massiR)
library(biomaRt)
library(org.Hs.eg.db)
library(readxl)


###############################################################################
# CARGA DE DATOS 
###############################################################################

load("../0_data/GSE28521.rda")
class(GSE28521)
dim(exprs(GSE28521))
annotation(GSE28521)
# Los datos no están normalizados con log

datos <- pData(GSE28521)
exp <- exprs(GSE28521)

# Se añade la información de sexo encontrada en el material suplementario de la 
# publicación
info <- read_xlsx("sondas.xlsx")
merged_data <- merge(info, datos, by.x = "GEO_ID", by.y = "geo_accession")
datos <- merged_data


# Se sustituyen los NA por 0
sum(is.na(exp))
exp[is.na(exp)] <- 0
sum(is.na(exp))

y <- exp
a <- y + abs(min(exp) + 1)
exp <- log2(a)
min(exp)
max(exp)


###############################################################################
# IMPUTACIÓN DEL SEXO
###############################################################################

## Se seleccionan solo las sondas del cromosoma Y (Bash)
# awk '{print $1}' sondas_GSE123302 > sondas_limpias_GSE123302

sondas_y <- read.table("sondas_y.txt", sep = "\t") # Nombres de las sondas
sondas_sexo <- data.frame(row.names = sondas_y$V1) # MassiR quiere un dataframe con los nombres de las filas


## 2. Imputación del sexo con massiR

# Se calcula la variacion de la expresion de las sondas
exprs(GSE28521) <- exp
massi_y_out <- massi_y(expression_data=GSE28521, y_probes=sondas_sexo)

# Se representa la variacion de las sondas del cromosoma Y
massi_y_plot(massi_y_out)

# Se extrae la informacion adecuada para realizar un clustering
massi_select_out <- massi_select(expression_data=GSE28521, y_probes=sondas_sexo, threshold=4)

# Se hace un clustering que predice el sexo de cada muestra
massi_cluster_out <- massi_cluster(massi_select_out)

# Se obtiene el sexo de todos los resultados
resultados <- data.frame(massi_cluster_out[[2]])

# Se adiciona el sexo imputado a los metadatos
datos$sexo_imputado <- resultados$sex
prueba <- data.frame(datos$SEX, datos$sexo_imputado)

prueba$sexo_detallado <- ifelse(prueba$datos.SEX == "M", "male", 
                               ifelse(prueba$datos.SEX == "F", "female", NA))

# Se calcula el número de veces que los valores son iguales
iguales <- sum(prueba$sexo_detallado == prueba$datos.sexo_imputado)

# Se calcula el número total de filas
total_filas <- nrow(datos)

# Se calcula el porcentaje
(iguales / total_filas) * 100


###############################################################################
# CREACIÓN DE LA COLUMNA SEXO DIAGNOSIS
###############################################################################

# Se extraen los metadatos y se crea una columna nueva con las condiciones de interes
datos$sex_diagnosis <- ifelse(
  grepl("M", datos$SEX) & grepl("autism", datos$`characteristics_ch1`), "male_ASD",
  ifelse(
    grepl("M", datos$SEX) & grepl("controls", datos$`characteristics_ch1`), "male_control",
    ifelse(
      grepl("F", datos$SEX) & grepl("autism", datos$`characteristics_ch1`), "female_ASD",
      ifelse(
        grepl("F", datos$SEX) & grepl("controls", datos$`characteristics_ch1`), "female_control",
        NA
      )
    )
  )
)

###############################################################################
# COMPARACIÓN ESTIMADORES DE DENSIDAD
###############################################################################

dir.create("images")
png(filename = "images/hist.png")
geneplotter::multidensity(exp)
dev.off()


###############################################################################
# BOXPLOT
###############################################################################

# Se convierte la matriz de expresión en un dataframe
df_boxplot = reshape2::melt(exp) 
colnames(df_boxplot) = c("Sonda", "Muestra", "Expresion")
head(df_boxplot)

# Se crea otro dataframe para colorear según la condición
grupos = data.frame(muestras = datos[, 1], 
                    grupos = datos[, "sex_diagnosis"])

conjunto = merge(df_boxplot, grupos, by.x = "Muestra", by.y = "muestras")

jpeg(filename = "images/boxplot.jpg",  width = 1600, height = 800, quality = 100)
ggplot(conjunto, aes(x=Muestra, y=Expresion, fill=grupos)) + 
  theme(axis.text.x = element_blank(), 
        plot.title = element_text(size = 26), 
        axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 24), 
        legend.title = element_text(size = 24),
        legend.text = element_text(size = 24)) +
  geom_boxplot() +
  labs(fill = "Grupos") +
  ggtitle("Comparación de la Expresión Génica entre los diferentes grupos experimentales") + 
  xlab("Muestras") + 
  ylab("Nivel de expresión normalizada")
dev.off()



###############################################################################
# PCA
###############################################################################

# En funcion del diagnostico
datos2 <- datos
exp2 <- exp
metadata2 <- datos2$sex_diagnosis
metadata2 <- factor(metadata2)
res.pca <- prcomp(t(exp2))

jpeg(filename = "images/PCA.jpg",  width = 800, height = 400, quality = 100)
result <- fviz_pca_ind(res.pca,
                       habillage = metadata2,
                       geom = "point",
                       pointsize = 4,
                       invisible = "quali")

# Personalizar el gráfico con ggplot2
result <- result +
  ggtitle("Análisis de componentes principales") +
  xlab("PC1") +
  ylab("PC2") +
  labs(color = "Grupos", shape = "Grupos") +
  theme(
    plot.title = element_text(size = 16),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12)
  ) +
  guides(shape = guide_legend(title = "Grupos"), color = guide_legend(title = "Grupos"))


result
dev.off()

# En funcion del estudio
datos2 <- datos
exp2 <- exp
metadata2 <- datos2$`tissue (brain region):ch1`
metadata2 <- factor(metadata2)
res.pca <- prcomp(t(exp2))

jpeg(filename = "images/PCA.jpg",  width = 800, height = 400, quality = 100)
result <- fviz_pca_ind(res.pca,
                       habillage = metadata2,
                       geom = "point",
                       pointsize = 4,
                       invisible = "quali")

# Personalizar el gráfico con ggplot2
result <- result +
  ggtitle("Análisis de componentes principales") +
  xlab("PC1") +
  ylab("PC2") +
  labs(color = "Grupos", shape = "Grupos") +
  theme(
    plot.title = element_text(size = 16),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12)
  ) +
  guides(shape = guide_legend(title = "Grupos"), color = guide_legend(title = "Grupos"))


result
dev.off()



###############################################################################
# CLUSTERING
###############################################################################

# Función para el clustering de genes
treeclust_from_ggplot <- function(matriz1, groups, title = "Clustering (distancia de correlación)", 
                                  subtitle = NULL, bottom = FALSE, dist = "cor", 
                                  save = NULL, width = 800, height = 600) {
  
  require(ggplot2)
  require(ggdendro)
  
  ## Create cluster
  if(dist == "cor"){
    correlacion <- cor(matriz1)
    distancia <- as.dist((1 - correlacion) / 2)
  }else if(dist == "euclid"){
    distancia <- dist(t(data), method = "euclidean") #dist trabaja por filas
  }else{
    stop("Please specify an accepted distance. Options are 'cor' or 'euclid'")
  }
  
  cluster <- hclust(distancia)
  cluster$clase <- groups
  dendr <- ggdendro::dendro_data(cluster, type = "rectangle")
  clases <- as.character(cluster$clase)
  clust.df <- data.frame(label = cluster$labels, Condition = factor(clases))
  dendr[["labels"]] <- merge(dendr[["labels"]], clust.df, by = "label")
  cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  
  g = ggplot() +
    geom_segment(data = ggdendro::segment(dendr), aes(x=x, y=y, xend=xend, yend=yend)) +
    geom_text(data = label(dendr), aes(x, y, label = label, hjust = 0, color = Condition), size = 4.5) +
    coord_flip() + scale_y_reverse(expand = c(0.2, 0)) +
    theme(axis.line.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.text.y=element_blank(),
          axis.title.y=element_blank(),
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          panel.background=element_rect(fill="white"),
          panel.grid=element_blank(),
          title = element_text(size = 10),
          text = element_text(size = 14),
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5)) +
    labs(color = "Grupos") +
    scale_colour_hue() +
    scale_fill_hue() 
  
  if(!is.null(subtitle)){
    # Añadir subtítulo
    g <- g + labs(title = title, subtitle = subtitle)
  }else{
    g <- g + ggtitle(title) + theme()
  }
  
  if (bottom == TRUE) {
    #leyenda bajo el gráfico
    g =  g + theme(legend.position ="bottom", legend.key.size = unit(0.5,"line"), 
                   legend.title=element_text(size=7))
  }
  
  if(!is.null(save)){
    png(filename = save, width = width, height = height)
    print(g)
    dev.off()
  }
  
  return (g)
}

# Visualización del clustering de genes
png(filename = "images/clustering.png")
metadata2 <- datos$sex_diagnosis
metadata2 <- factor(metadata2)
treeclust_from_ggplot(exp, metadata2, title = "Clustering (distancia de correlación)", 
                      subtitle = NULL, bottom = FALSE, dist = "cor", 
                      save = NULL, width = 800, height = 600)
dev.off()





###############################################################################
# SEPARACION DE LOS ESTUDIOS DE CEREBELO Y CORTEX
###############################################################################

# Se separa el estudio Frontal Cortex

datos2 <- subset(datos, characteristics_ch1.1 == "tissue (brain region): Frontal cortex")
rownames(datos2) <- datos2[, 1]
muestras <- data.frame(row.names = datos2$GEO_ID)
datos2 <- AnnotatedDataFrame(datos2)


exp2 <- subset(exp, select = rownames(muestras))
info <- experimentData(GSE28521)
anotacion <- fData(GSE28521)
anotacion_filtrada <- AnnotatedDataFrame(anotacion[rownames(fData(GSE28521)) %in% rownames(exp2), ])

a <- data.frame(colnames(exp2), rownames(datos2))

GSE28521 <- ExpressionSet(assayData = exp2,
                                   phenoData = datos2,
                                   experimentData = info,
                                   featureData = anotacion_filtrada)

save(GSE28521, file = "../0_data/GSE28521_frontal_cortex.rda")

# Se separa el estudio Temporal Cortex

datos2 <- subset(datos, characteristics_ch1.1 == "tissue (brain region): Temporal cortex")
rownames(datos2) <- datos2[, 1]
muestras <- data.frame(row.names = datos2$GEO_ID)
datos2 <- AnnotatedDataFrame(datos2)


exp2 <- subset(exp, select = rownames(muestras))
info <- experimentData(GSE28521)
anotacion <- fData(GSE28521)
anotacion_filtrada <- AnnotatedDataFrame(anotacion[rownames(fData(GSE28521)) %in% rownames(exp2), ])

a <- data.frame(colnames(exp2), rownames(datos2))

GSE28521 <- ExpressionSet(assayData = exp2,
                          phenoData = datos2,
                          experimentData = info,
                          featureData = anotacion_filtrada)

save(GSE28521, file = "../0_data/GSE28521_temporal_cortex.rda")


# Se separa el estudio Cerebellum

datos2 <- subset(datos, characteristics_ch1.1 == "tissue (brain region): Cerebellum")
rownames(datos2) <- datos2[, 1]
muestras <- data.frame(row.names = datos2$GEO_ID)
datos2 <- AnnotatedDataFrame(datos2)


exp2 <- subset(exp, select = rownames(muestras))
info <- experimentData(GSE28521)
anotacion <- fData(GSE28521)
anotacion_filtrada <- AnnotatedDataFrame(anotacion[rownames(fData(GSE28521)) %in% rownames(exp2), ])

a <- data.frame(colnames(exp2), rownames(datos2))

GSE28521 <- ExpressionSet(assayData = exp2,
                          phenoData = datos2,
                          experimentData = info,
                          featureData = anotacion_filtrada)

save(GSE28521, file = "../0_data/GSE28521_cerebellum.rda")




###############################################################################
# ANALISIS DE EXPRESION DIFERENCIAL
###############################################################################

load("../0_data/eset2_GSE6575.rda")
exp <- exprs(eset2)
datos <- pData(eset2)

# Se crea la matriz de diseño
condicion <- datos$sex_diagnosis
design <- model.matrix(~0 + condicion)
columnas <- levels(factor(datos$sex_diagnosis))
colnames(design)= columnas
head(design, n=200)

# Se ajusta la estimación de los niveles de expresión de cada gen a un modelo 
# lineal teniendo en cuenta el diseño experimental
fit <- limma::lmFit(eset2, design)


# Para determinar los genes expresados de forma diferencial debemos especificar 
# los contrastes que se van a considerar.
cont.matrix <- limma::makeContrasts( 
  dif1 = (female_ASD) - (female_control), # Contraste mujeres
  dif2 = (male_ASD) - (male_control),  # Contraste hombres
  dif3 = (male_ASD + female_ASD) - (male_control + female_control), # Enfermedad sin sexo
  dif4 = (female_ASD - female_control) - (male_ASD - male_control), # Diferencias de sexo
  
  levels=design
)

cont.matrix # 1 es caso -1 es control,


# Se realiza el ajuste de contrastes
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)


# La función limma::topTable() realiza un resumen de los ajustes, y ordena las 
# sondas por su p-valor ajustado.
resultados <- topTable(fit2,
                       coef=4,
                       adjust="BH",
                       number = nrow(exp))
head(resultados)
table(resultados[, "adj.P.Val"] < 0.05)
resultados_interes <- subset(resultados, adj.P.Val < 0.05)
datatable(resultados_interes)





