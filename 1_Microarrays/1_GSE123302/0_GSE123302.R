###############################################################################
# Autor: Gonzalo Anton Bernat
# Fecha: 03/09/2023
# Descripcion: 
#   Este programa permite realizar un análisis exploratorio de los
#   datos previamente descargados y el análisis de expresión diferencial
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


###############################################################################
# CARGA DE DATOS
###############################################################################

load("../../00_Datos/GSE123302.rda")
class(GSE123302)
dim(exprs(GSE123302))

datos <- pData(GSE123302)
exp <- exprs(GSE123302)
min(exp)
max(exp) # Los datos están normalizados


###############################################################################
# IMPUTACIÓN DEL SEXO
###############################################################################

## Se seleccionan solo las sondas del cromosoma Y (Bash)
# awk '{print $1}' sondas_GSE123302 > sondas_limpias_GSE123302

## Se cargan las sondas del cromosoma Y para imputar el sexo
sondas_y <- read.table("sondas_y.txt", sep = "\t") # Nombres de las sondas
sondas_sexo <- data.frame(row.names = sondas_y$V1) # MassiR quiere un dataframe con los nombres de las filas

# Se calcula la variacion de la expresion de las sondas

massi_y_out <- massi_y(expression_data=GSE123302, y_probes=sondas_sexo) # En la documentacion indica que acepta un expression set

# Se representa la variacion de las sondas del cromosoma Y
massi_y_plot(massi_y_out)

# Se extrae la informacion adecuada para realizar un clustering
massi_select_out <- massi_select(expression_data=GSE123302, y_probes=sondas_sexo, threshold=4)

# Se hace un clustering que predice el sexo de cada muestra
massi_cluster_out <- massi_cluster(massi_select_out)

# Se obtiene el sexo de todos los resultados
resultados <- data.frame(massi_cluster_out[[2]])

# Se adiciona el sexo imputado a los metadatos
datos$sexo_imputado <- resultados$sex



###############################################################################
# CREACIÓN DE LA COLUMNA SEXO DIAGNOSIS
###############################################################################

# Se reconoce un patron y se asigna el nombre correspondiente a cada muestra
datos$sex_diagnosis <- ifelse(
  grepl("^male", datos$`Sex:ch1`) & grepl("^ASD", datos$`diagnosis:ch1`), "male_ASD",
  ifelse(
    grepl("^female", datos$`Sex:ch1`) & grepl("^ASD", datos$`diagnosis:ch1`), "female_ASD",
    ifelse(
      grepl("^male", datos$`Sex:ch1`) & grepl("^TD", datos$`diagnosis:ch1`), "male_control",
      ifelse(
        grepl("^female", datos$`Sex:ch1`) & grepl("^TD", datos$`diagnosis:ch1`), "female_control",
        ifelse(
          grepl("^male", datos$`Sex:ch1`) & grepl("^Non-TD", datos$`diagnosis:ch1`), "male_pseudocontrol",
          ifelse(
            grepl("^female", datos$`Sex:ch1`) & grepl("^Non-TD", datos$`diagnosis:ch1`), "female_pseudocontrol",
            NA
          )
        )
      )
    )
  )
)


###############################################################################
# SELECCIÓN DE LOS DATOS DE INTERÉS (CASO VS CONTROL)
###############################################################################

# Se eliminan las filas con "pseudocontrol" en la columna "sex_diagnosis"
datos <- subset(datos, sex_diagnosis != "male_pseudocontrol" & sex_diagnosis != "female_pseudocontrol")
muestras <- rownames(datos)
exp <- subset(exp, select = muestras)


###############################################################################
# COMPARACIÓN ESTIMADORES DE DENSIDAD
###############################################################################

dir.create("images")
png(filename = "images/GSE123302_hist.png")
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
grupos = data.frame(muestras = rownames(datos[, 0]), 
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
metadata2 <- datos2$characteristics_ch1.2
metadata2 <- factor(metadata2)
res.pca <- prcomp(t(exp2))

jpeg(filename = "images/PCA_estudios.jpg",  width = 800, height = 400, quality = 100)
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
# SEPARACION DE LOS ESTUDIOS MARBLES Y EARLI
###############################################################################

# Se separa el estudio MARBLES
datos_marbles <- AnnotatedDataFrame(subset(datos, characteristics_ch1.2 == "study: MARBLES" ))
muestras_marbles <- rownames(datos_marbles)
exp_marbles <- subset(exp, select = muestras_marbles)
info_marbles <- experimentData(GSE123302)
anotacion <- fData(GSE123302)
anotacion_filtrada <- AnnotatedDataFrame(anotacion[rownames(fData(GSE123302)) %in% rownames(exp_marbles), ])

GSE123302_marbles <- ExpressionSet(assayData = exp_marbles,
                                   phenoData = datos_marbles,
                                   experimentData = info_marbles,
                                   featureData = anotacion_filtrada)

save(GSE123302_marbles, file = "../../00_Datos/GSE123302_marbles.rda")

# Se separa el estudio EARLI

datos_earli <- AnnotatedDataFrame(subset(datos, characteristics_ch1.2 == "study: EARLI" ))
muestras_earli <- rownames(datos_earli)
exp_earli <- subset(exp, select = muestras_earli)
info_earli <- experimentData(GSE123302)
anotacion <- fData(GSE123302)
anotacion_filtrada2 <- AnnotatedDataFrame(anotacion[rownames(fData(GSE123302)) %in% rownames(exp_earli), ])

GSE123302_earli <- ExpressionSet(assayData = exp_earli,
                                 phenoData = datos_earli,
                                 experimentData = info_earli,
                                 featureData = anotacion_filtrada2)

save(GSE123302_earli, file = "../../00_Datos/GSE123302_earli.rda")


