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

###############################################################################
# CARGA DE DATOS 
###############################################################################

load("../../00_Datos/GSE6575.rda")
class(GSE6575)
dim(exprs(GSE6575))

# Los datos no están normalizados con log

datos <- pData(GSE6575)
exp <- exprs(GSE6575)
min(exp)
max(exp)

# Se aplica el log2 a los datos de expresión
y <- exprs(GSE6575)
a <- y + abs(min(exprs(GSE6575))) + 1
exp <- log2(a)


###############################################################################
# IMPUTACIÓN DEL SEXO
###############################################################################

## Se seleccionan solo las sondas del cromosoma Y (Bash)
# awk '{print $1}' sondas_GSE123302 > sondas_limpias_GSE123302


## Se cargan las sondas del cromosoma Y para imputar el sexo
data("y.probes")
sondas_y <- y.probes$affy_hg_u133_plus_2


## 2. Imputación del sexo con massiR

# Se calcula la variacion de la expresion de las sondas
exprs(GSE6575) <- exp
massi_y_out <- massi_y(expression_data=GSE6575, y_probes=sondas_y)

# Se representa la variacion de las sondas del cromosoma Y
massi_y_plot(massi_y_out)

# Se extrae la informacion adecuada para realizar un clustering
massi_select_out <- massi_select(expression_data=GSE6575, y_probes=sondas_y, threshold=4)

# Se hace un clustering que predice el sexo de cada muestra
massi_cluster_out <- massi_cluster(massi_select_out)

# Se obtiene el sexo de todos los resultados
resultados <- data.frame(massi_cluster_out[[2]])

# Se adiciona el sexo imputado a los metadatos
datos$sexo_imputado <- resultados$sex
prueba <- data.frame(datos$`characteristics_ch1`, datos$sexo_imputado)

prueba$genero <- sub(".*, ", "", prueba$datos.characteristics_ch1)
prueba$genero <-tolower(prueba$genero)

# Se calcula el número de veces que los valores son iguales
iguales <- sum(prueba$genero == prueba$datos.sexo_imputado)

# Se calcula el número total de filas
total_filas <- nrow(datos)

# Se calcula el porcentaje
(iguales / total_filas) * 100

###############################################################################
# CREACIÓN DE LA COLUMNA SEXO DIAGNOSIS
###############################################################################

# Se extraen los metadatos y se crea una columna nueva con las condiciones de interes
datos$sex_diagnosis <- ifelse(
  grepl("Autism.*regression, Male", datos$`characteristics_ch1`), "male_ASD",
  ifelse(
    grepl("Autism.*regression, Female", datos$`characteristics_ch1`), "female_ASD",
    ifelse(
      grepl("General population, Male", datos$`characteristics_ch1`), "male_control",
      ifelse(
        grepl("General population, Female", datos$`characteristics_ch1`), "female_control",
        ifelse(
          grepl("mental retardation or developmental delay, Male", datos$`characteristics_ch1`), "male_pseudocontrol",
          ifelse(
            grepl("mental retardation or developmental delay, Female", datos$`characteristics_ch1`), "female_pseudocontrol",
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

# Se eliminan las filas con "WS" en la columna "sex_diagnosis"
datos <- subset(datos, sex_diagnosis != "male_pseudocontrol" & sex_diagnosis != "female_pseudocontrol")
muestras <- rownames(datos)
exp <- subset(exp, select = muestras)


# Se cuentan las diferentes categorías de sex_diagnosis
nrow(datos) # Totales
table(datos$sex_diagnosis) # Por categoria

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
# Construcción del nuevo Expression Set
###############################################################################

datos <- AnnotatedDataFrame(datos)
info <- experimentData(GSE6575)
anotacion <- fData(GSE6575)
anotacion <- data.frame(anotacion$`Gene Symbol`, anotacion$ID)
anotacion_filtrada <- AnnotatedDataFrame(anotacion[rownames(fData(GSE6575)) %in% rownames(exp), ])

eset2 <- ExpressionSet(assayData = exp,
                       phenoData = datos,
                       experimentData = info,
                       featureData = anotacion_filtrada)

save(eset2, file = "../../00_Datos/eset2_GSE6575.rda")


################################################################################
# Calcular medianas de genes repetidos
################################################################################

medianReps <- function(matriz){
  ID <- as.character(rownames(matriz))
  ID <- factor(ID, levels = unique(ID))
  df <- by(matriz, ID, function(x) apply(x, 2, stats::median))
  mat <- do.call("rbind", df)
  return(mat)
}

exp2 <- medianReps(exp)



###############################################################################
# ANALISIS DE EXPRESION DIFERENCIAL
###############################################################################

load("../../00_Datos/eset2_GSE6575.rda")
exp <- exprs(eset2)
datos <- pData(eset2)
info <- experimentData(eset2)

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
save(fit2, file = "../../3_Genera_todos_fit/fit_GSE6575.rda")


# La función limma::topTable() realiza un resumen de los ajustes, y ordena las 
# sondas por su p-valor ajustado.
resultados <- topTable(fit2,
                       coef=3,
                       adjust="BH",
                       number = nrow(exp))
head(resultados)
table(resultados[, "adj.P.Val"] < 0.05)
resultados_interes <- subset(resultados, adj.P.Val < 0.05)
ups <- subset(resultados_interes, logFC > 0)
dim(ups)
downs <- subset(resultados_interes, logFC < 0)
dim(downs)
datatable(resultados_interes)





