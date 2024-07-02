###############################################################################
# Autor: Gonzalo Anton Bernat
# Fecha: 03/09/2023
# Descripcion: 
#   Este programa permite realizar un análisis exploratorio de los
#   datos previamente descargados.
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
# CARGA DE DATOS Y ANÁLISIS EXPLORATORIO
###############################################################################

load("../0_data/GSE123302_marbles.rda")
eset <- GSE123302_marbles
class(eset)
dim(exprs(eset)) # Los datos están normalizados y con log

datos <- pData(eset)
exp <- exprs(eset)
dim(exp)
dim(datos)
min(exp)
max(exp)


##  Comparacion de los estimadores de densidad

dir.create("images_marbles")
png(filename = "images_marbles/GSE123302_hist.png")
geneplotter::multidensity(exp)
dev.off()



###############################################################################
# BOXPLOT
###############################################################################

df_boxplot = reshape2::melt(exp) # Se convierte la matriz de expresión en un dataframe
colnames(df_boxplot) = c("Sonda", "Muestra", "Expresion")
head(df_boxplot)

grupos = data.frame(muestras = rownames(datos[, 0]), # Se crea otro dataframe para colorear según la condición
                    grupos = datos[, "sex_diagnosis"])

conjunto = merge(df_boxplot, grupos, by.x = "Muestra", by.y = "muestras") # Se fusionan ambos dataframes

jpeg(filename = "images_marbles/boxplot.jpg",  width = 1600, height = 800, quality = 100)
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

jpeg(filename = "images_marbles/PCA.jpg",  width = 800, height = 400, quality = 100)
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
png(filename = "images_marbles/clustering.png")
metadata2 <- datos$sex_diagnosis
metadata2 <- factor(metadata2)
treeclust_from_ggplot(exp, metadata2, title = "Clustering (distancia de correlación)", 
                      subtitle = NULL, bottom = FALSE, dist = "cor", 
                      save = NULL, width = 800, height = 600)
dev.off()

###############################################################################
# Construcción del nuevo Expression Set
###############################################################################

condiciones <- colnames(exp)
muestras <- setdiff(condiciones, outliers)

exp2 <- subset(exp, select = muestras)
datos2 <- AnnotatedDataFrame(datos[!(rownames(datos) %in% outliers), ])
info <- experimentData(eset)
anotacion <- fData(eset)
anotacion_filtrada <- AnnotatedDataFrame(anotacion[rownames(fData(eset)) %in% rownames(exp2), ])

eset2 <- ExpressionSet(assayData = exp2,
                       phenoData = datos2,
                       experimentData = info,
                       featureData = anotacion_filtrada)

save(eset2, file = "../0_data/eset2_marbles.rda")


###############################################################################
# ANALISIS DE EXPRESION DIFERENCIAL
###############################################################################

load("../0_data/eset2_marbles.rda")
eset <- eset2
datos <- pData(eset)
exp <- exprs(eset)
anot <- fData(eset)
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
values = anot$GB_ACC
genes_buenos <- getBM(c("refseq_mrna", "external_gene_name"), "refseq_mrna", values,
                      mart=mart)
colnames(genes_buenos) <- c("GB_ACC", "Symbol")
save(genes_buenos, file = "../anotacion/marbles.rda")



exp <- exprs(eset2)
datos <- pData(eset2)
info <- experimentData(eset2)

genedata <- fData(eset2)
antes <- genedata
genes_buenos_unique <- unique(genes_buenos)
genes_buenos_unique <- unique(genes_buenos_unique)

# Realizar el merge manteniendo todas las filas de genedata
genedata <- merge(genedata, genes_buenos_unique, by = "GB_ACC", all.x = TRUE, no.dups = TRUE)

genedata <- genedata[,c("Symbol","ID")]
genedata$Symbol <- sub(" ///.*","",genedata$Symbol)
genedata[genedata == ""] <- "eliminar"

diff <- length(rownames(genedata)) - length(rownames(exp))
genedata <- genedata[!(is.na(genedata$Symbol) & seq_len(nrow(genedata)) <= diff), ]



medianReps <- function(matriz){
  ID <- as.character(rownames(matriz))
  ID <- factor(ID, levels = unique(ID))
  df <- by(matriz, ID, function(x) apply(x, 2, stats::median))
  mat <- do.call("rbind", df)
  return(mat)
}

exp2 <- exp
rownames(exp2) <- genedata$Symbol
exp2 <- medianReps(exp2)
exp3 <- exp2[rownames(exp2) != "eliminar",]
genedata <- genedata[match(rownames(exp3),genedata$Symbol),c("Symbol", "ID" )]
rownames(genedata) <- genedata$Symbol


eset2 <- ExpressionSet(assayData = exp3,
                       phenoData = AnnotatedDataFrame(datos),
                       experimentData = info,
                       featureData = AnnotatedDataFrame(data.frame(genedata)))


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
#save(fit2, file = "../fits/GSE1233022.rda")

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



