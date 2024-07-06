###############################################################################
# Autor: Gonzalo Anton Bernat
# Fecha: 04/09/2023
# Descripcion: 
#   Este programa permite descargar los datos de estudios de RNA-seq de interés
#   de forma programatica de la base de datos GEO.
###############################################################################


# 0. Se cargan las librerías

library(GEOquery)
library(affy)
library(Biobase)
library(geneplotter)
library(SummarizedExperiment)


###############################################################################
# Se prepara la ubicacion de los datos
###############################################################################

dir <- getwd()
dir_datos <- file.path((dirname(dir)), "00_datos")
setwd(dir_datos)


###############################################################################
# Se descargan los datos
###############################################################################

# Esta descarga puede fallar debido al tamaño de los ficheros o a que se ha excedido
# el tiempo de conexion a la base de datos. En este caso se tiene que descargar
# manualmente el ExpressionSet del estudio

# Si falla, se tiene que convertir la matriz de conteos en un SummarizedExperiemnt
# y guardarlo en la carpeta 00_datos

# Lista de estudios de GEO que se desean descargar
gse_list <- c("GSE212645", "GSE140702", "GSE102741")

# Descarga de los metadatos de cada estudio
metadata <- lapply(c(gse_list), function(x) getGEO(x, destdir = ".")[[1]])

# Descarga de la matriz de conteos de cada estudio
GSE <- lapply(c(gse_list), function(x) getGEOSuppFiles(x, makeDirectory = FALSE))

# Se asignan los nombres de cada elemento de la lista
names(GSE) <- gse_list 


################################# GSE212645 ####################################

# Descargar de https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE212645&format=file&file=GSE212645%5FCountsMatrix%2Etxt%2Egz
data_212645 <- read.delim("GSE212645_CountsMatrix.txt.gz") 
metadata_GSE21645 <- (metadata[1])

# Se crea la matriz de conteos en el formato adecuado
count_matrix <- data.frame(data_212645[, -1], row.names = data_212645$GeneID)

# Se crea el objeto SummarizedExperiment
seGSE212645 <- SummarizedExperiment(assays = list(counts = as.matrix(count_matrix)))

genes <- DataFrame(metadata_GSE21645)

# Asignar los metadatos al SummarizedExperiment
colData(seGSE212645) <- genes

################################# GSE140702 #####################################

eset <- getGEO("GSE140702")
gse <- eset$GSE140702_series_matrix.txt.gz
metadata <- pData(gse)
metadatos <- DataFrame(metadata)


# Se carga la matriz de conteos desde GEO
url <- "https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts"
path <- paste(url, "acc=GSE140702", "file=GSE140702_raw_counts_GRCh38.p13_NCBI.tsv.gz", sep="&");
conteos <- as.matrix(data.table::fread(path, header=T, colClasses="integer"), rownames=1)

# Se crea el objeto SummarizedExperiment
seGSE140702 <- SummarizedExperiment(assays = list(conteos),
                                    colData = metadatos)


################################# GSE102741 ####################################

# Descargar de https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE102741&format=file&file=GSE102741%5Flog2RPKMcounts%2Erda%2Egz
load("0_data/GSE102741_log2RPKMcounts.rda.gz")
class(geneRpkm2)
data_102741 <- geneRpkm2
a <- table(rowSums(data_102741))>1
table(a)
metadata <- getGEO("GSE102741")[1]
metadata <- metadata$GSE102741_series_matrix.txt.gz


# Se crea la amtriz de conteos en el formato adecuado
count_matrix <- data.frame(data_102741)

# Se crea el objeto SummarizedExperiment
seGSE102741 <- SummarizedExperiment(assays = list(counts = as.matrix(count_matrix)))

genes <- DataFrame(metadata)

# Asignar los metadatos al SummarizedExperiment
colData(seGSE102741) <- genes
colData(seGSE102741)


###############################################################################
# Se guardan los datos
###############################################################################

save(seGSE212645, file="GSE212645.rda")
save(seGSE140702, file="GSE140702.rda")
save(seGSE102741, file="GSE102741.rda")



