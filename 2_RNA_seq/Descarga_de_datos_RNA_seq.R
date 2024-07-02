###############################################################################
# Autor: Gonzalo Anton Bernat
# Fecha: 04/09/2023
# Descripcion: 
#   Este programa permite descargar los datos de estudios de RNA-seq de interés
#   de forma programatica de la base de datos GEO. Los dos primeros estudios se
#   han cargado a partir de tablas proporcionadas por los autores
###############################################################################


# 0. Se cargan las librerías

library(GEOquery)
library(affy)
library(Biobase)
library(geneplotter)
library(SummarizedExperiment)
library(data.table)

# 1. Se realiza la descarga de los datos

# Se establece el directorio donde se guardan los datos
current_dir <- getwd()
dir.create("0_data")
data_dir <- file.path(current_dir, "0_data")
setwd(data_dir)

# Lista de estudios de GEO que se desean descargar
gse_list <- c("GSE212645", "GSE140702", "GSE102741")

# Descarga de los metadatos de cada estudio
metadata <- lapply(c(gse_list), function(x) getGEO(x, destdir = ".")[[1]])

# Descarga de la matriz de conteos de cada estudio
GSE <- lapply(c(gse_list), function(x) getGEOSuppFiles(x, makeDirectory = FALSE))

# Se asignan los nombres de cada elemento de la lista
names(GSE) <- gse_list 


# 2. Se construye un Summarized experiment

################################# GSE212645 ####################################
data_212645 <- read.delim("GSE212645_CountsMatrix.txt.gz")
metadata_GSE21645 <- (metadata[1])

# Se crea la amtriz de conteos en el formato adecuado
count_matrix <- data.frame(data_212645[, -1], row.names = data_212645$GeneID)

# Se crea el objeto SummarizedExperiment
seGSE212645 <- SummarizedExperiment(assays = list(counts = as.matrix(count_matrix)))

genes <- DataFrame(metadata_GSE21645)

# Asignar los metadatos al SummarizedExperiment
colData(seGSE212645) <- genes


########################## GSE212645 sin normalizar ############################
data_212645 <- read.delim("0_data/counts_GSE212645.txt")
metadata_GSE21645 <- read.delim("0_data/metadata_GSE212645.txt")

# Se crea la amtriz de conteos en el formato adecuado
count_matrix <- data.frame(data_212645[, -1], row.names = data_212645$GeneID)

# Se crea el objeto SummarizedExperiment
seGSE212645_2 <- SummarizedExperiment(assays = list(counts = as.matrix(count_matrix)))

genes <- DataFrame(metadata_GSE21645)
class(metadata_GSE21645)

# Asignar los metadatos al SummarizedExperiment
colData(seGSE212645_2) <- genes


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

# 3. Se guardan los Summarized Experiments

setwd(current_dir)
setwd("../")

# Se guarda cada estudio descargado
save(seGSE212645, file="GSE212645.rda")
save(seGSE212645_2, file="0_data/GSE212645_2_raw.rda")
save(seGSE140702, file="0_data/GSE140702.rda")
save(seGSE102741, file="GSE102741.rda")



