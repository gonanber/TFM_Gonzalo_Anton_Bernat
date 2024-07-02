###############################################################################
# Autor: Gonzalo Anton Bernat
# Fecha: 04/09/2023
# Descripcion: 
#   Este programa permite descargar los datos de arrays de interés de forma
#   programatica de la base de datos GEO
###############################################################################


# 0. Se cargan las librerías

library(GEOquery)
library(affy)
library(Biobase)
library(geneplotter)

# 1. Se realiza la descarga de los datos

# Se establece el directorio donde se guardan los datos
current_dir <- getwd()
dir.create("0_data")
data_dir <- file.path(current_dir, "0_data")
setwd(data_dir)

# Lista de estudios de GEO que se desean descargar
gse_list <- c("GSE123302",
              "GSE89594",
              "GSE87847",
              "GSE18123",
              "GSE6575",
              "GSE28521",
              "GSE26415")

# Descarga de los estudios
GSE <- lapply(c(gse_list), function(x) getGEO(x, destdir = ".")[[1]])

# Se asignan los nombres de cada elemento de la lista
names(GSE) <- gse_list 


# 2. Se guardan los datos

setwd(data_dir)

# Se guarda cada estudio descargado
GSE123302 <- GSE$GSE123302
save(GSE123302, file="GSE123302.rda")

GSE89594 <- GSE$GSE89594
save(GSE89594, file="GSE89594.rda")

GSE87847 <- GSE$GSE87847
save(GSE87847, file="GSE87847.rda")

GSE18123 <- GSE$GSE18123 # Este estudio está formado por 2 plataformas
save(GSE18123, file="GSE18123_1.rda")

GSE6575 <- GSE$GSE6575
save(GSE6575, file="GSE6575.rda")

GSE28521 <- GSE$GSE28521
save(GSE28521, file="GSE28521.rda")

GSE26415 <- GSE$GSE26415
save(GSE26415, file="GSE26415.rda")

# 3. Se descargan los datos de la segunda plataforma de un estudio
gse <- getGEO("GSE18123")
GSE18123 <- gse[[2]]
save(GSE18123, file="GSE18123_2.rda")

