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


###############################################################################
# Se prepara la ubicacion de los datos
###############################################################################

dir <- getwd()
dir_datos <- file.path((dirname(dir)), "00_datos")
setwd(dir_datos)


###############################################################################
# Se descargan y guardan los datos
###############################################################################

# Esta descarga puede fallar debido al tamaño de los ficheros o a que se ha excedido
# el tiempo de conexion a la base de datos. En este caso se tiene que descargar
# manualmente el ExpressionSet del estudio

# Si falla, se tiene que convertir la matriz de expresión en un archivo.rda y guardarlo
# en la carpeta 00_datos

# Descargar de https://ftp.ncbi.nlm.nih.gov/geo/series/GSE123nnn/GSE123302/matrix/GSE123302_series_matrix.txt.gz
GSE123302 <- getGEO("GSE123302", GSEMatrix =TRUE, getGPL=FALSE)
GSE123302 <- GSE$GSE123302
save(GSE123302, file="GSE123302.rda")

# Descargar de https://ftp.ncbi.nlm.nih.gov/geo/series/GSE89nnn/GSE89594/matrix/GSE89594_series_matrix.txt.gz
GSE89594 <- getGEO("GSE89594", GSEMatrix =TRUE, getGPL=FALSE)
GSE89594 <- GSE$GSE89594
save(GSE89594, file="GSE89594.rda")

# Descargar de https://ftp.ncbi.nlm.nih.gov/geo/series/GSE87nnn/GSE87847/matrix/GSE87847_series_matrix.txt.gz
GSE87847 <- getGEO("GSE87847", GSEMatrix =TRUE, getGPL=FALSE)
GSE87847 <- GSE$GSE87847
save(GSE87847, file="GSE87847.rda")

# Descargar de https://ftp.ncbi.nlm.nih.gov/geo/series/GSE6nnn/GSE6575/matrix/GSE6575_series_matrix.txt.gz
GSE6575 <- getGEO("GSE6575", GSEMatrix =TRUE, getGPL=FALSE)
GSE6575 <- GSE$GSE6575
save(GSE6575, file="GSE6575.rda")

# Descargar de https://ftp.ncbi.nlm.nih.gov/geo/series/GSE28nnn/GSE28521/matrix/GSE28521_series_matrix.txt.gz
GSE28521 <- getGEO("GSE28521", GSEMatrix =TRUE, getGPL=FALSE)
GSE28521 <- GSE$GSE28521
save(GSE28521, file="GSE28521.rda")

# Descargar de https://ftp.ncbi.nlm.nih.gov/geo/series/GSE18nnn/GSE18123/matrix/GSE18123-GPL6244_series_matrix.txt.gz
gse <- getGEO("GSE18123", GSEMatrix =TRUE, getGPL=FALSE)
GSE18123 <- gse[[2]] # Solo nos quedamos con la segunda plataforma
save(GSE18123, file="GSE18123_2.rda")

