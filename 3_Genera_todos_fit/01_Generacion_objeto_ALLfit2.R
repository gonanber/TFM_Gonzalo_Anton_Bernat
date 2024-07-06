###############################################################################
# Autor: Gonzalo Anton Bernat
# Fecha: 20/03/2024
# Descripcion: 
#   Este programa permite generar el objeto ALLfit, necesario para realizar el
#   metaanálisis. En este objeto se añaden todos los fit2 de los diferentes
#   estudios
###############################################################################

dir <- getwd()

# Crear una lista para almacenar los objetos cargados
objetos_cargados <- list()

# Obtener la lista de archivos en la carpeta actual
archivos_todos <- list.files()

archivos_cerebro <- c("fit_GSE285211.rda", 
                      "fit_GSE285212.rda", 
                      "fit_GSE285213.rda",
                      "fit_GSE102741.rda"
                      )

archivos_sangre <- c("fit_GSE1233021.rda",
                     "fit_GSE1233022.rda",
                     "fit_GSE89594.rda",
                     "fit_GSE87847.rda",
                     "fit_GSE18123_2.rda",
                     "fit_GSE6575.rda",
                     "fit_GSE212645.rda",
                     "fit_GSE140702.rda")


###############################################################################
# TODOS LOS ESTUDIOS
###############################################################################
# Crear una lista para almacenar los objetos cargados
todos_fit2 <- list()


# Recorrer cada archivo
for (archivo in archivos_todos) {
  # Verificar si es un archivo de datos que se puede cargar
  if (endsWith(archivo, ".rda")) {
    # Cargar el archivo RDA
    objeto <- load(archivo)
    nombre <- gsub("\\.rda$", "", archivo)
    # Agregar el objeto a la lista
    todos_fit2[[nombre]] <- get(objeto)
  }
}

save(todos_fit2, file = "todos_fit2.RData")


###############################################################################
# ESTUDIOS SANGRE
###############################################################################

# Crear una lista para almacenar los objetos cargados
todos_fit2 <- list()


# Recorrer cada archivo
for (archivo in archivos_sangre) {
  # Verificar si es un archivo de datos que se puede cargar
  if (endsWith(archivo, ".rda")) {
    # Cargar el archivo RDA
    objeto <- load(archivo)
    nombre <- gsub("\\.rda$", "", archivo)
    # Agregar el objeto a la lista
    todos_fit2[[nombre]] <- get(objeto)
  }
}

save(todos_fit2, file = "sangre_fit2.RData")


###############################################################################
# ESTUDIOS CEREBRO
###############################################################################

# Crear una lista para almacenar los objetos cargados
todos_fit2 <- list()


# Recorrer cada archivo
for (archivo in archivos_cerebro) {
  # Verificar si es un archivo de datos que se puede cargar
  if (endsWith(archivo, ".rda")) {
    # Cargar el archivo RDA
    objeto <- load(archivo)
    nombre <- gsub("\\.rda$", "", archivo)
    # Agregar el objeto a la lista
    todos_fit2[[nombre]] <- get(objeto)
  }
}

save(todos_fit2, file = "cerebro_fit2.RData")

