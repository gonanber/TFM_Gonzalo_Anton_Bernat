###############################################################################
# Autor: Gonzalo Anton Bernat
# Fecha: 20/03/2024
# Descripcion: 
#   Este programa permite generar el objeto ALLfit, necesario para realizar el
#   metaanálisis. En este objeto se añaden todos los fit2 de los diferentes
#   estudios
###############################################################################

dir <- getwd()
dir2 <- paste0(dir, "/fits")
setwd(dir2)
archivos <- list.files(dir2)

# Crear una lista para almacenar los objetos cargados
objetos_cargados <- list()

# Obtener la lista de archivos en la carpeta actual
archivos <- list.files()

# Crear una lista para almacenar los objetos cargados
todos_fit2 <- list()


# Recorrer cada archivo
for (archivo in archivos) {
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

