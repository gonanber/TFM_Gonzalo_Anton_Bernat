###############################################################################
# Autor: Gonzalo Anton Bernat
# Fecha: 10/05/2024
# Descripcion: 
#   Este programa permite generar las listas que se utilizan en STRING para 
#   analizar las interacciones proteína-proteína. Los genes aparecen por consola
#   hay que copiarlos y pegarlos en STRING
###############################################################################


###############################################################################
# Carga de datos
###############################################################################


# Se leen los genes expresados diferencialmente  en todas las muestras de sangre y cerebro
all_females <- read.table(file = "Results/MA/Contrastes_alfa_0.05/1_MA_todos/1_Contraste_mujeres/sig.genes.df.txt",
                     header = T, 
                     sep = "\t")


all_ASD_control <- read.table(file = "Results/MA/Contrastes_alfa_0.05/1_MA_todos/3_Contraste_sin_sexo/sig.genes.df.txt",
                     header = T, 
                     sep = "\t")


all_ASD_sex <- read.table(file = "Results/MA/Contrastes_alfa_0.05/1_MA_todos/4_Contraste_con_sexo/sig.genes.df.txt",
                     header = T, 
                     sep = "\t")


# Se leen los genes expresados diferencialmente en la sangre
sangre_ASD_control <- read.table(file = "Results/MA/Contrastes_alfa_0.05/2_MA_sangre/3_Contraste_sin_sexo/sig.genes.df.txt",
                     header = T, 
                     sep = "\t")


# Se leen los genes expresados diferencialmente en el cerebro
cerebro_females <- read.table(file = "Results/MA/Contrastes_alfa_0.05/3_MA_cerebro/1_Contraste_mujeres/sig.genes.df.txt",
                     header = T, 
                     sep = "\t")


cerebro_ASD_control <- read.table(file = "Results/MA/Contrastes_alfa_0.05/3_MA_cerebro/3_Contraste_sin_sexo/sig.genes.df.txt",
                     header = T, 
                     sep = "\t")

cerebro_ASD_sex <- read.table(file = "Results/MA/Contrastes_alfa_0.05/3_MA_cerebro/4_Contraste_con_sexo/sig.genes.df.txt",
                     header = T, 
                     sep = "\t")

# combined_table_all <- rbind(all_females, all_ASD_control, all_ASD_sex)
# combined_table_all <- combined_table_all %>% distinct(ID, .keep_all = TRUE)

# combined_table_cerebro <- rbind(cerebro_females, cerebro_ASD_control, cerebro_ASD_sex)
# combined_table_cerebro <- combined_table_cerebro %>% distinct(ID, .keep_all = TRUE)


###############################################################################
# GENERACION DE LAS LISTAS DE GENES
###############################################################################

##### MA_todos #####
datos <- all_ASD_control
# datos <- subset(datos, logFC > 0)
datos <- datos$ID

# Crear una única cadena con cada ID seguido de un salto de línea
cadena_resultante <- paste(datos, collapse = "\n")

# Imprimir la cadena resultante
cat(cadena_resultante)

##### MA_cerebro #####
datos <- cerebro_ASD_control
# datos <- subset(datos, logFC > 0)
datos <- datos$ID

# Crear una única cadena con cada ID seguido de un salto de línea
cadena_resultante <- paste(datos, collapse = "\n")

# Imprimir la cadena resultante
cat(cadena_resultante)
