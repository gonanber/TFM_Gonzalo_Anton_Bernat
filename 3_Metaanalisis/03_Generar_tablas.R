###############################################################################
# Autor: Gonzalo Anton Bernat
# Fecha: 20/04/2024
# Descripcion: 
#   Este programa permite generar las tablas que resumen el análisis de 
#   expresión diferencial
###############################################################################

# Se cargan las librerías
library("openxlsx")


# Se cargan las funciones
generar_tabla <- function(tabla) {
  # Crear un nuevo data frame con las columnas seleccionadas
  datos <- data.frame(Symbol = tabla$ID,
                      p.adjust.fdr = tabla$p.adjust.fdr,
                      logFC = tabla$logFC,
                      lower_bound = tabla$lower_bound,
                      upper_bound = tabla$upper_bound,
                      QE = tabla$QE,
                      QEp = tabla$QEp,
                      SE = tabla$SE,
                      tau2 = tabla$tau2,
                      I2 = tabla$I2,
                      H2 = tabla$H2,
                      n.studies = tabla$n.studies)
  
  
  # Añadir la columna 'Regulation' basada en los valores de 'logFC'
  datos$Regulation <- ifelse(tabla$logFC > 0, "Up", "Down")
  
  datos2 <- data.frame(Symbol = datos$Symbol,
                       p.adjust.fdr = datos$p.adjust.fdr,
                       logFC = datos$logFC,
                       Regulation = datos$Regulation,
                       lower_bound = datos$lower_bound,
                       upper_bound = datos$upper_bound,
                       QE = datos$QE,
                       QEp = datos$QEp,
                       SE = datos$SE,
                       tau2 = datos$tau2,
                       I2 = datos$I2,
                       H2 = datos$H2,
                       n.studies = datos$n.studies)
  return(datos2)
}


################################################################################
# LECTURA DE LOS DATOS                                                         #
################################################################################

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


###############################################################################
# GENERACIÓN DE LAS TABLAS
###############################################################################

all_females <- generar_tabla(all_females)
all_ASD_control <- generar_tabla(all_ASD_control)
all_ASD_sex <- generar_tabla(all_ASD_sex)
sangre_ASD_control <- generar_tabla(sangre_ASD_control)
cerebro_females <- generar_tabla(cerebro_females)
cerebro_ASD_control <- generar_tabla(cerebro_ASD_control)
cerebro_ASD_sex <- generar_tabla(cerebro_ASD_sex)


###############################################################################
# GENERACIÓN DEL EXCEL
###############################################################################

# Crear un nuevo libro de Excel
wb <- createWorkbook()

# Función para añadir una pestaña por cada data frame
add_df_to_excel <- function(df, sheet_name) {
  addWorksheet(wb, sheet_name)
  writeDataTable(wb, sheet_name, df, tableStyle = "TableStyleMedium9")
}

# Aplicar la función a cada data frame
add_df_to_excel(all_females, "All Females")
add_df_to_excel(all_ASD_control, "All ASD vs Control")
add_df_to_excel(all_ASD_sex, "All ASD Sex")
add_df_to_excel(sangre_ASD_control, "Sangre ASD vs Control")
add_df_to_excel(cerebro_females, "Cerebro Females")
add_df_to_excel(cerebro_ASD_control, "Cerebro ASD vs Control")
add_df_to_excel(cerebro_ASD_sex, "Cerebro ASD Sex")

# Guardar el archivo Excel
saveWorkbook(wb, "Resultados_MA.xlsx", overwrite = TRUE)



