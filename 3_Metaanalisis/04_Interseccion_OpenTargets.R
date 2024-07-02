###############################################################################
# Autor: Gonzalo Anton Bernat
# Fecha: 02/05/2024
# Descripcion: 
#   Este programa permite realizar la intersección de los genes diferencialmente
#   expresados del MA contra los genes asociados al TEA de OpenTargets
###############################################################################

# Se cargan las librerías necesarias
library("openxlsx")
library("dplyr")

################################################################################
# LECTURA DE LOS DATOS                                                         #
################################################################################
ot <- read.table(file = "TEA_Open_Targets.tsv", sep = "\t", header = T)
ot <- data.frame(ot$symbol)
colnames(ot) <- "Symbol"

all_females <- read.table(file = "Results/MA/Contrastes_alfa_0.05/0_Contrastes todos/1_Females/sig.genes.df.txt",
                     header = T, 
                     sep = "\t")

all_ASD_control <- read.table(file = "Results/MA/Contrastes_alfa_0.05/0_Contrastes todos/3_ASD vs Control/sig.genes.df.txt",
                     header = T, 
                     sep = "\t")

all_ASD_sex <- read.table(file = "Results/MA/Contrastes_alfa_0.05/0_Contrastes todos/4_Female_ASD vs Male_ASD/sig.genes.df.txt",
                     header = T, 
                     sep = "\t")

sangre_ASD_control <- read.table(file = "Results/MA/Contrastes_alfa_0.05/3_Contraste_sangre_arrays+rnaseq/3_ASD vs Control/sig.genes.df.txt",
                     header = T, 
                     sep = "\t")

cerebro_males <- read.table(file = "Results/MA/Contrastes_alfa_0.05/4_Contraste_cerebro_arrays+rnaseq/1_Females/sig.genes.df.txt",
                     header = T, 
                     sep = "\t")

cerebro_ASD_control <- read.table(file = "Results/MA/Contrastes_alfa_0.05/4_Contraste_cerebro_arrays+rnaseq/3_ASD vs control/sig.genes.df.txt",
                     header = T, 
                     sep = "\t")

cerebro_ASD_sex <- read.table(file = "Results/MA/Contrastes_alfa_0.05/4_Contraste_cerebro_arrays+rnaseq/4_Female_ASD vs Male_ASD/sig.genes.df.txt",
                     header = T, 
                     sep = "\t")


################################################################################
# SE REALIZA LA INTERESECCIÓN                                                  #
################################################################################

datos2 <- merge(all_females, ot, by.x = "ID", by.y = "Symbol")
datos2 <- merge(all_ASD_control, ot, by.x = "ID", by.y = "Symbol")
datos2 <- merge(all_ASD_sex, ot, by.x = "ID", by.y = "Symbol")
datos2 <- merge(sangre_ASD_control, ot, by.x = "ID", by.y = "Symbol")
datos2 <- merge(cerebro_males, ot, by.x = "ID", by.y = "Symbol")
datos2 <- merge(cerebro_ASD_control, ot, by.x = "ID", by.y = "Symbol")
datos2 <- merge(cerebro_ASD_sex, ot, by.x = "ID", by.y = "Symbol")

datos2$Regulación <- ifelse(datos2$logFC > 0, "Up", "Down")
datos2 <- datos2 %>% rename(Symbol = ID)
datos3 <- data.frame(
  datos2$Symbol,
  datos2$p.adjust.fdr,
  datos2$logFC,
  datos2$Regulación,
  datos2$lower_bound,
  datos2$upper_bound,
  datos2$QE,
  datos2$QEp,
  datos2$SE,
  datos2$tau2,
  datos2$I2,
  datos2$H2,
  datos2$n.studies
  
)
colnames(datos3) <- sub("datos2\\.", "", colnames(datos3))

write.xlsx(datos3, file = "interseccion.xlsx")

ups <- subset(datos2, logFC >0)
dim(ups)
downs <- subset(datos2, logFC <0)
dim(downs)