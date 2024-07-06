###############################################################################
# Autor: Gonzalo Anton Bernat
# Fecha: 05/05/2024
# Descripcion: 
#   Este programa permite realizar el análisis funcional de los genes
#   diferencialmente expresados del metaanálisis
###############################################################################

# Se cargan las librerías
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(ggplot2)


# Se realiza la carga de los datos
all_females <- read.table(file = "2_Resultados_MA/1_MA_todos/1_Contraste_mujeres/sig.genes.df.txt",
                     header = T, 
                     sep = "\t")

all_ASD_control <- read.table(file = "2_Resultados_MA/1_MA_todos/3_Contraste_sin_sexo/sig.genes.df.txt",
                     header = T, 
                     sep = "\t")

all_ASD_sex <- read.table(file = "2_Resultados_MA/1_MA_todos/4_Contraste_con_sexo/sig.genes.df.txt",
                     header = T, 
                     sep = "\t")

sangre_ASD_control <- read.table(file = "2_Resultados_MA/2_MA_sangre/3_Contraste_sin_sexo/sig.genes.df.txt",
                     header = T, 
                     sep = "\t")

cerebro_females <- read.table(file = "2_Resultados_MA/3_MA_cerebro/1_Contraste_mujeres/sig.genes.df.txt",
                     header = T, 
                     sep = "\t")

cerebro_ASD_control <- read.table(file = "2_Resultados_MA/3_MA_cerebro/3_Contraste_sin_sexo/sig.genes.df.txt",
                     header = T, 
                     sep = "\t")

cerebro_ASD_sex <- read.table(file = "2_Resultados_MA/3_MA_cerebro/4_Contraste_con_sexo/sig.genes.df.txt",
                     header = T, 
                     sep = "\t")

#combined_table_all <- rbind(all_females, all_ASD_control, all_ASD_sex)
#combined_table_all <- combined_table_all %>% distinct(ID, .keep_all = TRUE)

#combined_table_cerebro <- rbind(cerebro_females, cerebro_males, cerebro_ASD_control, cerebro_ASD_sex)
#combined_table_cerebro <- combined_table_cerebro %>% distinct(ID, .keep_all = TRUE)


###############################################################################
# Preparacion de datos
###############################################################################

diffexp <- cerebro_ASD_control
head(diffexp)
dim(diffexp)

# differential expressed genes 
table(diffexp$adj.P.Val < 0.05)

sig.genes <- diffexp[diffexp$p.adjust.fdr < 0.05, ] 

###############################################################################
# Análisis de sobrerrepresentación (ORA)
###############################################################################

## Genes sobreexpresados
threshold <- 0
top <- sig.genes[sig.genes$logFC > threshold, ]
top_genes <- top$ID


oraTop <- enrichGO(gene        = top_genes,
                 OrgDb         = org.Hs.eg.db,
                 keyType       = 'SYMBOL',
                 ont           = "BP",
                 pAdjustMethod = "BH",
                 
                 qvalueCutoff  = 0.05)

dim(oraTop)
grafico <- barplot(oraTop)
grafico <- grafico +
  ggtitle("Análisis de sobrerrepresentación")
grafico


## Genes infraexpresados
threshold <- -0
bottom <- sig.genes[sig.genes$logFC < threshold, ]
bottom_genes = bottom$ID

oraBottom <- enrichGO(gene     = bottom_genes,
                 OrgDb         = org.Hs.eg.db,
                 keyType       = 'SYMBOL',
                 ont           = "BP",
                 pAdjustMethod = "BH"
                 )

dim(oraBottom)

grafico <- barplot(oraBottom, showCategory = 10)
grafico <- grafico +
  ggtitle("Análisis de sobrerrepresentación")
grafico

