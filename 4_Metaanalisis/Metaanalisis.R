###############################################################################
# Autor: Francisco García-García
# Modificado por: Gonzalo Antón Bernat
# Fecha: 14/04/2024
# Descripcion: 
#   Este programa permite realizar el metaanálisis. Requiere de forma previa la
#   creación de un objeto ALLfit y que ese objeto se encuentre en la carpeta de
#   este programa
###############################################################################

# Se cargan las librerías necesarias 
library(Biobase)
library(metafor)
library(stringr)
library(limma)

# Se cargan los datos
load("../3_Genera_todos_fit/todos_fit2.RData")
#load("../3_Genera_todos_fit/sangre_fit2.RData")
#load("../3_Genera_todos_fit/cerebro_fit2.RData")
ALLfit <- todos_fit2

# Se asignan los directorios
wd <- "./Results/MA/" #Output dir
metaanalysis_name <- "Metaanalisis"
wd <- paste0(wd,metaanalysis_name,"/")
dir.create(wd, recursive = TRUE)


###############################################################################
# PASO 0. PREPROCESADO DE LOS DATOS
###############################################################################

## Se calcula el SE
SE_array <- function(fit) {
  #The effect sizes are contained in fit$coefficients
  summary(fit$coefficients)
  head(fit$coefficients)
  SE.coef <- sqrt(fit$s2.post) * fit$stdev.unscaled[,3] # Se calcula el error estandar
  head(SE.coef)
  summary(SE.coef)
  
  # Se elige el contraste de interes
  allgenes <- topTable(fit, number = "all", confint=TRUE, adjust.method = "fdr",coef = 4)
  dim(allgenes)
  allgenes[, "SE"] <- (allgenes[, "CI.R"] - allgenes[, "CI.L"])/ 3.92
  head(allgenes)
  
  # Se procesan los resultados obtenidos
  table(rownames(SE.coef) == rownames(fit$coefficients))
  mat <- cbind(fit$coefficients[,3], SE.coef)
  colnames(mat) <- c("coef", "se.coef")
  head(mat)
  
  int <- intersect(rownames(allgenes), rownames(mat))
  length(int)
  res <- cbind(allgenes, mat[rownames(allgenes),])
  head(res)
  dim(res)
  return(res)
}

# Se utiliza la lista como input para el metaanálisis
EDs_sel <- lapply(ALLfit, SE_array)


###############################################################################
# PASO 1. PREPARACIÓN DEL INPUT DEL METAANÁLISIS: LOR y matriz SE
###############################################################################

# Se identifican los genes únicos de todos los estudios
genes <- NULL
for (fi in EDs_sel){
  genes <- c(genes, rownames(fi))
}

length (genes)
genes <- unique (genes)
length (genes)
genes[grepl("[0-9]+-[A-Z][a-z]{2}",genes)]
genes <- genes[!grepl("[0-9]+-[A-Z][a-z]{2}",genes)]
genes[grepl("[0-9]+-[A-Z][a-z]{2}",genes)]
genes <- base::sort(genes)


# Se genera una matriz con los logFC de todos los estudios
mat.logFC <- matrix (NA, nrow = length (genes), ncol = length(EDs_sel))
rownames (mat.logFC) <- genes
colnames (mat.logFC) <- names(EDs_sel)
head (mat.logFC)

for (i in 1:length(EDs_sel)){
  co <- names(EDs_sel[i])
  res <- EDs_sel[[i]]
  logFC <- res$logFC
  names (logFC) <- (rownames (res))
  mat.logFC[, co] <- logFC[rownames(mat.logFC)] 
}

head (mat.logFC)
tail(mat.logFC)
table (is.na(mat.logFC))
dim (mat.logFC)

# Los genes seleccionados deben estar en 2 o más estudios
mat.logFC.NA <- is.na(mat.logFC)
head(mat.logFC.NA)
sum.NA <-  apply(mat.logFC.NA, 1, sum)
table(sum.NA)
min.sum.NA <- sum.NA < ncol(mat.logFC) - 1
table(min.sum.NA)

# Se filtra por min.sum.NA
mat.logFC <- mat.logFC[min.sum.NA == T, ]
dim(mat.logFC)


# Se genera una matriz con todos los SE para todos los estudios
mat.SE <- matrix (NA, nrow = length (genes), ncol = length(EDs_sel))
rownames (mat.SE) <- genes
colnames (mat.SE) <- names(EDs_sel)
head (mat.SE)

## SE FROM GORDON: se.coef
for (i in 1:length(EDs_sel)){
  co <- gsub("_ED", "", names(EDs_sel[i]))
  res <- EDs_sel[[i]]
  SE <- res$se.coef
  names (SE) <- (rownames (res))
  mat.SE[, co] <- SE[rownames(mat.SE)] 
}

head (mat.SE)
tail(mat.SE)
table (is.na(mat.SE))
dim (mat.SE)

# Se filtra por min.sum.NA
mat.SE <- mat.SE[min.sum.NA == T, ]
dim(mat.SE)




###############################################################################
# PASO 2. METAANÁLISIS DE LOS GENES
###############################################################################

MA <- lapply(1:length(rownames(mat.logFC)),
             function(x){rma(yi = mat.logFC[x, ],
                             sei = mat.SE[x, ],
                             method = "DL")})

names (MA) <- rownames(mat.logFC)
class (MA)
length(MA)
head (MA)
MA[[1]]

# Se construye una tabla con todos los resultados
result_meta <- as.data.frame(do.call("rbind",
                                     lapply(MA,
                                            function(x){
                                              c(x$ci.lb, x$b, x$ci.ub, 
                                                x$pval, x$QE, x$QEp, x$se,
                                                x$tau2, x$I2, x$H2)
                                            })))

colnames(result_meta) <- c("lower_bound", "logFC", "upper_bound",
                           "pvalue", "QE", "QEp", "SE", "tau2", "I2", "H2")

p.adjust.fdr <- stats::p.adjust(result_meta[,4], method = "fdr")
p.adjust.BY  <- stats::p.adjust(result_meta[,4], method = "BY")
result_meta <- round(cbind(result_meta, p.adjust.fdr, p.adjust.BY), 3)
head(result_meta)


# Genes significativos
corte = 0.05
table(result_meta[, "pvalue"] < corte)
table(result_meta[, "p.adjust.fdr"] < corte)
table(result_meta[, "p.adjust.BY"] < corte)


# Se añade el número de estudios donde se evalua el gen
sum.NA <- sum.NA[sum.NA<ncol(mat.logFC)-1]
n.studies <-  ncol(mat.logFC) - sum.NA 
table(n.studies)
n.studies <- n.studies [rownames(mat.logFC)]
length(n.studies)
result_meta[, "n.studies"]  <- n.studies
head(result_meta)
summary(result_meta$p.adjust.fdr)


sig.genes.df = result_meta[result_meta$p.adjust.fdr < corte,] 
dim(sig.genes.df)

write.table(x = sig.genes.df[order(sig.genes.df$p.adjust.fdr),], file = paste0(wd,"sig.genes.tsv"), sep = "\t", quote = FALSE)
write.table(x = result_meta[order(result_meta$p.adjust.fdr),], file = paste0(wd,"all.genes.tsv"), sep = "\t", quote = FALSE)


###############################################################################
# PASO 3. ANÁLISIS DE INFLUENCIA Y SENSIBILIDAD
###############################################################################

# Se añaden nuevas variables  
for (i in rownames(sig.genes.df)){
  print(i)
  estudios <- colnames(mat.logFC)[!mat.logFC.NA[i,]]
  
  # number of studies where the sign of the logOR is the same  of the global logOR:
  sig.genes.df[i, "infl.same.sign.logFC"] <- sum(sign(MA[[i]]$yi)== rep(sign(MA[[i]]$b),length(estudios)))
  
  # how many studies could be influencers?
  inf <- influence(MA[[i]])
  res <- paste(estudios[inf$is.infl], collapse = ",")  
  sig.genes.df[i, "infl.nstudies"] <- ifelse(res =="", "non", res)
  
  # sensivity analysis
  l1 <-as.data.frame(leave1out(MA[[i]]))
  rownames(l1) <- estudios
  
  # p.value about differences between all estimates from leave one out
  #   and global estimate)
  sig.genes.df[i, "sensi.global"] <-t.test(x= l1$estimate,
                                           mu=as.numeric(MA[[i]]$b))$p.value
  # number of  studies where pvalue > 0.05 
  res2 <- paste(estudios[l1$pval > 0.05], collapse = ",")  
  sig.genes.df[i, "sensi.specific"] <- ifelse(res2 =="", "all.p.values < 0.05", res2)
}


#1. INFLUENCE STUDIES. How many logOR have the same sign to global logOR?
table(sig.genes.df$infl.same.sign.logFC)

#2. INFLUENCE STUDIES. How many functions including influence studies?
table(sig.genes.df$infl.nstudies=="non")

#3. SENSITIVITY. In global, are there many functions with differences in the estimate?
table(sig.genes.df$sensi.global < 0.05)

#4. SENSITIVITY.  How many functions including changes in the significance about 
# its new estimate  after leave1out? 
table(sig.genes.df$sensi.specific == "all.p.values < 0.05")


# Se guardan los resultados
cat ("ID\t", file = paste0(wd,"sig.genes.df.txt"))
write.table(sig.genes.df, file = paste0(wd,"sig.genes.df.txt"), sep ="\t", quote = F, 
            append = TRUE, row.names = T)



###############################################################################
# PASO 4. VISUALIZACIÓN DE LOS GENES SIGNIFICATIVOS
###############################################################################

# Se selecciona lo que se quiere visualizar
sig.results <- result_meta[result_meta[, "p.adjust.fdr"] < 0.05,]
sig.results
dim(sig.results)

dir.create(paste0(wd,"plots"), recursive = TRUE)

selMethod <- "DL"

for (i in 1:nrow(sig.results)){
  mygenes <- rownames(sig.results)[i]
  res <- rma(yi= mat.logFC[mygenes,], sei =mat.SE[mygenes,], method = "DL")
  
  # FOREST PLOT
  png (filename = paste0(wd,"plots/",gsub("-","_",mygenes),"_FOREST",".png"), width = 960 , 
       height = 960, res = 200) 
  forest(res, 
         slab = toupper(colnames(mat.logFC)), #Nombre de los estudios
         xlab="logFC", cex=0.7,
         mlab=paste(selMethod, "Model for All Studies", sep = " "), 
         border = "black", #Color del borde del rombo
         col = "red", #Color del rombo
         main = paste("\n", mygenes, sep=""))    
  text( 9,-3, "logFC [IC 95%]", pos=2, cex = 0.7)
  dev.off()
  
  # FUNNEL PLOT
  png (filename = paste0(wd,"plots/",gsub("-","_",mygenes),"_FUNNEL", ".png"), width = 960 , 
       height = 960, res = 200) 
  par(mfrow=c(2,2))
  funnel(res, main="Standard Error", back ="darkslategray1",
         xlab = paste("logFC (", mygenes, ")",sep =""))
  funnel(res, yaxis="vi", main="Sampling Variance", back ="darkslategray1",
         xlab = paste("logFC (", mygenes, ")",sep =""))
  funnel(res, yaxis="seinv", main="Inverse Standard Error",
         back ="darkslategray1", xlab = paste("logFC (", mygenes, ")",sep =""))
  funnel(res, yaxis="vinv", main="Inverse Sampling Variance", 
         back ="darkslategray1",  xlab = paste("logFC (", mygenes, ")",sep =""))
  par(mfrow=c(1,1))
  dev.off()
  
  # INFLUENCE PLOTS 
  png (filename = paste0(wd,"plots/",gsub("-","_",mygenes), "_INFLUENCE", ".png"), width = 960 , 
       height = 960, res = 200) ##CAMBIAR
  inf <- influence(res)
  plot(inf)
  dev.off()
  
}


###############################################################################
# PASO 5. GENERACIÓN DE UN REPORTE
###############################################################################

sig.genes.df <- sig.genes.df[order(sig.genes.df$p.adjust.fdr),]
save(sig.genes.df, result_meta, file = paste0(wd,metaanalysis_name, ".RData"))

# Función para crear múltiples tabs
make.tabs <- function(sig.genes.df){
  res <- NULL
  for(g in rownames(sig.genes.df)){
    file_name <- gsub("-","_", g)
    res <- c(res, '### ', g, '{-} \n',
             "**Statistics of ", g, " meta-analisys** \n",
             "```{r, fig.align='center'}", '\n',
             "kable(sig.genes.df['",g,"',1:11])", '\n',
             '```', '\n',
             "[Gene information](https://www.genecards.org/cgi-bin/carddisp.pl?gene=", g, ") \n\n",
             "**Forest plot** \n",
             "```{r, fig.align='center'}", '\n',
             'knitr::include_graphics("', 'plots/', file_name, '_FOREST.png")\n',
             '```', '\n',
             "**Funnel plot** \n",
             "```{r, fig.align='center'}", '\n',
             'knitr::include_graphics("', 'plots/', file_name, '_FUNNEL.png")\n',
             '```', '\n',
             "**Incluence plot** \n",
             "```{r, fig.align='center'}", '\n',
             'knitr::include_graphics("', 'plots/', file_name, '_INFLUENCE.png")\n',
             '```', '\n\n')
  }
  return(res)
}


# Se pasa de Rmd a knit
ns <- nrow(sig.genes.df)
cat(
  '---
title: "Meta-analysis of genes [DRAFT]"
output:
  html_document:
    toc: false
    toc_float: false
    code_folding: hide
    number_sections: true
    theme: spacelab
---
## ', metaanalysis_name, ' {.tabset .tabset-pills -}
  
```{r, warning=F, message=F}
library(dplyr)
library(knitr)
load("', metaanalysis_name, '.RData")
```  \n
### Significant results {-}  \n',
  "```{r, fig.align='center'}", '\n',
  "kable(sig.genes.df[,1:11], caption='Statistics of ", ns, " significant genes')", '\n',
  '```', '\n\n',
  make.tabs(sig.genes.df), "\n\n",
  '### sessionInfo {-}  \n',
  "```{r, fig.align='center'}", '\n',
  "date()", "\n",
  "sessionInfo()", '\n',
  '```', '\n\n',
  sep = "",
  file = paste0(wd,metaanalysis_name, ".Rmd"))

# Se genera el reporte HTML
rmarkdown::render(paste0(wd,metaanalysis_name, ".Rmd"))
