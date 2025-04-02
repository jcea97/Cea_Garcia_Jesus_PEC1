##################################################################
#Carga y preparación de los datos
##################################################################



# Cargo el csv 
cachexia <- read.csv("human_cachexia.csv", check.names = FALSE) 

# Creo la columna Cancer_Type
cachexia$Cancer_Type <- sub("_.*", "", cachexia$`Patient ID`) 

# Reordeno columnas
cachexia <- cachexia %>% relocate(Cancer_Type, .after = 2) 

# Renombramos pacientes
cachexia$`Patient ID`[cachexia$`Muscle loss` == "cachexic"] <- paste0("CA-", 1:sum(cachexia$`Muscle loss` == "cachexic")) 
cachexia$`Patient ID`[cachexia$`Muscle loss` == "control"] <- paste0("CT-", 1:sum(cachexia$`Muscle loss` == "control"))

# Preparo los metadatos
sample_metadata <- cachexia[, 1:3]
rownames(sample_metadata) <- cachexia$`Patient ID`

# Preparo la matriz de expresión
expr_matrix <- as.matrix(cachexia[, -(1:3)])
rownames(expr_matrix) <- sample_metadata$`Patient ID` 
expr_matrix <- t(expr_matrix)

# Creo el SummarizedExperiment
se <- SummarizedExperiment(
  assays = list(counts = expr_matrix),
  colData = sample_metadata
)


##################################################################
#Gráficos de distribución de los valores de concentración
##################################################################


# Accedo a la matriz de expresión
expr_unn <- assay(se)  

# Defino los grupos y sus colores
grupos <- colData(se)$`Muscle loss`  
colores <- ifelse(grupos == "cachexic", "red", "blue")

# Creo el boxplot de concentración por muestra sin normalizar
boxplot(expr_unn,
        las = 2, col = colores,
        cex.axis = 0.8,
        main = "Distribución de los valores de concentración")

#Normalización
expr_norm <- log2(assay(se))

# Creo el boxplot de concentración por muestra sin normalizar
boxplot(expr_norm,
        las = 2,
        col = colores,
        cex.axis = 0.8,
        main = "Distribución de los valores de concentración")


# Obtener expresión y grupos
grupo <- colData(se)$`Muscle loss`
# Calcular medianas por muestra
medianas <- apply(expr_norm, 2, median)

# Creo el boxplot de concentración por grupo sin normalizar
boxplot(medianas ~ grupo,
        col = c("red", "blue"),
        main = "Mediana de concentración por grupo",
        ylab = "Mediana de concentración")


##################################################################
# Análisis de componentes principales (PCA)
##################################################################


pc <-prcomp(t(expr_norm), scale=FALSE) #PCA sobre los datos normalizados
loads <- round(pc$sdev^2 / sum(pc$sdev^2) * 100, 1) #Varianza explicada por PC


# Representación de las muestras en el espacio definido por las dos primeras componentes principales (PCA), diferenciando por condición

pacientes_ID <- cachexia[ ,1] 

xlab<-c(paste("PC1",loads[1],"%"))
ylab<-c(paste("PC2",loads[2],"%"))

plot(pc$x[,1:2],xlab=xlab,ylab=ylab, col=colores, 
     main ="PCA")

legend("bottomleft",
       legend=c("cachexia", "control"), 
       col=c("red", "blue"), 
       pch=19, xpd=TRUE)

text(pc$x[,1],pc$x[,2],pacientes_ID, pos=3, cex=.6)

# Batch Effect

grupo_batch <- colData(se)$Cancer_Type

colores_batch <- ifelse(grupo_batch == "PIF", "red",
                        ifelse(grupo_batch == "NETCR", "blue",
                               "green"))

plot(pc$x[,1:2],xlab=xlab,ylab=ylab, col=colores_batch, 
     main ="PCA")

legend("bottomleft",  
       legend=c("PIF", "NETCR", "NETCL"), 
       col=c("red", "blue", "green"), 
       pch=19, xpd=TRUE)

text(pc$x[,1],pc$x[,2],pacientes_ID, pos=3, cex=.6)

##################################################################
# Agrupamiento jerárquico
##################################################################

d <- dist(t(expr_norm))  
hc <- hclust(d, method = "average")

# Dibujo el dendrograma
hc$labels <- gsub("_.*", "", hc$labels)  
plot(hc, main = "Dendrograma por muestra", cex = 0.7)
rect.hclust(hc, k = 2, border = c("red", "blue"))


##################################################################
#Análisis de expresión diferencial
##################################################################


ttest <- function(x) {
  tt = t.test(x[1:47], x[48:77]) # Aplico un test t entre los 47 caquéxicos y los 30 controles 
  return(c(
    t = tt$statistic,
    p.value = tt$p.value,
    log2FC = log2(mean(x[1:47]) / mean(x[48:77])) # Calculo el log2 Fold Change entre las medias de caquéxicos y controles
  ))
}

ans <- apply(expr_norm, 1, ttest) # Aplico la función ttest() a cada fila de la matriz expr_norm (cada metabolito)


# Convertimos la matriz a data.frame
ttest_df <- as.data.frame(t(ans))

# Añadir nombres de los metabolitos
ttest_df$metabolitos <- rownames(ttest_df)

EnhancedVolcano(ttest_df,
                lab = ttest_df$metabolitos,
                x = 'log2FC',
                y = 'p.value',
                pCutoff = 0.05,
                FCcutoff = 0.5, 
                pointSize = 2,
                labSize = 3,
                title = 'Volcano Plot - Metabolitos')

