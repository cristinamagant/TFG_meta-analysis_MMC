# Se cargan los diferentes paquetes que vamos a utilizar
library(GEOquery)
library(AnnotationDbi)
library(limma)
library(illuminaHumanv4.db)
library(lumi)
library(illuminaio)
library(geneplotter)
source("medianreps.R")

###### PARTE 1: ANÁLISIS EXPLORATORIO ######

parent_dir <- getwd()
dir.create("data")
dir.create("images")

setwd("data")
GEOquery::getGEOSuppFiles("GSE60464")
setwd("GSE60464")
system("gzip -d *.gz")

GSE60464 <- lumiR("GSE60464_non_normalized.txt")
save(GSE60464, file = "../../GSE60464_raw.RData")

# Con la función png y dev.off guardamos la figura, si queremos ejecutar sin
# guardar para no machacar siempre la imagen resultante, ejecutar solo la linea
# de código que crea la imagen, en este caso: lumi::MAplot(GSE60464, smoothScatter=TRUE)

# Margins too large
png("ma_plot2.png")
lumi::MAplot(GSE60464, smoothScatter=TRUE,logMode = TRUE)
dev.off()
lumi::MAplot(GSE60464, smoothScatter=TRUE)

# Realizamos un diagrama de cajas
# Tener en cuenta que por defecto el gráfico aplica un log2
png("images/boxplot_non_normalized")
lumi::boxplot(GSE60464, logMode = FALSE)
dev.off()

# Realizamos una normalización de los datos mediante logaritmo en base 2
expresion_original <- exprs(GSE60464)
expresion_log2 <- log2(expresion_original)
exprs(GSE60464) <- expresion_log2


# Realizamos un gráfico con los estimadores de densidades
png("images/hist_normalized")
lumi::hist(GSE60464, logMode = FALSE)
dev.off()

# Realizamos un diagrama de cajas
png("images/boxplot_normalized")
lumi::boxplot(GSE60464, logMode = FALSE)
dev.off()

###### PARTE 2: ADAPTACIÓN PARA METAFUN ######

#### Descarga de los datos de GEO ####

# Establecemos el directorio de trabajo
setwd("/home/biouser/Desktop/TFG_cristina/version_final/")
# dir.create("/home/biouser/Desktop/TFG_cristina/version_final/study_GSE60464") # creamos
setwd("/home/biouser/Desktop/TFG_cristina/version_final/study_GSE60464")

# Descarga
gset <- GEOquery::getGEO("GSE60464", GSEMatrix =TRUE, getGPL=FALSE)
GSE60464_raw <- gset$GSE60464_series_matrix.txt.gz

GSE60464 <- GSE60464_raw # para no sobreescribir

# Normalizamos por log2 y guardamos en el objeto no raw
exprs(GSE60464) <- log2(exprs(GSE60464_raw))

png("boxplot")
boxplot(exprs(GSE60464_raw))
dev.off()
boxplot(exprs(GSE60464))

### Anotación ###

#### Anotación ####

# obtenemos del paquete de anotación la relacion PROBEID - ENTREZID
probeid2eid <- select(illuminaHumanv4.db, keys = featureNames(GSE60464),
                      columns = c("ENTREZID"), keytype = "PROBEID")

# guardamos los índices donde coinciden los genes anotados
indices <- BiocGenerics::match(featureNames(GSE60464), probeid2eid$PROBEID)

fData(GSE60464) <- probeid2eid[indices,]

# pasamos a matriz
matriz <- exprs(GSE60464)

# las filas son el ENTREZID
rownames(matriz) <- fData(GSE60464)$ENTREZID
matriz <- matriz[!is.na(rownames(matriz)),]

# calculamos la mediana de los genes repetidos
matriz <- medianReps(matriz)

GSE60464_expr_final <- as.data.frame(matriz)

# Se guarda como CSV
write.table(
  file = "./GSE60464_expr_final.csv",
  GSE60464_expr_final,
  sep = ",",
  row.names = T
)

###### Obtención de los metadatos #####
GSE60464_metadata <- data.frame(GSE60464$source_name_ch1,GSE60464$geo_accession)

# Se cambia el nombre de los grupos de interés
GSE60464_metadata$groups[grepl("non-cerebrotropic", GSE60464_metadata$GSE60464.source_name_ch1)] <- "control"
GSE60464_metadata$groups[grepl("_cerebrotropic", GSE60464_metadata$GSE60464.source_name_ch1)] <- "case"

# Se omiten los valores NA
GSE60464_metadata <- na.omit(GSE60464_metadata)

# Se seleccionan todas las filas y solo las columnas 2 y 3
GSE60464_metadata <- GSE60464_metadata[,2:3]

# Se establecen los nombres de las columnas
colnames(GSE60464_metadata) <- c("sample","group")

# Se guarda en un fichero TSV para utilizarlo en el METAFUN
write.table(file="./GSE60464_metadata.tsv",GSE60464_metadata, sep = "\t",
            row.names = F, col.names = F)

