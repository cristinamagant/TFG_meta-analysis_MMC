# Se cargan los diferentes paquetes que vamos a utilizar
library(GEOquery)
library(AnnotationDbi)
library(limma)
library(hgu133a.db)
library(affy)
library(geneplotter)
source("medianreps.R")

#### PARTE 1: ANÁLISIS EXPLORATORIO ######

#### Descarga de los datos de GEO ####

parent_dir <- getwd()
dir.create("data")
dir.create("images")
# Establecemos el directorio de trabajo
setwd("data")

# Descargamos los datos de GEO
GEOquery::getGEOSuppFiles("GSE46517")

# Establecemos el nuevo directorio de trabajo donde están descargados los datos de GEO
setwd("GSE46517")

# Descomprimimos los datos
system("tar xvf GSE46517_RAW.tar")

# setwd("GSE46517")

#### Lectura y comprobación de los datos ####

# Leemos los datos del chip de Affimetrix
GSE46517_raw <- ReadAffy()  
# GSE46517_raw=x

# Comprobamos el tipo de datos que tenemos
class(GSE46517_raw)

#### Anotación de los datos ####

# Nos da información sobre con qué fichero se debe realizar la anotación
annotation(GSE46517_raw)

# Nos guardamos lo que tenemos hecho hasta el momento y lo cargamos para trabajar
# posteriormente con una mayor facilidad
save(GSE46517_raw, file = "../../GSE46517_raw.RData")

load(paste0(parent_dir,"/GSE46517_raw.RData"))

# Cambiamos el nombre para trabajar de una forma más sencilla
gse46517_raw <- GSE46517_raw

# Dimensiones de la matriz de expresión del estudio
dim(exprs(gse46517_raw))

## Análisis exploratorio

### Realización de MA-plots

# Nos guardamos las imágenes de los MA-plots
png(filename = 'images/gse46517_raw_MAplot.png')
affy::MAplot(gse46517_raw, pairs=F, plot.method = "smoothScatter")
dev.off()

# Vemos en R el resultado de los MA-plots
png(filename = 'images/gse46517_raw_MAplot.png')
MAplot = affy::MAplot(gse46517_raw, pairs = FALSE, which = 1, plot.method = "smoothScatter")

## Realización de los boxplots

# Nos guardamos las imágenes de los boxplots
png(filename = "images/gse46517_raw_boxplot.png")

affy::boxplot(gse46517_raw,
              main = "Boxplot sin normalizar",
              xlab = "Samples", names = FALSE)
dev.off()

# Lo visualizamos en R
boxplot <- affy::boxplot(gse46517_raw,
                         main = "Boxlot sin normalizar",
                         xlab = "Samples", names = FALSE)

# Comprobamos si los datos de la matriz de expresión están normalizados
head(exprs(gse46517_raw))

#### Normalización ####

# Realizamos una normalización de los datos mediante logaritmo en base 2
expresion_original <- exprs(gse46517_raw)
expresion_log2 <- log2(expresion_original)

# Sustituimos los valores de expresión sin normalizar por los valores normalizados
gse46517 <- gse46517_raw
exprs(gse46517) <- expresion_log2

#### Análisis exploratorio post-normalización ####
## Vemos los boxplots
png(filename = "images/gse46517_boxplot.png")
graphics::boxplot(exprs(gse46517))
dev.off()
graphics::boxplot(exprs(gse46517))

###### PARTE 2: ADAPTACIÓN PARA METAFUN ######

#### Descarga de los datos de GEO ####

# Establecemos el directorio de trabajo
setwd("/home/biouser/Desktop/TFG_cristina/version_final/")
# dir.create("/home/biouser/Desktop/TFG_cristina/version_final/study_GSE46517") # creamos
setwd("/home/biouser/Desktop/TFG_cristina/version_final/study_GSE46517")

# Se descargan los datos raw
gset <- GEOquery::getGEO("GSE46517", GSEMatrix =TRUE, getGPL=FALSE)
GSE46517_raw <- gset$GSE46517_series_matrix.txt.gz

GSE46517 <- GSE46517_raw # para o sobreescribir

exprs(GSE46517) <- log2(exprs(GSE46517_raw)) # normalizamos por log2 y guardamos en el objeto no raw

minimo <- min(exprs(GSE46517)) # el mínimo se encuentra entre 0-1, con log2 saldrán negativos

exprs(GSE46517) <- exprs(GSE46517) + abs(minimo)

boxplot(exprs(GSE46517_raw))
boxplot(exprs(GSE46517))

#### Anotación ####

# obtenemos del paquete de anotacion la relacion PROBEID - ENTREZID
probeid2eid <- select(hgu133a.db, keys = featureNames(GSE46517),
                      columns = c("ENTREZID"), keytype = "PROBEID")

# guardamos los índices donde coinciden los genes anotados
indices <- BiocGenerics::match(featureNames(GSE46517), probeid2eid$PROBEID)

fData(GSE46517) <- probeid2eid[indices,]

# pasamos a matriz
matriz <- exprs(GSE46517)

# las filas son el ENTREZID
rownames(matriz) <- fData(GSE46517)$ENTREZID
matriz <- matriz[!is.na(rownames(matriz)),]

# calculamos la media de los genes repetidos
matriz <- medianReps(matriz)

GSE46517_expr_final <- as.data.frame(matriz)

# guardamos
write.table(
  file = "./data2/GSE46517_expr_final.csv",
  GSE46517_expr_final,
  sep = ",",
  row.names = T
)

###### Obtención de los metadatos #####
GSE46517_metadata <- data.frame(GSE46517$source_name_ch1,GSE46517$geo_accession)

### Comparativa 1 ###
GSE46517_metadata_comp1 <- GSE46517_metadata

# Se cambia el nombre de los grupos de interés
GSE46517_metadata_comp1$groups[grepl("Primary Melanoma", GSE46517_metadata_comp1$GSE46517.source_name_ch1)] <- "control"
GSE46517_metadata_comp1$groups[grepl("Metastatic Melanoma", GSE46517_metadata_comp1$GSE46517.source_name_ch1)] <- "case"

# Se omiten los valores NA
GSE46517_metadata_comp1 <- na.omit(GSE46517_metadata_comp1)

# Se seleccionan todas las filas y solo las columnas 2 y 3
GSE46517_metadata_comp1 <- GSE46517_metadata_comp1[,2:3]

# Se establecen los nombres de las columnas
colnames(GSE46517_metadata_comp1) <- c("sample","group")

# Se guarda en un fichero TSV para utilizarlo en el METAFUN
write.table(file="./data2/GSE46517_metadata_comp1.tsv",GSE46517_metadata_comp1, sep = "\t",
            row.names = F, col.names = F)

### Comparativa 2 ###
GSE46517_metadata_comp2 <- GSE46517_metadata

# Se cambia el nombre de los grupos de interés
GSE46517_metadata_comp2$groups[grepl("Normal", GSE46517_metadata_comp2$GSE46517.source_name_ch1)] <- "control"
GSE46517_metadata_comp2$groups[grepl("Nevus", GSE46517_metadata_comp2$GSE46517.source_name_ch1)] <- "control"
GSE46517_metadata_comp2$groups[grepl("Metastatic Melanoma", GSE46517_metadata_comp2$GSE46517.source_name_ch1)] <- "case"

# Se omiten los valores NA
GSE46517_metadata_comp2 <- na.omit(GSE46517_metadata_comp2)
GSE46517_metadata_comp2 <- GSE46517_metadata_comp2[!GSE46517_metadata_comp2$GSE46517.geo_accession == "GSM1131685",]

# Se seleccionan todas las filas y solo las columnas 2 y 3
GSE46517_metadata_comp2 <- GSE46517_metadata_comp2[,2:3]

# Se establecen los nombres de las columnas
colnames(GSE46517_metadata_comp2) <- c("sample","group")

# Se guarda en un fichero TSV para utilizarlo en el METAFUN
write.table(file="./data2/GSE46517_metadata_comp2.tsv",GSE46517_metadata_comp2, sep = "\t",
            row.names = F, col.names = F)

###### adaptar los expr al metadata ###########

### Comparativa 1 ####
GSE46517_expr_comp1 <- GSE46517_expr_final[,GSE46517_metadata_comp1$sample]

# Se guarda en un fichero CSV para utilizarlo en el METAFUN
write.table(
  file = "./data2/GSE46517_expr_comp1.csv",
  GSE46517_expr_comp1,
  sep = ",",
  row.names = T
)

### Comparativa 2 ###
GSE46517_expr_comp2 <- GSE46517_expr_final[,GSE46517_metadata_comp2$sample]

# Se guarda en un fichero CSV para utilizarlo en el METAFUN
write.table(
  file = "./data2/GSE46517_expr_comp2.csv",
  GSE46517_expr_comp2,
  sep = ",",
  row.names = T
)
