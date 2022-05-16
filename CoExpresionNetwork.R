#Se cargan las librerías
library(GEOquery)
library(limma)
library(umap)
library(WGCNA)

#Opción para que permita descargarse los datos
Sys.setenv("VROOM_CONNECTION_SIZE" = 500000000)
#Se descargan los datos
gset <- getGEO("GSE106977", GSEMatrix =TRUE, AnnotGPL=FALSE)

if (length(gset) > 1) idx <- grep("GPL17586", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

#Se cambian los nombres de las caracteristicas para asegurarse de que sean compatibles con el formato
fvarLabels(gset) <- make.names(fvarLabels(gset))

#Grupos que se crearán en función de las muestras
gsms <- paste0("00010110000001000000100001010110001100010110010000",
               "01111101110000001100000001010010010100100001011010",
               "1100000110011110111")
sml <- strsplit(gsms, split="")[[1]]

#Vector que indica si tiene respuesta patológica completa o no para cada posición
gs <- factor(sml)
groups <- make.names(c("no pcr","pcr"))
levels(gs) <- groups

pcrTable <- data.frame(muestras = gset$geo_accession,ValorPCR = as.numeric(sml) , check.names = TRUE)
pcrTable[pcrTable$ValorPCR == 1,]$ValorPCR <- 2
pcrTable[pcrTable$ValorPCR == 0,]$ValorPCR <- 1
pcrTable[pcrTable$ValorPCR == 2,]$ValorPCR <- 0

write.table(pcrTable, file = "/Users/joseantoniomr/Desktop/TablaPCR.csv")

ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

MExpr <- t(ex)
colnames(MExpr) <- substr(colnames(MExpr),4, 10)
colnames(MExpr)[1] <- "0"

numGenes <- ncol(MExpr)
numSamples <- nrow(MExpr)
#Activo WGCNA
enableWGCNAThreads()
allowWGCNAThreads()

cor <- WGCNA::cor

#Creo la red de coexpresión génica 
Network = blockwiseModules(as.data.frame(MExpr), power = 7,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = FALSE, 
                       verbose = 3)



ModuleColors <- labels2colors(Network$colors)

 #Calculo los MEs, que son modulos eigengenes
SamplesMEs <- Network$MEs
SamplesMEs0 <- moduleEigengenes(MExpr,ModuleColors)$eigengenes
SamplesMEs0 <- orderMEs(SamplesMEs0)

#Correlación de cada módulo con la respuesta patológica completa y su p-valor
ModuleCorrelation <- cor(SamplesMEs0,pcrTable,use = "p")
ModulePValue <- corPvalueStudent(ModuleCorrelation,numSamples)
#Nombre de los módulos
ModuleNames <- substring(names(SamplesMEs0),3)
#Correlación de los módulos a los que puede pertenecer cada gen y el p-valor
geneModuleMembership <-  as.data.frame(cor(MExpr, SamplesMEs0, use = "p"))
MMPvalue <-  as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), numSamples))
names(geneModuleMembership) <-  paste("MM", ModuleNames, sep="")
names(MMPvalue) <-  paste("p.MM", ModuleNames, sep="")

#Correlación entre los genes expresados y la respuesta patológica completa y el p-valor
geneTraitSignificance = as.data.frame(cor(MExpr, pcrTable, use = "p"))[,-1]
GTSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), numSamples))

names(geneTraitSignificance) = paste("GS.", names(pcrTable[2]), sep="")
names(GTSPvalue) = paste("p.GS.", names(pcrTable), sep="")

#geneNames <- names(MExpr)
#Vector temporal que sustituye al nombre de los genes
geneNames <- 1:70523

#Módulos ordenados según la correlación que tengan con la variable de la respuesta patológica completa
modOrder = order(-abs(cor(SamplesMEs0, pcrTable, use = "p")[,-1]))

#Info recogida en un dataframe
geneInfo0 = data.frame(Genes = geneNames,
                       moduleColor = ModuleColors,
                       geneTraitSignificance,
                       GTSPvalue)

for (mod in 1:ncol(geneModuleMembership)){
  oldNames = names(geneInfo0)
  geneInfo <-  data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]], MMPvalue[, modOrder[mod]])
  names(geneInfo) <-  c(oldNames, paste("MM.", ModuleNames[modOrder[mod]], sep=""),
                         paste("p.MM.", ModuleNames[modOrder[mod]], sep=""))
}

write.table(geneInfo, "GeneInfo.csv")




