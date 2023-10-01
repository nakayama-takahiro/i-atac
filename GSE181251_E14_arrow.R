library(Signac)
library(Seurat)
library(ggplot2)
library(patchwork)
library(ArchR)
addArchRGenome("mm10")

E14 <- ArchRProject(
  ArrowFiles = "E14.arrow", 
  outputDirectory = "../",
  copyArrows = TRUE
)

df <- getCellColData(E14, select = "nFrags")

idxPass <- which(E14$TSSEnrichment >= 8 & df$nFrags >= 1000)
cellsPass <- E14$cellNames[idxPass]
E14[cellsPass, ]

E14 <- addIterativeLSI(
  ArchRProj = E14,
  useMatrix = "TileMatrix", 
  name = "IterativeLSI", 
  iterations = 2, 
  clusterParams = list(
    resolution = c(0.2), 
    sampleCells = 10000, 
    n.start = 10
  ), 
  varFeatures = 25000, 
  dimsToUse = 1:30
)

E14 <- addClusters(
  input = E14,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "Clusters",
  resolution = 0.8
)

E14 <- addUMAP(
  ArchRProj = E14, 
  reducedDims = "IterativeLSI", 
  name = "UMAP", 
  nNeighbors = 30, 
  minDist = 0.5, 
  metric = "cosine"
)

# saveArchRProject(ArchRProj = E14, outputDirectory = "../", load = FALSE)

E14_2 <- readRDS("GSE181251_E14_QC.rds")

p1 <- plotEmbedding(ArchRProj = E14_2, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
p1

E14_2 <- addImputeWeights(E14_2)

markerGenes  <- c(
  "Ccnd1",   #MPC
  "Plxna2",  #RPC
  "Mybl1",   #RPC
  "Notch1", "Notch2", "Notch3", "Notch4",
  "Jag1", "Jag2", "Dll1", "Dll3", "Dll4",
  "Hes1", "Hes5", "Hey1", "Hey2"
)

p <- plotEmbedding(
  ArchRProj = E14_2, 
  colorBy = "GeneScoreMatrix", 
  name = markerGenes, 
  embedding = "UMAP",
  imputeWeights = getImputeWeights(E14_2)
)

p$Ccnd1
p$Plxna2
p$Mybl1

p <- plotBrowserTrack(
  ArchRProj = E14_2, 
  groupBy = "Clusters", 
  geneSymbol = markerGenes, 
  upstream = 50000,
  downstream = 50000
)

grid::grid.newpage()
grid::grid.draw(p$Notch1)

p4 <- plotGroups(
  ArchRProj = E14_2, 
  groupBy = "Clusters", 
  colorBy = "GeneScoreMatrix", 
  name = markerGenes,
  plotAs = "violin",
  alpha = 0.4,
  addBoxPlot = TRUE
)
p4$Ccnd1
