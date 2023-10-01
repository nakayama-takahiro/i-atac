library(Signac)
library(Seurat)
library(ggplot2)
library(patchwork)
library(ArchR)
addArchRGenome("mm10")

E16 <- ArchRProject(
  ArrowFiles = "E16.arrow", 
  outputDirectory = "../",
  copyArrows = TRUE
)

df <- getCellColData(E16, select = "nFrags")

idxPass <- which(E16$TSSEnrichment >= 8 & df$nFrags >= 1000)
cellsPass <- E16$cellNames[idxPass]
E16[cellsPass, ]

E16 <- addIterativeLSI(
  ArchRProj = E16,
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

E16 <- addClusters(
  input = E16,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "Clusters",
  resolution = 0.8
)

E16 <- addUMAP(
  ArchRProj = E16, 
  reducedDims = "IterativeLSI", 
  name = "UMAP", 
  nNeighbors = 30, 
  minDist = 0.5, 
  metric = "cosine"
)

# saveArchRProject(ArchRProj = E16, outputDirectory = "../", load = FALSE)

E16_2 <- readRDS("GSE181251_E16_QC.rds")

p1 <- plotEmbedding(ArchRProj = E16_2, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
p1

E16_2 <- addImputeWeights(E16_2)

markerGenes  <- c(
  "Ccnd1",   #MPC
  "Plxna2",  #RPC
  "Mybl1",   #RPC
  "Notch1", "Notch2", "Notch3", "Notch4",
  "Jag1", "Jag2", "Dll1", "Dll3", "Dll4",
  "Hes1", "Hes5", "Hey1", "Hey2"
)

p <- plotEmbedding(
  ArchRProj = E16_2, 
  colorBy = "GeneScoreMatrix", 
  name = markerGenes, 
  embedding = "UMAP",
  imputeWeights = getImputeWeights(E16_2)
)

p$Ccnd1
p$Plxna2
p$Mybl1

p <- plotBrowserTrack(
  ArchRProj = E16_2, 
  groupBy = "Clusters", 
  geneSymbol = markerGenes, 
  upstream = 50000,
  downstream = 50000
)

grid::grid.newpage()
grid::grid.draw(p$Notch1)

p4 <- plotGroups(
  ArchRProj = E16_2, 
  groupBy = "Clusters", 
  colorBy = "GeneScoreMatrix", 
  name = markerGenes,
  plotAs = "violin",
  alpha = 0.4,
  addBoxPlot = TRUE
)
p4$Ccnd1
