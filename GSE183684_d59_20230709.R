library(dplyr)
library(Seurat)
library(patchwork)

rna <- ReadMtx(
  mtx = "GSM5567526_d59_matrix.mtx.gz", 
  features = "GSM5567526_d59_features.tsv.gz",
  cells = "GSM5567526_d59_barcodes.tsv.gz",
  feature.column = 2
)

pbmc <- CreateSeuratObject(counts = rna, min.cells = 3, min.features = 200)

pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

pbmc <- subset(pbmc, subset = nFeature_RNA > 1000 & nFeature_RNA < 4000 & percent.mt < 5)

pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

pbmc <- NormalizeData(pbmc)

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
ElbowPlot(pbmc)

pbmc <- FindNeighbors(pbmc, dims = 1:15)
pbmc <- FindClusters(pbmc, resolution = 0.10)

pbmc <- RunUMAP(pbmc, dims = 1:15)
DimPlot(pbmc, reduction = "umap")

saveRDS(pbmc, file = "GSE183684_d59_QC.rds")

FeaturePlot(pbmc, features = c("NOTCH1", "NOTCH2", "NOTCH3", "NOTCH4"))
DotPlot(pbmc, features = c("NOTCH1", "NOTCH2", "NOTCH3", "NOTCH4"))
DotPlot(pbmc, features = c("JAG1", "JAG2", "DLL1", "DLL3", "DLL4"))
DotPlot(pbmc, features = c("HES1", "HES5", "HEY1", "HEY2"))

# annotation according to Finkbeiner et al. Cell Rep. 2022
# this needs correction
FeaturePlot(pbmc, features = c("CCND1",   # MPC
                               "ATOH7",   # Npre
                               "POU4F2",  # RGC
                               "PTF1A",   # HOR
                               "CRX",     # CON
                               "ONECUT1", # AMA/HC
                               "VSX1",    # BIP 
                               "NR2E3"    # ROD 
))

DotPlot(pbmc, features = c("CCND1",   # MPC
                           "ELAVL4",  # RGC # PMID: 33674582
                           "PLXNA2",  # RPC # PMID: 28440030
                           "CPLX2",   # AMA # PMID: 19386896
                           "OTX2"     # CON/ROD
                           ))


cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, min.pct = 0.25)
head(cluster5.markers, n=50)
# write.table (cluster0.markers, file = "cluster0_markers.txt")

new.cluster.ids <- c("MPC", "RGC1", "RGC2","RPC", "AMA", "CON/ROD")
names(new.cluster.ids) <- levels(pbmc)
pbmc_whole <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc_whole, reduction = "umap", label = TRUE)

DotPlot(pbmc_whole, features = c("CCND1",   # MPC
                           "ELAVL4",  # RGC # PMID: 33674582
                           "PLXNA2",  # RPC # PMID: 28440030
                           "CPLX2",   # AMA # PMID: 19386896
                           "OTX2"     # CON/ROD
))


DotPlot(pbmc_whole, features = c("NOTCH1", "NOTCH2", "NOTCH3", "NOTCH4"))
DotPlot(pbmc_whole, features = c("JAG1", "JAG2", "DLL1", "DLL3", "DLL4"))
DotPlot(pbmc_whole, features = c("HES1", "HES5", "HEY1", "HEY2"))


###########################


library(ArchR)
#library(parallel)
#library(magick)
set.seed(1)

addArchRGenome("hg38")

genomeAnnotation <- createGenomeAnnotation(genome = "hg38")

inputFiles <- c("d59"="fragments.tsv.gz")

addArchRThreads(1)

ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  minTSS = 4,
  minFrags = 1000, 
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)

doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10, 
  knnMethod = "UMAP", 
  LSIMethod = 1
)

proj <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "d59",
  copyArrows = TRUE
)

getAvailableMatrices(proj)

proj <- filterDoublets(ArchRProj = proj)

proj <- saveArchRProject(ArchRProj = proj)

proj <- addIterativeLSI(ArchRProj = proj, useMatrix = "TileMatrix", name = "IterativeLSI")

proj <- addClusters(input = proj, reducedDims = "IterativeLSI")

proj <- addUMAP(ArchRProj = proj, reducedDims = "IterativeLSI")

p2 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
p2


genes_OF_interest <- c("CCND1","PLXNA2","CPLX2",
                       "NOTCH1", "NOTCH2", "NOTCH3", "NOTCH4",
                       "JAG1", "JAG2", "DLL1", "DLL3", "DLL4",
                       "HES1", "HES5", "HEY1", "HEY2")

p <- plotEmbedding(
  ArchRProj = proj, 
  colorBy = "GeneScoreMatrix", 
  name = genes_OF_interest, 
  embedding = "UMAP",
  imputeWeights = getImputeWeights(proj)
)

p$CCND1
p$NOTCH1
p$NOTCH2
p$NOTCH3
p$NOTCH4
p$JAG1
p$JAG2
p$DLL1
p$DLL3
p$DLL4
p$HES1
p$HES5
p$HEY1
p$HEY2

proj <- saveArchRProject(ArchRProj = proj)
# proj <- loadArchRProject(path = "./d59/", force = FALSE, showLogo = TRUE)

proj_gene <- getFeatures(proj, useMatrix = "GeneScoreMatrix")

p3 <- plotGroups(
  ArchRProj = proj, 
  groupBy = "Clusters", 
  colorBy = "GeneScoreMatrix", 
  name = c("CCND1","NOTCH1"),
  plotAs = "ridges",
  alpha = 0.4,
  addBoxPlot = TRUE
)
p3

p4 <- plotGroups(
  ArchRProj = proj, 
  groupBy = "Clusters", 
  colorBy = "GeneScoreMatrix", 
  name = c("CCND1", "PLXNA2", "CPLX2"),
  plotAs = "violin",
  alpha = 0.4,
  addBoxPlot = TRUE
)
p4


p4 <- plotGroups(
  ArchRProj = proj, 
  groupBy = "Clusters", 
  colorBy = "GeneScoreMatrix", 
  name = c("NOTCH1","NOTCH2","NOTCH3","NOTCH4"),
  plotAs = "violin",
  alpha = 0.4,
  addBoxPlot = TRUE
)
p4

p4 <- plotGroups(
  ArchRProj = proj, 
  groupBy = "Clusters", 
  colorBy = "GeneScoreMatrix", 
  name = c("JAG1","JAG2","DLL1","DLL3","DLL4"),
  plotAs = "violin",
  alpha = 0.4,
  addBoxPlot = TRUE
)
p4

p4 <- plotGroups(
  ArchRProj = proj, 
  groupBy = "Clusters", 
  colorBy = "GeneScoreMatrix", 
  name = c("HES1","HES5","HEY1","HEY2"),
  plotAs = "violin",
  alpha = 0.4,
  addBoxPlot = TRUE
)
p4





