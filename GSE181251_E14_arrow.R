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

# sessionInfo()
# R version 4.3.1 (2023-06-16)
# Platform: aarch64-apple-darwin20 (64-bit)
# Running under: macOS Sonoma 14.0
# 
# Matrix products: default
# BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
# LAPACK: /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.11.0
# 
# Random number generation:
#   RNG:     L'Ecuyer-CMRG 
#  Normal:  Inversion 
#  Sample:  Rejection 
#  
# locale:
# [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# time zone: Asia/Tokyo
# tzcode source: internal
# 
# attached base packages:
#  [1] parallel  stats4    grid      stats     graphics  grDevices utils     datasets  methods  
# [10] base     
# 
# other attached packages:
#  [1] ggrepel_0.9.3                      hexbin_1.28.3                     
#  [3] uwot_0.1.16                        BSgenome.Mmusculus.UCSC.mm10_1.4.3
#  [5] BSgenome_1.68.0                    rtracklayer_1.60.1                
#  [7] Biostrings_2.68.1                  XVector_0.40.0                    
#  [9] rhdf5_2.44.0                       SummarizedExperiment_1.30.2       
# [11] Biobase_2.60.0                     MatrixGenerics_1.12.3             
# [13] Rcpp_1.0.11                        Matrix_1.6-1.1                    
# [15] GenomicRanges_1.52.0               GenomeInfoDb_1.36.4               
# [17] IRanges_2.34.1                     S4Vectors_0.38.2                  
# [19] BiocGenerics_0.46.0                matrixStats_1.0.0                 
# [21] data.table_1.14.8                  stringr_1.5.0                     
# [23] plyr_1.8.9                         magrittr_2.0.3                    
# [25] gtable_0.3.4                       gtools_3.9.4                      
# [27] gridExtra_2.3                      ArchR_1.0.2                       
# [29] patchwork_1.1.3                    ggplot2_3.4.3                     
# [31] Seurat_4.9.9.9060                  SeuratObject_4.9.9.9091           
# [33] sp_2.1-0                           Signac_1.10.0                     
# 
# loaded via a namespace (and not attached):
#   [1] RcppAnnoy_0.0.21         splines_4.3.1            later_1.3.1             
#   [4] BiocIO_1.10.0            bitops_1.0-7             tibble_3.2.1            
#   [7] polyclip_1.10-6          XML_3.99-0.14            fastDummies_1.7.3       
#  [10] lifecycle_1.0.3          globals_0.16.2           lattice_0.21-9          
#  [13] MASS_7.3-60              plotly_4.10.2            yaml_2.3.7              
#  [16] httpuv_1.6.11            sctransform_0.4.0        spam_2.9-1              
#  [19] spatstat.sparse_3.0-2    reticulate_1.32.0        cowplot_1.1.1           
#  [22] pbapply_1.7-2            RColorBrewer_1.1-3       abind_1.4-5             
#  [25] zlibbioc_1.46.0          Rtsne_0.16               purrr_1.0.2             
#  [28] RCurl_1.98-1.12          GenomeInfoDbData_1.2.10  irlba_2.3.5.1           
#  [31] listenv_0.9.0            spatstat.utils_3.0-3     goftest_1.2-3           
#  [34] RSpectra_0.16-1          spatstat.random_3.1-6    fitdistrplus_1.1-11     
#  [37] parallelly_1.36.0        leiden_0.4.3             codetools_0.2-19        
#  [40] DelayedArray_0.26.7      RcppRoll_0.3.0           tidyselect_1.2.0        
#  [43] farver_2.1.1             spatstat.explore_3.2-3   GenomicAlignments_1.36.0
#  [46] jsonlite_1.8.7           ellipsis_0.3.2           progressr_0.14.0        
#  [49] ggridges_0.5.4           survival_3.5-7           tools_4.3.1             
#  [52] ica_1.0-3                glue_1.6.2               dplyr_1.1.3             
#  [55] withr_2.5.1              fastmap_1.1.1            rhdf5filters_1.12.1     
#  [58] fansi_1.0.4              digest_0.6.33            R6_2.5.1                
#  [61] mime_0.12                colorspace_2.1-0         scattermore_1.2         
#  [64] tensor_1.5               spatstat.data_3.0-1      utf8_1.2.3              
#  [67] tidyr_1.3.0              generics_0.1.3           httr_1.4.7              
#  [70] htmlwidgets_1.6.2        S4Arrays_1.0.6           pkgconfig_2.0.3         
#  [73] lmtest_0.9-40            htmltools_0.5.6.1        dotCall64_1.0-2         
#  [76] scales_1.2.1             png_0.1-8                rstudioapi_0.15.0       
#  [79] reshape2_1.4.4           rjson_0.2.21             nlme_3.1-163            
#  [82] zoo_1.8-12               KernSmooth_2.23-22       miniUI_0.1.1.1          
#  [85] restfulr_0.0.15          pillar_1.9.0             vctrs_0.6.3             
#  [88] RANN_2.6.1               promises_1.2.1           xtable_1.8-4            
#  [91] cluster_2.1.4            cli_3.6.1                compiler_4.3.1          
#  [94] Rsamtools_2.16.0         rlang_1.1.1              crayon_1.5.2            
#  [97] future.apply_1.11.0      labeling_0.4.3           stringi_1.7.12          
# [100] viridisLite_0.4.2        deldir_1.0-9             BiocParallel_1.34.2     
# [103] munsell_0.5.0            lazyeval_0.2.2           spatstat.geom_3.2-5     
# [106] RcppHNSW_0.5.0           future_1.33.0            Rhdf5lib_1.22.1         
# [109] shiny_1.7.5              ROCR_1.0-11              igraph_1.5.1            
# [112] fastmatch_1.1-4 
