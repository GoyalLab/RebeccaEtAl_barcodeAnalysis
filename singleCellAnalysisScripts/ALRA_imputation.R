#Code to perform imputation using ALRA on the single-cell RNAseq matrix

# Inputs to the code that might need to be modified
outputDirectory <- '~/CollabVito/temp_analysis/seuratOutputData/'
inputDirectory <-
  "~/CollabVito/inputData/cellRanger/pool1_filtered_feature_bc_matrix_GRCh38/"
input22Directory <- '/home/mzo5929/CollabVito/inputData/barcodes/'

# Importing the libraries
remove.packages('Matrix')
remotes::install_version("Matrix",version = "1.6.1.1")
library(Matrix)
library(Seurat)


library(SeuratObject)
library(Seurat)
library(ggplot2)
library(tidyverse)
library(patchwork)
BiocManager::install('glmGamPoi')
#devtools::install_github('satijalab/seurat-data')
library(SeuratData)
#remotes::install_github('satijalab/seurat-wrappers')
library(SeuratWrappers)


options(future.globals.maxSize = 40000 * 1024 ^ 2)
##################################################################################################

# Loading the sample11umPlx dataset
sample11umPlx.data <- Read10X(data.dir = inputDirectory)
sample11umPlx <-
  CreateSeuratObject(
    counts = sample11umPlx.data[["Gene Expression"]],
    project = "10X_Sample1_1uMPLX",
    min.cells = 3,
    min.features = 200
  )

#Calculating mitochondrial gene content
sample11umPlx[["percent.mt"]] <-
  PercentageFeatureSet(object = sample11umPlx, pattern = "^MT-")

#Quality Control and plots
VlnPlot(
  object = sample11umPlx,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3
)
plot1 <-
  FeatureScatter(object = sample11umPlx,
                 feature1 = "nCount_RNA",
                 feature2 = "percent.mt") + NoLegend()
plot2 <-
  FeatureScatter(object = sample11umPlx,
                 feature1 = "nCount_RNA",
                 feature2 = "nFeature_RNA") + NoLegend()

combined_plots <- plot1 + plot2
combined_plots
ggsave(
  combined_plots,
  file = paste0(outputDirectory, 'QC_plot_bothUMI.pdf'),
  width = 6,
  height = 4.321
)

sample11umPlx <-
  subset(x = sample11umPlx,
         subset = nFeature_RNA > 200 &
           nFeature_RNA < 7200 & percent.mt < 20)

#SCTransform to Seurat object
sample11umPlx <-
  SCTransform(object = sample11umPlx, variable.features.n = 5000, vst.flavor="v2", seed.use=1024)
#Imputation
sample11umPlx <- RunALRA(object = sample11umPlx, k=27)
DefaultAssay(sample11umPlx) <- "alra"
sample11umPlx <- ScaleData(sample11umPlx, assay = "alra")

sample11umPlx <- FindVariableFeatures(sample11umPlx, assay = 'alra')

sample11umPlx <-
  RunPCA(
    object = sample11umPlx,
    assay = "alra"
  )
sample11umPlx <-
  FindNeighbors(object = sample11umPlx,
                dims = 1:50,
                verbose = FALSE)
sample11umPlx <- FindClusters(object = sample11umPlx, verbose = FALSE)
sample11umPlx <- RunUMAP(object = sample11umPlx, dims = 1:50)
clusterUMAP = DimPlot(sample11umPlx, label = TRUE) + NoLegend()
clusterUMAP

#Saving the SCTransformed file
saveRDS(
  sample11umPlx,
  file = paste0(
    outputDirectory,
    'ScTransform_50pcs_filter_ALRA.rds'
  )
)
sample11umPlx <-
  readRDS(paste0(
    outputDirectory,
    'ScTransform_50pcs_filter_ALRA.rds'
  ))

#Reading the Sample information
samplesPDXSampleNum <-
  CreateSeuratObject(counts = sample11umPlx.data[["Antibody Capture"]])
samplesPDXSampleNum <-
  NormalizeData(samplesPDXSampleNum,
                assay = "RNA",
                normalization.method = "CLR")
samplesPDXSampleNum <-
  HTODemux(samplesPDXSampleNum,
           assay = "RNA",
           positive.quantile = 0.8)
samplesPDXSampleNum.singlet <-
  subset(samplesPDXSampleNum, idents = c("S1", "S2"))
FeatureScatter(samplesPDXSampleNum,
               feature1 = "S1",
               feature2 = "S2")
rm(samplesPDXSampleNum)
rm(sample11umPlx.data)
gc()

sampleAssignment <-
  as.matrix(samplesPDXSampleNum.singlet@active.ident)
cells_hash <-
  rownames(sampleAssignment) #CellIds with Sample number as prefix
sampleAssignment = as_tibble(sampleAssignment)
sampleAssignmentFinal = sampleAssignment %>% mutate(cellID = cells_hash)

#Aligning sample information with scRNAseq information using CellID
postFilter <- sampleAssignmentFinal$cellID
sample11umPlx_filtered <-
  sample11umPlx[, colnames(sample11umPlx) %in% postFilter]
cell_id_list <- colnames(sample11umPlx_filtered)
postFilter <- postFilter[postFilter %in% cell_id_list]

sampleAssignmentFinal = subset(sampleAssignmentFinal,
                               sampleAssignmentFinal$cellID %in% cell_id_list)
row.names(sampleAssignmentFinal) <- cell_id_list
colnames(sampleAssignmentFinal) <- c("Sample", "cellID")

select_temp <- (sampleAssignmentFinal %>% select(Sample))
row.names(select_temp) = cell_id_list

#Merging and creating a new Seurat Object with both information
sample11umPlx_Sample <- AddMetaData(object = sample11umPlx_filtered,
                                    metadata = as.vector(select_temp),
                                    col.name = "Sample")

# Creating Sample name labelled UMAP
sample_colors <- c("S1" = "turquoise3", "S2" = "hotpink3")
clusterUMAP <-
  DimPlot(
    sample11umPlx_Sample,
    reduction = 'umap',
    group.by = "Sample",
    cols = sample_colors,
    assay = "alra"
  )
ggsave(
  clusterUMAP,
  file = paste0(outputDirectory, 'Sample_1_2_umap.pdf'),
  width = 6,
  height = 4.321
)
clusterUMAP

#Saving the merged file and can be directly retrieved for future use.
saveRDS(sample11umPlx_Sample,
        file = paste0(outputDirectory, 'SCT_filtered_SampleID.rds'))
##################################################################################################
sample11umPlx_Sample <-
  readRDS(paste0(outputDirectory, 'SCT_filtered_SampleID.rds'))

#Adding Barcode Information
barcode40 = as_tibble(read.table(
  paste0(input22Directory, 'stepThreeStarcodeShavedReads.txt'),
  stringsAsFactors = F,
  header = T
))
barcode40 = barcode40  %>%
  select(cellID, BC40StarcodeD8, UMI) %>% unique() %>% select(-UMI)

umiCut = 3 # Minimum UMI cutoff for reliable analysis. The other option used is 2.
linCut = 1 # remove all the cases where one CellID is associated with X(one/two etc) barcodes.

barcode40_1 = barcode40 %>%
  group_by(cellID, BC40StarcodeD8) %>%
  summarise(nUMI = length(BC40StarcodeD8)) %>%
  filter(nUMI >= umiCut) %>%
  group_by(cellID) %>%
  mutate(nLineages = length(cellID))

###taking only single cellid-->barcode mappings
barcode40_2 = barcode40_1 %>%
  ungroup() %>%
  filter(nLineages <= linCut) %>%
  unique()

barcode40_2$cellID <- paste(barcode40_2$cellID, "-1", sep = "")

#Aligning barcode information with scRNAseq information using CellID
barcodeCellID <- barcode40_2$cellID

sample11umPlx_Sample <-
  sample11umPlx_Sample[, colnames(sample11umPlx_Sample) %in% barcodeCellID]
cell_id_list <- colnames(sample11umPlx_Sample)
barcodeCellID <- barcodeCellID[barcodeCellID %in% cell_id_list]

barcode40_fin = subset(barcode40_2, barcode40_2$cellID %in% barcodeCellID)
row.names(barcode40_fin) <- cell_id_list

select_col <- (barcode40_fin %>% select(BC40StarcodeD8))
row.names(select_col) = cell_id_list
names(select_col) <- NULL 

#Merging and creating a new Seurat Object with all information
sample11umPlx_new <- AddMetaData(object = sample11umPlx_Sample,
                                 metadata = as.vector(select_col),
                                 col.name = "Barcode")

Idents(object = sample11umPlx_new) <- sample11umPlx_new$Barcode

BarcodeSample <-
  as.data.frame(sample11umPlx_new@meta.data[, c("Barcode", "Sample")])
BarcodeSample <- BarcodeSample %>%
  rownames_to_column(var = "RowName")
BarcodeSample <- BarcodeSample %>%
  group_by(Sample, Barcode) %>%
  group_by(Sample) %>%
  mutate(ClusterIndex = dense_rank(Barcode)) %>%
  ungroup()
max_cluster_index_sample1 <-
  max(BarcodeSample$ClusterIndex[BarcodeSample$Sample == "S1"])
BarcodeSample <- BarcodeSample %>%
  mutate(
    ClusterIndex = ifelse(
      Sample == "S2",
      ClusterIndex + max_cluster_index_sample1,
      ClusterIndex
    )
  )

BarcodeSample <- BarcodeSample %>%
  column_to_rownames(var = "RowName")

scRNA_Barcode_Sample_Cluster <- AddMetaData(
  object = sample11umPlx_new,
  metadata = BarcodeSample$ClusterIndex,
  col.name = "Cluster"
)
saveRDS(
  scRNA_Barcode_Sample_Cluster,
  file = paste0(
    outputDirectory,
    'scRNA_Barcode_Sample_Cluster_ALRA_UMI_3.rds'
  )
)
