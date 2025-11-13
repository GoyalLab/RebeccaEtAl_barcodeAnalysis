#Code to compare clonal size and gene expression with and without Dasatinib treatment. 
#There are two plots - clone size distribution comparison and 

# Inputs to the code that might need to be modified

#Path to the ALRA imputed rds file - output of ALRA_imputation.R script
inputDirectory <- 'Z:/Basic_Sciences/CDB/GoyalLab/People/KeerthanaArun/vitoRebecca_project/plotData/'
#Path to save the plots/plot data
outputDirectory <-"Z:/Basic_Sciences/CDB/GoyalLab/People/KeerthanaArun/vitoRebecca_project/plotData/"

# Importing the libraries
library(Seurat)
BiocManager::install('EnhancedVolcano')
library(ggrepel)
library(EnhancedVolcano)
library(SeuratObject)
library(ggplot2)
library(tidyverse)
library(patchwork)
library(dplyr)
library(gridExtra)
library(ggplot2)
library(tidyverse)
library(patchwork)
library(grid)
library(cowplot)
library(flextable)
options(future.globals.maxSize = 4000 * 1024^2)

Dasitinib_gene<- c(
  "SRC", "YES", "LYN", "FYN", "KCK", "CSK", "BTK", "EGFR", "PTK9", "PCTK3",
  "PRKDC", "MAPKAPK2", "p38", "STK25", "RSK2", "eIF2A", "PIM3", "PKAC", 
  "PKC", "STK6", "CDK2", "PKN2", "Frk", "DDR1", "ABL2", "SIK2", "RIPK2", 
  "EPHA2", "EPHB2", "CRKL", "AK2", "MAPK26", "TP53RK", "RPS6KA3", "GSK3A", 
  "MLKL", "DCLK3", "NADK2", "DTYMK", "WNK4", "PGK1", "LAT1", "IRAK4", "AXL", 
  "MAP3K2", "JAK1", "NEK9", "LZTR1"
)

#Reading the input Seurat object which contains information about scRNAseq, Sample number, Barcode and cluster
scRNA_Barcode_Sample_Cluster <- readRDS(paste0(inputDirectory, 'scRNA_Barcode_Sample_Cluster_ALRA_UMI_2.rds'))

scRNA_rescaled <- ScaleData(scRNA_Barcode_Sample_Cluster, features = Dasitinib_gene, assay = "alra")
genes_interested = (intersect(rownames(scRNA_rescaled[["RNA"]]), Dasitinib_gene))

DoHeatmap(object = scRNA_rescaled, group.by = "Sample",assay = "alra",features= genes_interested) + theme(axis.text.y = element_text(size = 5))+ theme(plot.margin = margin(t = 20, r = 10, b = 10, l = 10))


# Get cluster sizes along with cluster number and Sample ID
cluster_sizes <- data.frame(Cluster = scRNA_Barcode_Sample_Cluster$Cluster, Sample = scRNA_Barcode_Sample_Cluster$Sample)
cluster_sizes$ClusterSize <- ave(cluster_sizes$Sample, cluster_sizes$Cluster, FUN = length)
cluster_sizes$ClusterSize <- (as.numeric(cluster_sizes$ClusterSize))
cluster_sizes_copy <- unique(cluster_sizes)

temp_S1 <- subset(cluster_sizes, Sample == "S1")
temp_S1 <- temp_S1[order(-temp_S1$ClusterSize), ]

temp_S2 <- subset(cluster_sizes, Sample == "S2")
temp_S2 <- temp_S2[order(-temp_S2$ClusterSize), ]


# Identify singlet clones in Sample 2
singlets_S2 <- cluster_sizes$Cluster[cluster_sizes$ClusterSize == 1 & cluster_sizes$Sample == "S2"]
singlets_S1 <- cluster_sizes$Cluster[cluster_sizes$ClusterSize == 1 & cluster_sizes$Sample == "S1"]

S1_df <- subset(cluster_sizes_copy, Sample == "S1")
S1_df <- S1_df[order(-S1_df$ClusterSize), ]

#Volume of Sample 1
V1 = 1889.05
#Volume of Sample 2
V2 = 749.96925

S1_df$ClusterSize <- S1_df$ClusterSize*V1/V2
S2_df <- subset(cluster_sizes_copy, Sample == "S2")
S2_df <- S2_df[order(-S2_df$ClusterSize), ]
res <- wilcox.test(S2_df$ClusterSize, S1_df$ClusterSize)
mean(S1_df$ClusterSize)
mean(S2_df$ClusterSize)

######PLOT to compare clone sizes in the treated with Dasatinib and control
ggplot(S1_df, aes(x = S1_df$ClusterSize,y = after_stat(count/sum(count)))) +
  geom_histogram(binwidth = 1, fill = "hotpink3", color = "black") +
  labs(x = "Clone Size", y = "Frequency", title = "Clone Size Distribution in Control\n(only BRAFi/MEKi )") +
  theme_minimal() +  theme(plot.title = element_text(hjust = 0.5)) + 
  scale_x_continuous(breaks = seq(0, max(S1_df$ClusterSize), by = 5)) 

ggplot(S2_df, aes(x = S2_df$ClusterSize,y = after_stat(count/sum(count)))) +
  geom_histogram(binwidth = 1, fill = "turquoise3", color = "black") +
  labs(x = "Clone Size", y = "Frequency", title = "Clone Size Distribution in \nBRAFi/MEKi + Dasatinib") +
  theme_minimal() +  theme(plot.title = element_text(hjust = 0.5)) + 
  scale_x_continuous(breaks = seq(0, 30, by = 5))  + xlim(0, max(S1_df$ClusterSize))


write.csv(S1_df, file=paste0(outputDirectory,'controlSampleCLoneSize.csv'))
write.csv(S2_df, file=paste0(outputDirectory, 'treatedSampleCLoneSize.csv'))

clonal_class_data <- data.frame(
  Cluster = c(S1_df$Cluster[1], S1_df$Cluster[2], S1_df$Cluster[3], S1_df$Cluster[4], S1_df$Cluster[5], S2_df$Cluster[1], S2_df$Cluster[2], S2_df$Cluster[3], S2_df$Cluster[4], S2_df$Cluster[5]),
  clonalClass = c("S1_R1", "S1_R2", "S1_R3", "S1_R4","S1_R5", "S2_R1", "S2_R2", "S2_R3", "S2_R4", "S2_R5")
)

# Set the "ClonalClass" column of the Seurat object based on the identified properties
scRNA_Barcode_Sample_Cluster$ClonalClass <- clonal_class_data$clonalClass[match(scRNA_Barcode_Sample_Cluster$Cluster, clonal_class_data$Cluster)]
scRNA_Barcode_Sample_Cluster$ClonalClass <- ifelse(scRNA_Barcode_Sample_Cluster$Cluster %in% singlets_S2, "singlets_S2", as.character(scRNA_Barcode_Sample_Cluster$ClonalClass))
scRNA_Barcode_Sample_Cluster$ClonalClass <- ifelse(scRNA_Barcode_Sample_Cluster$Cluster %in% singlets_S1, "singlets_S1", as.character(scRNA_Barcode_Sample_Cluster$ClonalClass))


# Setting Idents of the Seurat object
Idents(scRNA_Barcode_Sample_Cluster) <- scRNA_Barcode_Sample_Cluster$ClonalClass
print(Idents(scRNA_Barcode_Sample_Cluster))

##################
#Identifying Differentially Expressed genes between the largest clone of Sample 1 and the largest clone of Sample 2
deg_S1R1_S2R1 <- FindMarkers(scRNA_Barcode_Sample_Cluster, ident.1 = "S1_R1",ident.2 = "S2_R1", log = 0.1)
S1R1_S2R1 <- deg_S1R1_S2R1[intersect(rownames(deg_S1R1_S2R1),Dasitinib_gene), ]
S1R1_S2R1_FC <- FoldChange(scRNA_Barcode_Sample_Cluster, ident.1 = "S2_R1", ident.2 = "S1_R1", features=rownames(S1R1_S2R1 ))

FoldChange(scRNA_Barcode_Sample_Cluster, ident.1 = "S1_R3", ident.2 = "S2_R3", features=rownames(S1R1_S2R1 ))
FoldChange(scRNA_Barcode_Sample_Cluster, ident.1 = "S1_R4", ident.2 = "S2_R4", features=rownames(S1R1_S2R1 ))
FoldChange(scRNA_Barcode_Sample_Cluster, ident.1 = "S1_R5", ident.2 = "S2_R5", features=rownames(S1R1_S2R1 ))
VlnPlot(scRNA_Barcode_Sample_Cluster,features = rownames(S1R1_S2R1), flip = TRUE, stack = TRUE, fill.by = "ident", idents = c("S1_R1", "S2_R1"), cols=c("turquoise3", "hotpink3"))+ theme(legend.position = 'none') + theme(axis.text.x = element_text(angle = 0)) + theme(plot.margin = margin(t = 10, l = 5, r=5))

S1R1_S2R1_FC$Gene_Name <- rownames(S1R1_S2R1_FC)
ft_1 <- flextable(S1R1_S2R1_FC, col_keys = c("Gene_Name", "avg_log2FC" ,"pct.1","pct.2"))
ft_1 <- set_caption(ft_1, caption = "Log Fold change in genes between the largest clones in Sample 1 and Sample 2")
ft_1

Singlets_FC_same <- FoldChange(scRNA_Barcode_Sample_Cluster, ident.1 = "singlets_S2", ident.2 = "singlets_S1", features=rownames(S1R1_S2R1 ))
Singlets_FC_same$Gene_Name <- rownames(Singlets_FC_same)
flextable(Singlets_FC_same)

##############
#Identifying Differentially Expressed genes between the largest clone of Sample 2 and the singlets of Sample 2
deg_S2R1_S2singlets <- FindMarkers(scRNA_Barcode_Sample_Cluster, ident.1 = "S2_R1",ident.2 = "singlets_S2", test.use = "LR")
S2R1_S2singlets <- deg_S2R1_S2singlets[intersect(rownames(deg_S2R1_S2singlets),Dasitinib_gene), ]
FoldChange(scRNA_Barcode_Sample_Cluster, ident.1 = "S2_R1", ident.2 = "singlets_S2", features=rownames(S2R1_S2singlets ))
FoldChange(scRNA_Barcode_Sample_Cluster, ident.1 = "singlets_S1", ident.2 = "S1_R1", features=rownames(S2R1_S2singlets ))

VlnPlot(scRNA_Barcode_Sample_Cluster,features = rownames(S2R1_S2singlets), idents = c("S2_R1", "singlets_S2"))

##############
#Identifying Differentially Expressed genes between the second largest clones of Sample 1 and Sample 2
deg_S1R2_S2R2 <- FindMarkers(scRNA_Barcode_Sample_Cluster, ident.1 = "S1_R2",ident.2 =  "S2_R2", test.use = "LR")
S1R2_S2R2 <- deg_S1R2_S2R2[intersect(rownames(deg_S1R2_S2R2),Dasitinib_gene), ]
print(S1R2_S2R2)
FoldChange(scRNA_Barcode_Sample_Cluster, ident.1 = "S1_R2", ident.2 = "S2_R2", features=rownames(S1R1_S2R1 ))
VlnPlot(scRNA_Barcode_Sample_Cluster,features = rownames(S1R1_S2R1), idents = c("S1_R2", "S2_R2"))

################
set_1 <- c("S1_R1", "S1_R2", "S1_R3", "S1_R4", "S1_R5")
set_2 <- c("S2_R1", "S2_R2", "S2_R3", "S2_R4", "S2_R5")
scRNA_Top5 = scRNA_Barcode_Sample_Cluster

cells.use <- WhichCells(scRNA_Top5, idents = set_1)
scRNA_Top5 <- SetIdent(scRNA_Top5, cells = cells.use, value = "Sample 1")
cells.use <- WhichCells(scRNA_Top5, idents = set_2)
scRNA_Top5 <- SetIdent(scRNA_Top5, cells = cells.use, value = "Sample 2")
Idents(scRNA_Top5)
################################
deg_S1T5_S2T5 <- FindMarkers(scRNA_Top5, ident.1 = "Sample 2",ident.2 = "Sample 1", test.use = "LR", log = 0.1)
S1T5_S2T5 <- deg_S1T5_S2T5[intersect(rownames(deg_S1T5_S2T5),Dasitinib_gene), ]
S1T5_S2T5_FC <- FoldChange(scRNA_Top5, ident.1 = "Sample 2", ident.2 = "Sample 1", features=rownames(S1T5_S2T5))
S1T5_S2T5_FC$Gene_Name <- rownames(S1T5_S2T5_FC)
VlnPlot(scRNA_Top5, features = rownames(S1R1_S2R1), idents = c("Sample 2", "Sample 1"), flip = TRUE, stack = TRUE, fill.by = "ident", cols=c("turquoise3", "hotpink3"),pt.size = 0)+ theme(legend.position = 'none') + theme(axis.text.x = element_text(angle = 0))
VlnPlot(scRNA_Top5, features = rownames(S1T5_S2T5), idents = c("Sample 2", "Sample 1"), flip = TRUE, stack = TRUE, fill.by = "ident", cols=c("turquoise3", "hotpink3"),pt.size = 0)+ theme(legend.position = 'none') + theme(axis.text.x = element_text(angle = 0))

#################
#Identifying Differentially Expressed genes between the singlets of Sample 2 and the 5 largest clones of Sample 1
deg_S1T5_S2singlets <- FindMarkers(scRNA_Top5, ident.1 = "singlets_S2",ident.2 = "Sample 1", logfc.threshold = 0.1, assay = "alra")
write.csv(deg_S1T5_S2singlets,file='~/CollabVito/temp_analysis/newAnalysis/UMI_3/S1T5vsS2singlets_allDEG.csv')
S1T5_S2singlets_FC <- deg_S1T5_S2singlets[intersect(rownames(deg_S1T5_S2singlets),Dasitinib_gene), ]
S1T5_S2singlets_FC$Gene_Name <- rownames(S1T5_S2singlets_FC)
S1T5_S2singlets_FC <- S1T5_S2singlets_FC %>%  arrange(desc(abs(avg_log2FC)))

VlnPlot(scRNA_Top5, features = rownames(S1T5_S2singlets_FC), idents = c("singlets_S2", "Sample 1"), flip = TRUE, stack = TRUE, fill.by = "ident", cols=c("turquoise3", "hotpink3"),pt.size = 0)+ theme(legend.position = 'none') + theme(axis.text.x = element_text(angle = 0))
ft_maxdiff <- flextable(S1T5_S2singlets_FC, col_keys = c("Gene_Name", "avg_log2FC" ,"pct.1","pct.2"))
ft_maxdiff <- set_caption(ft_maxdiff, caption = "Log Fold change in genes between the 5 largest clones in Sample 1 and singlet clones of Sample 2")
ft_maxdiff
############################
#Identifying Differentially Expressed genes between the singlets of Sample 2 and the 5 largest clones of Sample 2
deg_S2T5_S2singlets <- FindMarkers(scRNA_Top5, ident.1 = "singlets_S2",ident.2 = "Sample 2",  logfc.threshold = 0.1)
S2T5_S2singlets_FC <- deg_S2T5_S2singlets[intersect(rownames(deg_S2T5_S2singlets),Dasitinib_gene), ]
FoldChange(scRNA_Top5, ident.1 = "singlets_S2", ident.2 = "Sample 2", features=rownames(S2T5_S2singlets))
S2T5_S2singlets_FC$Gene_Name <- rownames(S2T5_S2singlets_FC)
S2T5_S2singlets_FC <- S2T5_S2singlets_FC %>%  arrange(desc(abs(avg_log2FC)))
VlnPlot(scRNA_Top5, features = rownames(S2T5_S2singlets_FC), idents = c("singlets_S2", "Sample 2"), flip = TRUE, stack = TRUE, fill.by = "ident", cols=c("turquoise3", "hotpink3"),pt.size = 0)+ theme(legend.position = 'none') + theme(axis.text.x = element_text(angle = 0))
ft_subS2 <- flextable(S2T5_S2singlets_FC, col_keys = c("Gene_Name", "avg_log2FC" ,"pct.1","pct.2"))
ft_subS2 <- set_caption(ft_subS2, caption = "Log Fold change in genes between the 5 largest clones of Sample 2 and the singlet clones of Sample 2")
ft_subS2

##########################
#Identifying Differentially Expressed genes between the singlets of Sample 1 and Sample 2
deg_S1singlets_S2singlets <- FindMarkers(scRNA_Barcode_Sample_Cluster, ident.1 = "singlets_S1",ident.2 = "singlets_S2", test.use = "LR",  logfc.threshold = 0.1)
S1singlets_S2singlets <- deg_S1singlets_S2singlets[intersect(rownames(deg_S1singlets_S2singlets),Dasitinib_gene), ]
Singlets_FC <- FoldChange(scRNA_Barcode_Sample_Cluster, ident.1 = "singlets_S2", ident.2 = "singlets_S1", features=rownames(S1singlets_S2singlets))
Singlets_FC$Gene_Name <- rownames(Singlets_FC)
Singlets_FC <-  Singlets_FC %>%  arrange(desc(abs(avg_log2FC)))
ft_singlets <- flextable(Singlets_FC,col_keys = c("Gene_Name", "avg_log2FC" ,"pct.1","pct.2"))
ft_singlets <- set_caption(ft_singlets, caption = "Log Fold change in genes between the singlet clones in Sample 1 and Sample 2")
ft_singlets
VlnPlot(scRNA_Barcode_Sample_Cluster, features = Singlets_FC$Gene_Name, idents = c("singlets_S1", "singlets_S2"), flip = TRUE, stack = TRUE, fill.by = "ident", cols=c("turquoise3", "hotpink3"),pt.size = 0)+ theme(legend.position = 'none') + theme(axis.text.x = element_text(angle = 0))

#Identifying Differentially Expressed genes between the singlets of Sample 1 and the singlets of Sample 2
deg_S1singlets_S2singlets <- FindMarkers(scRNA_Barcode_Sample_Cluster, ident.1 = "singlets_S1",ident.2 = "singlets_S2", test.use = "bimod",  logfc.threshold = 0.01)
S1singlets_S2singlets <- deg_S1R1_S2singlets[intersect(rownames(deg_S1R1_S2singlets),Dasitinib_gene), ]
FoldChange(scRNA_Barcode_Sample_Cluster, ident.1 = "singlets_S1", ident.2 = "singlets_S2", features=rownames(S1R1_S2Singlets ))
VlnPlot(scRNA_Barcode_Sample_Cluster,features = rownames(S1singlets_S2singlets), idents = c("singlets_S1", "singlets_S2"), log = TRUE)



VlnPlot(scRNA_Barcode_Sample_Cluster,features = rownames(S2R1_S1singlets)[8:10], idents = c("singlets_S1", "S2_R2"), log = TRUE)
#################################
#Comparing the two samples overall
scRNA_overall = scRNA_Barcode_Sample_Cluster
Idents(scRNA_overall) <- scRNA_overall$Sample
print(Idents(scRNA_overall))
deg_S1_S2 <- FindMarkers(scRNA_overall, ident.1 = "S2",ident.2 = "S1", log = 0.1)
write.csv(deg_S1_S2,file='~/CollabVito/temp_analysis/newAnalysis/UMI_3/S1vsS2_allDEG.csv')
S1_S2_FC_new <- deg_S1_S2[intersect(rownames(deg_S1_S2),Dasitinib_gene),]
S1_S2_FC_new$Gene_Name <- rownames(S1_S2_FC_new)
S1_S2_FC_new <- S1_S2_FC_new %>% arrange(desc(abs(avg_log2FC)))

########################################################################
#Creating tables to show log fold change

plt_vln_plot <- function(seurat_object, gene_name, group_ident, log_scale = FALSE, ncol = 1) {
  # Check if the gene is present in the features of the Seurat object
  
  # Create the VlnPlot
  vln_plot <- VlnPlot(
    object = seurat_object,
    features = gene_name,
    idents = group_ident,
    log = log_scale,
    ncol = ncol,
    pt.size = 0.1
  ) + theme(legend.position = 'none') + theme(axis.text.x = element_text(angle = 0))
  return(vln_plot)
}

library(data.table)

z <- as.data.table(S2T5_S2singlets_FC)
z <- z[, .(Gene_Name, avg_log2FC, p_val_adj,pct.1,pct.2,violinPlot = lapply(Gene_Name, function(gene) plt_vln_plot(scRNA_Top5, gene, c("Sample 2", "singlets_S2"))))]

ft <- flextable(z)
ft <- mk_par(ft, j = "violinPlot",
             value = as_paragraph(
               gg_chunk(value = violinPlot, width = 3, height = 2)
             ))
ft <- fontsize(ft, size = 20, part = "body")
ft <- fontsize(ft, size = 15, part = "header")
ft <- set_caption(ft, caption = "Log Fold change in genes affected by cells of top 5 clones of DT+D vs all singlets of DT+D")
#ft <- bold(ft, ~ avg_log2FC > 0.1, bold = TRUE)
#ft <- bold(ft, ~ avg_log2FC < -0.1, bold = TRUE)
#ft <- add_footer_lines(ft, "The genes in bold are differentially expressed.")
ft

#################
#Volcano plots for genes for S1T5vsS2Singlets
deg_S1T5_S2singlets$geneName = rownames(deg_S1T5_S2singlets)
deg_S1T5_S2singlets$diffexpressed <- "Not significant or targets of Dasantinib"
deg_S1T5_S2singlets$diffexpressed[deg_S1T5_S2singlets$geneName %in% common_genes] <- "Downregulated Dasatinib Targets"
deg_S1T5_S2singlets$delabel <- NA
deg_S1T5_S2singlets$delabel[deg_S1T5_S2singlets$diffexpressed != "Not significant or targets of Dasantinib"] <- deg_S1T5_S2singlets$geneName[deg_S1T5_S2singlets$diffexpressed != "Not significant or targets of Dasantinib"]

S1T5_S2Singlets <- ggplot(data=deg_S1T5_S2singlets, aes(x=avg_log2FC, y=-log10(p_val_adj), col=diffexpressed, label=delabel)) +
  geom_point() + 
  geom_point(data = filter(deg_S1T5_S2singlets, diffexpressed == "Downregulated Dasatinib Targets"), alpha = 1) +
  geom_point(data = filter(deg_S1T5_S2singlets, diffexpressed == "Not significant or targets of Dasantinib"), alpha =0.001) +
  theme_minimal() +
  geom_text_repel(aes(color = diffexpressed),max.overlaps = 25) +   
  scale_color_manual(values = c("red", "lightgrey")) + NoLegend()
ggsave(S1T5_S2Singlets, file = paste0('~/CollabVito/temp_analysis/newAnalysis/UMI_2/S1T5_S2Singlets_VolcanoPlot.svg'), width = 6, height = 4)
svglite(filename = paste0('~/CollabVito/temp_analysis/newAnalysis/UMI_2/S1T5_S2Singlets_VolcanoPlot_6X4.svg'), width = 6, height = 4)
plot(S1T5_S2Singlets)
dev.off()
###############
#Volcano plots for genes for S1vsS2
deg_S1_S2$geneName = rownames(deg_S1_S2)
deg_S1_S2$diffexpressed <- "Not significant or targets of Dasantinib"
deg_S1_S2$diffexpressed[deg_S1_S2$geneName %in% common_genes] <- "Downregulated Dasatinib Targets"
deg_S1_S2$delabel <- NA
deg_S1_S2$delabel[deg_S1_S2$diffexpressed != "Not significant or targets of Dasantinib"] <- deg_S1_S2$geneName[deg_S1_S2$diffexpressed != "Not significant or targets of Dasantinib"]

S1_S2 <- ggplot(data=deg_S1_S2, aes(x=avg_log2FC, y=-log10(p_val_adj), col=diffexpressed, label=delabel)) +
  geom_point() + 
  geom_point(data = filter(deg_S1_S2, diffexpressed == "Downregulated Dasatinib Targets"), alpha = 1) +
  geom_point(data = filter(deg_S1_S2, diffexpressed == "Not significant or targets of Dasantinib"), alpha =0.001) +
  theme_minimal() +
  geom_text_repel(aes(color = diffexpressed),max.overlaps = 25) +   
  scale_color_manual(values = c("red", "lightgrey"))
svglite(filename = paste0('~/CollabVito/temp_analysis/newAnalysis/UMI_2/S1_S2_VolcanoPlot_6X6.svg'), width = 6, height = 6)
plot(S1_S2)
dev.off()










