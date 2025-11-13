library(ggplot2)
library(tidyverse)
library(dplyr)
library(svglite)
library(ggrepel)
options(future.globals.maxSize = 4000 * 1024^2)

inputFolder <- "Z:/Basic_Sciences/CDB/GoyalLab/People/KeerthanaArun/vitoRebecca_project/RebeccaEtAl_barcodeAnalysis/plotData/differentiallyExpressedGenes//"
outputFolder <- "Z:/Basic_Sciences/CDB/GoyalLab/People/KeerthanaArun/vitoRebecca_project/RebeccaEtAl_barcodeAnalysis/plots/differentiallyExpressedGenes/"

#Dasatinib target genes that were downregulated in all 4 cases: Control vs treatment and 
# 5 biggest clones of control vs singlet clones of Dasatinib treated sample
# for UMI cutoff 2 and 3
Dasitinib_gene<- c(
  "SRC", "YES", "LYN", "FYN", "KCK", "CSK", "BTK", "EGFR", "PTK9", "PCTK3",
  "PRKDC", "MAPKAPK2", "p38", "STK25", "RSK2", "eIF2A", "PIM3", "PKAC", 
  "PKC", "STK6", "CDK2", "PKN2", "Frk", "DDR1", "ABL2", "SIK2", "RIPK2", 
  "EPHA2", "EPHB2", "CRKL", "AK2", "MAPK26", "TP53RK", "RPS6KA3", "GSK3A", 
  "MLKL", "DCLK3", "NADK2", "DTYMK", "WNK4", "PGK1", "LAT1", "IRAK4", "AXL", 
  "MAP3K2", "JAK1", "NEK9", "LZTR1"
)

commonGeneList <- c("CDK2", "WNK4", "EPHB2", "PRKDC","EGFR", "DTYMK", "RPS6KA3", "STK25",
                    "DDR1", "TP53RK", "NEK9", "PKN2")

#UMI cutoff = 2
#compare control vs treatment - all clones
S1_S2_FC = read.csv(paste0(inputFolder, "UMI2_Imputed/", "S1vsS2_allDEG.csv"))
S1_S2_FC$geneName <- S1_S2_FC$X
S1_S2_FC <- S1_S2_FC %>% arrange(desc(abs(avg_log2FC)))
S1_S2_FC$diffexpressed <- "Not significant or targets of Dasantinib"
S1_S2_FC$diffexpressed[S1_S2_FC$geneName %in% commonGeneList] <- "Downregulated Dasatinib Targets"
S1_S2_FC$delabel <- NA
S1_S2_FC$delabel[S1_S2_FC$diffexpressed != "Not significant or targets of Dasantinib"] <- S1_S2_FC$geneName[S1_S2_FC$diffexpressed != "Not significant or targets of Dasantinib"]

S1_S2 <- ggplot(data=S1_S2_FC, aes(x=avg_log2FC, y=-log10(p_val_adj), col=diffexpressed, label=delabel)) +
  geom_point() + 
  geom_point(data = filter(S1_S2_FC, diffexpressed == "Downregulated Dasatinib Targets"), alpha = 1) +
  geom_point(data = filter(S1_S2_FC, diffexpressed == "Not significant or targets of Dasantinib"), alpha =0.001) +
  theme_minimal() +
  geom_text_repel(aes(color = diffexpressed),max.overlaps = 25) +   
  scale_color_manual(values = c("red", "lightgrey"))
plot(S1_S2)
svglite(filename = paste0(outputFolder, 'UMI2_Imputed/S1_S2_VolcanoPlot_6X6.svg'), width = 6, height = 6)
plot(S1_S2)
dev.off()

##Compare 5 biggest clones of control vs singlet clones of Dasatinib treated sample
S1T5vsS2singlets_allDEG = read.csv(paste0(inputFolder, "UMI2_Imputed/", "S1T5vsS2singlets_allDEG.csv"))
S1T5vsS2singlets_allDEG$geneName <- S1T5vsS2singlets_allDEG$X
S1T5vsS2singlets_allDEG <- S1T5vsS2singlets_allDEG %>% arrange(desc(abs(avg_log2FC)))
S1T5vsS2singlets_allDEG$diffexpressed <- "Not significant or targets of Dasantinib"
S1T5vsS2singlets_allDEG$diffexpressed[S1T5vsS2singlets_allDEG$geneName %in% commonGeneList] <- "Downregulated Dasatinib Targets"
S1T5vsS2singlets_allDEG$delabel <- NA
S1T5vsS2singlets_allDEG$delabel[S1T5vsS2singlets_allDEG$diffexpressed != "Not significant or targets of Dasantinib"] <- S1T5vsS2singlets_allDEG$geneName[S1T5vsS2singlets_allDEG$diffexpressed != "Not significant or targets of Dasantinib"]

S1T5vsS2singlets <- ggplot(data=S1T5vsS2singlets_allDEG, aes(x=avg_log2FC, y=-log10(p_val_adj), col=diffexpressed, label=delabel)) +
  geom_point() + 
  geom_point(data = filter(S1T5vsS2singlets_allDEG, diffexpressed == "Downregulated Dasatinib Targets"), alpha = 1) +
  geom_point(data = filter(S1T5vsS2singlets_allDEG, diffexpressed == "Not significant or targets of Dasantinib"), alpha =0.001) +
  theme_minimal() +
  geom_text_repel(aes(color = diffexpressed),max.overlaps = 25) +   
  scale_color_manual(values = c("red", "lightgrey"))
plot(S1T5vsS2singlets)
svglite(filename = paste0(outputFolder, 'UMI2_Imputed/S1T5vsS2singlets_VolcanoPlot_6X6.svg'), width = 6, height = 6)
plot(S1T5vsS2singlets)
dev.off()

#UMI cutoff = 3
#compare control vs treatment - all clones
S1_S2_FC = read.csv(paste0(inputFolder, "UMI3_Imputed/", "S1vsS2_allDEG.csv"))
S1_S2_FC$geneName <- S1_S2_FC$X
S1_S2_FC <- S1_S2_FC %>% arrange(desc(abs(avg_log2FC)))
S1_S2_FC$diffexpressed <- "Not significant or targets of Dasantinib"
S1_S2_FC$diffexpressed[S1_S2_FC$geneName %in% commonGeneList] <- "Downregulated Dasatinib Targets"
S1_S2_FC$delabel <- NA
S1_S2_FC$delabel[S1_S2_FC$diffexpressed != "Not significant or targets of Dasantinib"] <- S1_S2_FC$geneName[S1_S2_FC$diffexpressed != "Not significant or targets of Dasantinib"]

S1_S2 <- ggplot(data=S1_S2_FC, aes(x=avg_log2FC, y=-log10(p_val_adj), col=diffexpressed, label=delabel)) +
  geom_point() + 
  geom_point(data = filter(S1_S2_FC, diffexpressed == "Downregulated Dasatinib Targets"), alpha = 1) +
  geom_point(data = filter(S1_S2_FC, diffexpressed == "Not significant or targets of Dasantinib"), alpha =0.001) +
  theme_minimal() +
  geom_text_repel(aes(color = diffexpressed),max.overlaps = 25) +   
  scale_color_manual(values = c("red", "lightgrey"))
plot(S1_S2)
svglite(filename = paste0(outputFolder, 'UMI3_Imputed/S1_S2_VolcanoPlot_6X6.svg'), width = 6, height = 6)
plot(S1_S2)
dev.off()

##Compare 5 biggest clones of control vs singlet clones of Dasatinib treated sample
S1T5vsS2singlets_allDEG = read.csv(paste0(inputFolder, "UMI3_Imputed/", "S1T5vsS2singlets_allDEG.csv"))
S1T5vsS2singlets_allDEG$geneName <- S1T5vsS2singlets_allDEG$X
S1T5vsS2singlets_allDEG <- S1T5vsS2singlets_allDEG %>% arrange(desc(abs(avg_log2FC)))
S1T5vsS2singlets_allDEG$diffexpressed <- "Not significant or targets of Dasantinib"
S1T5vsS2singlets_allDEG$diffexpressed[S1T5vsS2singlets_allDEG$geneName %in% commonGeneList] <- "Downregulated Dasatinib Targets"
S1T5vsS2singlets_allDEG$delabel <- NA
S1T5vsS2singlets_allDEG$delabel[S1T5vsS2singlets_allDEG$diffexpressed != "Not significant or targets of Dasantinib"] <- S1T5vsS2singlets_allDEG$geneName[S1T5vsS2singlets_allDEG$diffexpressed != "Not significant or targets of Dasantinib"]

S1T5vsS2singlets <- ggplot(data=S1T5vsS2singlets_allDEG, aes(x=avg_log2FC, y=-log10(p_val_adj), col=diffexpressed, label=delabel)) +
  geom_point() + 
  geom_point(data = filter(S1T5vsS2singlets_allDEG, diffexpressed == "Downregulated Dasatinib Targets"), alpha = 1) +
  geom_point(data = filter(S1T5vsS2singlets_allDEG, diffexpressed == "Not significant or targets of Dasantinib"), alpha =0.001) +
  theme_minimal() +
  geom_text_repel(aes(color = diffexpressed),max.overlaps = 25) +   
  scale_color_manual(values = c("red", "lightgrey"))
plot(S1T5vsS2singlets)
svglite(filename = paste0(outputFolder, 'UMI3_Imputed/S1T5vsS2singlets_VolcanoPlot_6X6.svg'), width = 6, height = 6)
plot(S1T5vsS2singlets)
dev.off()
