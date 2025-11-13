library(ggplot2)
library(tidyverse)
library(svglite)
options(future.globals.maxSize = 4000 * 1024^2)

inputFolder <- "RebeccaEtAl_barcodeAnalysis/plotData/cloneSizePlot/"
outputFolder <- "RebeccaEtAl_barcodeAnalysis/plots/cloneSizePlot/"

controlSampleCloneData <- read.csv(paste0(inputFolder, "controlSampleCloneSize.csv"))
treatedSampleCloneData <- read.csv(paste0(inputFolder, "treatedSampleCloneSize.csv"))

plotControl <- ggplot(controlSampleCloneData, aes(x = controlSampleCloneData$ClusterSize,y = after_stat(count/sum(count)))) +
  geom_histogram(binwidth = 1, fill = "hotpink3", color = "black") +
  labs(x = "Clone Size", y = "Frequency", title = "Clone Size Distribution in Control\n(only BRAFi/MEKi )") +
  theme_minimal() +  theme(plot.title = element_text(hjust = 0.5)) + 
  scale_x_continuous(breaks = seq(0, max(controlSampleCloneData$ClusterSize), by = 5)) 

plotTreatment <-ggplot(treatedSampleCloneData, aes(x = treatedSampleCloneData$ClusterSize,y = after_stat(count/sum(count)))) +
  geom_histogram(binwidth = 1, fill = "turquoise3", color = "black") +
  labs(x = "Clone Size", y = "Frequency", title = "Clone Size Distribution in \nBRAFi/MEKi + Dasatinib") +
  theme_minimal() +  theme(plot.title = element_text(hjust = 0.5)) + 
  scale_x_continuous(breaks = seq(0, 30, by = 5))  + xlim(0, max(controlSampleCloneData$ClusterSize))

svglite(filename = paste0(outputFolder, 'controlSampleClonesize.svg'), width = 6, height = 4)
plot(plotControl)
dev.off()

svglite(filename = paste0(outputFolder, 'treatmentSampleClonesize.svg'), width = 6, height = 4)
plot(plotTreatment)
dev.off()