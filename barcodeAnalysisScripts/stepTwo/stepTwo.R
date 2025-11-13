#======================================================================================================================================
#Input Files: It takes inputs from files generated in stepOne and from 10XCellranger pipeline generated filteredMatrix -> barcode.tsv.gz
#Change these file PATHs based on your folder structure and where your datasets are stored. 
#$PATH needs to be changed/input at 3 places.
#======================================================================================================================================

#***************************************************************************************************************************************
#*****************************************************EDITED SPECIFICALLY FOR VITOs SAMPLES*****************************************************
#***************************************************************************************************************************************
input1Directory <- '/Volumes/YG_SCOPE_01/Temporary_VitoFiles/Analysis/pool1_filtered_feature_bc_matrix_GRCh38/'
input2Directory <- '/Volumes/YG_SCOPE_01/Temporary_VitoFiles/Analysis/stepOne/'
sampleFolders = c('L1/', 'L2/', 'L3/', 'L4/', 'L5/')
outputDirectory <- '/Volumes/YG_SCOPE_01/Temporary_VitoFiles/Analysis/stepTwo/'

#***************************************************************************************************************************************
#*****************************************************DO NOT EDIT BEYOND THIS POINT*****************************************************
#***************************************************************************************************************************************

library(stringdist)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(tidyverse)


numberSamples = length(sampleFolders);
data2file = list()

data1file = as_tibble(read.table(paste0(input1Directory,"barcodes.tsv.gz"), stringsAsFactors=F)) %>% dplyr::rename(cellID = V1) 
data1file = as_tibble(substring(data1file$cellID, 1,nchar(data1file[1,1])-2)) %>% dplyr::rename(cellID = value) 

for(i in 1:numberSamples){
  stepTwoCellIDUMIBarcodes = as_tibble(read.table(paste0(input2Directory, sampleFolders[i], 'uniqueShavedReads.txt'), stringsAsFactors=F)) %>% dplyr::rename(cellID = V1,UMI = V2, BC = V3) %>%
    mutate(BC50 = substring(BC,1,50),
           BC40 = substring(BC,1,40),
           BC30a = substring(BC,1,30),
           BC30b = substring(BC,1,30),
           sampleNum = 1)
  
  if(is.null(dim(data2file))){
    data2file = stepTwoCellIDUMIBarcodes
  } else {
    data2file = bind_rows(data2file, stepTwoCellIDUMIBarcodes)
  }
  
}

#data2file = as_tibble(read.table(paste0(input2Directory,"uniqueShavedReads.txt"), stringsAsFactors=F)) %>% dplyr::rename(cellID = V1,UMI = V2, BC = V3) %>%
  # mutate(BC50 = substring(BC,1,50),
  #        BC40 = substring(BC,1,40),
  #        BC30a = substring(BC,1,30),
  #        BC30b = substring(BC,1,30))
cellIDUMIBarcodes = inner_join(data1file, data2file, by = "cellID")
Barcodes = unique(cellIDUMIBarcodes$BC)
cellIDs = unique(cellIDUMIBarcodes$cellID)

set.seed(2059)
sampleSize = 2000
subsample1 = sample(Barcodes,sampleSize)
subsample2 = sample(Barcodes,sampleSize)
subsample3 = sample(Barcodes,sampleSize)
BarcodesLv1 = as.integer(stringdistmatrix(subsample1, method = "lv"))
BarcodesLv2 = as.integer(stringdistmatrix(subsample2, method = "lv"))
BarcodesLv3 = as.integer(stringdistmatrix(subsample3, method = "lv"))
lBarcodesLv = length(BarcodesLv1)

BarcodesLv = tibble(
  lvdist = c(BarcodesLv1, BarcodesLv2, BarcodesLv3),
  subsamNum = c(rep("subsamping1", lBarcodesLv), rep("subsamping2", lBarcodesLv), rep("subsamping3", lBarcodesLv)))

BarcodesLvHist <- BarcodesLv %>% group_by(subsamNum, lvdist) %>% summarise(length(lvdist)) %>%
  group_by(subsamNum) %>% mutate(totalNum = sum(`length(lvdist)`), fracLvDist = `length(lvdist)`/totalNum)

BarcodesLvHistPlot <- ggplot(BarcodesLvHist, aes(lvdist, fracLvDist)) +
  geom_bar(width = 0.5, stat = 'identity') +
  facet_wrap(facets = vars(subsamNum)) +
  theme_classic()

set.seed(2059)
cellIDsLv = tibble(lvdist = as.integer(stringdistmatrix(cellIDs, method = "lv")))
cellIDsHist <- cellIDsLv  %>% group_by(lvdist)%>% summarise(length(lvdist)) %>% mutate(totalNum = sum(`length(lvdist)`), fracLvDist = `length(lvdist)`/totalNum)

cellIDsHistPlot <- ggplot(cellIDsHist, aes(lvdist, fracLvDist)) +
  geom_bar(width = 0.5, stat = 'identity') +
  theme_classic()

#writing files
ggsave(BarcodesLvHistPlot,file=paste0(outputDirectory,'stepTwoBarcodesLvBeforeStarcode.pdf'))
ggsave(cellIDsHistPlot,file=paste0(outputDirectory,'stepTwoCellIdsLvBeforeStarcode.pdf'))
write.table(cellIDUMIBarcodes, file= paste0(outputDirectory,'stepTwoCellIDUMIBarcodes.txt'),row.names=F,col.names=T,quote=F,sep="\t")
write.table(cellIDUMIBarcodes[,4], file= paste0(outputDirectory,'stepTwoBarcodes50.txt'),row.names=F,col.names=F,quote=F,sep="\t")
write.table(cellIDUMIBarcodes[,5], file= paste0(outputDirectory,'stepTwoBarcodes40.txt'),row.names=F,col.names=F,quote=F,sep="\t")
write.table(cellIDUMIBarcodes[,6], file= paste0(outputDirectory,'stepTwoBarcodes30.txt'),row.names=F,col.names=F,quote=F,sep="\t")



