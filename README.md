This is the repository for code/data for the barcode-based single-cell analysis Rebecca et al.

The barcodeAnalysisScripts folder contains the steps that are needed to be performed to extract Fatemap barcodes associated with each cell from the side-reaction reads. More details about each step can be found in the methods section and in the FateMap paper (Goyal, et al, 2023 - [https://docs.google.com/document/d/1QVRW4Qn7IVNRu7fkkKeFV\\\_iPpvR7TQY98VShZUxzH1U/edit?tab=t.0](https://docs.google.com/document/d/1QVRW4Qn7IVNRu7fkkKeFV\_iPpvR7TQY98VShZUxzH1U/edit?tab=t.0){.uri})

In the singleCellAnalysis folder, there are scripts to process the single-cell RNA sequencing output from Cell Ranger pipeline and merge that with the barcode information. Then, compare clone size and differential gene expression between sample treated with Dasatinib + BRAF/MEKi vs sample with just BRAF/MEKi.

-   The ALRA_imputation.R is a script to impute the single-cell RNA sequencing matrix with ALRA and then merge with the barcodes associated with each cell ID.
-   The clonalAnalysis.R takes the output rds file of ALRA_imputation.R as input to make the plots in the paper. However, scAnalysisWithoutImputation.R can be used alternatively to ALRA_imputation.R and that Seurat object can also be used as input to clonalAnalysis.R (if interested in doing the analysis without any imputation.)

The data used to create the plots can also be found in the plotData folder and the plots are in plots folder.

The plots can be recreated from the plotData using scripts in plotScripts.
