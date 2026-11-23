---
title: "Untitled"
author: "Maurice C. Slaughter"
date: "`r Sys.Date()`"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

#################################################################################
#### Cancer data set ############################################################
#################################################################################

#######################
#### Load packages ####
#######################
# BiocManager::install("TCGAbiolinks")
# BiocManager::install("EDASeq") # required for sample normalization
library(TCGAbiolinks); library(EDASeq)
#### Vignette to TCGAbiolinks ####
browseVignettes("TCGAbiolinks")

###############################
#### Available target data ####
###############################
projects <- getGDCprojects()
projects[,c("id", "project_id", "released", "tumor")]

#####################################################
#### Identify data availability for each project ####
#####################################################
id <- projects$id
smpls <- list()
for(i in 1:length(id)){
  temp <- NULL
  query_Target <- NULL
  tryCatch({
    query_Target <- GDCquery(project = id[i], 
                             data.category = "Transcriptome Profiling", 
                             data.type = "Gene Expression Quantification",
                             workflow.type = "STAR - Counts")
  }, error= function(e){cat("ERROR", id[i], ":" ,conditionMessage(e), "\n")})
  if(!is.null(query_Target)){
  samplesDown_Target <- getResults(query_Target)#, cols = c("cases"))
  temp[[1]] <- table(samplesDown_Target$sample_type)
  names(temp) <- id[i]
  } else {
    temp[[1]] <- NA
    names(temp) <- id[i]
  }
  smpls <- c(smpls, temp)
}
smpls
posIDs <- c("CPTAC-3", "REBC-THYR", "TCGA-COAD", "TCGA-BRCA")
projects[projects$id %in% posIDs, ]
smpls[names(smpls) %in% posIDs[3:4]]

##########################################
#### Download and prepare target case ####
##########################################
TCGAbiolinks:::getProjectSummary("TCGA-BRCA")
#### Downloading and prepare TARGET CASE ####
TargetSamples <- GDCquery(project = "TCGA-BRCA", 
                         data.category = "Transcriptome Profiling", 
                         data.type = "Gene Expression Quantification",
                         workflow.type = "STAR - Counts")
#### obtain case information ####
CaseInfo <- getResults(TargetSamples)
head(CaseInfo)
#### subset samples so that there is an equal number of cancer and control samples ####
dataPrimary_Target <- TCGAquery_SampleTypes(barcode = CaseInfo$cases, typesample = "TP") # primary tumor
dataNormal_Target <- TCGAquery_SampleTypes(barcode = CaseInfo$cases, typesample = "NT") # normal tissue
dataPrimary_Target <- dataPrimary_Target[1:113]
dataNormal_Target <- dataNormal_Target[1:113]
#### downloaded samples of interest ####
TargetSamples <- GDCquery(project = "TCGA-BRCA",
                             data.category = "Transcriptome Profiling",
                             data.type = "Gene Expression Quantification",
                             workflow.type = "STAR - Counts",
                             barcode = c(dataPrimary_Target, dataNormal_Target))
GDCdownload(TargetSamples)

#############################################
#### Create summarized experiment object ####
#############################################
data <- GDCprepare(TargetSamples)
assays(data)
colData(data)
rowData(data)
table(rowData(data)$gene_type)

############################################
#### Subset summarizedExperiment object ####
############################################
SECoding <- data[rowData(data)$gene_type == "protein_coding", ]

#############################################################################################################
#### The following function will return the data from specified slots in the summarizedExperiment object ####
#############################################################################################################
dataPrep_Coding <- TCGAanalyze_Preprocessing(object = SECoding, cor.cut = 0.6,  datatype = "fpkm_unstrand")                      

########################################
#### quantile normalize the samples ####
########################################
dataNorm_Coding <- TCGAanalyze_Normalization(tabDF = dataPrep_Coding, geneInfo = geneInfoHT, method = "geneLength")
boxplot(dataPrep_Coding, outline = FALSE)
boxplot(dataNorm_Coding, outline = FALSE)

########################
#### Filter samples ####
########################
dataFilt_Coding <- TCGAanalyze_Filtering(tabDF = dataPrep_Coding, method = "quantile", qnt.cut =  0.25)   

##################################################
#### perform differential expression analysis ####
##################################################
DEGsCoding <- TCGAanalyze_DEA(mat1 = dataFilt_Coding[,dataNormal_Target],
                            mat2 = dataFilt_Coding[,dataPrimary_Target],
                            pipeline="limma",
                            Cond1type = "Normal",
                            Cond2type = "Tumor",
                            fdr.cut = 0.01 ,
                            logFC.cut = 1,
                            method = "glmLRT", ClinicalDF = data.frame())
head(DEGsCoding)
