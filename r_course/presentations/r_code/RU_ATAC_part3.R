params <-
list(isSlides = "no")

## ----include=FALSE------------------------------------------------------------
suppressPackageStartupMessages(require(knitr))
knitr::opts_chunk$set(echo = TRUE, tidy = T)


## ---- results='asis',include=TRUE,echo=FALSE----------------------------------
if(params$isSlides != "yes"){
  cat("# ATACseq (part 2)

---
"    
  )
  
}



## ----setup, include=FALSE-----------------------------------------------------
library(motifmatchr)
library(MotifDb)
library(JASPAR2020)



## ---- results='asis',include=TRUE,echo=FALSE----------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# Motifs in ATACseq

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("# Post-alignment processing

---
"    
  )
  
}


## ----eval=TRUE----------------------------------------------------------------
library(MotifDb)
MotifDb


## ----eval=TRUE----------------------------------------------------------------
class(MotifDb)


## ----eval=TRUE----------------------------------------------------------------
length(MotifDb)
MotifNames <- names(MotifDb)
MotifNames[1:10]


## ----eval=TRUE----------------------------------------------------------------
MotifDb[1]


## ----eval=TRUE----------------------------------------------------------------
MotifDb[[1]]
colSums(MotifDb[[1]])


## ----eval=TRUE----------------------------------------------------------------
values(MotifDb)[1:2,]


## ----eval=TRUE----------------------------------------------------------------
CTCFMotifs <- query(MotifDb,"CTCF")
CTCFMotifs


## ----eval=TRUE----------------------------------------------------------------
CTCFMotifs <- query(MotifDb,c("CTCF","hsapiens","jaspar2018"))
CTCFMotifs


## ----eval=TRUE----------------------------------------------------------------
library(JASPAR2020)
JASPAR2020


## ----eval=TRUE----------------------------------------------------------------
library(TFBSTools)


## ----eval=FALSE,echo=TRUE-----------------------------------------------------
## ?getMatrixByID


## ----eval=TRUE,echo=TRUE------------------------------------------------------
GATA2mat <- getMatrixByName(JASPAR2020,"GATA2")


## ----eval=TRUE,echo=TRUE------------------------------------------------------
GATA2mat <- getMatrixByID(JASPAR2020,"MA0036.2")


## ----eval=TRUE,echo=TRUE------------------------------------------------------
ID(GATA2mat)


## ----eval=TRUE----------------------------------------------------------------
myMatrix <- Matrix(GATA2mat)
myMatrixToo <- as.matrix(myMatrix)
myMatrix


## ----eval=TRUE,echo=FALSE,warning=FALSE,message=FALSE-------------------------
require(seqLogo)
CTCFMotifs <- query(MotifDb,"MYC")
seqLogo::seqLogo(CTCFMotifs[[1]],ic.scale = FALSE)


## ----eval=TRUE----------------------------------------------------------------
library(seqLogo)
seqLogo::seqLogo(CTCFMotifs[[1]],ic.scale = FALSE)


## ----eval=TRUE----------------------------------------------------------------
library(seqLogo)
seqLogo::seqLogo(CTCFMotifs[[1]])


## ----eval=TRUE----------------------------------------------------------------
myMatrix


## ----eval=TRUE,echo=TRUE------------------------------------------------------
ppm <- myMatrix/colSums(myMatrix)
ppm


## ----eval=TRUE,echo=TRUE------------------------------------------------------
seqLogo::seqLogo(ppm)


## ----eval=TRUE,echo=TRUE------------------------------------------------------
GATA2_IC <- toICM(GATA2mat)
TFBSTools::seqLogo(GATA2_IC)


## ----eval=TRUE,echo=TRUE------------------------------------------------------
library(ggseqlogo)
library(ggplot2)
ggseqlogo(myMatrix)+theme_minimal()


## ---- results='asis',include=TRUE,echo=FALSE----------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# Identifying Motifs in ATAC

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("# Identifying Motifs in ATAC

---
"    
  )
  
}


## ----eval=FALSE,echo=FALSE----------------------------------------------------
## require(tidyverse)
## require(DESeq2)
## load("/Volumes/TomBackup/Temp/rui_ATAC_IL12combined_hg38/salmon/Antibody_ATAC_fromMACSisBlacklisted/ATAC/dds.RData")
## ddsFilt <- dds[,colData(dds)$Group %in% c("ATAC_Control","ATAC_Tbet_mm")]
## colData(ddsFilt) <- colData(dds)[colData(dds)$Group %in% c("ATAC_Control","ATAC_Tbet_mm"),]
## colData(ddsFilt)$Group <- droplevels(colData(ddsFilt)$Group)
## dds <- DESeq(ddsFilt)
## myRes <- results(dds,contrast = c("Group","ATAC_Control","ATAC_Tbet_mm"))
## myRes <- myRes[order(myRes$padj),]
## upRegions <- rownames(myRes)[myRes$log2FoldChange > 0][1:500] %>% gsub("ID:","",.) %>% gsub(":","-",.) %>% sub("-",":",.) %>% GRanges
## # upRegions <- rownames(myRes)[myRes$log2FoldChange > 0 & myRes$padj < 0.05 & !is.na(myRes$padj)] %>% gsub("ID:","",.) %>% gsub(":","-",.) %>% sub("-",":",.) %>% GRanges
## downRegions <- rownames(myRes)[myRes$log2FoldChange < 0][1:500] %>% gsub("ID:","",.) %>% gsub(":","-",.) %>% sub("-",":",.) %>% GRanges
## 
## # DownRegions <- rownames(myRes)[myRes$log2FoldChange < 0 & myRes$padj < 0.05 & !is.na(myRes$padj)] %>% gsub("ID:","",.) %>% gsub(":","-",.) %>% sub("-",":",.) %>% GRanges
## 
## 
## library(BSgenome.Hsapiens.UCSC.hg38)
## upStrings <- getSeq(BSgenome.Hsapiens.UCSC.hg38,resize(upRegions,fix = "center",width = 100))
## downStrings <- getSeq(BSgenome.Hsapiens.UCSC.hg38,resize(downRegions,fix = "center",width = 100))
## names(upStrings) <- as.character(resize(upRegions,fix = "center",width = 100))
## names(downStrings) <- as.character(resize(downRegions,fix = "center",width = 100))
## writeXStringSet(upStrings,file="UpRegions.fa")
## writeXStringSet(downStrings,file="DownStrings.fa")
## 
## 
## require(motifmatchr)
## require(JASPAR2018)
## 
## 
## db <- file.path(system.file("extdata", package="JASPAR2018"),
##                 "JASPAR2018.sqlite")
## opts <- list()
## opts[["tax_group"]] <- "vertebrates"
## opts[["collection"]] <- "CORE"
## opts[["all_versions"]] <- FALSE
## require(TFBSTools)
## motifs <- getMatrixSet(db,opts)
## 
## 
## 
## motif_ixUp <- matchMotifs(motifs, upStrings)
## motif_ixDown <- matchMotifs(motifs, downStrings)
## totalsUp <- colSums(assay(motif_ixUp))
## totalsDow <- colSums(assay(motif_ixDown))
## total <- data.frame(totalsUp,totalsDow)
## total[total == 0] <- 1
## total$Diff <- total$totalsUp/total$totalsDow
## total[order(total$Diff,decreasing = TRUE),]
## 
## IDtoGRanges <- function(IDs,asIGV=FALSE){
##   if(!asIGV){
##     IDmat <- data.frame(matrix(unlist(strsplit(IDs,":")),ncol = 4,byrow = T)[,-1])
##     newSe <- GRanges(IDmat$X1,IRanges(as.numeric(as.vector(IDmat$X2)),as.numeric(as.vector(IDmat$X3))))
##   }
##   newSe
## }
## 
## load("data/myCounts.RData")
## # newDDS <- SummarizedExperiment(counts(dds,normalized=TRUE),
## #                                rowData=NULL,
## #                                rowRanges=IDtoGRanges(rownames(counts(dds))),
## #                                colData= colData(dds) %>%
## #                                  as.data.frame %>%
## #                                  cbind(data.frame(depth=colSums(counts(dds,normalized=FALSE)))) %>%
## #                                  dplyr::select(-sizeFactor)
## # )
## 
## newDDS <- myCounts
## Group <- factor(c("HindBrain","HindBrain","Kidney","Kidney",
##                   "Liver","Liver"))
## colData(newDDS) <- DataFrame(data.frame(Group,row.names=colnames(myCounts)))
## require(DESeq2)
## require(tidyverse)
## atacDDS <- DESeqDataSetFromMatrix(assay(myCounts),
##                                   colData(newDDS),
##                                   ~Group,
##                                   rowRanges=rowRanges(myCounts))
## atacDDS <- DESeq(atacDDS)
## myRes <- results(atacDDS,contrast = c("Group","Liver","Kidney"),format = "GRanges")
## myRes <- myRes[order(myRes$padj),]
## upRegions <- myRes[myRes$log2FoldChange > 0][1:1000]
## downRegions <- myRes[myRes$log2FoldChange < 0,][1:1000]
## 
## library(BSgenome.Mmusculus.UCSC.mm10)
## upStrings <- getSeq(BSgenome.Mmusculus.UCSC.mm10,resize(upRegions,fix = "center",width = 100))
## downStrings <- getSeq(BSgenome.Mmusculus.UCSC.mm10,resize(downRegions,fix = "center",width = 100))
## names(upStrings) <- as.character(resize(upRegions,fix = "center",width = 100))
## names(downStrings) <- as.character(resize(downRegions,fix = "center",width = 100))
## writeXStringSet(upStrings,file="UpRegions.fa")
## writeXStringSet(downStrings,file="DownStrings.fa")
## 
## # names(sxs) <- gsub("\\s.*","",names(sxs))
## assayNames(newDDS) <- "counts"
## require(chromVAR)
## newDDS <- newDDS[rowSums(assay(newDDS)) > 5,]
## require(BSgenome.Mmusculus.UCSC.mm10)
## newDDS <- addGCBias(newDDS,
##                     genome = BSgenome.Mmusculus.UCSC.mm10)
## 
## motif_ix <- matchMotifs(motifs, newDDS,
##                         genome = BSgenome.Mmusculus.UCSC.mm10)
## dev_Known <- computeDeviations(object = newDDS, annotations = motif_ix)
## 
## variability_Known <- computeVariability(dev_Known) %>% arrange(p_value)
## 
## plotVariability(variability_Known, use_plotly = FALSE)
## datatable(variability_Known)
## 
## temp <- deviationScores(dev_Known)
## colnames(temp) <- paste0("Var_",colnames(temp))
## mapOfIDs <- data.frame(id=names(rowData(dev_Known)$name),name=unname(rowData(dev_Known)$name))
## sevew <- merge(mapOfIDs,temp,by.x=1,by.y=0)
## sevew <- merge(variability_Known,sevew,by.x=1,by.y=2)
## sigZ <- sevew %>% arrange(p_value) %>%
##   head(n=20) %>%
##   mutate(newName=name) %>%
##   dplyr::select(newName,starts_with("Var_")) %>%
##   tibble::column_to_rownames("newName")
## pheatmap(sigZ)

