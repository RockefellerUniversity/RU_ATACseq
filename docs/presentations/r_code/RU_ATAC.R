params <-
list(isSlides = "no")

## ----include=FALSE------------------------------------------------------------
require(ShortRead)
library(ggplot2)
suppressPackageStartupMessages(require(knitr))
knitr::opts_chunk$set(echo = TRUE, tidy = T)


## ---- results='asis',include=TRUE,echo=FALSE----------------------------------
if(params$isSlides != "yes"){
  cat("# ATACseq (part 1)

---
"    
  )
  
}



## ---- results='asis',include=TRUE,echo=FALSE----------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# Alignment for ATACseq

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("# Alignment for ATACseq

---
"    
  )
  
}



## ----processData_BuildIndex, echo=TRUE,eval=FALSE,cache=FALSE,tidy=FALSE------
## library(BSgenome.Hsapiens.UCSC.hg19)
## mainChromosomes <- paste0("chr",c(1:21,"X","Y","M"))
## mainChrSeq <- lapply(mainChromosomes,
##                      function(x)BSgenome.Hsapiens.UCSC.hg19[[x]])
## names(mainChrSeq) <- mainChromosomes
## mainChrSeqSet <- DNAStringSet(mainChrSeq)
## writeXStringSet(mainChrSeqSet,
##                 "BSgenome.Hsapiens.UCSC.hg19.mainChrs.fa")


## ----index, echo=TRUE,eval=FALSE,tidy=FALSE-----------------------------------
## library(Rsubread)
## buildindex("BSgenome.Hsapiens.UCSC.hg19.mainChrs",
##            "BSgenome.Hsapiens.UCSC.hg19.mainChrs.fa",
##            indexSplit = TRUE,
##            memory = 1000)


## -----------------------------------------------------------------------------
read1 <- "ATAC_Data/ATAC_FQs/SRR891269_1.fastq.gz"
read2 <- "ATAC_Data/ATAC_FQs/SRR891269_2.fastq.gz"


## ----eval=FALSE,include=FALSE-------------------------------------------------
## require(ShortRead)
## read1 <- readFastq("~/Downloads/ENCFF175VOD.fastq.gz")
## read2 <- readFastq("~/Downloads/ENCFF447BGX.fastq.gz")
## writeFastq(read1[1:1000,],"ATACSample_r1.fastq.gz")
## writeFastq(read2[1:1000,],"ATACSample_r2.fastq.gz")
## # id(read2[1:1000,])
## # myRes <- bamQC("~/Downloads/Sorted_ATAC_50K_2.bam")


## ----eval=TRUE----------------------------------------------------------------
require(ShortRead)
read1 <- readFastq("data/ATACSample_r1.fastq.gz")
read2 <- readFastq("data/ATACSample_r2.fastq.gz")
id(read1)[1:2]
id(read2)[1:2]


## ----processData_acdc, include=FALSE, eval=F----------------------------------
## setwd("~/Projects/Results/chipseq/testRunforTalk/")


## ----processData_align, echo=TRUE,eval=FALSE,cache=FALSE,tidy=FALSE-----------
## 
## align("BSgenome.Hsapiens.UCSC.hg19.mainChrs",
##       readfile1=read1,readfile2=read2,
##       output_file = "ATAC_50K_2.bam",
##       nthreads=2,type=1,
##       unique=TRUE,maxFragLength = 2000)
## 


## ----indedsx, echo=TRUE,eval=FALSE,tidy=FALSE---------------------------------
## library(Rbowtie2)
## bowtie2_build(references="BSgenome.Hsapiens.UCSC.hg19.mainChrs.fa",
##               bt2Index="BSgenome.Hsapiens.UCSC.hg19.mainChrs_bowtie2")


## ----indeaax, echo=TRUE, eval=FALSE-------------------------------------------
## gunzip("ATAC_Data/ATAC_FQs/SRR891269_1.fastq.gz")
## gunzip("ATAC_Data/ATAC_FQs/SRR891269_2.fastq.gz")


## ----indexas, echo=TRUE,eval=FALSE,tidy=FALSE---------------------------------
## library(Rsamtools)
## bowtie2(bt2Index = "BSgenome.Hsapiens.UCSC.hg19.mainChrs_bowtie2",
##           samOutput = "ATAC_50K_2_bowtie2.sam",
##           seq1 = "ATAC_Data/ATAC_FQs/SRR891269_1.fastq",
##           seq1 = "ATAC_Data/ATAC_FQs/SRR891269_2.fastq"
##         )
## asBam("ATAC_50K_2_bowtie2.sam")


## ----processData_indexAndSort, echo=TRUE,eval=FALSE,cache=FALSE,tidy=FALSE----
## library(Rsamtools)
## sortedBAM <- file.path(dirname(outBAM),
##                        paste0("Sorted_",basename(outBAM))
##                        )
## 
## sortBam(outBAM,gsub("\\.bam","",basename(sortedBAM)))
## indexBam(sortedBAM)


## ---- eval=F, echo=F----------------------------------------------------------
## mappedReads <- idxstatsBam(sortedBAM)
## save(mappedReads,file="data/pres1_idxstats.RData")

## ---- eval=T, echo=F,cache=TRUE,dependson="processData_setBAM"----------------
load("data/pres1_idxstats.RData")


## ----quickMappingStatsPerChromosomea, echo=TRUE,eval=F------------------------
## library(Rsamtools)
## mappedReads <- idxstatsBam(sortedBAM)
## 


## ----quickMappingStatsPerChromosomes, echo=TRUE,eval=TRUE,cache=TRUE,dependson="processData_setBAM",fig.height=4,fig.width=6----
library(ggplot2)

ggplot(mappedReads,aes(seqnames,mapped,fill=seqnames))+
  geom_bar(stat="identity")+coord_flip()


## ---- results='asis',include=TRUE,echo=FALSE----------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# Post-alignment processing

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


## ----processData_readingInDatad, echo=TRUE,eval=TRUE,cache=TRUE,dependson="processData_setBAM"----
library(GenomicAlignments)
flags=scanBamFlag(isProperPair = TRUE)



## ----processData_readingInDatas, echo=TRUE,eval=TRUE,cache=TRUE,dependson="processData_setBAM"----
myParam=ScanBamParam(flag=flags,
                   what=c("qname","mapq","isize"),
                   which=GRanges("chr20", IRanges(1,63025520)))
myParam



## ---- eval=F, echo=F----------------------------------------------------------
## atacReads <- readGAlignmentPairs("~/Downloads/ATAC_Workshop/Sorted_ATAC_50K_2.bam",
##                                  param=myParam)
## save(atacReads, file="data/pres1_atacreads.RData")


## ----processData_readingInDataa, echo=TRUE,eval=F-----------------------------
## atacReads <- readGAlignmentPairs(sortedBAM,
##                                  param=myParam)
## class(atacReads)
## 


## ---- eval=T, echo=F----------------------------------------------------------

load("data/pres1_atacreads.RData")
class(atacReads)


## ----processData_readingInData, echo=TRUE,eval=TRUE,cache=TRUE,dependson="processData_setBAM"----
atacReads[1:2,]


## ----processData_readingInData2, echo=TRUE,eval=TRUE,cache=TRUE,dependson="processData_setBAM"----
read1 <- first(atacReads)
read2 <- second(atacReads)
read2[1,]


## ----processData_readingInData3, echo=TRUE,eval=TRUE,cache=TRUE,dependson="processData_readingInData2"----
read1MapQ <- mcols(read1)$mapq
read2MapQ <- mcols(read2)$mapq
read1MapQ[1:2]


## ----processData_readingInData4, echo=TRUE,eval=TRUE,cache=TRUE,dependson="processData_readingInData3"----
read1MapQFreqs <- table(read1MapQ)
read2MapQFreqs <- table(read2MapQ)
read1MapQFreqs
read2MapQFreqs


## ----processData_readingInData5,fig.width=9,fig.height=4,  echo=TRUE,eval=TRUE,cache=TRUE,dependson="processData_readingInData4"----
library(ggplot2)
toPlot <- data.frame(MapQ=c(names(read1MapQFreqs),names(read2MapQFreqs)),
           Frequency=c(read1MapQFreqs,read2MapQFreqs),
           Read=c(rep("Read1",length(read1MapQFreqs)),rep("Read2",length(read2MapQFreqs))))
toPlot$MapQ <- factor(toPlot$MapQ,levels = unique(sort(as.numeric(toPlot$MapQ))))
ggplot(toPlot,aes(x=MapQ,y=Frequency,fill=MapQ))+geom_bar(stat="identity")+facet_grid(~Read)


## ----processData_extractingReadsss1, echo=TRUE,eval=FALSE,cache=TRUE,dependson="processData_readingInData"----
## atacReads_read1 <- first(atacReads)
## insertSizes <- abs(elementMetadata(atacReads_read1)$isize)
## head(insertSizes)


## ----processData_extractingRead1, echo=FALSE,eval=TRUE,cache=TRUE,dependson="processData_readingInData"----
atacReads_read1 <- GenomicAlignments::first(atacReads)
insertSizes <- abs(elementMetadata(atacReads_read1)$isize)
head(insertSizes)


## ----processData_plottingFrffagmentLengths, echo=TRUE,eval=TRUE,cache=TRUE,dependson="processData_extractingRead1"----

fragLenSizes <- table(insertSizes)
fragLenSizes[1:5]



## ----processData_plottingFrdagmentLengths, echo=TRUE,eval=TRUE,cache=TRUE,dependson="processData_extractingRead1"----
library(ggplot2)
toPlot <- data.frame(InsertSize=as.numeric(names(fragLenSizes)),
                            Count=as.numeric(fragLenSizes))
fragLenPlot <- ggplot(toPlot,aes(x=InsertSize,y=Count))+geom_line()



## ----processData_plottfingFragmentLengths2, echo=TRUE,eval=TRUE,cache=TRUE,dependson="processData_extractingRead1"----

fragLenPlot+theme_bw()



## ----processData_plottingFragmentLengths3, echo=TRUE,eval=TRUE,cache=TRUE,dependson="processData_extractingRead1"----

fragLenPlot + scale_y_continuous(trans='log2')+theme_bw()


## ----processData_plottingFragmentLengths24, echo=TRUE,eval=TRUE,cache=TRUE,dependson="processData_plottingFragmentLengths",fig.width=6,fig.height=4----
fragLenPlot+ scale_y_continuous(trans='log2')+
  geom_vline(xintercept = c(180,247),colour="red")+
  geom_vline(xintercept = c(315,437),colour="darkblue")+
  geom_vline(xintercept = c(100),colour="darkgreen")+theme_bw()



## ----processData_createOpenRegionBAM, echo=TRUE,eval=TRUE,cache=TRUE,dependson=c("processData_extractingRead1","processData_readingInData")----
atacReads_NucFree <- atacReads[insertSizes < 100,]
atacReads_MonoNuc <- atacReads[insertSizes > 180 &
                                 insertSizes < 240,]
atacReads_diNuc <- atacReads[insertSizes > 315 &
                               insertSizes < 437,]


## ----processData_createOpenRegionBAM_2, echo=TRUE,eval=FALSE,cache=TRUE,dependson="processData_createOpenRegionBAM"----
## nucFreeRegionBam <- gsub("\\.bam","_nucFreeRegions\\.bam",sortedBAM)
## monoNucBam <- gsub("\\.bam","_monoNuc\\.bam",sortedBAM)
## diNucBam <- gsub("\\.bam","_diNuc\\.bam",sortedBAM)
## 
## library(rtracklayer)
## export(atacReads_NucFree,nucFreeRegionBam,format = "bam")
## export(atacReads_MonoNuc,monoNucBam,format = "bam")
## export(atacReads_diNuc,diNucBam,format = "bam")


## ----processData_createOpenRegionBigWig_2, echo=TRUE,eval=TRUE,cache=TRUE,dependson="processData_createOpenRegionBAM"----
atacReads[1,]
atacFragments<- granges(atacReads)
atacFragments[1,]


## ----processData_createOpenRegionBigWig_5, echo=TRUE,eval=TRUE,cache=TRUE,dependson="processData_createOpenRegionBigWig_2"----
duplicatedFragments <- sum(duplicated(atacFragments))
totalFragments <- length(atacFragments)
duplicateRate <- duplicatedFragments/totalFragments
nonRedundantFraction <-  1-duplicateRate
nonRedundantFraction


## ----processData_createOpenRegionBigWig_21, echo=TRUE,eval=FALSE,cache=TRUE,dependson="processData_createOpenRegionBigWig_2"----
## openRegionRPMBigWig <- gsub("\\.bam","_openRegionRPM\\.bw",sortedBAM)
## myCoverage <- coverage(atacFragments,
##                        weight = (10^6/length(atacFragments)))
## export.bw(myCoverage,openRegionRPMBigWig)


## ---- results='asis',include=TRUE,echo=FALSE----------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# ATACseqQC

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("# ATACseqQC

---
"    
  )
  
}


## ----eval=FALSE---------------------------------------------------------------
## BiocManager::install("ATACseqQC")


## ----eval=FALSE---------------------------------------------------------------
## library(ATACseqQC)
## ATACQC <- bamQC("~/Downloads/Sorted_ATAC_50K_2_ch17.bam")


## ----eval=TRUE,echo=FALSE,include=FALSE---------------------------------------
load("data/ATACQC.RData")


## -----------------------------------------------------------------------------
names(ATACQC)


## -----------------------------------------------------------------------------
ATACQC$PCRbottleneckCoefficient_1


## -----------------------------------------------------------------------------
ATACQC$PCRbottleneckCoefficient_2

