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
library(soGGi)
library(ChIPQC)
library(GenomicAlignments)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(DESeq2)
library(BSgenome.Mmusculus.UCSC.mm10)
library(tracktables)
library(clusterProfiler)
library(org.Mm.eg.db)

library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)



## ---- results='asis',include=TRUE,echo=FALSE----------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# Evaluating TSS signal

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("# Evaluating TSS signal

---
"    
  )
  
}


## ----processData_aligna, echo=TRUE,eval=FALSE,cache=FALSE---------------------
## BiocManager::install("soGGi")
## library(soGGi)


## ----processData_txdb, echo=TRUE,eval=TRUE,cache=FALSE------------------------
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
TxDb.Hsapiens.UCSC.hg19.knownGene


## ----processData_genes, echo=TRUE,eval=TRUE,cache=FALSE-----------------------
genesLocations <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
genesLocations


## ----processData_resize, echo=TRUE,eval=TRUE,cache=FALSE----------------------
tssLocations <- resize(genesLocations,fix="start",width = 1)
tssLocations



## ----processData_subset, echo=TRUE,eval=F-------------------------------------
## mainChromosomes <- paste0("chr",c(1:21,"X","Y","M"))
## 
## myindex<-(seqnames(tssLocations) %in% mainChromosomes)
## 
## tssLocations <- tssLocations[as.numeric(myindex)]
## 
## seqlevels(tssLocations) <- mainChromosomes
## 


## ----processData_subset_fix, echo=F,eval=F------------------------------------
## mainChromosomes <- paste0("chr",c(1:21,"X","Y","M"))
## 
## myindex<-(match(seqnames(tssLocations), mainChromosomes))
## 
## 
## tssLocations <- tssLocations[!is.na(myindex)]
## 
## seqlevels(tssLocations) <- mainChromosomes
## save(tssLocations, file="data/tss_pres2.Rmd")
## 


## ---- echo=F, eval=T----------------------------------------------------------
load("data/tss_pres2.Rmd")


## ----processData_soggi, echo=TRUE,eval=FALSE,cache=FALSE----------------------
## library(soGGi)
## sortedBAM <- "~/Downloads/ATAC_Workshop/Sorted_ATAC_50K_2.bam"
## 
## library(Rsamtools)
## # Nucleosome free
## allSignal <- regionPlot(bamFile = sortedBAM,
##                         testRanges = tssLocations)


## ----processData_soggia1, echo=F,eval=T---------------------------------------
load("data/nucFree_TSS.Rdata")


## ----processData_soggia, echo=TRUE,eval=FALSE,cache=FALSE---------------------
## nucFree <- regionPlot(bamFile = sortedBAM,
##                         testRanges = tssLocations,
##                         style = "point",
##                         format="bam",
##                         paired=TRUE,
##                         minFragmentLength = 0,
##                         maxFragmentLength = 100,
##                         forceFragment = 50)
## class(nucFree)


## ----processData_plot, echo=TRUE,eval=T,cache=TRUE,message=FALSE,warning=FALSE----
plotRegion(nucFree)


## ----processData_soggi4, echo=F,eval=FALSE,cache=FALSE,message=FALSE,warning=FALSE----
## monoNuc <- regionPlot(bamFile = sortedBAM,
##                         testRanges = tssLocations,
##                         style = "point",
##                         format="bam",
##                         paired=TRUE,
##                         minFragmentLength = 180,maxFragmentLength = 240,forceFragment = 80)
## save(monoNuc,file = "data/monoNuc_TSS.RData")


## ----processData_soggi2.5, echo=F,eval=T,cache=FALSE,message=FALSE,warning=FALSE----

load(file = "data/monoNuc_TSS.RData")


## ----processData_soggi3, echo=TRUE,eval=FALSE,cache=FALSE,message=FALSE,warning=FALSE----
## monoNuc <- regionPlot(bamFile = sortedBAM,
##                         testRanges = tssLocations,
##                         style = "point",
##                         format="bam",
##                         paired=TRUE,
##                         minFragmentLength = 180,maxFragmentLength = 240,forceFragment = 80)
## 


## ----processData_plot3, echo=TRUE,eval=T,cache=FALSE--------------------------
plotRegion(monoNuc)


## ---- results='asis',include=TRUE,echo=FALSE----------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# Peak calling

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("# Peak calling

---
"    
  )
  
}


## MACS2 callpeak -t singleEnd.bam --nomodel --shift -100

##                 --extsize 200 --format BAM -g MyGenome


## MACS2 callpeak -t singleEnd.bam --nomodel --shift 37

##                --extsize 73 --format BAM -g MyGenome


## MACS2 callpeak -t pairedEnd.bam -f BAMPE

##                --outdir path/to/output/

##                --name pairedEndPeakName -g MyGenome


## MACS2 callpeak  -t ~/Downloads/Sorted_ATAC_50K_2_openRegions.bam

##                 --outdir ATAC_Data/ATAC_Peaks/ATAC_50K_2

##                 --name Sorted_ATAC_50K_2_Small_Paired_peaks.narrowPeak

##                 -f BAMPE -g hs


## ----processData_callQC,messages=FALSE,warning=FALSE, echo=TRUE,eval=FALSE,cache=TRUE,dependson="readinPeakCalling"----
## library(ChIPQC)
## library(rtracklayer)
## library(DT)
## library(dplyr)
## library(tidyr)
## 
## blkList <- import.bed("~/Downloads/ENCFF001TDO.bed.gz")
## openRegionPeaks <- "~/Downloads/ATAC_Workshop/ATAC_Data/ATAC_Peaks/Sorted_ATAC_50K_2_Small_Paired_peaks.narrowPeak"
## 
## qcRes <- ChIPQCsample("~/Downloads/ATAC_Workshop/ATAC_Data/ATAC_BAM/Sorted_ATAC_50K_2_openRegions.bam",
##                       peaks = openRegionPeaks,
##                       annotation ="hg19",
##                       chromosomes = "chr20",
##                       blacklist = blkList,
##                       verboseT = FALSE)


## ----ded, include=FALSE,cache=TRUE--------------------------------------------
library(ChIPQC)
library(rtracklayer)
load("data/qcRes.RData")

blkList <- import.bed("data/ENCFF001TDO.bed.gz")
openRegionPeaks <- "data/Sorted_ATAC_50K_2_Small_Paired_peaks.narrowPeak"



## ----processData_callQC2,messages=FALSE,warning=FALSE, echo=TRUE,eval=TRUE,cache=TRUE,dependson="ded"----
myMetrics <- QCmetrics(qcRes)
myMetrics[c("RiBL%","RiP%")]
flgCounts <- flagtagcounts(qcRes)
DupRate <- flgCounts["DuplicateByChIPQC"]/flgCounts["Mapped"]
DupRate*100


## ----processData_filterBLKlist, echo=F,eval=F,cache=TRUE,dependson="processData_callQC"----
## 
## MacsCalls <- granges(qcRes[seqnames(qcRes) %in% "chr20"])
## 
## data.frame(Blacklisted=sum(MacsCalls %over% blkList),
##            Not_Blacklisted=sum(!MacsCalls %over% blkList))
## save(MacsCalls,file="data/MacsCall_pres2.RData")
## 


## ---- echo=T,eval=F-----------------------------------------------------------
## 
## MacsCalls <- granges(qcRes[seqnames(qcRes) %in% "chr20"])
## 
## data.frame(Blacklisted=sum(MacsCalls %over% blkList),
##            Not_Blacklisted=sum(!MacsCalls %over% blkList))
## 
## MacsCalls <- MacsCalls[!MacsCalls %over% blkList]
## 


## ---- echo=F, eval=T----------------------------------------------------------
load("data/MacsCall_pres2.RData")
MacsCalls <- MacsCalls[!MacsCalls %over% blkList]



## ---- results='asis',include=TRUE,echo=FALSE----------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# Annotating Peaks

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("# Annotating Peaks

---
"    
  )
  
}


## ----processData_annotatePeak, echo=TRUE,eval=TRUE,cache=TRUE,dependson="processData_filterBLKlist"----
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
MacsCalls_Anno <-  annotatePeak(MacsCalls,
                                TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene)
MacsCalls_Anno


## ----processData_Pie, echo=TRUE,eval=TRUE,cache=TRUE,dependson="processData_annotatePeak"----
plotAnnoPie(MacsCalls_Anno)


## ----processData_annotated, echo=TRUE,eval=TRUE,cache=TRUE,dependson="processData_annotatePeak"----
MacsGR_Anno <- as.GRanges(MacsCalls_Anno)
MacsGR_TSS <-   MacsGR_Anno[abs(MacsGR_Anno$distanceToTSS) < 500]
MacsGR_TSS[1,]


## ---- results='asis',include=TRUE,echo=FALSE----------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# Functional analysis of peaks

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("# Functional analysis of peaks

---
"    
  )
  
}


## ----processData_funAnalysise,echo=TRUE,eval=TRUE,cache=TRUE,dependson="processData_filterBLKlist"----
library(rGREAT)
great_Job <- submitGreatJob(MacsCalls, species = "hg19")
availableCategories(great_Job)


## ----processData_funAnalysis2,echo=TRUE,eval=TRUE,cache=TRUE,dependson="processData_funAnalysis2"----
great_ResultTable = getEnrichmentTables(great_Job, category = "GO")
names(great_ResultTable)          
great_ResultTable[["GO Biological Process"]][1:4, ]



## ---- results='asis',include=TRUE,echo=FALSE----------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# Differential ATACseq

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("# Differential ATACseq

---
"    
  )
  
}


## ----processData_consensus_read, echo=F,eval=T--------------------------------
peaks <- dir("data/ATAC_Peaks_forCounting/",
             pattern="*.narrowPeak",full.names=TRUE)


## ----processData_consensusa, echo=TRUE,eval=F---------------------------------
## peaks <- dir("~/Downloads/ATAC_Workshop/ATAC_Data/ATAC_Peaks_forCounting/",
##              pattern="*.narrowPeak",full.names=TRUE)


## ----processData_consensus_read2, echo=TRUE,eval=T----------------------------
myPeaks <- lapply(peaks,ChIPQC:::GetGRanges,simple=TRUE)
allPeaksSet_nR <- reduce(unlist(GRangesList(myPeaks)))
overlap <- list()
for(i in 1:length(myPeaks)){
  overlap[[i]] <- allPeaksSet_nR %over% myPeaks[[i]]
}
overlapMatrix <- do.call(cbind,overlap)
colnames(overlapMatrix) <- basename(peaks)
mcols(allPeaksSet_nR) <- overlapMatrix



## ----processData_consensusaa, echo=TRUE,eval=TRUE,cache=TRUE,dependson="processData_consensusa"----

allPeaksSet_nR[1:2,]


## ----processData_consensusCounting, echo=F,eval=F,cache=TRUE,dependson="processData_consensusa"----
## blklist <- import.bed("data/ENCFF547MET.bed.gz")
## nrToCount <- allPeaksSet_nR[!allPeaksSet_nR %over% blklist
##                             & !seqnames(allPeaksSet_nR) %in% "chrM"]
## save(nrToCount,file="data/nrToCount_pres2.RData")


## ---- eval=F, echo=T----------------------------------------------------------
## blklist <- import.bed("data/ENCFF547MET.bed.gz")
## nrToCount <- allPeaksSet_nR[!allPeaksSet_nR %over% blklist
##                             & !seqnames(allPeaksSet_nR) %in% "chrM"]
## nrToCount


## ---- eval=T, echo=F----------------------------------------------------------
load("data/nrToCount_pres2.RData")
nrToCount


## ----processData_consensusCountingf, echo=TRUE,eval=TRUE,cache=TRUE,dependson="processData_consensus"----
library(Rsubread)
occurrences <- rowSums(
  as.data.frame(elementMetadata(nrToCount)
                )
  )

nrToCount <- nrToCount[occurrences >= 2,]
nrToCount


## ----processData_consensusCounting2, echo=TRUE,eval=FALSE,cache=TRUE,dependson="processData_consensusCountingf"----
## library(GenomicAlignments)
## bamsToCount <- dir("~/Downloads/ATAC_Workshop/ATAC_Data/ATAC_BAM_forCounting/",
##                    full.names = TRUE,pattern = "*.\\.bam$")
## 
## myCounts <- summarizeOverlaps(consensusToCount,
##                               bamsToCount,singleEnd=FALSE)
## 
## colnames(myCounts) <- c("HindBrain_1","HindBrain_2","Kidney_1","Kidney_2",
##                         "Liver_1","Liver_2")


## ----processData_DEseq2_PCA, echo=TRUE,eval=TRUE,cache=TRUE-------------------
library(DESeq2)
load("data/myCounts.RData")
Group <- factor(c("HindBrain","HindBrain","Kidney","Kidney",
                  "Liver","Liver"))
metaData <- data.frame(Group,row.names=colnames(myCounts))

atacDDS <- DESeqDataSetFromMatrix(assay(myCounts),
                                  metaData,
                                  ~Group,
                                  rowRanges=rowRanges(myCounts))
atacDDS <- DESeq(atacDDS)



## ----processData_DEseq2_Results_ResultsTable, echo=TRUE,eval=TRUE,cache=TRUE,dependson="processData_DEseq2_PCA"----
KidneyMinusHindbrain <- results(atacDDS,
                                c("Group","Kidney","HindBrain"),
                               format="GRanges")
KidneyMinusHindbrain <- KidneyMinusHindbrain[
                                order(KidneyMinusHindbrain$pvalue)
                         ]
KidneyMinusHindbrain


## ----processData_DEseq2_ResultsToTSSregions,message=FALSE,warning=FALSE, echo=TRUE,eval=TRUE,cache=TRUE,dependson="processData_DEseq2_Results_ResultsTable"----
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
toOverLap <- promoters(TxDb.Mmusculus.UCSC.mm10.knownGene,
                       500,500)
KidneyMinusHindbrain <- KidneyMinusHindbrain[
                          (!is.na(KidneyMinusHindbrain$padj) 
                           & KidneyMinusHindbrain$padj < 0.05) 
                           & KidneyMinusHindbrain %over% toOverLap,]
makebedtable(KidneyMinusHindbrain,"KidneyMinusHindbrain.html",getwd())


## ----processData_DEseq2_functionalEnrichmentAnalysiss, echo=TRUE,eval=TRUE,cache=TRUE, dependson="processData_DEseq2_ResultsToTSSregions",message=FALSE,warning=FALSE----

anno_KidneyMinusHindbrain <- annotatePeak(KidneyMinusHindbrain,
                                          TxDb = TxDb.Mmusculus.UCSC.mm10.knownGene)
DB_ATAC <- as.data.frame(anno_KidneyMinusHindbrain)
DB_ATAC[1,]


## ----processData_DEseq2_functionalEnrichmentAnalysisd, echo=TRUE,eval=TRUE,cache=TRUE, dependson="processData_DEseq2_functionalEnrichmentAnalysiss",message=FALSE,warning=FALSE----
library(clusterProfiler)
go <- enrichGO(DB_ATAC$geneId, 
                OrgDb = "org.Mm.eg.db",ont = "BP",maxGSSize = 5000)
go[1:2,1:6]



## ---- results='asis',include=TRUE,echo=FALSE----------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# Finding Motifs in ATACseq data

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("# Finding Motifs in ATACseq data

---
"    
  )
  
}


## ----processData_motifCutsa, echo=TRUE,eval=TRUE,cache=TRUE,message=FALSE,warning=FALSE----
library(MotifDb)
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg19)
CTCF <- query(MotifDb, c("CTCF"))
CTCF


## ----processData_motifCutsb, echo=TRUE,eval=TRUE,cache=TRUE,dependson="processData_motifCutsa"----
names(CTCF)
ctcfMotif <- CTCF[[1]]
ctcfMotif[,1:4]


## ----processData_motifCutsc, echo=TRUE,eval=TRUE,cache=TRUE,dependson="processData_motifCutsb",fig.height=5,fig.width=7----
library(seqLogo)
seqLogo(ctcfMotif)


## ----processData_motifCutsd, echo=TRUE,eval=TRUE,cache=TRUE,dependson="processData_motifCutsc",warning=FALSE,message=FALSE----

myRes <- matchPWM(ctcfMotif,BSgenome.Hsapiens.UCSC.hg19[["chr20"]])
myRes


## ----processData_motifCutse, echo=TRUE,eval=TRUE,cache=TRUE,dependson="processData_motifCutsd"----

toCompare <- GRanges("chr20",ranges(myRes))
toCompare


## ----processData_motifCutsdwdw, echo=TRUE,eval=FALSE--------------------------
## BAM <- "~/Downloads/ATAC_Workshop/ATAC_Data/ATAC_BAM/Sorted_ATAC_50K_2_openRegions.bam"
## atacReads_Open <- readGAlignmentPairs(BAM)
## read1 <- first(atacReads_Open)
## read2 <- second(atacReads_Open)
## read2[1,]


## ---- echo=F, eval=F----------------------------------------------------------
## save(read1, read2, file="data/pres2_reads1and2.RData")
## 
## 


## ----processData_motifCutsfa, echo=FALSE,eval=TRUE----------------------------
load(file="data/pres2_reads1and2.RData")
read2[1,]


## ----processData_motifCutsf, echo=TRUE,eval=TRUE,cache=TRUE,dependson="processData_motifCutsfa"----

Firsts <- resize(granges(read1),fix="start",1)
First_Pos_toCut <- shift(granges(Firsts[strand(read1) == "+"]),
                                         4)
First_Neg_toCut <- shift(granges(Firsts[strand(read1) == "-"]),
                                         -5)

Seconds <- resize(granges(read2),fix="start",1)
Second_Pos_toCut <- shift(granges(Seconds[strand(read2) == "+"]),
                                4)
Second_Neg_toCut <- shift(granges(Seconds[strand(read2) == "-"]),
                                -5)

test_toCut <- c(First_Pos_toCut,First_Neg_toCut,
                Second_Pos_toCut,Second_Neg_toCut)
test_toCut[1:2,]


## ----processData_motifCutsg, echo=TRUE,eval=TRUE,cache=TRUE,dependson="processData_motifCutsf"----

cutsCoverage <- coverage(test_toCut)
cutsCoverage20 <- cutsCoverage["chr20"]
cutsCoverage20[[1]]


## ----processData_motifCutsh, echo=TRUE,eval=TRUE,cache=TRUE,dependson="processData_motifCutsg",message=FALSE,warning=FALSE,fig.height=5,fig.width=7----


CTCF_Cuts_open <- regionPlot(cutsCoverage20,
                         testRanges = toCompare,
                         style = "point",
                         format="rlelist",distanceAround = 500)



## ----processData_motifCutsi, echo=TRUE,eval=TRUE,cache=TRUE,dependson="processData_motifCutsh",message=FALSE,warning=FALSE,fig.height=5,fig.width=7----

plotRegion(CTCF_Cuts_open,outliers = 0.001)+
  ggtitle("NucFree Cuts Centred on CTCF")+theme_bw()


