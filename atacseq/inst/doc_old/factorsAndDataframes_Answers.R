## ----message=T--------------------------------------------------------------------------
ordinalFactor <- factor(c("high", "low", "mid", "low", "mid", "low", "mid", "high", "mid", "high"),ordered=T,levels=c("low", "mid", "high"))
ordinalFactor


## ----message=T--------------------------------------------------------------------------
filteredFactor <- ordinalFactor[ordinalFactor > "low"]
filteredFactor


## ----message=T--------------------------------------------------------------------------
newFactor <- factor(ordinalFactor,ordered=T,levels=c("low", "mid", "high","veryHigh"))
newFactor[length(newFactor)] <- "veryHigh"
newFactor[newFactor > "mid"]


## ----message=T--------------------------------------------------------------------------
Annotation <- data.frame(geneNames=c("Gene_1", "Gene_2", "Gene_3","Gene_4","Gene_5"), ensembl=c("Ens001", "Ens003", "Ens006", "Ens007", "Ens010"),pathway=c("Glycolysis", "TGFb", "Glycolysis", "TGFb", "Glycolysis"),geneLengths=c(100, 3000, 200, 1000,1200))


## ----message=T--------------------------------------------------------------------------
Sample1 <- data.frame(ensembl=c("Ens001", "Ens003", "Ens006","Ens010"),expression=c(1000, 3000, 10000,5000))


## ----message=T--------------------------------------------------------------------------
Sample2 <- data.frame(ensembl=c("Ens001", "Ens003", "Ens006","Ens007","Ens010"),expression=c(1500, 1500, 17000,500,10000))


## ----message=T--------------------------------------------------------------------------
AnnoSample1 <- merge(Annotation,Sample1,by.x=2,by.y=1,all=F) 
AnnoSample1And2 <- merge(AnnoSample1,Sample2,by=1,all=F) 
AnnoSample1And2


## ----message=T--------------------------------------------------------------------------
AnnoSample1And2$Sample1_lne <- AnnoSample1And2$expression.x/AnnoSample1And2$geneLengths
AnnoSample1And2$Sample2_lne <- AnnoSample1And2$expression.y/AnnoSample1And2$geneLengths
AnnoSample1And2


## ----message=T--------------------------------------------------------------------------
rownames(AnnoSample1And2) <- AnnoSample1And2$ensembl
mean(c(AnnoSample1And2["Ens006","Sample1_lne"],AnnoSample1And2["Ens006","Sample2_lne"]))


## ----message=T--------------------------------------------------------------------------
log2FoldChange <- log2(AnnoSample1And2$Sample2_lne) - log2(AnnoSample1And2$Sample1_lne)
names(log2FoldChange) <- AnnoSample1And2$geneNames
log2FoldChange


## ----message=T--------------------------------------------------------------------------
sum(AnnoSample1And2[AnnoSample1And2$pathway == "Glycolysis","geneLengths"])

