## ----message=T--------------------------------------------------------------------------
geneNames <- c("Gene_1", "Gene_2", "Gene_3","Gene_4")
expression <- c(1000, 3000, 10000, 12000)
geneLengths <- c(100, 3000, 200, 1000)
geneMatrix <- cbind(expression,geneLengths)
rownames(geneMatrix) <- geneNames
geneMatrix


## ----message=T--------------------------------------------------------------------------
lne <- geneMatrix[,"expression"]/geneMatrix[,"geneLengths"]
geneMatrix <- cbind(geneMatrix,lne)
geneMatrix


## ----message=T--------------------------------------------------------------------------
smallGeneMatrix <- geneMatrix[geneMatrix[,"geneLengths"] > 200,]
smallGeneMatrix


## ----message=T--------------------------------------------------------------------------
smallGeneMatrix <- geneMatrix[geneMatrix[,"geneLengths"] > 200 & geneMatrix[,"expression"] > 300,c("expression","lne")]
smallGeneMatrix


## ----message=T--------------------------------------------------------------------------
expressionSum <- sum(geneMatrix[geneMatrix[,"geneLengths"] > 100,"expression"])
geneLengthSum <- sum(geneMatrix[geneMatrix[,"geneLengths"] > 100,"geneLengths"])
#OR
expressionAndGeneLengthSum <- colSums(geneMatrix[geneMatrix[,"geneLengths"] > 100,c("expression","geneLengths")])
expressionAndGeneLengthSum

