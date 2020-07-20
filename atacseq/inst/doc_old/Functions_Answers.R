## ----message=T--------------------------------------------------------------------------
findSmallestFactorial <- function(x){
  if(!is.numeric(x)){
    message("Please provide a numeric argument!")
  }else{
    factorialAnswer <- 0
    count <- 0
    while(factorialAnswer <= x){
      count <- count+1
      if(count == 1){
        factorialAnswer <- 1
      }else{
        factorialAnswer <- factorialAnswer * count 
      }
    }
    return(count)
  }
}

findSmallestFactorial(3000)
findSmallestFactorial(10^100)


## ----message=T--------------------------------------------------------------------------
findPerfectSquares <- function(x){
  vecOfSquares <- c()
  evenSquares <- 0
  for(i in x){
    if(sqrt(i) %% 1 == 0){
      vecOfSquares <- c(vecOfSquares,i)
      if(sqrt(i) %% 1 == 0 & i %% 2 == 0){
        evenSquares <- evenSquares+1
      }
    }
  }
  return(list(vectorOfSquares=vecOfSquares,nOfEvenSquares=evenSquares))
}

findPerfectSquares(1:100)

#Or

altFindPerfectSquares <- function(x){
  perfectSquare <- sqrt(x) %% 1 == 0
  vecOfSquares <- x[perfectSquare]
  evenSquares <- table(perfectSquare & x %% 2 == 0)[T]
  return(list(vectorOfSquares=vecOfSquares,nOfEvenSquares=evenSquares))
}

altFindPerfectSquares(1:100)

#Res <- findPerfectSquares(1:10^8)
#altRes <- altFindPerfectSquares(1:10^8)


## ----message=T--------------------------------------------------------------------------
summariseResults <- function(dirName=getwd(),annotation){
  if(missing(annotation)){
    message("Annotation file must be provided")
    break()
  }
  filesToRead <- dir(dirName,pattern = "*\\.txt",full.names=T)
  if(length(filesToRead) == 0){
    message("No expression files found in ",dirName)
    break()
  }
  
  fileRead <- vector("list",length=length(filesToRead))
  for(i in 1:length(filesToRead)){
    fileRead[[i]] <- read.delim(filesToRead[i],header=F,sep="\t")
    colnames(fileRead[[i]]) <- c("GeneNames",basename(filesToRead[i]))
  }
  mergedTable <- NULL
  for(i in fileRead){
    if(is.null(mergedTable)){
      mergedTable <- i
    }else{
      mergedTable <- merge(mergedTable,i,by=1,all=T)
    }
  }
   if(!file.exists(annotation)){
    message("No annotation file found:",annotation)
    break()     
   }
   Annotation <- read.table(annotation,sep="\t",h=T)  
  annotatedExpression <- merge(Annotation,mergedTable,by=1,all.x=F,all.y=T)
  
  outName <- paste0(basename((normalizePath(dirName))),"_Summarised.csv")

  write.table(annotatedExpression,file=outName,sep=",",col.names=T,row.names=F,quote=F)
  
  return(outName)
}

summariseResults("ExpressionResults/","ExpressionResults/Annotation.ann")


## ----message=T--------------------------------------------------------------------------
allSamples <- c("sample1.txt","sample2.txt","sample3.txt","sample4.txt","sample5.txt","sample5.txt")
testSamples <- c("sample1.txt","sample5.txt")
match(testSamples,allSamples)
allSamples[match(testSamples,allSamples)]


## ----message=T--------------------------------------------------------------------------

allSamples <- c("sample1.txt","sample2.txt","sample3.txt","sample4.txt","sample5.txt")
testSamples <- c("sample1.txt","sample5.txt")
allSamples %in% testSamples
allSamples[allSamples %in% testSamples]


## ----message=T--------------------------------------------------------------------------
summariseResults2 <- function(dirName=getwd(),annotation,sampleGroups=NULL){
  if(missing(annotation)){
    message("Annotation file must be provided")
    break()
  }
  filesToRead <- dir(dirName,pattern = "*\\.txt",full.names=T)
  if(length(filesToRead) == 0){
    message("No expression files found in ",dirName)
    break()
  }
  
  fileRead <- vector("list",length=length(filesToRead))
  for(i in 1:length(filesToRead)){
    fileRead[[i]] <- read.delim(filesToRead[i],header=F,sep="\t")
    colnames(fileRead[[i]]) <- c("GeneNames",basename(filesToRead[i]))
  }
  mergedTable <- NULL
  for(i in fileRead){
    if(is.null(mergedTable)){
      mergedTable <- i
    }else{
      mergedTable <- merge(mergedTable,i,by=1,all=T)
    }
  }
  
  if(!file.exists(annotation)){
    message("No annotation file found:",annotation)
    break()     
  }  
  Annotation <- read.table(annotation,sep="\t",h=T)  
  annotatedExpression <- merge(Annotation,mergedTable,by=1,all.x=F,all.y=T)
  
  outAnnotatedExpression <- paste0(basename((normalizePath(dirName))),"_Summarised.csv")

  write.table(annotatedExpression,file=outAnnotatedExpression,sep=",",col.names=T,row.names=F,quote=F)
  
  
  ### New section
  
  if(!is.null(sampleGroups)){
      groupAsamples <- unique(sampleGroups[[1]])
      groupBsamples <- unique(sampleGroups[[2]])
      groupA <- colnames(annotatedExpression) %in% groupAsamples
      groupB <- colnames(annotatedExpression) %in% groupBsamples

      indexGroupOne <- groupA
      indexGroupTwo <- groupB

      ttestResults <- apply(annotatedExpression,1,function(x)
        t.test(as.numeric(x[indexGroupOne]),as.numeric(x[indexGroupTwo]))
        )
      
      testResult <- sapply(ttestResults,function(x)
        c(log2(x$estimate[2]) - log2(x$estimate[1]), x$statistic,x$p.value)
        )
      testResult <- t(testResult)
      colnames(testResult) <- c("logFC","tStatistic","pValue")
      annotatedResult <- cbind(annotatedExpression[,1:3],testResult)
      annotatedResult <- annotatedResult[order(annotatedResult$tStatistic),]
      outAnnotatedTestResults <- paste0(basename((normalizePath(dirName))),"_ttestResults.csv")
      write.table(annotatedResult,file=outAnnotatedTestResults,sep=",",col.names=T,row.names=F,quote=F)    
  }
  
  
  
}
expressionFiles <- dir("ExpressionResults/",pattern="*.txt")
groupA <- expressionFiles[grep("[1-5].txt",expressionFiles)]
groupB <- expressionFiles[grep("[6-9,0].txt",expressionFiles)]
myGroups <- list(groupA,groupB)
summariseResults2("ExpressionResults/","ExpressionResults/Annotation.ann",myGroups)


