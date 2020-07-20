## ----message=T--------------------------------------------------------------------------

for(x in 1:10){
  if(x == 1){
    factorialAnswer <- 1
  }else{
    factorialAnswer <- factorialAnswer * x 
  }
}
factorialAnswer


## ----message=T--------------------------------------------------------------------------

factorialAnswer <- 0
count <- 0

while(factorialAnswer <= 1000){
  count <- count+1
  if(count == 1){
    factorialAnswer <- 1
  }else{
    factorialAnswer <- factorialAnswer * count 
  }
}
count


## ----message=T--------------------------------------------------------------------------
condExercise <- 1:40
condExercise
vectorResult <- ifelse(condExercise<10,"small",ifelse(condExercise < 31,"mid","big"))
temp <- factor(vectorResult,levels=c("small","mid","big"),order=T)
temp


## ----message=T--------------------------------------------------------------------------

filesToRead <- dir("ExpressionResults/",pattern = "*\\.txt",full.names=T)
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
  
  print(nrow(mergedTable))
}

mergedTable[1:3,] 



## ----message=T--------------------------------------------------------------------------
Annotation <- read.table("ExpressionResults/Annotation.ann",sep="\t",h=T)
annotatedExpression <- merge(Annotation,mergedTable,by=1,all.x=F,all.y=T)
annotatedExpression[1:2,]

summary(annotatedExpression$Pathway)
summary(Annotation$Pathway)

