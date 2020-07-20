## ----message=T--------------------------------------------------------------------------
 catAndExpr <- read.delim("data/categoriesAndExpression.txt")


## ----message=T--------------------------------------------------------------------------
 summary(catAndExpr)


## ----message=T--------------------------------------------------------------------------
 summary(catAndExpr)
 summary(catAndExpr[catAndExpr$pathway == "Glycolysis",])
 summary(catAndExpr[catAndExpr$pathway == "TGFb",])


## ----message=T--------------------------------------------------------------------------
contingency <- ftable(catAndExpr[,c("ofInterest","pathway")])
contingency
#fisher.test(contingency)



## ----message=T--------------------------------------------------------------------------
meanExpression <- mean(catAndExpr$Expr)
sdExpression <- sd(catAndExpr$Expr)
Gene13Expression <- catAndExpr[catAndExpr$geneName == "Gene13",]$Expression
1-pnorm(Gene13Expression,meanExpression,sdExpression)


## ----message=T--------------------------------------------------------------------------
glycolysisExpression <- catAndExpr[catAndExpr$pathway == "Glycolysis",]$Expression
tgfbExpression <- catAndExpr[catAndExpr$pathway == "TGFb",]$Expression
var.test(glycolysisExpression,tgfbExpression)
t.test(glycolysisExpression,tgfbExpression,var.equal = FALSE)


## ----message=T--------------------------------------------------------------------------
 lmExercise <- read.delim("data/lmExercise.txt")


## ----message=T--------------------------------------------------------------------------
 lmXY <- lm(x~y,lmExercise)
 lmXZ <- lm(x~z,lmExercise)
lmXY
 plot(x~y,lmExercise)
 abline(lmXY,col="red")
lmXZ
plot(x~z,lmExercise)
 abline(lmXZ,col="red")


## ----message=T--------------------------------------------------------------------------
str(summary(lmXY))
summary(lmXY)$r.squared
summary(lmXZ)$r.squared


## ----message=T--------------------------------------------------------------------------
str(lmXY)
predictedXfromY <- unname(lmXY$coefficients[1] + lmXY$coefficients[2]*100)
predictedXfromZ <- unname(lmXZ$coefficients[1] + lmXZ$coefficients[2]*100)
predictedXfromY
predictedXfromZ

