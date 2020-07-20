params <-
list(toMessage = FALSE)

## ---- setup-------------------------------------------------------------------------
toMessage <- ifelse(!is.null(params$toMessage),params$toMessage,FALSE)


## ----message=toMessage--------------------------------------------------------------
x <- c(1,2,3,4,5)
x


## ----message=toMessage--------------------------------------------------------------
x <- 1:100
x


## ----message=toMessage--------------------------------------------------------------
x <- seq(0,20,5)
x


## ----message=toMessage--------------------------------------------------------------
x <- c(rep(1,2),rep(2,2),rep(3,2))
x


## ----message=toMessage--------------------------------------------------------------
x <- 1:10
x


## ----message=toMessage--------------------------------------------------------------
y <- x[-c(1,length(x))]
y


## ----message=toMessage--------------------------------------------------------------
y <- x[-c(1,5)]
y


## ----message=toMessage--------------------------------------------------------------
y <- sqrt(x[6:7])
y


## ----message=toMessage--------------------------------------------------------------
y <- x[seq(1,10,2)]
y


## ----message=toMessage--------------------------------------------------------------
# Depend on current working directory..
length(dir())


## ----message=toMessage--------------------------------------------------------------
x <- dir()
x[1]


## ----message=toMessage--------------------------------------------------------------
geneNames <- c("Gene_1", "Gene_2", "Gene_3","Gene_4")
expression <- c(1000, 3000, 10000, 12000)
geneLengths <- c(100, 3000, 200, 1000)
names(expression) <- geneNames
names(geneLengths) <- geneNames
expression
geneLengths


## ----message=toMessage--------------------------------------------------------------
names(geneLengths[geneLengths == max(geneLengths)])
names(geneLengths[which.max(geneLengths)])


## ----message=toMessage--------------------------------------------------------------
geneNames[geneLengths > 100 & expression < 10000]



## ----message=toMessage--------------------------------------------------------------
lne <- expression/geneLengths
lne


## ----message=toMessage--------------------------------------------------------------
geneNames[lne > mean(lne)]


