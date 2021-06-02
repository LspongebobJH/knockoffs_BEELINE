unloadNamespace("rlookc")
library("rlookc")
args <- commandArgs(trailingOnly = T)
inFile <- args[1]
outFile <-  args[2]

# input expression data
inputExpr <- read.table(inFile, sep=",", header = 1, row.names = 1)
inputExpr = matrix(rnorm(1e3), nrow = 10, ncol = 100)
rownames(inputExpr) = LETTERS[1:10]
geneNames <- rownames(inputExpr)
rownames(inputExpr) <- c(geneNames)

knockoffResults = rlookc::generateLooks(t(inputExpr), mu = 0, Sigma = cor(t(inputExpr)), 
                                        statistic = knockoff::stat.lasso_lambdasmax, 
                                        output_type = "statistics")
DF = list()
for(i in seq(nrow(inputExpr))){
  DF[[i]] = data.frame(
    Gene1 = geneNames[ i],
    Gene2 = geneNames[-i],
    knockoff_stat = knockoffResults[[i]], 
    q_value = rlookc::knockoff.qvals(knockoffResults[[i]])
  )
}
DF = data.table::rbindlist(DF)
# Write output to a file
# https://stackoverflow.com/questions/38664241/ranking-and-counting-matrix-elements-in-r
outDF <- DF[order(DF$q_value, decreasing=FALSE), ]
write.table(outDF, outFile, sep = "\t", quote = FALSE, row.names = FALSE)
