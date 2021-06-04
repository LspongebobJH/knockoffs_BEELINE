unloadNamespace("rlookc")
library("rlookc")
args <- commandArgs(trailingOnly = T)
inFile <- args[1]
outFile <-  args[2]

# input expression data
inputExpr <- read.table(inFile, sep=",", header = 1, row.names = 1)
geneNames <- rownames(inputExpr)
rownames(inputExpr) <- c(geneNames)
for(i in seq(nrow(inputExpr))){
  inputExpr[i,] = inputExpr[i,] - mean(inputExpr[i,])
  inputExpr[i,] = inputExpr[i,] / (1e-8 + sd(inputExpr[i,]))
}
knockoffResults = rlookc::generateLooks(
  t(inputExpr), 
  mu = 0,
  Sigma = diag(rep(1, nrow(inputExpr))), #  0.95*cor(t(inputExpr)) + 0.05*
  statistic = rlookc::marginal_screen,
  output_type = "statistics"
)

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
