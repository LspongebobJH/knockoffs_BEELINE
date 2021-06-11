args <- commandArgs(trailingOnly = T)
inFile <- args[1]
outFile <-  args[2]

# Input expression data
inputExpr <- read.table(inFile, sep=",", header = 1, row.names = 1)
inputExpr = as.matrix(inputExpr)
geneNames <- rownames(inputExpr)

rownames(inputExpr) <- c(geneNames)
# Optional smoothing
# neighbors = FNN::get.knn(t(inputExpr), k = 20)
# inputExpr = sapply(seq(nrow(neighbors$nn.dist)),
#                    function(i) rowMeans(inputExpr[,c(i, neighbors$nn.index[i,])]) )

# Optional transformation to gaussian marginals
# div_by_max = function(x) x/max(x)
# for(i in seq(nrow(inputExpr))){
#   inputExpr[i,] = rank( inputExpr[i,] , ties.method = "random" ) %>% div_by_max %>% qnorm
# }

# Input must be standardized
for(i in seq(nrow(inputExpr))){
  inputExpr[i,] = inputExpr[i,] - mean(inputExpr[i,])
  inputExpr[i,] = inputExpr[i,] / (1e-8 + sd(inputExpr[i,]))
}

# Optional calibration check
knockoffs = rlookc::computeGaussianKnockoffs(X = t(inputExpr), mu = 0, Sigma = cor(t(inputExpr)), num_realizations = 100)
calibration_results = rlookc::simulateY(
  X = t(inputExpr), 
  knockoffs = knockoffs,
  plot_savepath = paste0(outFile, "_calibration.pdf")
)
saveRDS(calibration_results, paste0(outFile, "_calibration.Rda"))

# Core functionality: GRN inference via LOOKs
knockoffResults = rlookc::generateLooks(
  t(inputExpr), 
  mu = 0,
  Sigma = cor(t(inputExpr)), 
  statistic = knockoff::stat.lasso_lambdasmax,
  output_type = "statistics"
)

DF = list()
for(i in seq(nrow(inputExpr))){
  DF[[i]] = data.frame(
    Gene1 = geneNames[ i],
    Gene2 = geneNames[-i],
    knockoff_stat = knockoffResults[[i]], 
    q_value = rlookc::knockoffQvals(knockoffResults[[i]], offset = 0)
  )
}
DF = data.table::rbindlist(DF)
# print( 
#   ggplot(DF %>% subset(Gene1=="g15")) + 
#     geom_point(aes(x = knockoff_stat, y = q_value, color = Gene1))
# )
# Write output to a file
outDF <- DF[order(DF$q_value, decreasing=FALSE), ]
write.table(outDF, outFile, sep = "\t", quote = FALSE, row.names = FALSE)
warnings()