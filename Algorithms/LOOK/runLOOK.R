suppressPackageStartupMessages(library (magrittr, warn.conflicts = FALSE, quietly = TRUE))
suppressPackageStartupMessages(library (optparse, warn.conflicts = FALSE, quietly = TRUE))
option_list <- list (
  make_option(c("-e","--expressionFile"), type = 'character',
              help= "Path to comma separated file containing gene-by-cell matrix with
              cell names as the first row and gene names as 
              the first column. Required."),
  make_option(c("-o","--outFile"), , type = 'character',
              help= "outFile name to write the output ranked edges. Required."),
  make_option(c("-c","--calibrate"), action = 'store_true', default = FALSE,
              type = 'character',
              help= "")
)
parser <- OptionParser(option_list = option_list)
arguments <- parse_args(parser, positional_arguments = FALSE)


# Input expression data
inputExpr <- read.table(arguments$expressionFile, sep=",", header = 1, row.names = 1)
inputPT = 
  arguments$expressionFile %>% 
  dirname %>% 
  dirname %>%
  file.path("PseudoTime.csv") %>% 
  read.table(sep = ",", header = 1, row.names = 1)
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
if(arguments$calibrate){
  knockoffs = rlookc::computeGaussianKnockoffs(X = t(inputExpr),
                                               mu = 0,
                                               Sigma = cor(t(inputExpr)), 
                                               num_realizations = 100)
  calibration_results = rlookc::simulateY(
    X = t(inputExpr), 
    knockoffs = knockoffs,
    plot_savepath = paste0(outFile, "_calibration.pdf"), 
    active_set_size = 2, 
    FUN = function(x) all(x>0) + rbinom(n = 1, size = 1, prob = 0.5)
  )
  saveRDS(calibration_results, paste0(outFile, "_calibration.Rda"))
}

# Core functionality: GRN inference via LOOKs
arguments$method = "time_series"
if( arguments$method == "look" ){
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
} else if(arguments$method == "time_series" ){
  # We just need one set of knockoffs
  knockoffs = knockoff::create.gaussian(
    t(inputExpr), 
    mu = 0,
    Sigma = cor(t(inputExpr))
  )
  # For multi-branching datasets, do each branch separately
  # and combine between calculation of symmetric stats and q-values.
  knockoffResults = list()
  for(branch in seq(ncol(inputPT))){
    # Select and order the right cells
    # train on the next closest cell ahead of you, regardless of how far away it is...
    cells_to_use = data.frame( pt = inputPT[,branch] )
    cells_to_use[["original_position"]] = seq(nrow(cells_to_use))
    cells_to_use = cells_to_use[!is.na(cells_to_use[["pt"]]), ]
    cells_to_use = cells_to_use[order( cells_to_use[["pt"]]), ]
    cells_to_use = cells_to_use[["original_position"]]
    ncell = length(cells_to_use)
    # Do each gene separately
    for(k in seq_along(geneNames)){
      if(length(knockoffResults) < k){
        knockoffResults[[k]] = rep(0, length(geneNames))
      }
      X = t(inputExpr)[cells_to_use[1:(ncell-1)], ]
      X_k = knockoffs [cells_to_use[1:(ncell-1)], ]
      y = t(inputExpr)[cells_to_use[2: ncell   ],k]
      knockoffResults[[k]] = knockoffResults[[k]] +
        knockoff::stat.glmnet_lambdasmax(X, X_k, y)
    }  
  }
  # Assemble results
  DF = list()
  for(k in seq(nrow(inputExpr))){
    w = knockoffResults[[k]][-k] # Disallow autoregulation
    DF[[k]] = data.frame(
      Gene1 = geneNames[ k],
      Gene2 = geneNames[-k],
      knockoff_stat = w, 
      q_value = rlookc::knockoffQvals(w, offset = 0)
    )
  }
  DF = data.table::rbindlist(DF)
}


# Write output to a file
outDF <- DF[order(DF$q_value, decreasing=FALSE), ]
write.table(outDF, arguments$outFile, sep = "\t", quote = FALSE, row.names = FALSE)
warnings()