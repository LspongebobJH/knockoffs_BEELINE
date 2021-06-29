suppressPackageStartupMessages(library (ggplot2, warn.conflicts = FALSE, quietly = TRUE))
suppressPackageStartupMessages(library (magrittr, warn.conflicts = FALSE, quietly = TRUE))
suppressPackageStartupMessages(library (optparse, warn.conflicts = FALSE, quietly = TRUE))
option_list <- list (
  make_option(c("-e","--expressionFile"), type = 'character',
              help= "Path to comma separated file containing gene-by-cell matrix with
              cell names as the first row and gene names as 
              the first column. Required."),
  make_option(c("-o","--outFile"), , type = 'character',
              help= "outFile name to write the output ranked edges. Required."),
  make_option(c("-c","--calibrate"), action = 'store', default = FALSE,
              type = 'logical',
              help= "")
)
parser <- OptionParser(option_list = option_list)
arguments <- parse_args(parser, positional_arguments = FALSE)

nonparametricMarginalScreen = function(X, knockoffs, y){
  sapply(1:ncol(X), function(k) loess(y ~ knockoffs[,k])$s - loess(y ~ X[,k])$s )
}

knockoffEmpiricalCorrection = function(w){
  inflation = max( median(w), 0 )
  w - inflation
}
# For manual inspection of the process
# arguments = list(
#   expressionFile = "~/Desktop/jhu/research/projects/Beeline/inputs/Synthetic_with_protein_and_velocity/dyn-LL/dyn-LL-500-10/LOOK/ExpressionData.csv",
#   calibrate = T
# )

standardize = function(inputExpr){
  for(i in seq(nrow(inputExpr))){
    inputExpr[i,] = inputExpr[i,] - mean(inputExpr[i,])
    inputExpr[i,] = inputExpr[i,] / (1e-8 + sd(inputExpr[i,]))
  }
  inputExpr
}
# Input expression data
inputPT = 
  arguments$expressionFile %>% 
  dirname %>% 
  dirname %>% 
  file.path("PseudoTime.csv") %>% 
  read.table(sep = ",", header = 1, row.names = 1)
inputExpr <- read.table(arguments$expressionFile, sep=",", header = 1, row.names = 1)
inputExpr = inputExpr[,order(inputPT[[1]])]
inputPT   = inputPT[order(inputPT[[1]]),]
stopifnot("Pseudotime and expression don't have the same number of cells.\n"=
            nrow(inputPT)==ncol(inputExpr))
# Separate different types of measurements
inputProtein     = inputExpr[grepl("^p_", rownames(inputExpr)),]
inputRNA         = inputExpr[grepl("^x_|^g", rownames(inputExpr)),]
inputRNAvelocity = inputExpr[grepl("^velocity_x_", rownames(inputExpr)),]
inputRNA = as.matrix(inputRNA) %>% standardize
# Gene name handling:
# Clean
geneNames_more_like_cleanNames = function(x){
  x %>% gsub("^velocity_", "", .) %>% gsub("^(p|x)_", "", .) 
}
geneNames <- rownames(inputRNA) %>% geneNames_more_like_cleanNames 
rownames(inputRNA) <- geneNames
# Sort
geneNames %<>% gtools::mixedsort()
inputRNA = inputRNA[geneNames, ]
# Clean and sort other types of measurements
if(nrow(inputProtein)>0){
  inputProtein = as.matrix(inputProtein) %>% sqrt %>% standardize
  rownames(inputProtein) %<>% geneNames_more_like_cleanNames
  inputProtein = inputProtein[geneNames, ]
} 
if(nrow(inputRNAvelocity)>0){
  inputRNAvelocity = as.matrix(inputRNAvelocity) %>% standardize
  rownames(inputRNAvelocity) %<>% geneNames_more_like_cleanNames
  inputRNAvelocity = inputRNAvelocity[geneNames, ]
} 
rm(inputExpr)

# Optional smoothing
# neighbors = FNN::get.knn(t(inputExpr), k = 20)
# inputExpr = sapply(seq(nrow(neighbors$nn.dist)),
#                    function(i) rowMeans(inputExpr[,c(i, neighbors$nn.index[i,])]) )

# Optional transformation to gaussian marginals
# div_by_max = function(x) x/max(x)
# for(i in seq(nrow(inputExpr))){
#   inputExpr[i,] = rank( inputExpr[i,] , ties.method = "random" ) %>% div_by_max %>% qnorm
# }


runCalibrationCheck = function(X, noiselevel = 1){
  if(arguments$calibrate){
    knockoffs = rlookc::computeGaussianKnockoffs(X = X,
                                                 mu = 0,
                                                 Sigma = cor(X), 
                                                 num_realizations = 100)
    message("Checking calibration with a univariate sigmoid activation function.\n")
    calibration_results = rlookc::simulateY(
      X = X, 
      knockoffs = knockoffs,
      plot_savepath = paste0(arguments$outFile, "_calibration.pdf"), 
      # Univariate sigmoidal
      active_set_size = 1,
      FUN = function(x) 1/(1+exp(-10*(x-1.5)))
      
      # Bivariate bool-ish
      # active_set_size = 2, 
      # FUN = function(x) all(x>0) + rbinom(n = 1, size = noiselevel, prob = 0.5)
    )
    saveRDS(calibration_results, paste0(arguments$outFile, "_calibration.Rda"))
    return(invisible(calibration_results))
  }
}


# Core functionality: GRN inference via knockoff-based tests
# of carefully constructed null hypotheses
arguments$method = "rna_production_protein_predictor" # "steady_state"
{
  if( arguments$method == "steady_state" )
  {
    # Optional calibration check
    runCalibrationCheck(X = t(inputRNA))
    # Use leave-one-out knockoffs
    knockoffResults = rlookc::generateLooks(
      t(inputRNA), 
      mu = 0,
      Sigma = cor(t(inputRNA)), 
      statistic = knockoff::stat.lasso_lambdasmax,
      output_type = "statistics"
    ) %>% lapply(knockoffEmpiricalCorrection)
    
    DF = list()
    for(i in seq(nrow(inputRNA))){
      DF[[i]] = data.frame(
        Gene1 = geneNames[ i],
        Gene2 = geneNames[-i],
        knockoff_stat = knockoffResults[[i]], 
        q_value = rlookc::knockoffQvals(knockoffResults[[i]], offset = 0)
      )
    }
    DF = data.table::rbindlist(DF)
  } 
  else if(arguments$method == "return_to_average" )
  {
    # Test E[X(t+1, k)] indep. X(t, j) given X(t, -j)
    # (Averaging is used for the future state, but not the current state)
    # We just need one set of knockoffs
    knockoffs = knockoff::create.gaussian(
      t(inputExpr), 
      mu = 0,
      Sigma = cor(t(inputExpr))
    )
    # For multi-branching datasets, do each branch separately
    # and combine symmetric stats before getting q-values.
    knockoffResults = list()
    for(branch in seq(ncol(inputPT))){
      # Select and order the right cells
      cells_to_use = data.frame( pt = inputPT[,branch] )
      cells_to_use[["original_position"]] = seq(nrow(cells_to_use))
      cells_to_use = cells_to_use[!is.na(cells_to_use[["pt"]]), ]
      cells_to_use = cells_to_use[order( cells_to_use[["pt"]]), ]
      ncell = nrow(cells_to_use)
      # Estimate E[X(t)]
      num_abcissae = 1e3
      abcissae = seq(min(cells_to_use[["pt"]]), 
                     max(cells_to_use[["pt"]]), 
                     length.out = num_abcissae)
      smooth_one_gene =  function(gene, evaluate_at = abcissae) {
        to_smooth = cells_to_use
        to_smooth[["y"]] = inputExpr[gene,cells_to_use[["original_position"]]] 
        predict(
          loess(y ~ pt, data = to_smooth), 
          newdata = data.frame(pt = evaluate_at)
        )
      }
      # Do each gene separately
      for(k in seq_along(geneNames)){
        if(length(knockoffResults) < k){
          knockoffResults[[k]] = rep(0, length(geneNames))
        }
        # Use observed current state
        X = t(inputExpr)[cells_to_use[["original_position"]],]
        X_k = knockoffs [cells_to_use[["original_position"]],]
        # Use expected future state assuming a rapid return to average
        timestep_length = diff(abcissae)[[1]]
        future_timepoints = cells_to_use[["pt"]] + timestep_length
        Y =  smooth_one_gene(gene = k, evaluate_at = future_timepoints)
        # Don't run off the end of the trajectory and get garbage or missing future estimates
        keep = !is.na(Y)
        X   = X[keep, ]
        X_k = X_k[keep, ]
        Y   = Y[keep]
        # Different branches are combined by adding the symmetric stats
        knockoffResults[[k]] = knockoffResults[[k]] +
          knockoff::stat.glmnet_lambdasmax(X, X_k, Y)
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
  else if(arguments$method == "average_along_trajectory" )
  {
    # For multi-branching datasets, do each branch separately
    # and combine between calculation of symmetric stats and q-values.
    knockoffResults = list()
    for(branch in seq(ncol(inputPT))){
      # Select and order the right cells
      cells_to_use = data.frame( pt = inputPT[,branch] )
      cells_to_use[["original_position"]] = seq(nrow(cells_to_use))
      cells_to_use = cells_to_use[!is.na(cells_to_use[["pt"]]), ]
      cells_to_use = cells_to_use[order( cells_to_use[["pt"]]), ]
      ncell = nrow(cells_to_use)
      # Estimate E[X(t)]
      num_abcissae = 1e3
      abcissae = seq(min(cells_to_use[["pt"]]), 
                     max(cells_to_use[["pt"]]), 
                     length.out = num_abcissae)
      smooth_one_gene =  function(i) {
        to_smooth = cells_to_use
        to_smooth[["y"]] = inputExpr[i,cells_to_use[["original_position"]]] 
        predict(
          loess(y ~ pt, data = to_smooth), 
          newdata = data.frame(pt = abcissae)
        )
      }
      smoothed_expression = sapply( 1:nrow(inputExpr), smooth_one_gene )
      stopifnot(ncol(smoothed_expression)==nrow(inputExpr))
      runCalibrationCheck(X = smoothed_expression)
      
      # Test E[X(t+1, k)] indep. E[X(t, j)] given E[X(t, -j)]
      # (Averaging is used for both the future state and the current state)
      # Make a set of knockoffs
      knockoffs = knockoff::create.gaussian(
        smoothed_expression, 
        mu = 0,
        Sigma = cor(smoothed_expression)
      )
      
      # Do each gene separately
      for(k in seq_along(geneNames)){
        if(length(knockoffResults) < k){
          knockoffResults[[k]] = rep(0, length(geneNames))
        }
        knockoffResults[[k]] = knockoffResults[[k]] +
          knockoff::stat.glmnet_lambdasmax(smoothed_expression[1:(num_abcissae-1), ],
                                           knockoffs          [1:(num_abcissae-1), ], 
                                           smoothed_expression[2:num_abcissae, k])
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
  else if(arguments$method == "next_cell" )
  {
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
  else if(arguments$method == "steady_state_protein" )
  { 
    stopifnot("Protein levels must be provided with prefix 'p_'.\n"=nrow(inputProtein)>0)
    # Input must be standardized
    for(i in seq(nrow(inputProtein))){
      inputProtein[i,] = inputProtein[i,] - mean(inputProtein[i,])
      inputProtein[i,] = inputProtein[i,] / (1e-8 + sd(inputProtein[i,]))
    }
    # Optional calibration check
    # runCalibrationCheck(X = t(inputProtein))
    # Generate knockoffs for the protein levels
    knockoffs = knockoff::create.gaussian(
      t(inputProtein), 
      mu = 0,
      Sigma = cor(t(inputProtein))
    )
    knockoffResults = list()
    # Do each gene separately
    for(k in seq_along(geneNames)){
      X = t(inputProtein)
      X_k = knockoffs 
      y = t(inputRNA)[,k]
      knockoffResults[[k]] = knockoff::stat.glmnet_lambdasmax(X, X_k, y)
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
  else if(arguments$method == "rna_production_protein_predictor" )
  { 
    stopifnot("Protein levels must be provided with prefix 'p_'.          \n"=nrow(inputProtein)>0)
    stopifnot("Velocity levels must be provided with prefix 'velocity_x_'.\n"=nrow(inputRNAvelocity)>0)
    # Generate knockoffs for combined RNA and protein levels
    X = t(inputProtein) 
    knockoffs = knockoff::create.gaussian(
      X, 
      mu = 0,
      Sigma = cor(X)
    )
    
    # Optional calibration check
    x = runCalibrationCheck(X, noiselevel = 1)
    
    # Do each gene separately
    q = w = list()
    for(k in seq_along(geneNames)){
      y = t(inputRNAvelocity)[,k]
      # Infer the decay rate in a robust way (piecewise linear)
      concentration = inputRNA[k,]
      plot(concentration, y, pch = ".")
      concentration_bins = cut(concentration, breaks = 5)
      decay_rate = list()
      for(bin in levels(concentration_bins)){
        idx = concentration_bins==bin
        if(sum(idx)<10){next}
        decay_rate[[bin]] = coef( quantreg::rq(y[idx] ~ concentration[idx] ) )
        clip(min(concentration[idx]), max(concentration[idx]), y1 = -100, y2 = 100)
        abline(decay_rate[[bin]][[1]], decay_rate[[bin]][[2]])
      }
      nona = function(x) x[!is.na(x)]
      negative_only = function(x) x[x<0]
      decay_rate %<>% sapply(extract2, "concentration[idx]") %>% nona %>% negative_only %>% median
      clip(min(concentration), max(concentration[idx]), y1 = -100, y2 = 100)
      abline(a = 0, b = decay_rate, col = "red")
    
      # subtract off decay rate; only production rate remains to be modeled
      y = y - concentration*decay_rate
      # w[[k]] = nonparametricMarginalScreen(X, knockoffs, y)
      w[[k]] = knockoff::stat.glmnet_lambdasmax(X, knockoffs, y)
      
      # For interactive use
      # data.frame(
      #   production = y,
      #   protein_regulator = X[,k-1],
      #   protein_product = X[,k],
      #   protein_regulator_knockoff = knockoffs[,k-1],
      #   protein_product_knockoff = knockoffs[,k]
      # ) %>%
      #   tidyr::pivot_longer(cols = !production) %>%
      #   ggplot() +
      #   geom_point(aes(x = value, y = production, colour = name, shape = name)) +
      #   ggtitle(paste0("Candidate regulators and their knockoffs versus gene", k, " production rate"))
        

    }  
    # Assemble results
    DF = list()
    for(k in seq_along(geneNames)){
      keep = seq_along(geneNames) #[-k] allow autoregulation
      DF[[k]] = data.frame(
        Gene1 = geneNames[ k],
        Gene2 = geneNames[keep],
        knockoff_stat = w[[k]][keep], 
        q_value = rlookc::knockoffQvals(w[[k]][keep], offset = 1)
      )
    }
    DF = data.table::rbindlist(DF)
  }
  else if(arguments$method == "rna_velocity_protein_predictor" )
  { 
    stopifnot("Protein levels must be provided with prefix 'p_'.          \n"=nrow(inputProtein)>0)
    stopifnot("Velocity levels must be provided with prefix 'velocity_x_'.\n"=nrow(inputRNAvelocity)>0)
    # Input must be standardized
    for(i in seq(nrow(inputProtein))){
      inputProtein[i,] = inputProtein[i,] - mean(inputProtein[i,])
      inputProtein[i,] = inputProtein[i,] / (1e-8 + sd(inputProtein[i,]))
    }
    # Optional calibration check
    # runCalibrationCheck(X = t(inputProtein))
    # Generate knockoffs for the protein levels
    knockoffs = knockoff::create.gaussian(
      t(inputProtein), 
      mu = 0,
      Sigma = cor(t(inputProtein))
    )
    knockoffResults = list()
    # Do each gene separately
    for(k in seq_along(geneNames)){
      X = t(inputProtein)
      X_k = knockoffs 
      y = t(inputRNAvelocity)[,k]
      knockoffResults[[k]] = knockoff::stat.glmnet_lambdasmax(X, X_k, y)
    }  
    # Assemble results
    DF = list()
    for(k in seq(nrow(inputRNA))){
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
  else if(arguments$method == "rna_velocity_rna_predictor" )
  { 
    stopifnot("Velocity levels must be provided with prefix 'velocity_x_'.\n"=nrow(inputRNAvelocity)>0)
    # Generate knockoffs for the rna levels
    knockoffs = knockoff::create.gaussian(
      t(inputRNA), 
      mu = 0,
      Sigma = cor(t(inputRNA))
    )
    knockoffResults = list()
    # Do each gene separately
    for(k in seq_along(geneNames)){
      X = t(inputRNA)
      X_k = knockoffs 
      y = t(inputRNAvelocity)[,k]
      knockoffResults[[k]] = knockoff::stat.glmnet_lambdasmax(X, X_k, y)
    }  
    # Assemble results
    DF = list()
    for(k in seq(nrow(inputRNA))){
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
}

# What's the pattern of discoveries?
try({
  p = ggplot(DF) + 
    geom_tile(
      aes(
        x = Gene1 %>% gsub("^g", "", .) %>% as.numeric,
        y = Gene2 %>% gsub("^g", "", .) %>% as.numeric,
        fill = knockoff_stat
      ) 
    ) + 
    geom_point(
      aes(
        x = Gene1 %>% gsub("^g", "", .) %>% as.numeric,
        y = Gene2 %>% gsub("^g", "", .) %>% as.numeric,
        alpha = q_value < 0.25 
      ), 
      colour = "red"
    ) + 
    xlab("Target") + ylab("TF") + 
    ggtitle("Pattern of discoveries") 
  print(p)
  ggsave(paste0(arguments$outFile, "_discoveries.pdf"), p, width = 6, height = 6)
})

# Write output to a file
outDF <- DF[order(DF$q_value, decreasing=FALSE), ]
write.table(outDF, arguments$outFile, sep = "\t", quote = FALSE, row.names = FALSE)
warnings()
