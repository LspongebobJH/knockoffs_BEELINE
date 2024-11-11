suppressPackageStartupMessages({
  library(ggplot2)
  library(magrittr)
  library(optparse)
  library(mclust)
  library(GeneNet)
})
option_list <- list (
  make_option(c("-e","--expressionFile"), type = 'character',
              help= "Path to comma separated file containing gene-by-cell matrix with
              cell names as the first row and gene names as
              the first column. Required."),
  make_option(c("-o","--outFile"), type = 'character',
              help= "outFile name to write the output ranked edges. Required."),
  make_option(c("-c","--calibrate"), action = 'store', default = FALSE,
              type = 'logical',
              help= "Run a simple check with additional simulated target genes?"),
  make_option(c("-d","--data_mode"),
              type = 'character',
              help= "Data to make available. Options are 'easy' (protein concentration and RNA production revealed) and  'rna_only'.")

  )
parser <- OptionParser(option_list = option_list)
arguments <- parse_args(parser, positional_arguments = FALSE)

# # For manual testing of the script
# arguments = list(
#   expressionFile = "~/Desktop/jhu/research/projects/knockoffs/applications/Beeline/inputs/Synthetic/dyn-LL/dyn-LL-500-1/GENENET/ExpressionData.csv",
#   outFile = "~/Desktop/jhu/research/projects/knockoffs/applications/Beeline/temp/outputs.txt",
#   calibrate = T,
#   data_mode = "rna_only"
# )
# dir.create(dirname(arguments$outFile), showWarnings = F, recursive = T)

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
inputPT   = inputPT[order(inputPT[[1]]),1]
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

# Robust heuristic to separate velocity into decay + production; based on
# piecewise quantile regression of RNA velocity against RNA concentration.
getProductionRate = function(inputRNAvelocity, inputRNA){
  production = inputRNAvelocity * 0
  for(k in seq_along(geneNames)){
    y = t(inputRNAvelocity)[,k]
    # Diagnostic plots will be saved; we add each quant reg separately so the pdf stays open. Ugly; sorry.
    dir.create(file.path(dirname(arguments$outFile), "decay_estimation"), recursive = T, showWarnings = F)
    pdf(       file.path(dirname(arguments$outFile), "decay_estimation", paste0("g", k, ".pdf")))
    {
      concentration = inputRNA[k,]
      plot(concentration, y, pch = ".")
      concentration_bins = cut(concentration, breaks = 10)
      decay_rate_constant = list()
      for(bin in levels(concentration_bins)){
        idx = concentration_bins==bin
        if(sum(idx)<10){next}
        decay_rate_constant[[bin]] = coef( quantreg::rq(y[idx] ~ concentration[idx] ) )
        clip(min(concentration[idx]), max(concentration[idx]), y1 = -100, y2 = 100)
        abline(decay_rate_constant[[bin]][[1]], decay_rate_constant[[bin]][[2]])
      }
      nona = function(x) x[!is.na(x)]
      negative_only = function(x) x[x<0]
      decay_rate_constant %<>% sapply(extract2, "concentration[idx]") %>% nona %>% negative_only %>% median
      clip(min(concentration), max(concentration[idx]), y1 = -100, y2 = 100)
      abline(a = 0, b = decay_rate_constant, col = "red")
    }
    dev.off()
    production[k, ] = y - concentration*decay_rate_constant
  }
  return(production)
}

# Simulate additional target genes with simple dependence on X to check subset selection fdr control
runCalibrationCheck = function(X, noiselevel = 1){
  if(arguments$calibrate){
    diverse_y = rlookc::calibrate__chooseDiverseY(X)
    y = diverse_y$y %>% Reduce(cbind, .)
    Xy = cbind(X,y)
    m.pcor = GeneNet::ggm.estimate.pcor(Xy)
    fdr = GeneNet::network.test.edges(m.pcor)
    qvals = list()
    for(y_idx in seq_along(diverse_y$ground_truth)){
      cat(".")
      y_idx_shifted = y_idx + ncol(X)
      qvals[[y_idx]] =
        subset(fdr, node1==y_idx_shifted | node2==y_idx_shifted ) %>%
        dplyr::mutate(other_node = ifelse(node1==y_idx_shifted, node2, node1)) %>%
        subset(other_node <= ncol(X)) %>%
        dplyr::arrange(other_node) %>%
        extract2("qval")
    }
    calibration_results = rlookc::checkCalibration(ground_truth = diverse_y$ground_truth, 
                                                   qvals = qvals)
    saveRDS(calibration_results, file.path(dirname(arguments$outFile), "calibration.Rda"))
    return(invisible(calibration_results))
  }
}

# Core functionality: GRN inference via tests of nonzero partial correlation
applyGeneNet = function(data_mode = arguments$data_mode){
  if( data_mode == "easy" ){
    stopifnot("Protein levels must be provided with prefix 'p_'.          \n"=nrow(inputProtein)>0)
    stopifnot("Velocity levels must be provided with prefix 'velocity_x_'.\n"=nrow(inputRNAvelocity)>0)
    runCalibrationCheck(X = t(inputProtein))
    inputRNAproduction = getProductionRate(inputRNAvelocity = inputRNAvelocity, inputRNA = inputRNA)
    Xy = cbind(t(inputProtein),t(inputRNAproduction))
    m.pcor = GeneNet::ggm.estimate.pcor(Xy)
    fdr = GeneNet::network.test.edges(m.pcor)
  } else if (data_mode == "rna_only") {
    Xy = cbind(t(inputRNA))
    m.pcor = GeneNet::ggm.estimate.pcor(Xy)
    fdr = GeneNet::network.test.edges(m.pcor)
    runCalibrationCheck(X = t(inputRNA))
  } else {
    stop("data_mode not recognized")
  }
  DF = fdr
  colnames(DF) = c("pcor", "Gene1", "Gene2", "pval", "q_value", "prob")
  DF[["Gene1"]] = geneNames[DF[["Gene1"]]]
  DF[["Gene2"]] = geneNames[DF[["Gene2"]]]
  return(DF[c("Gene1", "Gene2", "pcor", "q_value")])
}

# DO IT
DF = applyGeneNet()

# Write output to a file
outDF <- DF[order(DF$q_value, decreasing=FALSE), ]
write.table(outDF, arguments$outFile, sep = "\t", quote = FALSE, row.names = FALSE)
warnings()
