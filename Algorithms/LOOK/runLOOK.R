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

# For manual inspection of the process
# arguments = list(
#   expressionFile = "~/Desktop/jhu/research/projects/Beeline/inputs/Synthetic_with_protein_and_velocity/dyn-LL/dyn-LL-500-1/LOOK/ExpressionData.csv",
#   calibrate = T
# )

# Optional transformation to gaussian marginals
makeMarginalsGaussian = function(X){
  div_by_max = function(x) (x-0.5)/max(x)
  for(i in seq(nrow(X))){
    X[i,] = rank( X[i,] , ties.method = "random" ) %>% div_by_max %>% qnorm
  }
  X
}
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

# Optional smoothing
# neighbors = FNN::get.knn(t(inputExpr), k = 20)
# inputExpr = sapply(seq(nrow(neighbors$nn.dist)),
#                    function(i) rowMeans(inputExpr[,c(i, neighbors$nn.index[i,])]) )

# Apply "genomic control" if input options say so
arguments$center_knockoff_stats = T
knockoffEmpiricalCorrection = function(w){
  if(arguments$center_knockoff_stats){
    inflation = max( median(w), 0 )
  } else {
    inflation = 0
  }
  w - inflation
}

runCalibrationCheck = function(X, noiselevel = 1){
  if(arguments$calibrate){
    X_k = rlookc::computeGaussianKnockoffs(X = X,
                                           mu = 0,
                                           Sigma = cor(X),
                                           num_realizations = 1)
    diverse_y = rlookc::chooseDiverseY(X, n_quantiles = 20)
    calibration_results = rlookc::findWorstY(
      X,
      X_k,
      y = diverse_y$y,
      ground_truth = diverse_y$ground_truth
    )

    # rlookc::simulateY(
    #   X = X,
    #   knockoffs = knockoffs,
    #   plot_savepath = paste0(arguments$outFile, "_calibration.pdf"),
    #   # Heteroskedasticity Hunter
    #   FUN = "adversarial",
    #   kmeans_centers =5
    #
    #   # Univariate sigmoidal
    #   # active_set_size = 1,
    #   # FUN = function(x) 1/(1+exp(-10*(x-1.5)))
    #
    #   # Bivariate bool-ish
    #   # active_set_size = 2,
    #   # FUN = function(x) all(x>0) + rbinom(n = 1, size = noiselevel, prob = 0.5)
    # )
    saveRDS(calibration_results, paste0(arguments$outFile, "_calibration.Rda"))
    return(invisible(calibration_results))
  }
}


# Core functionality: GRN inference via knockoff-based tests
# of carefully constructed null hypotheses
arguments$method =  "steady_state" # "rna_production_protein_predictor_mixture" #
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
        Gene1 = geneNames[-i],
        Gene2 = geneNames[ i],
        knockoff_stat = knockoffResults[[i]]
      )
    }
    DF = data.table::rbindlist(DF)
    DF[["q_value"]] = rlookc::knockoffQvals(DF[["knockoff_stat"]], offset = 0)
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
    calib = runCalibrationCheck(X)
    plot(calib$calibration$targeted_fdrs, colMeans(calib$calibration$fdr))
    abline(a = 0, b = 1)

    # Do each gene separately
    q = w = list()
    for(k in seq_along(geneNames)){
      y = t(inputRNAvelocity)[,k]
      # Infer the decay rate in a robust way (piecewise linear)
      dir.create(file.path(dirname(arguments$outFile), "decay_estimation"), recursive = T, showWarnings = F)
      pdf(file.path(dirname(arguments$outFile), "decay_estimation", paste0("g", k, ".pdf")))
      {
        concentration = inputRNA[k,]
        plot(concentration, y, pch = ".")
        concentration_bins = cut(concentration, breaks = 10)
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
        decay_rate %<>% sapply(extract2, "concentration[idx]") %>% nona %>% negative_only %>% quantile(0.2)
        clip(min(concentration), max(concentration[idx]), y1 = -100, y2 = 100)
        abline(a = 0, b = decay_rate, col = "red")
      }
      dev.off()

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

    # w %<>% lapply(knockoffEmpiricalCorrection)
    # Assemble results
    DF = list()
    for(k in seq_along(geneNames)){
      keep = seq_along(geneNames)[-k] # disallow autoregulation
      DF[[k]] = data.frame(
        Gene1 = geneNames[ keep],
        Gene2 = geneNames[k],
        knockoff_stat = w[[k]][keep]
      )
    }
    DF = data.table::rbindlist(DF)
    DF[["q_value"]] = rlookc::knockoffQvals(DF[["knockoff_stat"]], offset = 0)

  }
  else if(arguments$method == "rna_production_protein_predictor_mixture" )
  {
    stopifnot("Protein levels must be provided with prefix 'p_'.          \n"=nrow(inputProtein)>0)
    stopifnot("Velocity levels must be provided with prefix 'velocity_x_'.\n"=nrow(inputRNAvelocity)>0)
    # Generate knockoffs for protein levels
    X = t(inputProtein)
    # Fit mixture model knockoffs using the BIC criterion to select how many clusters
    library(mclust)
    mixtureModel = mclust::Mclust(X, G = 100, modelNames = "EII")
    mixtureModel$z %>% image
    mus    = lapply( 1:mixtureModel$G, function(g) mixtureModel$parameters[["mean"    ]]           [,g] )
    sigmas = lapply( 1:mixtureModel$G, function(g) mixtureModel$parameters[["variance"]][["sigma"]][,,g] )
    crap_knockoffs = rlookc::computeGaussianKnockoffs(X, output_type = "knockoffs")
    knockoffs = rlookc::computeGaussianMixtureKnockoffs(X, mus, sigmas, posterior_probs = mixtureModel$z, output_type = "knockoffs")
    # Plotting code: assess mixture model fit
    try(silent = T, {
      exprByCluster = inputProtein %>%
        t %>%
        as.data.frame() %>%
        cbind( cluster = apply( mixtureModel$z, 1, which.max ) ) %>%
        cbind( time = inputPT ) %>%
        tidyr::pivot_longer(seq(ncol(X)), names_to = "gene", values_to = "expression")
      ggplot(exprByCluster) +
        geom_point(aes(x = time, y = expression, color = cluster)) +
        facet_wrap(~gene) +
        ggtitle("Gaussian mixture model fit")
        ggsave(file.path(dirname(arguments$outFile), "mixture_model_fit.pdf"), width = 8, height = 8)

      g1 = inputProtein %>% rownames         %>% extract2(1)
      other_genes =  inputProtein %>% rownames %>% extract(c(2, floor( nrow(inputProtein)/2 ), nrow(inputProtein)))
      for( g2 in other_genes){
        title = paste("Genes", g1, "and", g2)
        inputProtein %>%
          t %>%
          as.data.frame() %>%
          cbind( cluster = apply( mixtureModel$z, 1, which.max ) %>% paste0("cluster", .) ) %>%
          cbind( time = inputPT ) %>%
          ggplot() +
          geom_point(aes_string(x = g1, y = g2, color = "time")) +
          facet_wrap(~cluster) +
          ggtitle(title)
          ggsave(file.path(dirname(arguments$outFile), paste0("mixture_model_fit_", title, ".pdf")),
                 width = 8, height = 8)
      }
    })

    # Optional calibration check
    calib = runCalibrationCheck(X)
    plot(calib$calibration$targeted_fdrs, colMeans(calib$calibration$fdr))
    abline(a = 0, b = 1)

    # Do each gene separately
    q = w = list()
    for(k in seq_along(geneNames)){
      y = t(inputRNAvelocity)[,k]
      # Infer the decay rate in a robust way (piecewise linear)
      dir.create(file.path(dirname(arguments$outFile), "decay_estimation"), recursive = T, showWarnings = F)
      pdf(       file.path(dirname(arguments$outFile), "decay_estimation", paste0("g", k, ".pdf")))
      {
        concentration = inputRNA[k,]
        plot(concentration, y, pch = ".")
        concentration_bins = cut(concentration, breaks = 10)
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
        decay_rate %<>% sapply(extract2, "concentration[idx]") %>% nona %>% negative_only %>% quantile(0.2)
        clip(min(concentration), max(concentration[idx]), y1 = -100, y2 = 100)
        abline(a = 0, b = decay_rate, col = "red")
      }
      dev.off()

      # subtract off decay rate; only production rate remains to be modeled
      y = y - concentration*decay_rate
      w[[k]] = knockoff::stat.glmnet_lambdasmax(X, knockoffs, y)

      # Plotting code for ad hoc use: knockoffs vs X
      # This will fail often for boring reasons but it also produces a lot of useful plots.
      goodness_of_fit_plots = file.path(dirname(arguments$outFile), "knockoff_goodness_of_fit_plots")
      dir.create(goodness_of_fit_plots, recursive = T, showWarnings = F)
      try(silent = T, {
        plot_data = data.frame(
          production = y,
          protein_regulator = X[,k-1],
          protein_product = X[,k],
          protein_regulator_mixture_knockoff = knockoffs[,k-1],
          protein_product_mixture_knockoff   = knockoffs[,k],
          protein_regulator_gaussian_knockoff = crap_knockoffs[,k-1],
          protein_product_gaussian_knockoff   = crap_knockoffs[,k],
          time = inputPT
        )
        plot_data %<>%
          tidyr::pivot_longer(cols = protein_regulator:protein_product_gaussian_knockoff,
                              names_to = "feature",
                              values_to = "expression")
        # A null variable vs Y
        ggplot(plot_data %>% subset(!grepl("regulator", feature))) +
          geom_point(aes(x = expression, y = production, colour = feature, shape = feature)) +
          ggtitle(paste0("Candidate regulators and their knockoffs versus gene", k, " production rate"),
                  "Regulator of gene k is assumed to be k-1.\nThis works only for linear long (and even then, not for gene 1.)\nSorry.")
          ggsave(goodness_of_fit_plots %>% file.path(paste0("association_", k, ".pdf")),
                 width = 6, height = 4)
        # A variable over time
        ggplot(plot_data %>% subset(!grepl("product", feature))) +
          geom_point(aes(x = time, y = expression, colour = feature, shape = feature)) +
          geom_smooth(aes(x = time, y = expression, colour = feature, shape = feature), se = F) +
          ggtitle(paste0("Protein ", k, " expression and knockoffs over time")) +
          scale_color_manual(values = c("protein_regulator" ="black",
                                        "protein_regulator_gaussian_knockoff" ="red",
                                        "protein_regulator_mixture_knockoff" = "blue"))
          ggsave(goodness_of_fit_plots %>% file.path(paste0("timeseries_", k, ".pdf")),
                 width = 6, height = 4)


      } )
    }

    # w %<>% lapply(knockoffEmpiricalCorrection)
    # Assemble results
    DF = list()
    for(k in seq_along(geneNames)){
      keep = seq_along(geneNames)[-k] # disallow autoregulation
      DF[[k]] = data.frame(
        Gene1 = geneNames[keep],
        Gene2 = geneNames[ k],
        knockoff_stat = w[[k]][keep],
        q_value = rlookc::knockoffQvals(w[[k]][keep], offset = 0)
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
