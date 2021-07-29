# Create data for each synthetic network from a sparse Gaussian MRF.
library("magrittr")
setwd("~/Desktop/jhu/research/projects/Beeline/inputs/Synthetic_gaussian")
toy_networks = list.files(pattern = "^dyn")
for(toy in toy_networks){
  # read base network
  network_structure = read.csv(file.path(toy, paste0(toy, "-500-1"), "refNetwork.csv"))
  network_structure$Gene1 %<>% gsub("^g", "", .) %>% as.numeric
  network_structure$Gene2 %<>% gsub("^g", "", .) %>% as.numeric
  n_genes = pmax(
    max(network_structure$Gene1),
    max(network_structure$Gene2)
  )
  # fill in precision matrix with mostly zeroes
  precision = diag(n_genes)
  for(edge in seq_along(network_structure$Gene1)){
    if(network_structure[edge,"Gene1"] != network_structure[edge,"Gene2"]){
      precision[network_structure[edge,"Gene1"],
                network_structure[edge,"Gene2"]] = -0.2
      precision[network_structure[edge,"Gene2"],
                network_structure[edge,"Gene1"]] = -0.2
    }
  }
  # Simulate data with this precision matrix
  set.seed(0)
  M = chol(solve(precision))
  for(rep in 1:10){
    for(n_cell in c(100, 200, 500, 2000, 5000)){
      X = t(M) %*% matrix(rnorm(n_cell*n_genes, mean = 0, sd = 1), nrow = n_genes)
      rownames(X) = paste0("g", 1:nrow(X))
      # BEELINE violates DRY pretty flagrantly here and now it's my problem.
      write.csv(X, file.path(toy, paste(toy, n_cell, rep, sep = "-"), "ExpressionData.csv"))
      write.csv(X, file.path(toy, paste(toy, n_cell, rep, sep = "-"), "PPCOR", "ExpressionData.csv"))
      write.csv(X, file.path(toy, paste(toy, n_cell, rep, sep = "-"), "LOOK", "ExpressionData.csv"))
      # Visual check for correctness of precision matrix
      remove_diagonal_that_otherwise_ruins_colorscale = function(X){
        diag(X) = mean(X)
        X
      }
      image( main = paste(toy, n_cell, rep, sep = "-"),
             rbind( solve(cov(t(X))) %>% remove_diagonal_that_otherwise_ruins_colorscale,
                    precision %>% remove_diagonal_that_otherwise_ruins_colorscale ) )
    }
  }
}
