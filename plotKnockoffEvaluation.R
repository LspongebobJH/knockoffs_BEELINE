setwd("~/Desktop/jhu/research/projects/knockoffs/applications/Beeline")
library(magrittr)
library(ggplot2)
ggplot2::theme_update(text = element_text(family = "ArialMT"))
methods_include_figure_a = c(
  #"GENENET", 
  "LOOK")
knockoff_type_color_mapping =
  scale_color_manual(limits = force, 
                     values = c("y=x" = "black",
                                "original data" = "black",
                                "geneNet" = "forestgreen",
                                "knockoffs (permuted)" = "blue",
                                "knockoffs (gaussian)" ="red",
                                "knockoffs (mixture)" = "orange"))
grabResults = function(pattern,
                       reader = read.csv,
                       base_dir = "outputs/Synthetic_with_protein_and_velocity",
                       ...){
  x = list.files(base_dir,
                 pattern = paste0("(-|_)", pattern),
                 ignore.case = T,
                 full.names = T,
                 recursive = T)
  lapply(x, reader, ...) %>% setNames(x)
}

plot_protein_vs_rna = function(other_quantity = "protein"){
  prefix = switch(other_quantity, 
                  protein="^p_",
                  velocity="^velocity_x_")
  expressionData =
    "inputs/Synthetic_with_protein_and_velocity" %>%
    list.files(full.names = T) %>%
    list.files(pattern = "dyn-", full.names = T) %>%
    list.files(pattern = "ExpressionData", full.names = T) 
  expressionData = as.list(expressionData) %>% setNames(basename(dirname(expressionData)))
  expressionData %<>% lapply(read.csv, row.names = 1)
  get_rna_protein_correlation = function(X){
    index_protein = grep(prefix, rownames(X))
    index_rna     = grep("^x_", rownames(X))
    data.frame(
      Pearson_correlation = diag( cor( 
        t( X[index_protein,] ),
        t( X[index_rna,] )
      ) ),
      gene = rownames(X[index_rna,]) %>% gsub("^x_g", "", .) %>% as.numeric
      )
  }
  rna_protein_correlation = lapply( expressionData, get_rna_protein_correlation )
  rna_protein_correlation = 
    mapply( 
      SIMPLIFY = F,
      function(X, dataset){ X[["dataset"]] = dataset; X },
      X = rna_protein_correlation,
      dataset = names(rna_protein_correlation) 
    ) %>%
    data.table::rbindlist()
  rna_protein_correlation %<>% tidyr::separate(dataset, into = c("prefix", "network", "cellcount", "replicate"))
  ggplot(rna_protein_correlation) + 
    geom_jitter(aes(x = gene, group = gene, y = Pearson_correlation)) + 
    facet_wrap(~network, scale = "free_x") +
    ylab("Pearson correlation") +
    ggtitle(paste0("RNA versus ", other_quantity, " levels"))
}
plot_protein_vs_rna(other_quantity = "protein")
ggsave("rna_vs_protein.pdf", height = 4, width = 6)
plot_protein_vs_rna(other_quantity = "velocity")
ggsave("rna_vs_velocity.pdf", height = 4, width = 6)

# Calibration checks based on BEELINE ground truth
for(metric in c("FDR", "undirectedFDR")){
  plot_data = grabResults(pattern = metric) %>%
    data.table::rbindlist() %>%
    tidyr::pivot_longer(cols = X0:X9, values_to = "observed_fdr", names_to = "expected_fdr") %>%
    dplyr::mutate(expected_fdr = gsub("^X", "", expected_fdr) %>% as.numeric %>% divide_by(10)) %>%
    tidyr::separate(X, into = c(NA, "network", "cellcount", "replicate")) %>%
    tidyr::separate(method, into = c("method", "knockoff_type", "data_mode"), sep = "_") %>%
    subset(method %in% methods_include_figure_a)  %>%
    dplyr::mutate(is_knockoff_based = method=="LOOK") %>%
    dplyr::mutate(data_mode =  ifelse(data_mode=="easy", "RNA+protein", data_mode)) %>%
    dplyr::mutate(data_mode =  ifelse(!is_knockoff_based, knockoff_type, data_mode)) %>%
    dplyr::mutate(knockoff_type =  ifelse(!is_knockoff_based, "", knockoff_type)) %>%
    dplyr::mutate(method_summary =  paste0(method, " (", knockoff_type, ")") ) %>%
    dplyr::mutate(method_summary = gsub("GENENET ()", "geneNet", method_summary, fixed = T)) %>%
    dplyr::mutate(method_summary = gsub("LOOK", "knockoffs", method_summary, fixed = T))  %>%
    dplyr::mutate(knockoff_type = gsub("naive", "permuted", knockoff_type, fixed = T)) %>%
    dplyr::mutate(method_summary = gsub("naive", "permuted", method_summary, fixed = T)) 
  ggplot(plot_data) +
    geom_smooth(aes(x = expected_fdr, y = observed_fdr, colour = method_summary, group = method_summary), se = F) +
    facet_grid(data_mode ~ network) +
    ylab(metric) +
    ggtitle("Calibration on BEELINE simple network simulations") +
    geom_abline(aes(slope = 1, intercept = 0)) +
    knockoff_type_color_mapping + 
    scale_y_continuous(breaks = setNames((0:2)/2, c("0", "0.5", "1")), limits = 0:1) + 
    scale_x_continuous(breaks = setNames((0:2)/2, c("0", "0.5", "1")), limits = 0:1) + 
    coord_fixed()
  ggsave(paste0(metric, ".pdf"), height = 4, width = 8)
}

# Example knockoffs + original data tsne
knockoffs_and_orig_data = Reduce(rbind, list(
  read.csv("outputs/Synthetic_with_protein_and_velocity/dyn-BFC/dyn-BFC-500-1/LOOK_gaussian_easy/knockoffs.csv", row.names = 1),
  read.csv("outputs/Synthetic_with_protein_and_velocity/dyn-BFC/dyn-BFC-500-1/LOOK_mixture_easy/knockoffs.csv", row.names = 1),
  read.csv("outputs/Synthetic_with_protein_and_velocity/dyn-BFC/dyn-BFC-500-1/LOOK_naive_easy/knockoffs.csv", row.names = 1),
  read.csv("outputs/Synthetic_with_protein_and_velocity/dyn-BFC/dyn-BFC-500-1/LOOK_naive_easy/data.csv", row.names = 1)
))
embedding = tsne::tsne(knockoffs_and_orig_data, 2, max_iter = 800) %>% as.data.frame %>% set_colnames(c("tsne1", "tsne2"))
metadata = data.frame(
  type = rep(c("knockoffs (gaussian)", 
                            "knockoffs (mixture)", 
                            "knockoffs (permuted)", 
                            "original data"), each=500),
  "time" = read.csv("inputs/Synthetic_with_protein_and_velocity/dyn-BFC/dyn-BFC-500-1/PseudoTime.csv")[[2]] %>% rep(times=4)
)  
ggplot(cbind(embedding, metadata)) + 
  geom_point(aes(x = tsne1, y = tsne2, color = type, shape = type)) + 
  knockoff_type_color_mapping +
  coord_fixed() + 
  ggtitle("Input data and all knockoffs", "BFC network protein concentration")
ggsave("tsne.pdf", width = 6, height= 6)
# Example knockoffs + original data timeseries
knockoff_types = c(
  "knockoffs (gaussian)", 
  "knockoffs (mixture)", 
  "knockoffs (permuted)", 
  "original data")
knockoffs_and_orig_data = data.table::rbindlist( list(
  read.csv("outputs/Synthetic_with_protein_and_velocity/dyn-CY/dyn-CY-500-1/LOOK_gaussian_easy/knockoffs.csv"),
  read.csv("outputs/Synthetic_with_protein_and_velocity/dyn-CY/dyn-CY-500-1/LOOK_mixture_easy/knockoffs.csv"),
  read.csv("outputs/Synthetic_with_protein_and_velocity/dyn-CY/dyn-CY-500-1/LOOK_naive_easy/knockoffs.csv"),
  read.csv("outputs/Synthetic_with_protein_and_velocity/dyn-CY/dyn-CY-500-1/LOOK_naive_easy/data.csv")
)) %>% 
  set_colnames(c("cell", colnames(.)[-1])) %>%
  dplyr::mutate(
    type=knockoff_types %>% rep(each=500)
  )
metadata = data.frame(
  type=knockoff_types %>% rep(each=500),
  time = read.csv("inputs/Synthetic_with_protein_and_velocity/dyn-CY/dyn-CY-500-1/PseudoTime.csv")[[2]] %>%
    rep(times=4),
  cell = read.csv("inputs/Synthetic_with_protein_and_velocity/dyn-CY/dyn-CY-500-1/PseudoTime.csv")[[1]] %>% 
    rep(times=4)
) 
ggplot(merge(knockoffs_and_orig_data, metadata, by = c("cell", "type"))) + 
  geom_point(aes(x = time, y = g1, color = type, shape = type)) +
  geom_smooth(aes(x = time, y = g1, color = type, group = type), se = F) + 
  knockoff_type_color_mapping +
  ggtitle("Input data and all knockoffs", "CY network protein concentration")
ggsave("time.pdf", width = 8, height= 6)
