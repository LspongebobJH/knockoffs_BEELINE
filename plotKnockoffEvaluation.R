setwd("~/Desktop/jhu/research/projects/knockoffs/applications/Beeline")
library(magrittr)
library(ggplot2)
ggplot2::theme_update(text = element_text(family = "ArialMT"))
knockoff_type_color_mapping =
  scale_color_manual(values = c("y=x" ="black",
                                "geneNet" = "forestgreen",
                                "knockoffs (naive)" = "blue",
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
    tidyr::pivot_longer(cols = X0:X9, values_to = "empirical_fdr", names_to = "targeted_fdr") %>%
    dplyr::mutate(targeted_fdr = gsub("^X", "", targeted_fdr) %>% as.numeric %>% divide_by(10)) %>%
    tidyr::separate(X, into = c(NA, "network", "cellcount", "replicate")) %>%
    tidyr::separate(method, into = c("method", "knockoff_type", "data_mode"), sep = "_") %>%
    subset(method %in% c("GENENET", "LOOK"))  %>%
    dplyr::mutate(is_knockoff_based = method=="LOOK") %>%
    dplyr::mutate(data_mode =  ifelse(!is_knockoff_based, knockoff_type, data_mode)) %>%
    dplyr::mutate(knockoff_type =  ifelse(!is_knockoff_based, "", knockoff_type)) %>%
    dplyr::mutate(method_summary =  paste0(method, " (", knockoff_type, ")") ) %>%
    dplyr::mutate(method_summary = gsub("GENENET ()", "geneNet", method_summary, fixed = T)) %>%
    dplyr::mutate(method_summary = gsub("LOOK", "knockoffs", method_summary, fixed = T)) 
  ggplot(plot_data) +
    geom_smooth(aes(x = targeted_fdr, y = empirical_fdr, colour = method_summary, group = method_summary), se = F) +
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

# BEELINE metrics
metric = "aupr"
aupr = grabResults(metric)
plot_data = list()
for(fname in names(aupr)){
  aupr[[fname]] %<>%
    tidyr::pivot_longer(cols = !X, values_to = "metric") %>%
    tidyr::separate(name, into = c(NA, "network", "cellcount", "replicate")) %>%
    dplyr::rename(method=X) %>%
    dplyr::mutate(cellcount = factor(cellcount, levels = gtools::mixedsort(unique(cellcount)))) %>%
    tidyr::separate(method, into = c("method", "knockoff_type", "data_mode"), sep = "_") %>%
    subset(!is.na(knockoff_type))
}
aupr %<>% data.table::rbindlist()
ggplot(aupr) +
  geom_violin(aes(x = knockoff_type, y = metric, group = knockoff_type)) +
  geom_point(aes(x = knockoff_type, y = metric, group = knockoff_type)) +
  facet_grid(data_mode~network) +
  ylab(metric) +
  ggtitle(paste0(metric, " on BEELINE simple simulations"))
ggsave(paste0(metric, ".pdf"), height = 4, width = 8)


# Calibration checks based on simulated Y
calibration_checks = grabResults(pattern = "calibration.Rda",
                                 reader = readRDS)
plot_data = list()
for(fname in names(calibration_checks)){
  plot_data[[fname]] = data.frame(
    empirical_fdr = calibration_checks[[fname]]$calibration$fdr %>% colMeans,
    targeted_fdr  = calibration_checks[[fname]]$calibration$targeted_fdrs,
    network = basename(dirname(dirname(dirname(fname)))),
    replicate = basename(dirname(dirname(fname))) %>% strsplit("-") %>% extract2(1) %>% extract2(4),
    cellcount = basename(dirname(dirname(fname))) %>% strsplit("-") %>% extract2(1) %>% extract2(3)
  )
}
plot_data %<>% data.table::rbindlist()
ggplot(plot_data) +
  geom_point(aes(x = targeted_fdr, y = empirical_fdr, colour = cellcount, shape = cellcount)) +
  geom_smooth(aes(x = targeted_fdr, y = empirical_fdr, colour = cellcount, shape = cellcount)) +
  facet_wrap(~network) +
  ggtitle("Calibration on BEELINE simple network simulations") +
  geom_abline(slope=1, intercept = 0)

