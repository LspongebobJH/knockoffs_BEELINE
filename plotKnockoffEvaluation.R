setwd("~/Desktop/jhu/research/projects/knockoffs/applications/Beeline")
library(magrittr)
library(ggplot2)
ggplot2::theme_update(text = element_text(family = "ArialMT"))
knockoff_type_color_mapping =
  scale_color_manual(values = c("identical" ="black",
                                "naive" = "blue",
                                "gaussian" ="red",
                                "mixture" = "orange"))
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

# Calibration checks based on BEELINE ground truth
for(metric in c("FDR", "undirectedFDR")){
  plot_data = grabResults(pattern = metric) %>%
    data.table::rbindlist() %>%
    tidyr::pivot_longer(cols = X0:X9, values_to = "empirical_fdr", names_to = "targeted_fdr") %>%
    dplyr::mutate(targeted_fdr = gsub("^X", "", targeted_fdr) %>% as.numeric %>% divide_by(10)) %>%
    tidyr::separate(X, into = c(NA, "network", "cellcount", "replicate")) %>%
    tidyr::separate(method, into = c("method", "knockoff_type", "data_mode"), sep = "_") %>%
    subset(!is.na(knockoff_type))
  ggplot(plot_data) +
    geom_smooth(aes(x = targeted_fdr, y = empirical_fdr, colour = knockoff_type, group = knockoff_type), se = F) +
    facet_grid(data_mode ~ network) +
    ylab(metric) +
    ggtitle("Calibration on BEELINE simple network simulations") +
    geom_abline(aes(slope = 1, intercept = 0)) +
    knockoff_type_color_mapping + 
    scale_x_continuous(breaks = setNames((0:2)/2, c("0", "0.5", "1")), limits = 0:1) 
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

