library(tidyverse)
setwd("projects/def-wainberg/karbabi/maze-neptune")

de_a = read_csv("output/DE_Additive_Genotype_table.csv") %>%
    mutate(model = "Additive_Genotype")
de_d = read_csv("output/DE_Dominant_Genotype_table.csv") %>%
    mutate(model = "Dominant_Genotype")
de_r = read_csv("output/DE_Recessive_Genotype_table.csv") %>%
    mutate(model = "Recessive_Genotype")

de_table = rbind(de_a, de_d, de_r)

de_table %>%  
    mutate(Dir = factor(sign(logFC)), cell_type = as.factor(cell_type)) %>%
    mutate(Dir = recode_factor(Dir, `-1` = "Down", `1` = "Up")) %>%
    group_by(model, cell_type, Dir) %>% 
    dplyr::summarize(
        "pt01" = sum(FDR <= 0.01),
        "pt05" = sum(FDR <= 0.05),
        "pt10" = sum(FDR <= 0.10)) %>%
    ungroup() %>%
    pivot_longer(cols = c(pt01, pt05, pt10),
        values_to = "Freq",
        names_to = "FDR_thresh") %>%
    mutate(Freq = if_else(Dir == "Down", -Freq, Freq)) %>%
    mutate(Freq = if_else(Freq == 0, NA, Freq)) %>%
    filter(!is.na(Freq)) %>%
    mutate(model = str_replace_all(model, "_", " ")) %>%
ggplot(., aes(x = reorder(cell_type, abs(Freq)), y = Freq, fill = Dir)) +
    geom_bar(data = . %>% filter(FDR_thresh == "pt01"),
        stat = "identity", width = 0.8, alpha = 1) +
    geom_bar(data = . %>% filter(FDR_thresh == "pt05"), 
        stat = "identity", width = 0.8, alpha = 0.7) +
    geom_bar(data = . %>% filter(FDR_thresh == "pt10"),
        stat = "identity", width = 0.8, alpha = 0.5) +
    geom_text(data = . %>% filter(FDR_thresh == "pt10"), 
            aes(label = abs(Freq), y = Freq-5), 
            hjust = -0.1, size = 3) +
    geom_tile(aes(y = NA_integer_, alpha = factor(FDR_thresh))) + 
    scale_fill_manual(values = c("#D62728FF", "#1F77B4FF")) +
    scale_alpha_manual("FDR Threshold", values = c(1, 0.7, 0.5), 
        labels = c("0.01","0.05","0.10")) + 
    labs(y = "Number of DE genes", x = "", fill = "Direction") +
    coord_flip() +
    facet_wrap(~ model, ncol = 1, drop=F) +
    theme_classic() +
      theme(strip.text = element_text(size = 12, color = "black"),
        axis.text = element_text(size = 12, color = "black"))

ggsave("figures/number_of_degs.png", height = 8, width = 7)


