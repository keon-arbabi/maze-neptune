library(tidyverse)
setwd("projects/def-wainberg/karbabi/maze-neptune")

de_a = read_csv("output/DE_Additive_Genotype/table.csv") %>%
    mutate(model = "Additive_Genotype")
de_d = read_csv("output/DE_Dominant_Genotype/table.csv") %>%
    mutate(model = "Dominant_Genotype")
de_r = read_csv("output/DE_Recessive_Genotype/table.csv") %>%
    mutate(model = "Recessive_Genotype")

de_table = rbind(de_a, de_d, de_r)
