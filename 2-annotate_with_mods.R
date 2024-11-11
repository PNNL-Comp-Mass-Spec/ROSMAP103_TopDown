
library(MSnSet.utils)
library(tidyverse)
library(TopPICR)



load("./output_data/shotgun_topdown_int_20240730.RData")
m <- annotate_modifications(m, mod_file = "./source_data/TopPIC_Dynamic_Mods.txt")
save(m, file="./output_data/shotgun_topdown_int_20240730_modann.RData")



load("./output_data/shotgun_topdown_sc_20240730.RData")
m <- annotate_modifications(m, mod_file = "./source_data/TopPIC_Dynamic_Mods.txt")
save(m, file="./output_data/shotgun_topdown_sc_20240730_modann.RData")




