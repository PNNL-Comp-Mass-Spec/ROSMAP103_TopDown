

library(tidyverse)
library(MSnSet.utils)
library(purrr)
library(writexl)



cols <- c("Gene", 
          "UniProtAcc",
          "AnnType",
          "proteoform_id",
          "Proteoform",
          "firstAA",
          "lastAA",
          "ModsString",
          "#unexpected modifications",
          "mass",
          "rt",
          "collision",
          "count",
          "gpath.apv",
          "cogng_demog_slope.apv")




# intensity data
load("./output_data/shotgun_topdown_int_20240730_modann_cnt_wres.RData")
mi <- m
x <- fData(m)
x$ModsString <- map_chr(fData(m)$mods, \(x) x[['mods_str']])
x <- x[,cols]
x <- cbind(x, signif(exprs(m), digits = 3))
write_xlsx(x, path="./figures_tables/Supplementary Table S1.xlsx")

# key nums
nrow(m)
length(unique(fData(m)$Gene))



# spectral counts
load("./output_data/shotgun_topdown_sc_20240730_modann_wres.RData")
mc <- m
x <- fData(m)
x$ModsString <- map_chr(fData(m)$mods, \(x) x[['mods_str']])
x <- x[,cols]
x <- cbind(x, signif(exprs(m), digits = 3))
write_xlsx(x, path="./figures_tables/Supplementary Table S2.xlsx")

# key nums
nrow(m)
length(unique(fData(m)$Gene))
