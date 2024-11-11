
library(MSnSet.utils)


load("./output_data/shotgun_topdown_int_20240730_modann.RData")
mi <- m
load("./output_data/shotgun_topdown_sc_20240730_modann.RData")

fData(mi)$spectralCount <- NULL
fData(mi)$count <- fData(m[featureNames(mi),])$count
m <- mi

save(m, file="./output_data/shotgun_topdown_int_20240730_modann_cnt.RData") # overwrite


