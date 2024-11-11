
library(MSnSet.utils)
library(tidyverse)
library(TopPICR)


# almost certainly to be artifacts due to long tailing of the other, more abundant species
to_remove <- c("APP_31", "APP_55", "APP_58", "APP_57")


load("./output_data/shotgun_topdown_int_20240730.RData")
dim(m)
m <- m[!(featureNames(m) %in% to_remove),]
dim(m)
save(m, file="./output_data/shotgun_topdown_int_20240730.RData")




load("./output_data/shotgun_topdown_sc_20240730.RData")
dim(m)
m <- m[!(featureNames(m) %in% to_remove),]
dim(m)
save(m, file="./output_data/shotgun_topdown_sc_20240730.RData")


