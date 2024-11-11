
library(knitr)
library(TopPICR)
library(tidyverse)
library(MSnSet.utils)
library(ggplot2)


load("./output_data/wgcna_clusters.RData")
load("./output_data/shotgun_topdown_int_20240730_modann_cnt_wres.RData")



app_clusters <- c("coral", "darkolivegreen4", "mediumpurple1")
x <- filter(clstrs, cluster %in% app_clusters)
x <- inner_join(x, fData(m))


x$cluster <- ordered(x$cluster, levels = app_clusters)
ggplot(x) +
   aes(y = count, x = cluster, fill = cluster) +
   geom_boxplot() +
   ylab(bquote(log[10](spectral~count))) +
   xlab(NULL) +
   scale_y_log10() +
   scale_fill_manual(values = levels(x$cluster),
                     labels=c('APP-1', 'APP-2', "APP-3")) +
   theme_bw()


ggsave(path = "figures_tables", filename = "FigS9_WGCNA_APP_abundances.png", width = 9, height = 6, scale = 0.55)




vgf_clusters <- c("yellowgreen", "lightcoral")
x <- filter(clstrs, cluster %in% vgf_clusters)
x <- inner_join(x, fData(m))


x$cluster <- ordered(x$cluster, levels = vgf_clusters)

ggplot(x) +
   aes(y = count, x = cluster, fill = cluster) +
   geom_boxplot() +
   ylab(bquote(log[10](spectral~count))) +
   xlab(NULL) +
   scale_y_log10() +
   scale_fill_manual(values = levels(x$cluster),
                     labels=c('VGF-1', 'VGF-2')) +
   theme_bw()

ggsave(path = "figures_tables", filename = "FigS9_WGCNA_VGF_abundances.png", width = 7.5, height = 6, scale = 0.55)


