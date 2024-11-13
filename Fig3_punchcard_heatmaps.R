
# scale factor for tuning the size of the heatmaps


#original scripts are 6a/6b 
library(tidyverse)
library(MSnSet.utils)
library(dplyr)
library(latex2exp)

###############

#######Intensity

###############

load("./output_data/shotgun_topdown_int_20240730_modann_cnt_wres.RData")

# apply > 50 filter
m <- m[rowSums(!is.na(exprs(m))) > 50, ] # 70 out of 103


#' # Overall test results across clinical and neuropath traits
get_sig_res <- function(m, pheno) {
   modl <- sprintf("~ %s + pmi", pheno)
   
   limma_gen(m, model.str = modl, coef.str = pheno) %>%
      rownames_to_column("proteoform") %>%
      dplyr::select(proteoform, logFC, P.Value, adj.P.Val)
}

vls <- c("sqrt_tangles",
         "sqrt_amyloid",
         "lb_neo",
         "hip_scl_3reg_yn",
         "cogn_global_lv",
         "cogng_demog_slope",
         "ci_num2_gct",
         "caa_4gp",
         "arteriol_scler",
         "anye4")

results <- lapply(vls, function(pheno_i) {
   res_i <- get_sig_res(m, pheno = pheno_i)
   res_i$variable <- pheno_i
   
   return(res_i)
})
results <- data.table::rbindlist(results)

results %>%
   filter(adj.P.Val < 0.05) %>%
   distinct(proteoform) %>%
   nrow() 
# 42 proteoforms show a statistically significant difference in at least one
# variable

results %>%
   filter(adj.P.Val < 0.05) %>%
   distinct(variable, proteoform) %>%
   group_by(variable) %>%
   tally() %>%
   arrange(desc(n))
# Most of these differences are related to sqrt_amyloid

# Filter to those 42 proteoforms
filtered_results <- results %>%
   filter(adj.P.Val < 0.05) %>%
   distinct(proteoform) %>%
   inner_join(results)


filtered_results2 <- results %>%
   group_by(variable) %>%
   summarise(logFC_mean = mean(logFC, na.rm = TRUE),
             logFC_sd = sd(logFC, na.rm = TRUE)) %>%
   inner_join(filtered_results) %>%
   mutate(Zscore = (logFC - logFC_mean) / logFC_sd,
          Zscore2 = case_when(Zscore > 3 ~ 3,
                              Zscore < -3 ~ -3,
                              TRUE ~ Zscore),
          is_sig = ifelse(adj.P.Val < 0.05, 1, NA_real_))

# check collisions
collision_pairs <- fData(m)[unique(filtered_results2$proteoform), ] %>%
   dplyr::filter(collision != "-") %>%
   rownames_to_column(var="feature_name") %>%
   select(feature_name, collision)

# and let's remove those that have collisions
filtered_results2 <- filtered_results2 %>%
   filter(!proteoform %in% collision_pairs$feature_name) %>%
   mutate(variable = gsub("sqrt_amyloid", "amyloid", variable)) %>%
   mutate(variable = gsub("sqrt_tangles", "tangles", variable))
 #changes name of variable   

## Plotting --------------------------------------------------------------------

# Arrange proteoforms by gene and number
df <- filtered_results2 %>% 
   mutate(gene = sub("^(.*)_(\\d+)$", "\\1", proteoform),
          number =  as.numeric(sub("^(.*)_(\\d+)$", "\\2", proteoform))) %>% 
   arrange(gene, number) %>% 
   # Convert to factor to preserve order in the plot
   mutate(proteoform = factor(proteoform, levels = unique(proteoform)))


df <- df %>%
   mutate(variable = ordered(variable, 
                             levels = c("amyloid", "tangles", "anye4", 
                                        "caa_4gp", "ci_num2_gct", 
                                        "hip_scl_3reg_yn", "cogn_global_lv",
                                        "cogng_demog_slope",
                                        "arteriol_scler",
                                        "lb_neo")))

df_int <- df
#geom tile with border
int <- ggplot(data = df, 
              aes(y = proteoform, x = variable, 
                  fill = Zscore2, size = is_sig)) + 
   geom_tile(aes(fill = Zscore2, color=is_sig, width=0.9, height=0.9), size=1) +
   scale_colour_gradient( high = "black")+
   scale_fill_gradientn(colours = colorspace::lighten(gplots::bluered(100), amount = 0.5), 
                        # values = c(0, seq(qn01[1], qn01[2], length.out = 98),1), 
                        limits = c(-3, +3),
                        name = "Z-score") +
   guides(color = FALSE, size = FALSE) +
   geom_point( aes( y = proteoform, x = variable, size=ifelse(is_sig, "dot", "no_dot"))) +
   scale_size_manual(values=c(dot=1.5, no_dot=NA), guide="none") +
   theme_classic() +
   scale_y_discrete(limits=rev) +
   scale_x_discrete(position="top",labels=c(
      "amyloid" = "amyloid",
      "tangles" = "tangles",
      "anye4" = expression(italic("APOE")~ε4), 
      "caa_4gp" = "cerebral amyloid angiopathy",
      "ci_num2_gct" = "macroscopic infarcts",
      "hip_scl_3reg_yn" = "hippocampal sclerosis",
      "cogn_global_lv" = "global cognition", 
      "cogng_demog_slope" = "slope of cognitive decline", 
      "arteriol_scler" = "arteriolosclerosis",
      "lb_neo" = "neocortical Lewy bodies")) +
   theme(aspect.ratio = 2.6*(39/44),
         axis.text.x = element_text(angle = 35, vjust = 0.5, hjust = 0),
         legend.box = "horizontal", 
         # legend.position = "none",
         legend.key = element_rect(color = "black"),
         axis.text = element_text(color = "black"),
         axis.ticks = element_line(color = "black"),
         text=element_text(family="Helvetica"),
         plot.background = element_rect(fill = "transparent"))
int







#number of pfrs to report in manuscript 
length(unique(df$proteoform))


#Vlad's version of data
load("./output_data/spectral_counting_test_results.RData")

#Vlad's version of data
load("./output_data/shotgun_topdown_sc_20240730_modann_wres.RData")

# significant traits

filtered_results <- results %>%
   filter(adj.P.Val < 0.05) %>%
   distinct(proteoform) %>%
   inner_join(results)

nrow(filtered_results) # 819 # 378
length(unique(filtered_results$variable)) # 7 # 7

results %>%
   filter(adj.P.Val < 0.05) %>%
   group_by(variable) %>%
   tally() %>%
   arrange(-n)



filtered_results2 <- results %>%
   group_by(variable) %>%
   summarise(logFC_mean = mean(logFC, na.rm = TRUE),
             logFC_sd = sd(logFC, na.rm = TRUE)) %>%
   inner_join(filtered_results) %>%
   mutate(Zscore = (logFC - logFC_mean)/logFC_sd) %>%
   mutate(Zscore2 = case_when(Zscore > 3 ~ 3,
                              Zscore < -3 ~ -3,
                              TRUE ~ Zscore)) %>%
   mutate(is_sig = case_when(adj.P.Val < 0.05 ~ 1, TRUE ~ NA_real_))


# check collisions
collision_pairs <- fData(m)[unique(filtered_results2$proteoform),] %>%
   dplyr::filter(collision != "-") %>%
   select(proteoform_id, collision)

# and let's remove those that have collisions
filtered_results2 <- filtered_results2 %>%
   filter(!proteoform %in% collision_pairs$proteoform_id)



## Plotting --------------------------------------------------------------------

# Arrange proteoforms by gene and number
df <- filtered_results2 %>% 
   mutate(gene = sub("^(.*)_(\\d+)$", "\\1", proteoform),
          number =  as.numeric(sub("^(.*)_(\\d+)$", "\\2", proteoform))) %>% 
   arrange(gene, number) %>% 
   # Convert to factor to preserve order in the plot
   mutate(proteoform = factor(proteoform, levels = unique(proteoform)))


library(extrafont)
loadfonts()

df2 <- df

rankdf <- data.frame(proteoform = unique(df2$proteoform)) %>%
   separate(proteoform, c("Gene", "feature"), remove = FALSE) %>%
   arrange(Gene) %>%
   #dont do gene, feature for arrange
   mutate(rank = 1:length(unique(df2$proteoform))) %>%
   mutate(group = case_when(rank <= 44 ~ "low", 
                            rank > 88 ~ "high",
                            rank > 44 | rank < 88 ~ "mid"))

df3 <- merge(df2, rankdf, by = "proteoform")


df3 <- df3 %>%
   mutate(variable = gsub("sqrt_amyloid", "amyloid", variable)) %>%
   mutate(variable = gsub("sqrt_tangles", "tangles", variable))


df3 <- df3 %>%
   mutate(variable = ordered(variable, 
                             levels = c("amyloid", "tangles", "anye4", 
                                        "caa_4gp", "ci_num2_gct", 
                                        "hip_scl_3reg_yn", "cogn_global_lv",
                                        "cogng_demog_slope",
                                        "arteriol_scler",
                                        "lb_neo")))





aspect_ratio_sc <- length(unique(df2$variable)) / length(unique(df2$proteoform))

df_sc <- df3
#to make plot add int and + sc using patchwork 
sc <- ggplot(data = df3, 
             aes(y = proteoform, x = variable, 
                 fill = Zscore2, size = is_sig)) + 
   geom_tile(aes(fill = Zscore2, color=is_sig, width=0.9, height=0.9), size=1)+
   scale_colour_gradient( high = "black")+
   scale_fill_gradientn(colours = colorspace::lighten(gplots::bluered(100), amount = 0.5), 
                        # values = c(0, seq(qn01[1], qn01[2], length.out = 98),1), 
                        limits = c(-3, +3),
                        name = "Z-score") +
   guides(color = FALSE, size = FALSE)+
   geom_point( aes( y = proteoform, x = variable, size=ifelse(is_sig, "dot", "no_dot")))+
   scale_size_manual(values=c(dot=1.5, no_dot=NA), guide="none") +
   theme_classic()+
   scale_y_discrete(limits=rev)+
   scale_x_discrete(position="top",labels=c(
      "amyloid" = "amyloid",
      "tangles" = "tangles",
      "anye4" = expression(italic("APOE")~ε4), 
      "caa_4gp" = "cerebral amyloid angiopathy",
      "ci_num2_gct" = "macroscopic infarcts",
      "hip_scl_3reg_yn" = "hippocampal sclerosis",
      "cogn_global_lv" = "global cognition", 
      "cogng_demog_slope" = "slope of cognitive decline", 
      "arteriol_scler" = "arteriolosclerosis",
      "lb_neo" = "neocortical Lewy bodies")) +
   # coord_fixed(ratio = aspect_ratio_sc)+
   theme(aspect.ratio = 2.6,
         axis.text.x = element_text(angle = 35, vjust = 0.5, hjust = 0),
         legend.box = "horizontal", 
         legend.key = element_rect(color = "black"),
         axis.text = element_text(color = "black"),
         axis.ticks = element_line(color = "black"),
         text=element_text(family="Helvetica"),
         strip.text = element_blank(),
         plot.background = element_rect(fill = "transparent"))+
   #facet_wrap(~group, scales="free")+
   facet_wrap(~factor(group, levels=c('low', 'mid', 'high')), scales="free")


sc




sc2 <- ggplot(data = df3, 
             aes(y = proteoform, x = variable, 
                 fill = Zscore2, size = is_sig)) + 
   geom_tile(aes(fill = Zscore2, color=is_sig, width=0.9, height=0.9), size=1)+
   scale_colour_gradient( high = "black")+
   scale_fill_gradientn(colours = colorspace::lighten(gplots::bluered(100), amount = 0.5), 
                        # values = c(0, seq(qn01[1], qn01[2], length.out = 98),1), 
                        limits = c(-3, +3),
                        name = "Z-score") +
   guides(color = FALSE, size = FALSE)+
   geom_point( aes( y = proteoform, x = variable, size=ifelse(is_sig, "dot", "no_dot")))+
   scale_size_manual(values=c(dot=1.5, no_dot=NA), guide="none") +
   theme_classic()+
   scale_y_discrete(limits=rev)+
   scale_x_discrete(position="top",labels=c(
      "amyloid" = "amyloid",
      "tangles" = "tangles",
      "anye4" = expression(italic("APOE")~ε4), 
      "caa_4gp" = "cerebral amyloid angiopathy",
      "ci_num2_gct" = "macroscopic infarcts",
      "hip_scl_3reg_yn" = "hippocampal sclerosis",
      "cogn_global_lv" = "global cognition", 
      "cogng_demog_slope" = "slope of cognitive decline", 
      "arteriol_scler" = "arteriolosclerosis",
      "lb_neo" = "neocortical Lewy bodies")) +
   # coord_fixed(ratio = aspect_ratio_sc)+
   theme(aspect.ratio = 2.6*(43/44),
         axis.text.x = element_text(angle = 35, vjust = 0.5, hjust = 0),
         legend.box = "horizontal", 
         legend.key = element_rect(color = "black"),
         axis.text = element_text(color = "black"),
         axis.ticks = element_line(color = "black"),
         text=element_text(family="Helvetica"),
         strip.text = element_blank(),
         plot.background = element_rect(fill = "transparent"))+
   #facet_wrap(~group, scales="free")+
   facet_wrap(~factor(group, levels=c('low', 'mid', 'high')), scales="free")


sc2



length(unique(df3$proteoform))




# the plot below is just for a reference
library(patchwork)

#changes ratio of panels but at some limit doesn't change exact aspect ratio of individual plots
panelratio <- 4

#guides = "collect" removes extra legend

plot <- int + sc +
   plot_layout(guides = "collect", ncol = 2, nrow = 1, heights = c(1, 1, 1, 1), widths = c(1, panelratio)) &
   plot_annotation(tag_levels = "a") &
   theme(plot.tag = element_text(size=24,face = 'bold'))
plot









# saving pieces for assembing in affinity designer
SCALE = 0.75

ggsave(plot = int,
       filename = "./figures_tables/Fig3_intensity_punchcard.png",
       units = c("px"),
       bg = 'white',
       dpi = 300,
       scale = SCALE,
       width = 4388, # 4951 * 39/44,
       height = 2897 * 39/44) # 2567) # (39/44)

ggsave(plot = int,
       filename = "./figures_tables/Fig3_intensity_punchcard.pdf",
       units = c("px"),
       bg = 'white',
       dpi = 300,
       scale = SCALE,
       width = 4388, # 4951 * 39/44,
       height = 2897 * 39/44) # 2567) # (39/44)

#------------------------------------

ggsave(plot = sc,
       filename = "./figures_tables/Fig3_spectral_counts_punchcard.png",
       units = c("px"),
       bg = 'white',
       dpi = 300,
       scale = SCALE,
       width = 4951,
       height = 2897) # 44/44

ggsave(plot = sc,
       filename = "./figures_tables/Fig3_spectral_counts_punchcard.pdf",
       units = c("px"),
       bg = 'white',
       dpi = 300,
       scale = SCALE,
       width = 4951,
       height = 2897) # 44/44


#------------------------------------

# the last panel is shorter (43 vs 44 rows), thus plotting in different scale
ggsave(plot = sc2 ,
       filename = "./figures_tables/Fig3_spectral_counts_punchcard_43.png",
       units = c("px"),
       bg = 'white',
       dpi = 300,
       scale = SCALE,
       width = 4951, # * 43/44,
       height = 2897 * 43/44) # 43/44

ggsave(plot = sc2 ,
       filename = "./figures_tables/Fig3_spectral_counts_punchcard_43.pdf",
       units = c("px"),
       bg = 'white',
       dpi = 300,
       scale = SCALE,
       width = 4951, # * 43/44,
       height = 2897 * 43/44) # 43/44







