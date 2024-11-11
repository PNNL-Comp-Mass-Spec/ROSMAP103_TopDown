

library(MSnSet.utils)
library(tidyverse)
library(ggplot2)
library(ggsci)


load("./output_data/shotgun_topdown_int_20240730_modann_cnt_wres.RData")
td103 <- m

load("./source_data/processed_protein_level_msnset.RData")
emory <- m1

load("./source_data/final_r1r2r3_peptide_srm_data.RData")
srm <- m

# load("../ROSMAP103_4kDa_20231023/Emory_SRM_corrs/all_Emory_SRM_cors.RData")



# the confusing renaming of the MSnSet objects is because this code is a patchwork
# of pieces taken from multiple places
if(!file.exists("./output_data/TMT_SRM_cors.RData")){
   
   m1 <- emory
   m2 <- srm
   
   ### SRM data ###################################################################
   srm_representative_proteoform_ids <- fData(m2) %>%
      filter(!grepl("80", Peptide.Modified.Sequence)) %>% # remove phosphopeptides
      group_by(Gene) %>%
      slice_max(s2n) %>%
      pull(specie)
   m2 <- m2[srm_representative_proteoform_ids,]
   featureNames(m2) <- fData(m2)$Gene
   ################################################################################
   
   # intersecting samples
   common_samples <- intersect(sampleNames(m1), sampleNames(m2)) # 55
   common_features <- intersect(featureNames(m1), featureNames(m2)) # 55
   m1 <- m1[common_features, common_samples]
   m2 <- m2[common_features, common_samples]
   
   cors <- c()
   for(i in 1:nrow(m1)){
      cor_i <- cor(as.numeric(exprs(m1)[i,]), 
                   as.numeric(exprs(m2)[i,]),
                   use = "pairwise.complete.obs")
      cors <- c(cors, cor_i)
   }
   
   dafr <- data.frame(Gene = featureNames(m1), cor = cors)
   
   # appending variances to dafr
   dafr$var_tmt <- apply(exprs(m1), 1, var, na.rm = T)
   dafr$var_srm  <- apply(exprs(m2), 1, var, na.rm = T)
   
   save(dafr, file="./output_data/TMT_SRM_cors.RData")
}
(load("./output_data/TMT_SRM_cors.RData")) # dafr















# restrict to 50% presence
td103 <- td103[rowSums(!is.na(exprs(td103))) > 50,]

td103 <- remove_covariate(td103, "pmi")


selected_genes <- fData(td103) %>%
   inner_join(dafr) %>%
   group_by(Gene) %>%
   summarise(number_of_pforms = n(),
             cor = unique(cor)) %>%
   filter(cor > 0.35,                 # 0.35
          number_of_pforms > 1) %>%  # 3
   distinct(Gene) %>%
   pull()




# deriving cors for Emory
################################################################

m1 <- td103
m2 <- emory

# intersecting samples
common_samples <- intersect(sampleNames(m1), sampleNames(m2)) # 55
m1 <- m1[, common_samples]
m2 <- m2[, common_samples]

cors <- c()
pform <- c()
for(s_i in selected_genes){
   m1_s_i <- m1[fData(m1)$Gene == s_i,]
   for(i in 1:nrow(m1_s_i)){
      cor_i <- cor(as.numeric(exprs(m1_s_i)[i,]), 
                   as.numeric(exprs(m2)[s_i,]),
                   use = "pairwise.complete.obs")
      pform <- c(pform, featureNames(m1_s_i)[i])
      names(cor_i) <- s_i
      cors <- c(cors, cor_i)
   }
}

x_emory <- data.frame(Gene = names(cors), cor = as.numeric(cors), pform = pform)
################################################################



# deriving corsf for SRM
################################################################

m1 <- td103
# srm
srm_representative_proteoform_ids <- fData(srm) %>%
   filter(!grepl("80", Peptide.Modified.Sequence)) %>%
   group_by(Gene) %>%
   slice_max(s2n) %>%
   pull(specie)
srm <- srm[srm_representative_proteoform_ids,]
featureNames(srm) <- fData(srm)$Gene
m2 <- srm


# intersecting samples
common_samples <- intersect(sampleNames(m1), sampleNames(m2)) # 55
m1 <- m1[, common_samples]
m2 <- m2[, common_samples]

cors <- c()
pform <- c()
for(s_i in selected_genes){
   m1_s_i <- m1[fData(m1)$Gene == s_i,]
   for(i in 1:nrow(m1_s_i)){
      cor_i <- cor(as.numeric(exprs(m1_s_i)[i,]), 
                   as.numeric(exprs(m2)[s_i,]),
                   use = "pairwise.complete.obs")
      pform <- c(pform, featureNames(m1_s_i)[i])
      names(cor_i) <- s_i
      cors <- c(cors, cor_i)
   }
}

x_srm <- data.frame(Gene = names(cors), cor = as.numeric(cors), pform = pform)
################################################################



# prep data
################################################################
y <- dafr %>%
   filter(Gene %in% selected_genes) %>%
   arrange(-cor) %>%
   mutate(Gene = ordered(Gene, levels = Gene))

x_emory$data <- "TMT"
x_srm$data <- "SRM"
x <- bind_rows(x_emory, x_srm)

x <- group_by(x, pform) %>%
   summarise(mean_cor = mean(cor, na.rm = TRUE)) %>%
   ungroup() %>%
   inner_join(x)

x <- arrange(x, mean_cor)
x$pform <- ordered(x$pform, levels = unique(x$pform))
x$Gene <- ordered(x$Gene, levels = levels(y$Gene))
################################################################




# x <- filter(x, Gene == "SNCA")
# x$Gene <- as.character(x$Gene)
xw <- pivot_wider(x, names_from = "data", values_from = "cor")
p <- ggplot(x) +
   aes(x = cor, y = pform, color=data) +
   # scale_color_npg() +
   scale_color_manual(values = c("TMT" = "#00BCD4", "SRM" = "#E53935"), name = "Data Type") +
   geom_segment(aes(x=TMT, xend = SRM, y=pform, yend=pform), color="lightgrey", data=xw) +
   geom_point(size = 2) +
   # geom_point(aes(x = mean_cor), color="darkgrey") +
   geom_vline(xintercept = 0, color = "#78909C33", linewidth = 1) +
   geom_vline(aes(xintercept = cor), color = "#78909C55", linewidth = 1.2, data = y, linetype = "dashed") +
   facet_wrap(~Gene, scales = "free_y", nrow = 2) +
   xlim(-1,+1) +
   theme_bw(base_size = 18) +
   theme(axis.text.y = element_blank(),
         axis.ticks.y = element_blank(),
         panel.grid = element_blank(),
         strip.background = element_rect(fill = "#ECEFF1",
                                         color = "#78909C"),
         panel.border = element_rect(fill = NA, color = "#78909C")) +
   xlab("Correlation Coefficient") +
   ylab("Rank of Correlation")

plot(p)

ggsave(plot = p, path = "figures_tables", filename = "Fig2_correlations_with_prior.png", width = 12, height = 6, scale = 1.2)





















