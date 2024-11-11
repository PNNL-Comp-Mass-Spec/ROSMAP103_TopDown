


library(tidyverse)
library(MSnSet.utils)
library(msmsTests)

load("./output_data/shotgun_topdown_sc_20240730_modann.RData")

# apply > 50 filter out of 103
# m <- m[rowSums(!is.na(exprs(m))) > 50,]

# selecting features
# subset to at least 3 batches with 2 samples
selected_features <- m %>%
   exprs() %>%
   as.data.frame() %>%
   rownames_to_column("feature_name") %>%
   pivot_longer(cols = -feature_name, names_to = "ProjID", values_to = "spectral_counts") %>%
   inner_join(pData(m)) %>%
   select(feature_name, ProjID, spectral_counts, batch) %>%
   filter(spectral_counts > 0) %>%
   group_by(feature_name, batch) %>%
   tally() %>%
   filter(n >= 2) %>%
   group_by(feature_name) %>%
   tally() %>%
   filter(n >= 3) %>%
   pull(feature_name)


m <- m[selected_features,]
m$batch <- as.factor(m$batch)



#' # Overall test results across clinical and neuropath traits
get_sig_res_for_counts <- function(m, pheno){
   m <- m[,!is.na(m[[pheno]])]
   alt.f <- sprintf("y ~ %s + batch + pmi + 1", pheno)
   null.f <- "y ~ 1 + batch + pmi"
   div <- colSums(exprs(m)) # normalization factor
   res <- msms.glm.qlll(m, alt.f, null.f, div=div, facs = pData(m))
   res$p.val.adj <- p.adjust(res$p.value, "BH")
   res <- res %>%
      rownames_to_column("proteoform") %>%
      dplyr::rename(logFC = LogFC) %>%
      dplyr::rename(P.Value = p.value) %>%
      dplyr::rename(adj.P.Val = p.val.adj) %>%
      select(proteoform, logFC, P.Value, adj.P.Val)
   return(res)
}

# vl_2_remove <- c("ProjID", "study", "scaled_to", "amyloid", "tangles", 
#                  "tangles_mf", "pmi", "race", "nlw", "spanish", "batch")
# vls <- setdiff(varLabels(m), vl_2_remove)
vls <- c("gpath", "cogng_demog_slope", "anye2", "anye4")
results <- list()
for(vl in vls){
   res_i <- get_sig_res_for_counts(m, vl)
   res_i$variable <- vl
   results <- c(results, list(res_i))
}
results <- bind_rows(results)



sr_long_pv <- pivot_wider(results, id_cols = proteoform, names_from = variable, values_from = P.Value)
sr_long_apv <- pivot_wider(results, id_cols = proteoform, names_from = variable, values_from = adj.P.Val)
sr_long_eff_size <- pivot_wider(results, id_cols = proteoform, names_from = variable, values_from = logFC)

# constrain to tested proteins only
m <- m[sr_long_pv$proteoform,]

# add significant results to featureData
names(sr_long_pv) <- paste(names(sr_long_pv), ".pv", sep = "")
names(sr_long_apv) <- paste(names(sr_long_apv), ".apv", sep = "")
names(sr_long_eff_size) <- paste(names(sr_long_eff_size), ".eff_size", sep = "")
sr_long <- inner_join(sr_long_pv, sr_long_apv, by=c("proteoform.pv" = "proteoform.apv")) %>%
   inner_join(sr_long_eff_size, by=c("proteoform.pv" = "proteoform.eff_size"))
fData(m) <- cbind(fData(m), sr_long[,-1]) # no need for another proteoform ID

save(m, file="./output_data/shotgun_topdown_sc_20240730_modann_wres.RData")


