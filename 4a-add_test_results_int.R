
library(tidyverse)
library(MSnSet.utils)


load("./output_data/shotgun_topdown_int_20240730_modann_cnt.RData")

# apply > 50 filter out of 103
m <- m[rowSums(!is.na(exprs(m))) > 50,]

#' # Overall test results across clinical and neuropath traits
get_sig_res <- function(m, pheno){
   modl <- sprintf("~ %s + pmi", pheno)
   res <- limma_gen(m, modl, pheno) %>%
      rownames_to_column("proteoform") %>%
      select(proteoform, logFC, P.Value, adj.P.Val, starts_with("SE"))
   colnames(res)[grepl("^SE", colnames(res))] <- "SE"
   return(res)
}

vl_2_remove <- c("ProjID", "study", "scaled_to", "amyloid", "tangles", 
                 "tangles_mf", "pmi", "race", "nlw", "spanish", "batch")
vls <- setdiff(varLabels(m), vl_2_remove)
results <- list()
for(vl in vls){
   res_i <- get_sig_res(m, vl)
   res_i$variable <- vl
   results <- c(results, list(res_i))
}
results <- bind_rows(results)




# Linking with stat results
results <- filter(results, variable %in% c("gpath", "cogng_demog_slope", "anye2", "anye4", "sqrt_amyloid"))
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

save(m, file="./output_data/shotgun_topdown_int_20240730_modann_cnt_wres.RData")


