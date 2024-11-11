
library(tidyverse)
library(MSnSet.utils)

load("./output_data/shotgun_topdown_int_20240730_modann_cnt_wres.RData")

# apply > 50 filter
m <- m[rowSums(!is.na(exprs(m))) > 50,] # 70 out of 103

#' # Overall test results across clinical and neuropath traits
get_sig_res <- function(m, pheno){
   modl <- sprintf("~ %s + pmi", pheno)
   res <- limma_gen(m, modl, pheno) %>%
      rownames_to_column("proteoform") %>%
      dplyr::select(proteoform, logFC, P.Value, adj.P.Val)
   return(res)
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


results <- list()
for(vl in vls){
   res_i <- get_sig_res(m, vl)
   res_i$variable <- vl
   results <- c(results, list(res_i))
}
results <- bind_rows(results)

save(results, file = "./output_data/intensity_test_results.RData")

