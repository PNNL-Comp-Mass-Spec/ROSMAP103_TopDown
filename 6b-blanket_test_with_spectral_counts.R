#' ---
#' title: "Spectral Counting Statistical Analysis"
#' output: 
#'      BiocStyle::html_document:
#'          toc: true
#'          number_sections: true
#'          highlight: haddock
#' ---

#+ echo=F
knitr::opts_chunk$set(echo=F, message=F, warning=F, fig.align='center')
#, out.width='10cm')



#+ libraries
library(tidyverse)
library(MSnSet.utils)

#+ loading MSnSet
load("./output_data/shotgun_topdown_sc_20240730_modann_wres.RData")


### Filter works, but please be careful. The data is non-imputed and thus sparse.

#' # Original numbers Proteoforms/Genes
nrow(m)
featureNames(m) %>% sub("([^_]*)_\\d+","\\1",.) %>% unique() %>% length()

# apply > 50 filter
# m <- m[rowSums(exprs(m) > 0) > 10,] # 70 out of 103

# alternative filter. Better suited is the data is non-batch-corrected and the 
# batch is used as a covariate.
# 
# subset to at least 2 batches with 2 samples
selected_features <- m %>%
   exprs() %>%
   as.data.frame() %>%
   rownames_to_column("feature_name") %>%
   pivot_longer(cols = -feature_name, names_to = "ProjID", values_to = "spectral_counts") %>%
   inner_join(pData(m)) %>%
   dplyr::select(feature_name, ProjID, spectral_counts, batch) %>%
   filter(spectral_counts > 0) %>%
   group_by(feature_name, batch) %>%  # filter below retains features/batches with >= given number of counts
   tally() %>%
   filter(n >= 2) %>%
   group_by(feature_name) %>%  # filter below retains features present across >= given number of batches
   tally() %>%
   filter(n >= 3) %>%
   pull(feature_name)


m <- m[selected_features,]




library(msmsTests)

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
   res_i <- get_sig_res_for_counts(m, vl)
   res_i$variable <- vl
   results <- c(results, list(res_i))
}
results <- bind_rows(results)

save(results, file = "./output_data/spectral_counting_test_results.RData")

