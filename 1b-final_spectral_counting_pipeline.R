#' ---
#' title: "Preprocessing of Top-Down Data for 103 Human Subjects"
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
library(TopPICR) # remotes::install_github("evanamartin/TopPICR", ref="0.0.3")
library(memoise)
library(MSnbase)
library(MSnSet.utils)
library(PNNL.DMS.utils)


load("./output_data/pre_identifications.RData")
load("./output_data/pre_features.RData")


################################################################################################################################
# Create cross-tab with counts

x <- x_grp %>% ungroup()
pttrn <- ".*Alz_TD_(\\d{8})_([ABC])_(\\d+)_(AB|WO)_.*"
x <- x %>%
   mutate(sample_name = sub(pttrn,"X\\1",Dataset),
          batch = sub(pttrn,"B\\3",Dataset),
          fraction = sub(pttrn,"\\4",Dataset)) %>%
   mutate(feature_name = paste(Gene, pcGroup, sep = "_")) %>%
   select(-Dataset)


x_expr <- x %>%
   select(feature_name, sample_name, `Scan(s)`) %>%
   pivot_wider(values_from = `Scan(s)`,
               names_from = sample_name,
               values_fn = length,
               values_fill = 0) %>%
   as.data.frame() %>%
   {rownames(.) <- .$feature_name;.} %>%
   select(-feature_name) %>%
   as.matrix()


# features
x_feat <- x %>%
   group_by(Gene, pcGroup, feature_name) %>%
   dplyr::summarize(count = n(), .groups = "keep") %>%
   ungroup() %>%
   left_join(x_meta, by=c("Gene","pcGroup")) %>%
   mutate(proteoform_id = paste(Gene, pcGroup, sep="_")) %>%
   as.data.frame() %>%
   {rownames(.) <- .$feature_name;.}


# phenodata
x_pheno <- x %>%
   distinct(sample_name, batch) %>%
   as.data.frame() %>%
   {rownames(.) <- .$sample_name;.}

m <- MSnSet(x_expr, x_feat[rownames(x_expr),], x_pheno[colnames(x_expr),]) #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
################################################################################################################################





#' # Linking with clinical data
#'
#+ linking with clinical data, echo=TRUE
load("./source_data/clinical_data.RData")
m <- MSnSet(exprs = exprs(m),
            fData = fData(m),
            pData = pall[sampleNames(m),])

fData(m)$spectralCount <- NULL

#' # Saving MSnSet
save(m, file = "./output_data/shotgun_topdown_sc_20240730.RData")







