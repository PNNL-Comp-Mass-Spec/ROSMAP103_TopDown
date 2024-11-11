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


# memoisation of the slow steps
#check if memoised is.memoised(match_features)
# match_features <- memoise(match_features, cache = cache_filesystem("./caching"))
# correct_batch_effect_NA <- memoise(correct_batch_effect_NA, cache = cache_filesystem("./caching"))
# rrollup <- memoise(rrollup, cache = cache_filesystem("./caching"))


load("./output_data/pre_identifications.RData")
load("./output_data/pre_features.RData")

#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!
#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!
# Unidentified feature steps
#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!
#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!

# Align unidentified features ---------------

# Align the unidentified feature retention times with the model created by the
# identified feature retention times.
feat_art <- align_rt(
   x = feat,
   model = the_model,
   var_name = "Time_apex")

# Recalibrate unidentified feature mass ---------------

# Recalibrate the unidentified feature masses with the model created by the
# identified feature masses.
feat_rcm <- recalibrate_mass(
  x = feat_art,
  errors = the_error,
  var_name = "Mass")




#



# Retrieve unidentified features ---------------
#takes forever 
feat_rtv <- match_features(ms2 = x_grp,
                           ms1 = feat_rcm,
                           errors = the_error,
                           n_mme_sd = 4/the_error$ppm_sd,
                           n_rt_sd = 150/the_error$rt_sd, 
                           cores = 60
                           # ppm_cutoff = 4,
                           # rt_cutoff = 150
                           )



#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!
#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!
# Creating MSnSet
#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!
#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!


#' # Making MSnSet out of fused

#+ converting to MSnSet
x <- feat_rtv

pttrn <- ".*Alz_TD_(\\d{8})_([ABC])_(\\d+)_(AB|WO)_.*"
x <- x %>%
  mutate(sample_name = sub(pttrn,"X\\1",Dataset),
         batch = sub(pttrn,"B\\3",Dataset),
         fraction = sub(pttrn,"\\4",Dataset)) %>%
  mutate(feature_name = paste(Gene, pcGroup, CV, fraction, sep = "_")) %>%
  # select(-Dataset, -RTapex)
  select(-Dataset)

# expression
x_expr <- x %>%
  pivot_wider(id_cols = "feature_name",
              names_from = "sample_name",
              values_from = "Intensity") %>%
  as.data.frame() %>%
  {rownames(.) <- .$feature_name;.} %>%
  select(-feature_name) %>%
  as.matrix()



# features
x_feat <- x %>%
  group_by(Gene, pcGroup, CV, fraction, feature_name) %>%
   dplyr::summarize(median_intensity = median(Intensity),
            count = n(),
            .groups = "keep") %>%
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

m <- MSnSet(x_expr, x_feat[rownames(x_expr),], x_pheno[colnames(x_expr),])


# used for FigS2
save(m, file="./output_data/shotgun_topdown_int_20240730_beforecorrections_mid1a.RData")



#' # Preprocessing of MSnSet

#+ transform
m <- log2_zero_center(m)

#+ normalization
#' # Normalization
#'
m <- normalizeByGlob(m)


#' ## Trends before
#+ echo=T
evaluate_parameter_effect(m,
                          sample = "X04879843",
                          property = "rt",
                          partition_by = "fraction",
                          partition_value = "WO",
                          span = 1)

evaluate_parameter_effect(m,
                          sample = "X04879843",
                          property = "rt",
                          partition_by = "fraction",
                          partition_value = "AB",
                          span = 1)

#' ## Normalization by trend
m <- normalize_by_feature_property_partitioned(m,
                                               property = "rt",
                                               partition_by = "fraction",
                                               method="loess", span=1)

#' ## Trends after
#+ echo=T
evaluate_parameter_effect(m,
                          sample = "X04879843",
                          property = "rt",
                          partition_by = "fraction",
                          partition_value = "WO",
                          span = 1)

evaluate_parameter_effect(m,
                          sample = "X04879843",
                          property = "rt",
                          partition_by = "fraction",
                          partition_value = "AB",
                          span = 1)










#'
#' Batch effect after non-parametric ComBat. 13 minutes
#'
#+ batch correction, echo = TRUE, cache = TRUE
t1 <- Sys.time()
m <- correct_batch_effect_NA(m, "batch", par.prior = FALSE)
t2 <- Sys.time()
print(t2 - t1)




#'
#' # Roll-Up from Individual Measurement to Proteoforms
#'
#'
#' ## What did we measure? What are the features in the cross-tab? What are proteoforms?
#'
#' We tracked the intensities of gene/cluster groups within each fraction and CV value.
#' What is the intensity of a cluster group needs to be explained.
#' Features in the cross-tab are the these Gene/pcGroup/fraction/CV combinations.
#' Proteoforms are Gene/pcGroup combinations.
#'
#' We need to work on explanation of the terminology.
#'
#'
#+ rollup, echo=TRUE, cache=TRUE
m0 <- m
m <- rrollup(m, "proteoform_id", rollFun = "-", verbose = FALSE)

#' recover proteoform info
x_meta <- x_meta %>%
   mutate(proteoform_id = paste(Gene, pcGroup, sep="_")) %>%
   as.data.frame() %>%
   {rownames(.) <- .$proteoform_id;.}
fData(m) <- x_meta[featureNames(m),]



#' # Remove Proteoforms not present in at least 10 samples
nrow(m)
m <- m[rowSums(!is.na(exprs(m))) >= 10,]

nrow(m)

featureNames(m) %>% sub("([^_]*)_\\d+","\\1",.) %>% unique() %>% length()


#' # Linking with clinical data
#'
#+ linking with clinical data, echo=TRUE
load("./source_data/clinical_data.RData")
m <- MSnSet(exprs = exprs(m),
            fData = fData(m),
            pData = pall[sampleNames(m),])


#' # Saving MSnSet
save(m, file = "./output_data/shotgun_topdown_int_20240730.RData")






