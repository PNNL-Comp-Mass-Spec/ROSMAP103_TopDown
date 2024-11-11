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


# just in case
if(!dir.exists("output_data"))
   dir.create("output_data")



# # memoisation of the slow steps
# augment_annotation <- memoise(augment_annotation, cache = cache_filesystem("./caching"))
# form_model <- memoise(form_model, cache = cache_filesystem("./caching"))
# align_rt <- memoise(align_rt, cache = cache_filesystem("./caching"))
# calc_error <- memoise(calc_error, cache = cache_filesystem("./caching"))
# cluster <- memoise(cluster, cache = cache_filesystem("./caching"))
# create_pcg <- memoise(create_pcg, cache = cache_filesystem("./caching"))
# create_mdata <- memoise(create_mdata, cache = cache_filesystem("./caching"))


if("toppic_output.RData" %in% list.files(path = "./source_data/")){
   load("./source_data/toppic_output.RData")
}else{
   #data set with max mod = 4000 Da
   toppic_output <- read_TopPIC_DMS(data_package_num = 5149)
   save(toppic_output, file = "./source_data/toppic_output.RData")
}


ids <- toppic_output$ms2identifications
feat <- toppic_output$ms1features


# remove NA annotations. We can't handle them at the FDR filter tuning step
# Those are non-uniprot IDs, like short ORFs and contaminants
ids <- ids %>%
  dplyr::filter(!is.na(AnnType))


# add CVs
# Note, CVs are encoded in the dataset by consecutive numbering.
# I believe there is no table converting them from consecutive number to CV.
# In this case it has been sorted out by mass ranges. Also, the conversion was
# consistent across the datasets. In the future, we need a non-manual approach
# for mapping the suffixes to CVs.
cv_codes <- c("0" = "-50", "1" = "-40", "2" = "-30")
ids <- ids %>%
   mutate(CV = as.character(cv_codes[sub(".*(\\d)$","\\1",Dataset)])) %>%
   mutate(Dataset = sub("(.*)_\\d$","\\1",Dataset))
feat <- feat %>%
   mutate(CV = as.character(cv_codes[sub(".*(\\d)$","\\1",Dataset)])) %>%
   mutate(Dataset = sub("(.*)_\\d$","\\1",Dataset))


#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!
#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!
# Identified feature steps
#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!
#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!




# Remove erroneous genes ---------------
# What does this do? In case a proteoform assigned to multiple genes,
# this function selects the gene with the lowest E-value.
x_err <- rm_false_gene(ids) #rename !!!


# getting protLength (parent protein length)
# No need to match. Just need to link accessions.
# 1. read fasta and 
# 2. link by `Protein accession` or UniProtAcc

if(!file.exists("./source_data/ID_008032_8627C6BD.fasta")){
   file.copy(file.path(r"(\\gigasax\DMS_FASTA_File_Archive\Dynamic\Forward\)",
                       "ID_008032_8627C6BD.fasta"), "./source_data/")
}
fst <- Biostrings::readAAStringSet("./source_data/ID_008032_8627C6BD.fasta")



add_protein_length <- function(x, fst_obj){
  fst_df <- data.frame(`Protein accession` = names(fst_obj),
                       protLength = BiocGenerics::width(fst_obj), 
                       row.names = NULL, 
                       check.names = FALSE) %>%
    mutate(`Protein accession` = word(`Protein accession`))
  x <- inner_join(x, fst_df)
  return(x)
}


x_err <- add_protein_length(x_err, fst)

compute_fdr(x_err)


# Filter by cleanSeq count ---------------

# really minimal filter
x_fbc <- filter_by_count(x = x_err,
                         count_within = c("Dataset", "Scan(s)"),
                         count = "cleanSeq",
                         threshold = 2) # better to switch to 3 and drop the E-value filter altogether.


compute_fdr(x_fbc)

# FDR control ---------------

# First find the E-value threshold for each of the three annotation types.
# FDR is calculated at Gene level.
the_cutoff <- find_evalue_cutoff(x = x_fbc,
                                 fdr_threshold = 0.01)

# Apply the actual FDR control/filter with the threshold values from above.
x_fdr <- apply_evalue_cutoff(x = x_fbc,
                             e_vals = the_cutoff)

compute_fdr(x_fdr)


# Proteoform inference ---------------
x_ipf <- infer_prot(x_fdr)

compute_fdr(x_ipf)


# I think this concludes polishing the MS/MS identifications



# Align retention time ---------------

# Create the model that will be used to align each retention time. The model is
# created between a reference data set and all other data sets.
the_model <- form_model(
  x = x_ipf,
  ref_ds = find_ref_ds(x = x_ipf),
  control = loess.control(surface = "direct"), # Use direct to avoid NAs.
  span = 0.5,
  family = "symmetric")


# Align retention times according to the model created previously.
x_art <- align_rt(
  x = x_ipf,
  model = the_model,
  var_name = "Feature apex")

# Recalibrate the mass ---------------

# Calculate the error between the `Precursor mass` and the `Adjusted precursor
# mass`. This acts as the model for recalibrating the mass. A reference data set
# is not used when calculating the mass error model.
the_error <- calc_error(x = x_art,
                        ref_ds = find_ref_ds(x = x_ipf))


# Recalibrate the mass according the the errors computed previously.
x_rcm <- recalibrate_mass(x = x_art,
                          errors = the_error,
                          var_name = "Precursor mass")



# Cluster ---------------
# Key step. But doesn't take too long.
x_clu <- cluster(x = x_rcm,
                 errors = the_error,
                 method = "single",
                 height = 1.5,
                 min_size = 3)


# Group clusters ---------------
x_grp <- create_pcg(x = x_clu,
                      errors = the_error,
                      n_mme_sd = 4/the_error$ppm_sd, # 5.4 st.devs
                      # ppm_cutoff = 4,
                      n_Da = 4,
                      n_rt_sd = 1.5)


# Create metadata ---------------
# Summary at Gene/pcGroup combo level
# NOTE: The number of distinct Gene_pcGroups can be different between the x_meta
# object and the fused object. This happens because there are some Gene_clusters
# that are present in the identified data but do not exist in the unidentified
# data.
# The object created here is not used in this script. This object will become
# part of the MSnSet created in another script. This object is also not needed
# for any of the functions that manipulate the feature data.
x_meta <- create_mdata(x = x_grp,
                       errors = the_error,
                       n_mme_sd = 4/the_error$ppm_sd,
                       n_rt_sd = 150/the_error$rt_sd
                       # ppm_cutoff = 4,
                       # rt_cutoff = 150
                       )


# Maybe this where the first break supposed to be.
# Then we'll split the pipeline into two: intensity and spectral counts-based.
# What do I need to save/retain?
# SPC: x_grp, x_meta
# INT: feat, the_model, the_error


save(x_grp, x_meta, file="./output_data/pre_identifications.RData")
save(feat, the_model, the_error, file="./output_data/pre_features.RData")


